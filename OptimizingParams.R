library(tidyverse)
library(readxl)
library(pracma)
library(NMOF)
library(doParallel)
setwd("~/Desktop/FFx for Henry copy/")
cores = registerDoParallel(cores = 4)

all_paths = list.files(pattern = '[SW].xlsx')
all_paths = sort(all_paths); all_paths = all_paths[-c(1, 2)] #august doesnt have reliable scores

get_cor = function(x1, x2){
  return(cor(x1, x2, "complete.obs"))
}
get_results = function(x, my_data, labels) {
  
  # method for doing the analysis on an xlsx file
  time_cutoff = unlist(x[1])
  thresh = unlist(x[2])
  Time = unname(unlist(my_data[,1]))
  cells = select(my_data, -1)
  #cells = as.data.frame(apply(cells, 2, as.numeric))
  #cells = drop_na(cells)
  #cells = cells[, -rmcell]
  get_area = function(col){
    #my_min = min(col, na.rm = T)
    #col = col - my_min
    #my_max = max(col, na.rm = T)
    #max_amp = stats[5]
    baseline = median(unlist(col))   # get median of all cell data
    variance = mad(unlist(col))  
    threshold_line = baseline + thresh*variance
    col = col - threshold_line
    #col = col / my_max
    
    all_spike_locations = (col > 0)
    ## method from https://masterr.org/r/how-to-find-consecutive-repeats-in-r/
    rle.all_spike_locations = rle(as.numeric(all_spike_locations))
    myruns = which(rle.all_spike_locations$values == TRUE & rle.all_spike_locations$lengths >= time_cutoff)
    runs.lengths.cumsum = cumsum(rle.all_spike_locations$lengths)
    ends = runs.lengths.cumsum[myruns]
    newindex = ifelse(myruns>1, myruns-1, 0)
    starts = runs.lengths.cumsum[newindex] + 1   #starts indices of runs > 10 in length
    
    if (0 %in% newindex) starts = c(1,starts)
    results = length(ends)                                 # number of spikes

      return(results)
    
  }
  results = sapply(as.data.frame(cells), get_area )

  results = as.numeric(unname(unlist(results)))
  #View(results)
  #View(labels)
  x = get_cor(labels, results)+1 #add one to keep numbers (0, 1)
  return(-1*x)
}


labels = lapply(all_paths, function(str){
  labels <- read_excel(str, sheet = 2)
  labels = as.data.frame(apply(labels, 2, as.numeric))
  #rmcell = which(is.na(labels))
  
  #labels = labels[-rmcell,]

})
all_data = lapply(all_paths, function(str){
   return(as_tibble(read_excel(str, skip = 1)))
})

dummy = function(x){
  score = foreach(i=1:length(all_paths), .combine = '+') %dopar% {
    score = get_results(x, all_data[[i]], labels[[i]] )
    if(!is.na(score)){
      return(score)
      }
    else{
          return(0)
    }
  }
  return(score/length(all_paths))
}
levels <- list(time_cutoff = seq(15, 32, by=1), thresh = seq(2, 4, by= 0.5))
res <- gridSearch(dummy, levels)





