library(tidyverse)
library(readxl)
library(pracma)
library(NMOF)
library(doParallel)
setwd("~/Desktop/FFx for Henry/")
cores = registerDoParallel(cores = 4)
get_data = function(file){
  labels <- read_excel(file, sheet = 2)
  labels = as.data.frame(apply(labels, 2, as.numeric))
  rmcell = which(is.na(labels))
  labels = labels[-rmcell,]
  data = read_excel(file, skip = 1 )
  return(list(data, labels, rmcell))
}
all_paths = c("July1_4302FFx_115944W.xlsx",
  "June3_4231FFx_150939W.xlsx",	"June3_4231FFx_160955S.xlsx",
  "Oct13_4226FFx_115758W.xlsx",
  "July16_4223FFx_130633W.xlsx"	,"Oct13_4226FFx_122902S.xlsx",
  "July16_4223FFx_145147S.xlsx"	,"Oct7_4254FFx_101225W.xlsx",
  "Oct7_4254FFx_122012S.xlsx",
  "July1_4302FFx_132412S.xlsx"	,"Sept1_4255FFx_112451W.xlsx",
  "July2_4300FFx_140029W.xlsx"	,"Sept1_4255FFx_133441S.xlsx",
  "July2_4300FFx_150500S.xlsx")

# paths = all_paths[grep("S", all_paths, value=F)]
# euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2, na.rm = T))
get_cor = function(x1, x2){
  return(cor(x1, x2))
}
get_results = function(x, data, labels, rmcell) {
  
  # method for doing the analysis on an xlsx file
  time_cutoff = unlist(x[1])
  thresh = unlist(x[2])
  Time = unname(unlist(data[,1]))
  cells = select(data, -1)
  cells = as.data.frame(apply(cells, 2, as.numeric))
  cells = drop_na(cells)
  cells = cells[, -rmcell]
  get_area = function(col){
    my_min = min(col, na.rm = T)
    col = col - my_min
    my_max = max(col, na.rm = T)
    #max_amp = stats[5]
    col = col / my_max
    baseline = median(unlist(col))   # get median of all cell data
    variance = mad(unlist(col))  
    all_spike_locations = (col > thresh*variance + baseline)
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
  x = get_cor(labels, results)+1
  return(-1*x)
}

dummy = function(x){
  score = foreach(i=1:length(all_paths), .combine = '+') %dopar% {
    labels <- read_excel(all_paths[i], sheet = 2)
    labels = as.data.frame(apply(labels, 2, as.numeric))
    rmcell = which(is.na(labels))

    labels = labels[-rmcell,]
    data = read_excel(all_paths[i], skip = 1 )
    get_results(x, data, labels, rmcell)
  }
  return(score/length(all_paths))
}
levels <- list(time_cutoff = seq(15, 32, by=1), thresh = seq(2, 4, by= 0.5))
res <- gridSearch(dummy, levels)
