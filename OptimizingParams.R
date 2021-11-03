library(tidyverse)
library(readxl)
library(pracma)
library(NMOF)
library(doParallel)
setwd("~/Desktop/FFx for Henry/")
cores = registerDoParallel(cores = 4)
paths = c("July1_4302FFx_115944W.xlsx", "Aug2_4225FFx_120742W.xlsx",
  "June3_4231FFx_150939W.xlsx",
  "Aug2_4225FFx_123257S.xlsx",	"June3_4231FFx_160955S.xlsx",
  "Oct13_4226FFx_115758W.xlsx",
  "July16_4223FFx_130633W.xlsx"	,"Oct13_4226FFx_122902S.xlsx",
  "July16_4223FFx_145147S.xlsx"	,"Oct7_4254FFx_101225W.xlsx",
  "Oct7_4254FFx_122012S.xlsx",
  "July1_4302FFx_132412S.xlsx"	,"Sept1_4255FFx_112451W.xlsx",
  "July2_4300FFx_140029W.xlsx"	,"Sept1_4255FFx_133441S.xlsx",
  "July2_4300FFx_150500S.xlsx")


euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2, na.rm = T))
get_results = function(x, data, labels, rmcell) {
  # method for doing the analysis on an xlsx file
  time_cutoff = x[1L]
  thresh = x[2L]
  Time = unname(unlist(data[,1]))
  cells = select(data, -1)
  cells = cells[, -rmcell]
  cells = as.data.frame(apply(cells, 2, as.numeric))
  get_area = function(col){
    my_min = min(col)
    col = col - my_min
    my_max = max(col)
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
    
    
    if (length(ends)==0) {
      return(c(0.0))
      
    }else{
      
      results = c(  
        #cum_area/length(ends),                        # average spike area
        length(ends)                                 # number of spikes
        # max_spike_area,                               # max_spike_area
        # longest_time,                                 # longest time
        # mean_rising_spike_duration,
        # mean_falling_spike_duration,
        # thresh*variance+baseline,      # threshold
        # mean_exp_coeff                               # average coeff from exp linear regression
        # # mean_lin_coeff                                # avg coeff from simple linear regression
        
      )
      return(results)
    }
  }
  results = mapply(get_area, as.data.frame(cells))
  x = euc.dist(labels, unname(unlist(results)))
  return(x)
}

dummy = function(x){
  score = foreach(i=1:length(paths), .combine = '+') %dopar% {
    labels <- read_excel(paths[i], sheet = 2)
    labels = as.data.frame(apply(labels, 2, as.numeric))
    rmcell = which(is.na(labels))

    labels = labels[-rmcell]
    data = read_excel(paths[i], skip = 1 )
    
    get_results(x, data, labels, rmcell)
  }
  return(score)
}
levels <- list(time_cutoff = seq(5, 35, by=.5), thresh = seq(1, 4, by= 0.1))
res <- gridSearch(dummy, levels)
