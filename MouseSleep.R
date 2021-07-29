setwd("~/Documents/R/MouseSleep/Analysis")
library(tidyverse)
library(readxl)
library(pracma)
data = read_excel("../rested.xlsx", skip = 1 )
data = drop_na(data)
data[,-1] = data[,-1]-min(data[,-1])
Time = unname(unlist(data[,1]))
thresholdFactor = 1


cells = select(data, -1)

cell = cells$undecided...6

get_area = function(col){
  baseline = median(unlist(col))   # get median of all cell data
  variance = mad(unlist(col))      # measure of variance (adds 1.4826 factor)
  #col_minus_baseline = (col - baseline)
  all_spike_locations = (col > thresholdFactor*variance + baseline)
  
  ## method from https://masterr.org/r/how-to-find-consecutive-repeats-in-r/
  rle.all_spike_locations = rle(as.numeric(all_spike_locations))
  myruns = which(rle.all_spike_locations$values == TRUE & rle.all_spike_locations$lengths >= 10)  #dont hardcode
  runs.lengths.cumsum = cumsum(rle.all_spike_locations$lengths)
  ends = runs.lengths.cumsum[myruns]
  newindex = ifelse(myruns>1, myruns-1, 0)
  starts = runs.lengths.cumsum[newindex] + 1   #starts indices of runs > 10 in length
  if (0 %in% newindex) starts = c(1,starts)
 
  
  cum_area = 0
  if (length(ends)==0) {
    return(c(0.0,0.0,0.0,0.0,3*baseline, 0.0, 0.0))
    
  }else{
    max_spike_area = 0.0
    longest_time = 0.0
    lin_coeff = 0.0
    exp_coeff = 0.0
    for(i in 1:length(ends)){                              # for every spike...
      time = Time[starts[i]:ends[i]]                       # time interval for spike_i
      if (Time[ends[i]]-Time[starts[i]] > longest_time){   # identify longest spike-time interval
        longest_time = Time[ends[i]]-Time[starts[i]]
      }
      
      spike_points = col[starts[i]:ends[i]]       # y values of spikes
      
      max_spike_val = max(spike_points) 
      #print(spike_points)
      max_spike_idx = which(spike_points == max_spike_val) + starts[i] -1  # includes the maximum value
      falling_spike_times = Time[max_spike_idx:ends[i]]         # spike time-interval, this is point in time of maximal value
      falling_spike_vals = col[max_spike_idx:ends[i]]
      lm_result = lm(falling_spike_vals ~ falling_spike_times)             # run linear regression

      exp_result = lm(log(falling_spike_vals) ~ falling_spike_times)
      lin_coeff = lin_coeff + unname(lm_result$coefficients[2])  
      exp_coeff = exp_coeff + unname(exp_result$coefficients[2])
      # print(summary(lm_result))
      spike_area = trapz(time, spike_points)
      
      if (spike_area > max_spike_area){          #identify the max spike area
        max_spike_area = spike_area
      }
      cum_area = cum_area+spike_area
    }
    mean_lin_coeff =lin_coeff/length(ends)
    mean_exp_coeff = exp_coeff/length(ends)
    
    results = c(  
      cum_area/length(ends), # average spike area
      length(ends),          # number of spikes
      max_spike_area,        # max_spike_area
      longest_time,          # longest time
      3*baseline,             # threshold
      mean_lin_coeff,
      mean_exp_coeff
      
    )
    return(results)
  }
}

ans = apply(cells, 2, get_area)
#results = get_area(cell)
ans = as.data.frame(do.call(cbind, ans))
ans = t(ans)
colnames(ans) = c("avg_spike_area", "num_spikes", "max_spike", "longest_spike", "threshold", "mean_lin_coeff", "mean_exp_coeff")


get_spike_inx = function(col){
  baseline = median(unlist(col))   # get median of all cell data
  variance = mad(unlist(col))      # measure of variance (adds 1.4826 factor)
  #col_minus_baseline = (col - baseline)
  all_spike_locations = (col > thresholdFactor*variance + baseline)
  
  ## method from https://masterr.org/r/how-to-find-consecutive-repeats-in-r/
  rle.all_spike_locations = rle(as.numeric(all_spike_locations))          
  for (i in 1:length(rle.all_spike_locations$lengths)){
    if(rle.all_spike_locations$lengths[i] < 10){          #TODO FIXME
      rle.all_spike_locations$values[i] = 0
    }
  }
  
  return(inverse.rle(rle.all_spike_locations))
  
}


logical = as.logical(get_spike_inx(cell))
# print(which(logical != 0 ))
 
# ggplot(data, aes(x = as.vector(as.matrix(data[,1])), y = cell))+
#   geom_point(aes(color = logical))+
#   scale_colour_manual(labels = c("<baseline", ">baseline"), values=c('blue', 'red')) +
#   labs(color = "Spikes")


library(gridExtra)

p1 = ggplot(as_tibble(ans))+
   geom_histogram(mapping = aes(mean_lin_coeff), alpha = 0.5, fill = "red") +
   geom_histogram(data = ans_frag, mapping = aes(mean_lin_coeff), alpha = .5, fill = "green")


p2 =  ggplot(as_tibble(ans))+
  geom_histogram(mapping = aes(mean_exp_coeff), fill = "red", alpha = 0.5) +
  geom_histogram(data = ans_frag, mapping = aes(mean_exp_coeff), fill = "green", alpha = .5)
grid.arrange(p1, p2)
