setwd("~/Documents/R/MouseSleepApp/analysis/")
library(tidyverse)
library(readxl)
library(pracma)
data_rested = read.csv("rested_spike_stats.csv")
data_frag = read.csv("sleep__spike_stats.csv")
deltaF_reted = read.csv("rested_deltaF.csv")
deltaF_frag = read.csv("sleep__deltaF.csv")

library(gridExtra)

p1 = ggplot(data_rested)+
  geom_histogram(mapping = aes(mean_lin_coeff, y  = ..density..), alpha = 0.3, fill = "blue") +
  geom_histogram(data = data_frag, mapping = aes(mean_lin_coeff, y  = ..density..), alpha = .5, fill = "orange")+
  ggtitle("linear coeff")


p3 =  ggplot(data_rested)+
  geom_histogram(mapping = aes(x = log(-mean_exp_coeff), y  = ..density..), fill = "blue", alpha = 0.3) +
  geom_histogram(data = data_frag, mapping = aes(x = log(-mean_exp_coeff),y = ..density..), fill = "orange", alpha = .5)+
  ggtitle("density of log exp coeff ")

p2 =  ggplot(data_rested)+
  geom_histogram(mapping = aes(x = mean_exp_coeff, y  = ..density..), fill = "blue", alpha = 0.3) +
  geom_histogram(data = data_frag, mapping = aes(x = mean_exp_coeff, y  = ..density..), fill = "orange", alpha = .5)+
  ggtitle("exp coeff")
grid.arrange(p1, p2, p3)

# ggplot(deltaF_frag)+
#   geom_point(mapping = aes(x = seq(1:nrow(deltaF_frag))), y = select(deltaF_frag, 2))
