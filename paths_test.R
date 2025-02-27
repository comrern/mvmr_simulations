
set.seed(1234)

setwd("/user/work/kb22541/simulations")
output_path <- "./results"
.libPaths("/user/work/kb22541/rlib")



library(dplyr)
library(tidyverse)
library(MASS)
library(TwoSampleMR)
library(MVMR)
library(truncnorm)


source('modes_sims.R')
source('functions_sims.R')






