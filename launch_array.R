args <- as.numeric(commandArgs(T))
set.seed((args[1]*100000))
job_id <- ((args[1]))
message("job number ", job_id)


setwd("/user/work/kb22541/simulations")
output_path <- "./results"
.libPaths("/user/work/kb22541/rlib")



library(dplyr)
library(MASS)
library(TwoSampleMR)
library(MVMR)
library(truncnorm)
library(tidyverse)

source("run_sims.R")
source('modes_sims.R')
source('functions_sims.R')


output <- 


save(output, file=sprintf("rmvmr_sims_%s.Rda", job_id))

