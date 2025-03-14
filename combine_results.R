
setwd("C:/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/raw")


files <- list.files(pattern = "^results_.*\\.csv$", full.names = TRUE)

data_list <- lapply(files, function(f) {
  env <- new.env()
  load(f, envir = env)
  get(ls(env)[1], envir = env)  # Extract the first object
})


final_data <- do.call(rbind, data_list)  # Merge all into one data frame


write.table(final_data, "../results_full_sims.csv")