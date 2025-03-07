

files <- list.files(path = "/user/work/kb22541/simulations/results", pattern = "results_.*\\.csv", full.names = TRUE)
files <- list.files(path = "/user/work/kb22541/simulations/results", pattern = "^results_.*\\.csv$", full.names = TRUE)
data_list <- lapply(files, function(f) {
  load(f)   # This loads the object into the environment
  get(ls()[1])  # Retrieve the first object in the environment
})

final_data <- do.call(rbind, data_list)  # Combine all into one data frame
