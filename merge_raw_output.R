

files <- list.files(path = "/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/real_data_sims/", pattern = "results_.*\\.csv", full.names = TRUE)
# files <- list.files(path = "/user/work/kb22541/simulations/results", pattern = "^results_.*\\.csv$", full.names = TRUE)
data_list <- lapply(files, function(f) {
  e <- new.env()             # Create a temporary environment
  load(f, envir = e)         # Load the .RData file into that environment
  obj_names <- ls(envir = e)
  
  if (length(obj_names) != 1) {
    stop(paste("File", f, "contains more than one object."))
  }
  
  get(obj_names[[1]], envir = e)  # Safely retrieve the one object
})

final_data <- do.call(rbind, data_list)


write.table(final_data, "/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/combined_results10k.csv", row.names=F, quote=F)