library(ggplot2)
library(TwoSampleMR)

data <- read.table("/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/avergaed_results_fullsims.csv", sep=",", header=T)

data <- data[!(data$Exposure==2),]

data$key <- paste0(data$Model,"-", data$Method ,"-",data$Exposure)

data$labels <- c("IVW - Exposure 1", "MVMR - Exposure 1","IVW - Exposure 1", "MVMR - Exposure 1", "IVW - Exposure 1", "MVMR - Exposure 1", "IVW - Exposure 1", "MVMR - Exposure 1")

data$lci95 <- data$beta - (data$se * 1.96)
data$uci95 <- data$beta + (data$se * 1.96)


### subplots


for (mode in c(1,2,3,4)) {
  
  mode_dat <- data[data$setup.mode == mode, ]
  
  # Create a dataframe to hold vertical lines for Models B and D
  vline_data <- data.frame(
    xintercept = 0.4,
    Model = c("B", "D")
  )
  
  combined_plot <- ggplot(mode_dat, aes(y = key, x = beta, xmin = lci95, xmax = uci95, group = model)) +
    geom_point() +
    xlab("Beta") + 
    theme_bw() +
    geom_errorbarh(height = 0.1) +
    geom_vline(xintercept = 0) +
    geom_vline(data = vline_data, aes(xintercept = xintercept), 
               linetype = "dashed", color = "red") +
    scale_y_discrete(drop = TRUE, labels = data$labels) +  
    facet_wrap(~ Model, scales = "free_y") +
    ggtitle(paste("Setup mode", mode)) +
    theme(legend.position = "none")
  
  print(combined_plot)
}


