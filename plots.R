library(ggplot2)
library(TwoSampleMR)

data <- read.table("/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/avg_test_ld_exp.csv")
data <- generate_odds_ratios(data)

data <- data[!(data$exposure==2),]

data$key <- paste0(data$model,"-", data$method ,"-",data$exposure)

data$labels <- c("IVW - exposure 1", "MVMR - exposure 1","IVW - exposure 1", "MVMR - exposure 1", "IVW - exposure 1", "MVMR - exposure 1", "IVW - exposure 1", "MVMR - exposure 1")

### subplots


for (mode in c("TRUE","FALSE")){
  
  mode_dat <- data[data$LD_mod == mode,]
  
  
  # Create a dataframe to hold vertical lines for Models B and D
  vline_data <- data.frame(
    xintercept = 0.4,
    model = c("B", "D")
  )
  
  combined_plot <- ggplot(mode_dat, aes(y = key, x = b, xmin = lo_ci, xmax = up_ci, group = model)) +
    geom_point() +
    xlab("Beta") + 
    theme_bw() +
    geom_errorbarh(height = 0.1) +
    geom_vline(xintercept = 0) +
    geom_vline(data = vline_data, aes(xintercept = xintercept), 
               linetype = "dashed", color = "red") +
    scale_y_discrete(drop = TRUE, labels = data$labels) +  
    facet_wrap(~ model, scales = "free_y") +
    ggtitle(paste("Setup mode", mode)) +
    theme(legend.position = "none")
  
  print(combined_plot)
}

