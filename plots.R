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
  
  combined_plot <- ggplot(mode_dat, aes(y=(key), x=(b),xmin=(lo_ci), xmax=(up_ci),group=model))
  
  
  print(combined_plot + geom_point() +
    xlab("beta") + 
    theme_bw() +
    geom_errorbarh(height=.1) +
    geom_vline(xintercept=0) +
    scale_y_discrete(drop=T, labels=data$labels) +  
    facet_wrap(~ model, scale= "free_y") +
    theme(legend.position="none"))
   
  
}

