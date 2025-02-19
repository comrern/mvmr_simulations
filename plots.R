library(ggplot2)
library(TwoSampleMR)

data <- read.csv("./results/results_averaged_ld_mod.csv")
data <- generate_odds_ratios(data)

data <- data[!(data$method== "IVW" & data$exposure==2),]

data$key <- paste0(data$model,"-", data$method ,"-",data$exposure)

data$labels <- c("IVW - exposure 1", "MVMR - exposure 1", "MVMR - ancestry","IVW - exposure 1", "MVMR - exposure 1", "MVMR - ancestry","IVW - exposure 1", "MVMR - exposure 1", "MVMR - ancestry","IVW - exposure 1", "MVMR - exposure 1", "MVMR - ancestry")

### subplots


for (mode in c(1,2,3,4)){
  
  mode_dat <- data[data$setup_mode == mode,]
  
  combined_plot <- ggplot(mode_dat, aes(y=(key), x=(or),xmin=(or_lci95), xmax=(or_uci95),group=model))
  
  
  print(combined_plot + geom_point() +
    xlab("Odds ratio") + 
    geom_errorbarh(height=.1) +
    geom_vline(xintercept=1) +
    scale_y_discrete(drop=T, labels=data$labels) +  
    facet_wrap(~ model, scale= "free_y") +
    ggtitle(paste0("Setup mode"), mode))
  
  
}

