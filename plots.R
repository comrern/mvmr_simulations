library(ggplot2)
library(TwoSampleMR)

data <- read.table("/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/real_data_10k.csv")
data <- generate_odds_ratios(data)

# data <- data[!(data$method== "IVW" & data$exposure==2),]

data$key <- paste0(data$model,"-", data$method ,"-",data$exposure)

data$labels <- c("IVW - exposure 1", "MVMR - exposure 1", "MVMR - ancestry","IVW - exposure 1", "MVMR - exposure 1", "MVMR - ancestry","IVW - exposure 1", "MVMR - exposure 1", "MVMR - ancestry","IVW - exposure 1", "MVMR - exposure 1", "MVMR - ancestry")

### subplots




  combined_plot <- ggplot(data, aes(y=(key), x=(or),xmin=(or_lci95), xmax=(or_uci95),group=model))
  
  
  print(combined_plot + geom_point() +
    xlab("Odds ratio") + 
    geom_errorbarh(height=.1) +
    geom_vline(xintercept=1) +
    scale_y_discrete(drop=T, labels=data$labels) +  
    facet_wrap(~ model, scale= "free_y") +
    ggtitle(paste0("Setup mode"), model))
  
  
  
  ### run without ancestry estimate
  
  data_noanc <- data[!(data$method== "MVMR" & data$exposure==2),]
  combined_plot <- ggplot(data_noanc, aes(y=(key), x=(or),xmin=(or_lci95), xmax=(or_uci95),group=model))
  
  
  print(combined_plot + geom_point() +
          xlab("Odds ratio") + 
          geom_errorbarh(height=.1) +
          geom_vline(xintercept=1) +
          scale_y_discrete(drop=T, labels=data$labels) +  
          facet_wrap(~ model, scale= "free_y") +
          ggtitle(paste0("Setup mode"), model))
  
  

