library(ggplot2)
library(TwoSampleMR)

data <- read.table("/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/real_data_sims_processed_260925")

# data <- data[!(data$method== "IVW" & data$exposure==2),]

data$key <- paste0(data$model,"-", data$method ,"-",data$exposure)

data$labels <- c("IVW - exposure 1", "MVMR - exposure 1", "MVMR - ancestry","IVW - exposure 1", "MVMR - exposure 1", "MVMR - ancestry","IVW - exposure 1", "MVMR - exposure 1", "MVMR - ancestry","IVW - exposure 1", "MVMR - exposure 1", "MVMR - ancestry")

### subplots

data$lci95 <- data$b - (data$se * 1.96)
data$uci95 <- data$b + (data$se * 1.96)



  combined_plot <- ggplot(data, aes(y=(key), x=(b),xmin=(lci95), xmax=(uci95),group=model))

  vline_data <- data.frame(
    xintercept = 0.4,
    model = c("B", "D")   # lowercase to match facet variable
  )
  
  print(combined_plot + 
          geom_point() +
          xlab("Effect size (Beta)") + 
          geom_errorbarh(height=.1) +
          geom_vline(xintercept=0) +
          geom_vline(
            data = vline_data, 
            aes(xintercept = xintercept), 
            linetype = "dashed", color = "red"
          ) +
          scale_y_discrete(drop=TRUE, labels = data$labels) +  
          facet_wrap(~ model, scale = "free_y") +
          ggtitle("Setup mode"))
  
  
  
  
  ### run without ancestry estimate
  
  data_noanc <- data[!(data$method== "MVMR" & data$exposure==2),]
  combined_plot <- ggplot(data_noanc, aes(y=(key), x=(b),xmin=(lci95), xmax=(uci95),group=model))
  
  
  print(combined_plot + 
          geom_point() +
          xlab("Effect size (Beta)") + 
          geom_errorbarh(height=.1) +
          geom_vline(xintercept=0) +
          geom_vline(
            data = vline_data, 
            aes(xintercept = xintercept), 
            linetype = "dashed", color = "red"
          ) +
          scale_y_discrete(drop=TRUE, labels = data$labels) +  
          facet_wrap(~ model, scale = "free_y") )
  

