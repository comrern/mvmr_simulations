library(ggplot2)
library(TwoSampleMR)

data <- read.table("C:/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/effect_size_sims/avergaed_effect_size.csv")

data <- cbind.data.frame(data, (data$b - (1.96 * data$se)),((data$b + (1.96 * data$se)))) 

# data <- data[!(data$exposure==2),]

colnames(data)[13:14] <- c("lci","uci")

data$key <- paste0(data$model,"-", data$method ,"-",data$exposure)

data$labels <- c("IVW - exposure 1", "MVMR - exposure 1","IVW - exposure 1", "MVMR - exposure 1", "IVW - exposure 1", "MVMR - exposure 1", "IVW - exposure 1", "MVMR - exposure 1")

### subplots


for (mode in c(1,2,3,4)){
  
  mode_dat <- data[data$setup_mode == mode,]
  
  combined_plot <- ggplot(mode_dat, aes(y=(key), x=(b),xmin=(lci), xmax=(uci),group=model))
  
  
  print(combined_plot + geom_point() +
    xlab("effect size") + 
    theme_bw() +
    geom_errorbarh(height=.1) +
    geom_vline(xintercept=0) +
    scale_y_discrete(drop=T, labels=data$labels) +  
    facet_wrap(~ model, scale= "free_y") +
    ggtitle(paste("Setup mode", mode))) +
    theme(legend.position="none")
   
  
}


##### plot increasing effect size for model B

data_eff_size_plot <- data[data$model == "B" & data$method == "MVMR",]

ggplot(data_eff_size_plot, aes(x = setup_mode, y = b)) +
  geom_point(size = 3) +  # dots for point estimates
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.2) +  # CI bars
  geom_hline(yintercept = 0.4, linetype = "dashed", color = "red") +  # intercept line
  scale_x_continuous(breaks = 1:4) +
  labs(x = "setup mode", y = "Effect size (beta)", title = "MVMR causal effect estimate for exposure 1 with increasing effect size of SNPs on expousre 1") +
  geom_hline(yintercept = 0) +
  theme_minimal()


data_eff_size_plot <- data[data$model == "D" & data$method == "MVMR",] ## re-run on model D

ggplot(data_eff_size_plot, aes(x = setup_mode, y = b)) +
  geom_point(size = 3) +  # dots for point estimates
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.2) +  # CI bars
  geom_hline(yintercept = 0.4, linetype = "dashed", color = "red") +  # intercept line
  scale_x_continuous(breaks = 1:4) +
  labs(x = "setup mode", y = "Effect size (beta)") +
  geom_hline(yintercept = 0) +
  theme_minimal()




