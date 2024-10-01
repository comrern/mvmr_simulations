library(ggplot2)
library(TwoSampleMR)


data <- read.csv("./results/results_averaged.csv")
data <- generate_odds_ratios(data)

data <- data[!(data$method== "IVW" & data$exposure==2),]

data$key <- paste0(data$model,"-", data$method ,"-",data$exposure)

### megaplot ###


combined_plot <- ggplot(data, aes(y=(key), x=(or),xmin=(or_lci95), xmax=(or_uci95),color=(model)))

combined_plot + geom_point() + xlab("Odds ratio") + geom_errorbarh(height=.1) + geom_vline(xintercept=1)