library(ggplot2)
library(TwoSampleMR)

data <- read.table("/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/LD_sims_MC_CIS.csv", header=T)

# data <- data[!(data$exposure==2),]

data$key <- paste0(data$model,"-",data$LD_mod)


data <- data %>%
  mutate(
    method_group = sub("_.*", "", key)  # e.g., 'Method1' from 'Method1_A'
  )

# Create spaced key_id with gaps between method groups
data <- data %>%
  group_by(method_group, key) %>%
  summarise(.groups = "drop") %>%
  arrange(method_group, key) %>%
  mutate(key_id = row_number() + cumsum(c(0, diff(as.numeric(factor(method_group))) != 0))) %>%
  right_join(data, by = c("key", "method_group"))

# Offset y by setup.mode to avoid overlap
data <- data %>%
  mutate(
    setup_mode_offset = (as.numeric(factor(LD_mod)) - 2.5) * 0.25,
    y_pos = key_id + setup_mode_offset
  )

# Y-axis labels at central position per key
y_labels <- data %>%
  distinct(key, key_id)

# Red line only in models B and D

mc_offset <- 0.8   # controls vertical separation (tune if needed)

data <- data %>%
  mutate(
    y_pos_mc = y_pos - mc_offset
  )



combined_plot <- ggplot(
  data,
  aes(
    y = y_pos,
    x = b,
    xmin = lci,
    xmax = uci,
    color = factor(model)
  )
) +
  ## Main point estimate
  geom_point(size = 2) +
  
  ## Main (analytic / bootstrap) CI
  geom_errorbar(
    height = 1,
    linewidth = 0.8,
    alpha = 0.6
  ) +
  
  ## Monte Carlo CI — jittered BELOW  🔽
  geom_errorbar(
    aes(
      y = y_pos_mc,
      xmin = mc_lci,
      xmax = mc_uci
    ),
    height = 0.8,
    linewidth = 0.8,
    color = "black",
  ) +
  
  
  ## Reference line
  geom_vline(
    xintercept = 0.4,
    linetype = "dashed",
    color = "red"
  ) +
  
  ## Axes & facets
  scale_y_continuous(
    breaks = y_labels$key_id,
    labels = y_labels$key
  ) +
  facet_wrap(~ LD_mod, scales = "free_y") +
  xlab("Beta") +
  ylab("Model / LD mode") +
  theme_bw() +
  labs(color = "Model") +
  xlim(0.2, 0.6) + theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    
  )

print(combined_plot)