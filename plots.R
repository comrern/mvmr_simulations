library(ggplot2)
library(TwoSampleMR)

data <- read.table("/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/mainsims_MCCIs.csv", header=T, sep= ",")


data <- data[1:48,]
data <- data[!(data$exposure==2),]

data$key <- paste0(data$model,"-", data$method ,"-",data$exposure)

data$labels <- c("IVW - Exposure 1", "MVMR - Exposure 1","IVW - Exposure 1", "MVMR - Exposure 1", "IVW - Exposure 1", "MVMR - Exposure 1", "IVW - Exposure 1", "MVMR - Exposure 1")

data$lci95 <- data$b - (data$se * 1.96)
data$uci95 <- data$b + (data$se * 1.96)




### subplots


# for (mode in c(1,2,3,4)) {
#   
#   mode_dat <- data[data$setup.mode == mode, ]
#   
#   # Create a dataframe to hold vertical lines for Models B and D
#   vline_data <- data.frame(
#     xintercept = 0.4,
#     Model = c("B", "D")
#   )
#   
#   combined_plot <- ggplot(mode_dat, aes(y = key, x = beta, xmin = lci95, xmax = uci95, group = model)) +
#     geom_point() +
#     xlab("Beta") + 
#     theme_bw() +
#     geom_errorbarh(height = 0.1) +
#     geom_vline(xintercept = 0) +
#     geom_vline(data = vline_data, aes(xintercept = xintercept), 
#                linetype = "dashed", color = "red") +
#     scale_y_discrete(drop = TRUE, labels = data$labels) +  
#     facet_wrap(~ Model, scales = "free_y") +
#     ggtitle(paste("Setup mode", mode)) +
#     theme(legend.position = "none")
#   
#   print(combined_plot)
# }


# Assume 'key' is something like 'Method1_A', 'Method2_A', etc.
# First, extract 'method' grouping for spacing
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
    setup_mode_offset = (as.numeric(factor(setup_mode)) - 2.5) * 0.5,
    y_pos = key_id + setup_mode_offset
  )

# Y-axis labels at central position per key
y_labels <- data %>%
  distinct(key, key_id)

# Red line only in models B and D
vline_data <- data.frame(
  xintercept = 0.4,
  model = c("B", "D")
)

# Plot
combined_plot <- ggplot(
  data,
  aes(
    y = y_pos,
    x = b,
    xmin = lci95,
    xmax = uci95,
    color = factor(setup_mode)
  )
) +
  ## Point estimate (keep centered)
  geom_point(size = 2) +
  
  ## Analytic CI — nudged slightly UP
  geom_errorbarh(
    height = 0.18,
    linewidth = 0.8,
    alpha = 0.6
  ) +
  
  ## Monte Carlo CI — nudged slightly DOWN
  geom_errorbarh(
    aes(xmin = mc_lci, xmax = mc_uci),
    height = 0.08,
    linewidth = 0.8,
    color = "black",
    inherit.aes = TRUE,
    position = position_nudge(y = -0.08)
  ) +
  
  ## Reference lines
  geom_vline(xintercept = 0, color = "black") +
  geom_vline(
    data = vline_data,
    aes(xintercept = xintercept),
    linetype = "dashed",
    color = "red"
  ) +
  
  ## Axes & facets
  scale_y_continuous(
    breaks = y_labels$key_id,
    labels = y_labels$key
  ) +
  facet_wrap(~ model, scales = "free_y") +
  xlab("Beta") +
  ylab("Key") +
  ggtitle("Combined Results by Setup Mode") +
  theme_bw() +
  labs(color = "Setup Mode")

print(combined_plot)