library(ggplot2)
library(TwoSampleMR)

data <- read.table("/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/averaged_ld_sims_results.txt", header=T)

# data <- data[!(data$exposure==2),]

data$key <- paste0(data$model,"-", data$method ,"-",data$LD_mod)

data$labels <- c("IVW - exposure 1", "MVMR - exposure 1","IVW - exposure 1", "MVMR - exposure 1", "IVW - exposure 1", "MVMR - exposure 1", "IVW - exposure 1", "MVMR - exposure 1")

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


# Plot
combined_plot <- ggplot(data, aes(y = y_pos, x = b, xmin = lci, xmax = uci, color = factor(model))) +
  geom_point() +
  geom_errorbar(height = 0.15) +
  geom_vline(xintercept = 0.4, 
             linetype = "dashed", color = "red") +
  scale_y_continuous(
    breaks = y_labels$LD_mod,
    labels = y_labels$LD_mod
  ) +
  facet_wrap(~ LD_mod, scales = "free_y") +
  xlab("Beta") +
  ylab("Model / Mode") +
  theme_bw() +
  labs(color = "Model") +
  xlim(0.2, 0.6)

print(combined_plot)
