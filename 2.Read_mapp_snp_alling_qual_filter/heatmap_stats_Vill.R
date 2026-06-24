library(ggplot2)
library(tidyr)

# Vill data
df <- read.csv("/shared/home/toufiktakhi/neurosporagbs/mapping_stats/mapping_stats_vill.csv", stringsAsFactors = FALSE)

df_clean <- separate(df, isolate, into = c("Sample", "Reference"), sep = " vs ")
df_clean$alignment_rate_pct <- as.numeric(df_clean$alignment_rate_pct)

# plot
plot_vill <- ggplot(df_clean, aes(x = Sample, y = Reference, fill = alignment_rate_pct)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradient(low = "#f7fbff", high = "#084594", na.value = "grey80") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 20)),
    axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1, color = "black"), 
    axis.text.y = element_text(face = "bold", size = 11, color = "black"), 
    panel.grid = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Mapping Alignment Rates (%) - Villeveyrac",
    subtitle = "Isolates mapped against Vill-A1-3 reference",
    x = "Isolate ID",
    y = "Reference Genome",
    fill = "Alignment %"
  )

ggsave("alignment_heatmap_vill.pdf", plot = plot_vill, width = 12, height = 4, limitsize = FALSE)

cat("Success: 'alignment_heatmap_vill.pdf' has been generated.\n")
