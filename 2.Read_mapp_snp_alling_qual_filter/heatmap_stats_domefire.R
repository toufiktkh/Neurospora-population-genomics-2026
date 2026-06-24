library(ggplot2)
library(tidyr)

# open the data
df <- read.csv("mapping_stats.csv", stringsAsFactors = FALSE)


df_clean <- separate(df, isolate, into = c("Sample", "Reference"), sep = " vs ")

df_clean$alignment_rate_pct <- as.numeric(df_clean$alignment_rate_pct)

# plot
plot_full <- ggplot(df_clean, aes(x = Sample, y = Reference, fill = alignment_rate_pct)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#f7fbff", high = "#084594", na.value = "grey80") + # Unified Palette
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1, color = "black"), 
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 20)),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Mapping Alignment Rates",
    subtitle = "All isolates against E-DF1 and E-DF3 references",
    x = "Isolate ID",
    y = "Reference Genome",
    fill = "Alignment %"
  )

ggsave("alignment_heatmap_domefire.pdf", plot = plot_full, width = 16, height = 5)

cat("Success: Heatmap saved as 'alignment_heatmap_domefire.pdf'\n")
