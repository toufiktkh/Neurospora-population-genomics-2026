library(ggplot2)
library(tidyr)

df <- read.csv("mapping_stats_gm_neutre-22.csv", stringsAsFactors = FALSE)
df_clean <- separate(df, isolate, into = c("Sample", "Reference"), sep = " vs ")
df_clean$alignment_rate_pct <- as.numeric(df_clean$alignment_rate_pct)

# plot
plot_gm <- ggplot(df_clean, aes(x = Sample, y = Reference, fill = alignment_rate_pct)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradient(low = "#f7fbff", high = "#084594", na.value = "grey80") + # Unified Palette
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
    title = "Mapping Alignment Rates (%) - La Grande Motte",
    subtitle = "Isolates mapped against neutre-22 (Neurospora tetrasperma) reference",
    x = "Isolate ID",
    y = "Reference Genome",
    fill = "Alignment %"
  )

ggsave("alignment_heatmap_LGM.pdf", plot = plot_gm, width = 12, height = 4, limitsize = FALSE)

cat("Success: Heatmap saved as 'alignment_heatmap_LGM.pdf'\n")
