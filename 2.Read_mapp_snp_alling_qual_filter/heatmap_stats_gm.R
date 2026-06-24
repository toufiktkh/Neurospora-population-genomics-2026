library(ggplot2)
library(tidyr)
library(RColorBrewer)

# load data
df <- read.csv("mapping_stats_gm_neutre-22.csv", stringsAsFactors = FALSE)

# split the 'isolate' column into 'Sample' and 'Reference'
df_clean <- separate(df, isolate, into = c("Sample", "Reference"), sep = " vs ")
df_clean$alignment_rate_pct <- as.numeric(df_clean$alignment_rate_pct)

# plot
plot_vill <- ggplot(df_clean, aes(x = Reference, y = Sample, fill = alignment_rate_pct)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1, na.value = "grey90") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 20)),
    axis.text.x = element_text(face = "bold", size = 10, color = "black"),
    axis.text.y = element_text(size = 7, color = "black"), # Slightly larger for Vill samples
    
    panel.grid = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Mapping Alignment Rates (%) - La Grande Motte",
    subtitle = "Isolates mapped against neutre-22 (Neurospora tetrasperma) reference",
    x = "Reference Genome",
    y = "Isolate ID",
    fill = "Alignment %"
  )

ggsave("heatmap_gm_neutre-22.pdf", plot = plot_vill, width = 8, height = 12, limitsize = FALSE)

cat("Success: 'heatmap_gm_neutre-22.pdf' has been generated.\n")
