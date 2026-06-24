library(ggplot2)
library(ggrepel)

pca <- read.table("/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_pca.eigenvec", header = FALSE)
eval <- read.table("/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_pca.eigenval", header = FALSE)

# format columns and calculate Percent Variance Explained (PVE)
colnames(pca)[1:2] <- c("FID", "IID")
colnames(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))
pve <- round(eval$V1 / sum(eval$V1) * 100, 1)

# Extract Group/Site from Isolate ID
pca$group <- sub("\\.[0-9]+.*", "", pca$IID)

# Generate Plot
ggplot(pca, aes(x = PC1, y = PC2, color = group)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray80") +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray80") +
  
  # 95% confidence ellipses (Note: may not appear for groups with < 4 points)
  stat_ellipse(aes(fill = group), type = "t", geom = "polygon", alpha = 0.1, show.legend = FALSE) +
  
  geom_point(size = 4, alpha = 0.8) +
  
  geom_text_repel(aes(label = IID), size = 2.5, color = "black",
                  show.legend = FALSE, max.overlaps = Inf, box.padding = 0.5) +
  
  labs(title = "PCA: La Grande Motte Neurospora spp. population isolates",
       subtitle = paste0("Total Variation - PC1: ", pve[1], "% | PC2: ", pve[2], "%"),
       x = paste0("PC1 (", pve[1], "%)"),
       y = paste0("PC2 (", pve[2], "%)"),
       color = "Sample tree") + 

  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_neutre-22_PCA_snps_only_V2.pdf", 
       width = 10, height = 7, dpi = 300)
