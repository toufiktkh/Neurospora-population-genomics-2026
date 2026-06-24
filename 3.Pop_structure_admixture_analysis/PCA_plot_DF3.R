library(ggplot2)
library(ggrepel)

pca_df1 <- read.table("/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_pca.eigenvec", header=FALSE)
eval_df1 <- read.table("/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_pca.eigenval", header=FALSE)

#Format columns and calculate Percent Variance Explained
colnames(pca_df1)[1:2] <- c("FID", "IID")
colnames(pca_df1)[3:ncol(pca_df1)] <- paste0("PC", 1:(ncol(pca_df1)-2))
pve_df1 <- round(eval_df1$V1 / sum(eval_df1$V1) * 100, 1)

# Extract Transect grouping from Sample ID
pca_df1$transect <- substr(pca_df1$IID, 5, 6)

ggplot(pca_df1, aes(x = PC1, y = PC2, color = transect)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  stat_ellipse(aes(fill = transect), type = "t", geom = "polygon", alpha = 0.1, show.legend = FALSE) +
  geom_point(size = 3.5, alpha = 0.8) +

# ggrepel was used here to prevent label overlapping
  geom_text_repel(aes(label = IID), size = 2.5, color = "black", 
                  show.legend = FALSE, max.overlaps = Inf, box.padding = 0.5) +
  labs(title = "Intraspecific PCA of Neurospora discreta PS4-001B populations mapped against E-DF3",
       subtitle = paste0("Total Variation - PC1: ", pve_df1[1], "% | PC2: ", pve_df1[2], "%"),
       x = paste0("PC1 (", pve_df1[1], "%)"),
       y = paste0("PC2 (", pve_df1[2], "%)"),
       color = "Sample tree") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

ggsave("/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_PCA_snps_only.pdf", width = 10, height = 7, dpi = 300)
