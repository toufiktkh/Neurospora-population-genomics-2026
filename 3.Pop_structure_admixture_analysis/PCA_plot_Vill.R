library(ggplot2)
library(ggrepel)

pca <- read.table("/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_vill/all_sites/Vill_pca.eigenvec", header = FALSE)
eval <- read.table("/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_vill/all_sites/Vill_pca.eigenval", header = FALSE)

colnames(pca)[1:2] <- c("FID", "IID")
colnames(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))
pve <- round(eval$V1 / sum(eval$V1) * 100, 1)

# Extract tree ID (Keeping your precise filtering rules)
pca$tree_id <- pca$IID
pca$tree_id <- gsub("^Vill\\.([A-C][0-9]+)[BIS.0-9]*.*", "\\1", pca$tree_id)
pca$tree_id <- gsub("^Vill-([A-C][0-9]+)[-_].*", "\\1", pca$tree_id)

pca$site <- substr(pca$tree_id, 1, 1)

n_trees <- length(unique(pca$tree_id))
tree_colors <- setNames(
    colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
                       "#FF7F00","#A65628","#F781BF","#999999"))(n_trees),
    sort(unique(pca$tree_id))
)

# plot
ggplot(pca, aes(x = PC1, y = PC2, color = tree_id)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +

  stat_ellipse(aes(fill = tree_id), type = "t", geom = "polygon", alpha = 0.1, show.legend = FALSE) +

  geom_point(aes(shape = site), size = 3.5, alpha = 0.8) +

  geom_text_repel(aes(label = IID), size = 2.5, color = "black",
                  show.legend = FALSE, max.overlaps = Inf, box.padding = 0.5) +

  scale_color_manual(values = tree_colors) +
  scale_fill_manual(values = tree_colors) + 
  scale_shape_manual(values = c("A" = 15, "B" = 16, "C" = 17)) + # 17 is a clean triangle for site C

  labs(title = "Intraspecific PCA: Neurospora crassa populations from Villeveyrac",
       subtitle = paste0("Total Variation - PC1: ", pve[1], "% | PC2: ", pve[2], "%"),
       x = paste0("PC1 (", pve[1], "%)"),
       y = paste0("PC2 (", pve[2], "%)"),
       color = "Sample tree",
       shape = "Sampling site") +

  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

ggsave("/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_vill/all_sites/Vill_PCA_snps_only.pdf", width = 10, height = 7, dpi = 300)
