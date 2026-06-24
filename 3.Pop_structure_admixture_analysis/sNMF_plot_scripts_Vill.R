library(LEA)
library(ggplot2)
library(tidyr)

setwd("/shared/home/toufiktakhi/neurosporagbs/")

cat("Launching the Villeveyrac project...\n")
project_vill <- load.snmfProject("/shared/projects/neurosporagbs/admixture_analysis/Vill_snps_only_clone_corr_snmf.snmfProject")

vill_samples <- system("bcftools query -l ~/neurosporagbs/pca_splitstree_vill/all_sites/Vill_snps_invar_clone_corr_snmf.vcf.gz", intern = TRUE)

vill_tree <- gsub("^Vill[.-]([A-Z][0-9]+).*", "\\1", vill_samples)
vill_tree[vill_tree == vill_samples] <- "Other"

#Extracting Q matrix for K=3
cat("Extraction for K=3...\n")
k <- 3

best_run <- which.min(cross.entropy(project_vill, K = k))
Q_matrix <- Q(project_vill, K = k, run = best_run)
colnames(Q_matrix) <- paste0("Pop", 1:k)

df_k <- as.data.frame(Q_matrix)
df_k$sample <- vill_samples
df_k$tree <- vill_tree

vill_multi_K_df <- pivot_longer(df_k,
                                cols = starts_with("Pop"),
                                names_to = "population",
                                values_to = "proportion")

vill_multi_K_df <- vill_multi_K_df[!is.na(vill_multi_K_df$proportion), ]

#Harmonizing and selecting factors
vill_multi_K_df$population <- factor(vill_multi_K_df$population, levels = paste0("Pop", 1:k))

sample_order <- vill_samples[order(vill_tree)]
vill_multi_K_df$sample <- factor(vill_multi_K_df$sample, levels = unique(sample_order))

# Plotting
cat("Generating the graphic...\n")

p <- ggplot(vill_multi_K_df, aes(x = sample, y = proportion, fill = population)) +
    geom_bar(stat = "identity", width = 1) +
    facet_grid(. ~ tree, scales = "free_x", space = "free_x") +
    scale_fill_brewer(palette = "Set1", name = "Ancestral\nClusters", na.value = "transparent") +
    labs(title = "sNMF Ancestry Proportions of N. crassa clade B from Villeveyrac site at K = 3",
         x = "Isolates (grouped by sample shrubs)", y = "Ancestry Proportion") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(size = 9, face = "bold", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(face = "bold", size = 9),
        panel.spacing.x = unit(0.1, "lines"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )

#Saving
output_file <- "/shared/home/toufiktakhi/neurosporagbs/snmf/Vill_sNMF_bestK.pdf"
ggsave(output_file, p, width = 16, height = 6, dpi = 300)

cat("Success! The graph has been generated at:", output_file, "\n")
