if (!require(ggplot2)) install.packages("ggplot2", repos="https://cloud.r-project.org")
if (!require(dplyr)) install.packages("dplyr", repos="https://cloud.r-project.org")
library(ggplot2)
library(dplyr)

out_dir <- "~/neurosporagbs/pixy_plots/"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

pop_colors <- c(
  "Domefire_DF1" = "#008cff",
  "Domefire_DF3" = "#1e4a66",
  "Villeveyrac"  = "#16beab",
  "La_Grande_Motte"          = "#4f8b00"
)

load_pi <- function(path, pop_label, min_sites = 0) {
  if (!file.exists(path)) return(data.frame())
  df <- read.table(path, header=TRUE, stringsAsFactors=FALSE, sep="\t")
  df$site <- pop_label
  df$pop <- paste(pop_label, df$pop, sep="_")
  df$window_mid <- (df$window_pos_1 + df$window_pos_2) / 2
  
  # Eliminating NAs and applying the filters on the "no_sites" values
  df <- df[!is.na(df$avg_pi), ]
  df <- df[df$no_sites >= min_sites, ] 
  
  return(df)
}
# PI calculations
SEUIL_SITES <- 1000 

pi_df1  <- load_pi("~/neurosporagbs/pca_splitstree_domefire/all_sites/pixy_output_DF1_90/pixy_pi.txt", "Domefire_DF1", min_sites = SEUIL_SITES)
pi_df3  <- load_pi("~/neurosporagbs/pca_splitstree_domefire/all_sites/pixy_output_DF3_90/pixy_pi.txt", "Domefire_DF3", min_sites = SEUIL_SITES)
pi_vill <- load_pi("~/neurosporagbs/pca_splitstree_vill/all_sites/pixy_output_90/pixy_pi.txt", "Villeveyrac", min_sites = SEUIL_SITES)
pi_lgm  <- load_pi("~/neurosporagbs/pca_splitstree_LGM/all_sites/pixy_output/pixy_pi.txt", "La_Grande_Motte", min_sites = SEUIL_SITES)

pi_all <- rbind(pi_df1, pi_df3, pi_vill, pi_lgm)
pi_all$site <- factor(pi_all$site, levels=c("Domefire_DF1", "Domefire_DF3", "Villeveyrac", "La_Grande_Motte"))

pi_all <- pi_all %>%
  arrange(site, chromosome, window_pos_1) %>%
  group_by(pop) %>%
  mutate(window_index = row_number()) %>%
  ungroup()

load_tajd <- function(path, pop_label, min_sites = 0) {
  if (!file.exists(path)) return(data.frame())
  df <- read.table(path, header=TRUE, stringsAsFactors=FALSE, sep="\t")
  df$site <- pop_label
  df$pop <- paste(pop_label, df$pop, sep="_")
  df$window_mid <- (df$window_pos_1 + df$window_pos_2) / 2
  
  df <- df[!is.na(df$tajima_d), ]
  df <- df[df$no_sites >= min_sites, ] # Filtre de Pierre
  return(df)
}

# Calculating Tajima's D

taj_df1  <- load_tajd("~/neurosporagbs/pca_splitstree_domefire/all_sites/pixy_output_DF1_90/pixy_tajima_d.txt", "Domefire_DF1", min_sites = SEUIL_SITES)
taj_df3  <- load_tajd("~/neurosporagbs/pca_splitstree_domefire/all_sites/pixy_output_DF3_90/pixy_tajima_d.txt", "Domefire_DF3", min_sites = SEUIL_SITES)
taj_vill <- load_tajd("~/neurosporagbs/pca_splitstree_vill/all_sites/pixy_output_90/pixy_tajima_d.txt", "Villeveyrac", min_sites = SEUIL_SITES)
taj_lgm  <- load_tajd("~/neurosporagbs/pca_splitstree_LGM/all_sites/pixy_output/pixy_tajima_d.txt", "La_Grande_Motte", min_sites = SEUIL_SITES)

taj_all <- rbind(taj_df1, taj_df3, taj_vill, taj_lgm)
taj_all$site <- factor(taj_all$site, levels=c("Domefire_DF1", "Domefire_DF3", "Villeveyrac", "La_Grande_Motte"))

taj_all <- taj_all %>%
  arrange(site, chromosome, window_pos_1) %>%
  group_by(pop) %>%
  mutate(window_index = row_number()) %>%
  ungroup()

load_dxy <- function(path, pop_label, min_sites = 0) {
  if (!file.exists(path)) return(data.frame())
  df <- read.table(path, header=TRUE, stringsAsFactors=FALSE, sep="\t")
  if (nrow(df) == 0) return(data.frame())
  df$site <- pop_label
  df$pop <- paste(df$pop1, "vs", df$pop2)
  df$window_mid <- (df$window_pos_1 + df$window_pos_2) / 2
  
  df <- df[!is.na(df$avg_dxy), ]
  df <- df[df$no_sites >= min_sites, ] # Filtre de Pierre
  return(df)
}

# Calculating Dxy

dxy_df1  <- load_dxy("~/neurosporagbs/pca_splitstree_domefire/all_sites/pixy_output_DF1_90/pixy_dxy.txt", "Domefire_DF1", min_sites = SEUIL_SITES)
dxy_df3  <- load_dxy("~/neurosporagbs/pca_splitstree_domefire/all_sites/pixy_output_DF3_90/pixy_dxy.txt", "Domefire_DF3", min_sites = SEUIL_SITES)
dxy_vill <- load_dxy("~/neurosporagbs/pca_splitstree_vill/all_sites/pixy_output_90/pixy_dxy.txt", "Villeveyrac", min_sites = SEUIL_SITES)
dxy_lgm  <- load_dxy("~/neurosporagbs/pca_splitstree_LGM/all_sites/pixy_output/pixy_dxy.txt", "La_Grande_Motte", min_sites = SEUIL_SITES)

dxy_all <- rbind(dxy_df1, dxy_df3, dxy_vill, dxy_lgm)
dxy_all$site <- factor(dxy_all$site, levels=c("Domefire_DF1", "Domefire_DF3", "Villeveyrac", "La_Grande_Motte"))

dxy_all <- dxy_all %>%
  arrange(site, chromosome, window_pos_1) %>%
  group_by(pop) %>%
  mutate(window_index = row_number()) %>%
  ungroup()

# Plotting

p1 <- ggplot(pi_all, aes(x=pop, y=avg_pi, fill=site)) +
  geom_violin(alpha=0.65, trim=FALSE, linewidth=0.3) +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA, linewidth=0.5) +
  scale_fill_manual(values=pop_colors) +
    theme_minimal(base_size=14) +
  labs(title="Nucleotide Diversity (Pi)", x="", y="Average Pi per window") +
  theme(legend.position="none", panel.grid.major.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~site, scales="free")

ggsave(file.path(out_dir, "Plot_1_Pi_Violin_90.pdf"), plot=p1, width=12, height=6)

p2 <- ggplot(taj_all, aes(x=pop, y=tajima_d, fill=site)) +
  geom_violin(alpha=0.65, trim=FALSE, linewidth=0.3) +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA, linewidth=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=0.8) +
  scale_fill_manual(values=pop_colors) +
  theme_minimal(base_size=14) +
  labs(title="Tajima's D", x="", y="Tajima's D") +
  theme(legend.position="none", panel.grid.major.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~site, scales="free")

ggsave(file.path(out_dir, "Plot_2_TajimasD_Violin_90.pdf"), plot=p2, width=12, height=6)

p3 <- ggplot(taj_all, aes(x=window_index, y=tajima_d, color=pop)) +
  geom_point(alpha=0.35, size=0.8) +
  geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=0.7) +
  geom_smooth(method="loess", span=0.15, se=FALSE, linewidth=0.8, color="black") +
  facet_wrap(~pop, ncol=1, scales="free_x") +
  theme_minimal(base_size=13) +
  labs(title="Tajima's D - Sliding Windows (50kb)",
       x="Window index (ordered by contig/position)",
       y="Tajima's D") +
  theme(legend.position="none", strip.text=element_text(size=13, face="bold"))

ggsave(file.path(out_dir, "Plot_3_TajimasD_Scatter_90.pdf"), plot=p3, width=10, height=12)

p4 <- ggplot(pi_all, aes(x=window_index, y=avg_pi, color=pop)) +
  geom_point(alpha=0.35, size=0.8) +
  geom_smooth(method="loess", span=0.15, se=FALSE, linewidth=0.8, color="black") +
  facet_wrap(~pop, ncol=1, scales="free_x") +
  theme_minimal(base_size=13) +
  labs(title="Nucleotide Diversity (Pi) - Sliding Windows (50kb)",
       x="Window index (ordered by contig/position)",
       y="Average Pi per window") +
  theme(legend.position="none", strip.text=element_text(size=13, face="bold"))

ggsave(file.path(out_dir, "Plot_4_Pi_Scatter_90.pdf"), plot=p4, width=10, height=12)

p5 <- ggplot(dxy_all, aes(x=pop, y=avg_dxy, fill=site)) +
  geom_violin(alpha=0.65, trim=FALSE, linewidth=0.3) +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA, linewidth=0.5) +
  scale_fill_manual(values=pop_colors) +
    theme_minimal(base_size=14) +
  labs(title="Absolute Divergence (Dxy)", x="", y="Average Dxy per window") +
  theme(legend.position="none", panel.grid.major.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~site, scales="free")

ggsave(file.path(out_dir, "Plot_5_Dxy_Violin_90.pdf"), plot=p5, width=12, height=6)

p6 <- ggplot(dxy_all, aes(x=window_index, y=avg_dxy, color=pop)) +
  geom_point(alpha=0.35, size=0.8) +
  geom_smooth(method="loess", span=0.15, se=FALSE, linewidth=0.8, color="black") +
  facet_wrap(~pop, ncol=1, scales="free_x") +
  theme_minimal(base_size=13) +
  labs(title="Absolute Divergence (Dxy) - Sliding Windows (50kb)",
       x="Window index (ordered by contig/position)",
       y="Average Dxy per window") +
  theme(legend.position="none", strip.text=element_text(size=13, face="bold"))

ggsave(file.path(out_dir, "Plot_6_Dxy_Scatter_90.pdf"), plot=p6, width=10, height=12)

message("Done. Plots saved to: ", out_dir)
