if (!require(ggplot2)) install.packages("ggplot2", repos="https://cloud.r-project.org")
if (!require(dplyr)) install.packages("dplyr", repos="https://cloud.r-project.org")
library(ggplot2)
library(dplyr)

out_dir <- "~/neurosporagbs/pixy_plots/"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

pop_colors <- c(
  "Domefire_DF1" = "#E63946",
  "Domefire_DF3" = "#457B9D",
  "Villeveyrac"  = "#2A9D8F",
  "La_Grande_Motte"          = "#E9C46A"
)

load_pi <- function(path, pop_label) {
  if (!file.exists(path)) return(data.frame())
  df <- read.table(path, header=TRUE, stringsAsFactors=FALSE, sep="\t")
  df$pop <- pop_label
  df$window_mid <- (df$window_pos_1 + df$window_pos_2) / 2
  df <- df[!is.na(df$avg_pi), ]
  return(df)
}

# PI calculations

pi_df1  <- load_pi("~/neurosporagbs/pca_splitstree_domefire/all_sites/pixy_output_DF1_90/pixy_pi.txt", "Domefire_DF1")
pi_df3  <- load_pi("~/neurosporagbs/pca_splitstree_domefire/all_sites/pixy_output_DF3_90/pixy_pi.txt", "Domefire_DF3")
pi_vill <- load_pi("~/neurosporagbs/pca_splitstree_vill/all_sites/pixy_output_90/pixy_pi.txt", "Villeveyrac")
pi_lgm  <- load_pi("~/neurosporagbs/pca_splitstree_LGM/all_sites/pixy_output/pixy_pi.txt", "La_Grande_Motte")

pi_all <- rbind(pi_df1, pi_df3, pi_vill, pi_lgm)
pi_all$pop <- factor(pi_all$pop, levels=c("Domefire_DF1", "Domefire_DF3", "Villeveyrac", "La_Grande_Motte"))

pi_all <- pi_all %>%
  arrange(pop, chromosome, window_pos_1) %>%
  group_by(pop) %>%
  mutate(window_index = row_number()) %>%
  ungroup()

load_tajd <- function(path, pop_label) {
  if (!file.exists(path)) return(data.frame())
  df <- read.table(path, header=TRUE, stringsAsFactors=FALSE, sep="\t")
  df$pop <- pop_label
  df$window_mid <- (df$window_pos_1 + df$window_pos_2) / 2
  df <- df[!is.na(df$tajima_d), ]
  return(df)
}

# Calculating Tajima's D

taj_df1  <- load_tajd("~/neurosporagbs/pca_splitstree_domefire/all_sites/pixy_output_DF1_90/pixy_tajima_d.txt", "Domefire_DF1")
taj_df3  <- load_tajd("~/neurosporagbs/pca_splitstree_domefire/all_sites/pixy_output_DF3_90/pixy_tajima_d.txt", "Domefire_DF3")
taj_vill <- load_tajd("~/neurosporagbs/pca_splitstree_vill/all_sites/pixy_output_90/pixy_tajima_d.txt", "Villeveyrac")
taj_lgm  <- load_tajd("~/neurosporagbs/pca_splitstree_LGM/all_sites/pixy_output/pixy_tajima_d.txt", "La_Grande_Motte")

taj_all <- rbind(taj_df1, taj_df3, taj_vill, taj_lgm)
taj_all$pop <- factor(taj_all$pop, levels=c("Domefire_DF1", "Domefire_DF3", "Villeveyrac", "La_Grande_Motte"))

taj_all <- taj_all %>%
  arrange(pop, chromosome, window_pos_1) %>%
  group_by(pop) %>%
  mutate(window_index = row_number()) %>%
  ungroup()

load_dxy <- function(path, pop_label) {
  if (!file.exists(path)) return(data.frame())
  df <- read.table(path, header=TRUE, stringsAsFactors=FALSE, sep="\t")
  if (nrow(df) == 0) return(data.frame())
  df$pop <- pop_label
  df$window_mid <- (df$window_pos_1 + df$window_pos_2) / 2
  df <- df[!is.na(df$avg_dxy), ]
  return(df)
}

# Calculating Dxy

dxy_df1  <- load_dxy("~/neurosporagbs/pca_splitstree_domefire/all_sites/pixy_output_DF1_90/pixy_dxy.txt", "Domefire_DF1")
dxy_df3  <- load_dxy("~/neurosporagbs/pca_splitstree_domefire/all_sites/pixy_output_DF3_90/pixy_dxy.txt", "Domefire_DF3")
dxy_vill <- load_dxy("~/neurosporagbs/pca_splitstree_vill/all_sites/pixy_output_90/pixy_dxy.txt", "Villeveyrac")
dxy_lgm  <- load_dxy("~/neurosporagbs/pca_splitstree_LGM/all_sites/pixy_output/pixy_dxy.txt", "La_Grande_Motte")

dxy_all <- rbind(dxy_df1, dxy_df3, dxy_vill, dxy_lgm)
dxy_all$pop <- factor(dxy_all$pop, levels=c("Domefire_DF1", "Domefire_DF3", "Villeveyrac", "La_Grande_Motte"))

dxy_all <- dxy_all %>%
  arrange(pop, chromosome, window_pos_1) %>%
  group_by(pop) %>%
  mutate(window_index = row_number()) %>%
  ungroup()

# Plotting

p1 <- ggplot(pi_all, aes(x=pop, y=avg_pi, fill=pop)) +
  geom_violin(alpha=0.65, trim=FALSE, linewidth=0.3) +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA, linewidth=0.5) +
  scale_fill_manual(values=pop_colors) +
  theme_minimal(base_size=14) +
  labs(title="Nucleotide Diversity (Pi)", x="", y="Average Pi per window") +
  theme(legend.position="none", panel.grid.major.x=element_blank())

ggsave(file.path(out_dir, "Plot_1_Pi_Violin_90.pdf"), plot=p1, width=10, height=6)

p2 <- ggplot(taj_all, aes(x=pop, y=tajima_d, fill=pop)) +
  geom_violin(alpha=0.65, trim=FALSE, linewidth=0.3) +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA, linewidth=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=0.8) +
  scale_fill_manual(values=pop_colors) +
  theme_minimal(base_size=14) +
  labs(title="Tajima's D", x="", y="Tajima's D") +
  theme(legend.position="none", panel.grid.major.x=element_blank())

ggsave(file.path(out_dir, "Plot_2_TajimasD_Violin_90.pdf"), plot=p2, width=10, height=6)

p3 <- ggplot(taj_all, aes(x=window_index, y=tajima_d, color=pop)) +
  geom_point(alpha=0.35, size=0.8) +
  geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=0.7) +
  geom_smooth(method="loess", span=0.15, se=FALSE, linewidth=0.8, color="black") +
  scale_color_manual(values=pop_colors) +
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
  scale_color_manual(values=pop_colors) +
  facet_wrap(~pop, ncol=1, scales="free_x") +
  theme_minimal(base_size=13) +
  labs(title="Nucleotide Diversity (Pi) - Sliding Windows (50kb)",
       x="Window index (ordered by contig/position)",
       y="Average Pi per window") +
  theme(legend.position="none", strip.text=element_text(size=13, face="bold"))

ggsave(file.path(out_dir, "Plot_4_Pi_Scatter_90.pdf"), plot=p4, width=10, height=12)

p5 <- ggplot(dxy_all, aes(x=pop, y=avg_dxy, fill=pop)) +
  geom_violin(alpha=0.65, trim=FALSE, linewidth=0.3) +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA, linewidth=0.5) +
  scale_fill_manual(values=pop_colors) +
  theme_minimal(base_size=14) +
  labs(title="Absolute Divergence (Dxy)", x="", y="Average Dxy per window") +
  theme(legend.position="none", panel.grid.major.x=element_blank())

ggsave(file.path(out_dir, "Plot_5_Dxy_Violin_90.pdf"), plot=p5, width=10, height=6)

p6 <- ggplot(dxy_all, aes(x=window_index, y=avg_dxy, color=pop)) +
  geom_point(alpha=0.35, size=0.8) +
  geom_smooth(method="loess", span=0.15, se=FALSE, linewidth=0.8, color="black") +
  scale_color_manual(values=pop_colors) +
  facet_wrap(~pop, ncol=1, scales="free_x") +
  theme_minimal(base_size=13) +
  labs(title="Absolute Divergence (Dxy) - Sliding Windows (50kb)",
       x="Window index (ordered by contig/position)",
       y="Average Dxy per window") +
  theme(legend.position="none", strip.text=element_text(size=13, face="bold"))

ggsave(file.path(out_dir, "Plot_6_Dxy_Scatter_90.pdf"), plot=p6, width=10, height=12)

message("Done. Plots saved to: ", out_dir)
