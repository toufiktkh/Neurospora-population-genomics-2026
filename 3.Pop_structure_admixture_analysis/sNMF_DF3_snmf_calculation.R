
if (!requireNamespace("LEA", quietly = TRUE)) {
    user_lib <- Sys.getenv("R_LIBS_USER")
    if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
    .libPaths(c(user_lib, .libPaths()))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", lib = user_lib, repos = "https://cloud.r-project.org")
    }
    BiocManager::install("LEA", lib = user_lib, update = FALSE, ask = FALSE)
}
library(LEA)

vcf_path  <- "/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_snps_only_clone_corr.vcf"
geno_path <- "/shared/home/toufiktakhi/neurosporagbs/admixture_analysis/DF3_snps_only_clone_corr.geno"
out_dir   <- "/shared/home/toufiktakhi/neurosporagbs/admixture_analysis/"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

#Convert VCF to GENO
vcf2geno(vcf_path, geno_path)

# CLAUDE (Sonnet 4.3 Fix) PLINK outputs diploid format, since sNMF haploid mode rejects '2'
geno_lines <- readLines(geno_path)
geno_lines <- gsub("2", "1", geno_lines)
writeLines(geno_lines, geno_path)
cat("Geno file fixed: 2 → 1 for haploid compatibility\n")

if (any(grepl("2", geno_lines))) {
    stop("Geno file still contains '2' values — fix failed.")
} else {
    cat("Verification passed: no '2' values remaining\n")
}

# Remove old project if it exists to avoid conflict
if (file.exists("DF3_snps_only_snmf.snmfProject")) {
    remove.snmfProject("DF3_snps_only_snmf.snmfProject")
    cat("Old snmfProject removed\n")
}

#Run sNMF
project_df3 <- snmf(geno_path,
                    K = 1:5,
                    ploidy = 1,
                    entropy = TRUE,
                    repetitions = 10,
                    project = "new",
                    CPU = 10,
                    seed = 42)

# Cross-entropy plot
pdf(file.path(out_dir, "cross_entropy_DF3_snps_only.pdf"), width = 8, height = 6)
plot(project_df3, col = "steelblue", pch = 19, cex = 1.2, main = "Cross-Entropy DF3")
dev.off()

cat("Done. Plot saved to: admixture_analysis/cross_entropy_DF3_snps_only.pdf\n")
