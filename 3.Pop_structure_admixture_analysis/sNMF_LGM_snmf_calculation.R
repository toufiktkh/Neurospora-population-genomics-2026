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

# since LGM isolates are most probably  N. tetrasperma, which is DIPLOID (need to add ploidy = 2)

vcf_path  <- "/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_snps_only_clone_corr.vcf"
geno_path <- "/shared/home/toufiktakhi/neurosporagbs/admixture_analysis/LGM_snps_only_clone_corr.geno"
out_dir   <- "/shared/home/toufiktakhi/neurosporagbs/admixture_analysis/"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

#Convert VCF to GENO
vcf2geno(vcf_path, geno_path)

if (!file.exists(geno_path)) {
    stop("vcf2geno failed — .geno file not found. Check VCF path.")
}
cat("Geno file created successfully\n")

if (file.exists("LGM_snps_only_snmf.snmfProject")) {
    remove.snmfProject("LGM_snps_only_snmf.snmfProject")
    cat("Old snmfProject removed\n")
}

# Run snmf analysis
project_gm <- snmf(geno_path,
                   K = 1:10,
                   ploidy = 2,
                   entropy = TRUE,
                   repetitions = 10,
                   project = "new",
                   CPU = 10,
                   seed = 42)

# Cross-entropy plot
pdf(file.path(out_dir, "cross_entropy_LGM_snps_only_clone_corr_V1.pdf"), width = 8, height = 6)
plot(project_gm, col = "steelblue", pch = 19, cex = 1.2, main = "Cross-Entropy LGM")
dev.off()

cat("Done. Plot saved to: admixture_analysis/cross_entropy_LGM_snps_only_clone_corr_V1.pdf\n")
