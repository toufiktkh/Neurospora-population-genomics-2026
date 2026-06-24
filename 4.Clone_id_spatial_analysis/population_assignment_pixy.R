library(LEA)

setwd("/shared/home/toufiktakhi/neurosporagbs/")

process_snmf_popmap <- function(project_path, K_value, sample_id_file, output_txt_name) {
    
    cat("\n------------------------------------------------------------\n")
    cat("Processing:", basename(project_path), "at K =", K_value, "\n")
    
# Extract properly the file and the name of the file
    target_dir <- dirname(output_txt_name)   # ex: "pca_splitstree_domefire/all_sites"
    file_base <- basename(output_txt_name)   # ex: "pixy_popmap_DF1.txt"

# Create the file if it doesn't exsit on the cluster
    if (!dir.exists(target_dir)) {
        dir.create(target_dir, recursive = TRUE)
    }
    
    diagnostic_path <- file.path(target_dir, paste0("diagnostic_", file_base))
    
# Loading the snmf project
    project <- load.snmfProject(project_path)
    
# Find the best run based on the lowest cross-entropy value
    best_run <- which.min(cross.entropy(project, K = K_value))
    
# Extraction the Q matrix
    q_matrix <- Q(project, K = K_value, run = best_run)
    
# Changing sample IDs to the corresponding isolate ID extracted from the clone-corrected VCF in a .txt file
    if (!file.exists(sample_id_file)) {
        stop(paste("ERROR: Sample ID file not found:", sample_id_file))
    }
    sample_ids <- read.table(sample_id_file, header = FALSE)[, 1]
    
# CLAUDE solution to err in names 
    if (length(sample_ids) != nrow(q_matrix)) {
        stop(paste("ROW MISMATCH! The sNMF matrix has", nrow(q_matrix), 
                   "samples, but your ID file has", length(sample_ids), "samples."))
    }
    
# Constructing the dataframe
    ancestry_df <- data.frame(Sample_ID = sample_ids)
    
    for (i in 1:K_value) {
        ancestry_df[paste0("Cluster_", i)] <- q_matrix[, i]
    }
    
# The most important part, assigning population based on the 90% strict threshold
    if (K_value == 1) {
        ancestry_df$Final_Assignment <- "Population_1"
    } else {
        max_val <- apply(q_matrix, 1, max)
        max_cluster <- max.col(q_matrix)
        
        ancestry_df$Final_Assignment <- ifelse(max_val >= 0.90, 
                                               paste0("Population_", max_cluster), 
                                               "Admixed/Exclude")
    }
    
    write.table(ancestry_df, diagnostic_path, 
                sep="\t", row.names=FALSE, quote=FALSE)
    
#Export the popmap files for PIXY (Exclusion of isolates using based on the ancestry proportions with the 90% threshold)
    pixy_popmap <- ancestry_df[ancestry_df$Final_Assignment != "Admixed/Exclude", c("Sample_ID", "Final_Assignment")]
    
    write.table(pixy_popmap, output_txt_name, 
                sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
    
    cat("Success! Generated files:\n")
    cat("  - Diagnostic Table: ", diagnostic_path, "\n", sep="")
    cat("  - Pixy Popmap:      ", output_txt_name, " (", nrow(pixy_popmap), " samples kept)\n", sep="")
}

# 1. DOMEFIRE DF1 (Clone-Corrected, K=3)
process_snmf_popmap(
    project_path = "admixture_analysis/DF1_snps_only_clone_corr.snmfProject",
    K_value = 3,
    sample_id_file = "admixture_analysis/sample_ID_DF1.txt", 
    output_txt_name = "pca_splitstree_domefire/all_sites/pixy_popmap_DF1_90_.txt"
)

# 2. DOMEFIRE DF3 (Not Clone-Corrected, since there are no clones, K=1)
process_snmf_popmap(
    project_path = "admixture_analysis/DF3_snps_only_clone_corr.snmfProject",
    K_value = 1,
    sample_id_file = "admixture_analysis/sample_ID_DF3.txt",
    output_txt_name = "pca_splitstree_domefire/all_sites/pixy_popmap_DF3_90_.txt"
)

# 3. VILLEVEYRAC (Clone-Corrected, K=3) 
process_snmf_popmap(
    project_path = "admixture_analysis/Vill_snps_only_clone_corr.snmfProject",
    K_value = 3,
    sample_id_file = "admixture_analysis/sample_ID_Vill.txt",
    output_txt_name = "pca_splitstree_vill/all_sites/pixy_popmap_Vill_90_.txt"
)

# 4. LA GRANDE MOTTE (Not Clone-Corrected, since there are no clones, K=3)
process_snmf_popmap(
    project_path = "admixture_analysis/LGM_snps_only_clone_corr.snmfProject",
    K_value = 3,
    sample_id_file = "admixture_analysis/sample_ID_LGM.txt",
    output_txt_name = "pca_splitstree_LGM/all_sites/pixy_popmap_LGM_90_.txt"
)

cat("\nAll populations processed successfully!\n")
