###########################################################################################
# Script to create ObtainTIDE function:
# Compute TIDE using published model (tidepy)
# Input: log2(TPM + 1) + quantile normalization + mean centralization (author's suggestions)
###########################################################################################

ObtainTIDE <- function(gene.tpm, Cancer){
  
  # ************
  # Packages
  library(preprocessCore)
  
  # ************
  # Pre-process 
  
  ## log2(TPM + 1)
  log2.gene.tpm <- as.data.frame(log2(gene.tpm + 1))
  
  ## Quantile normalization
  log2.gene.tpm.norm <- normalize.quantiles(as.matrix(log2.gene.tpm))
  dimnames(log2.gene.tpm.norm) <- dimnames(log2.gene.tpm)
  
  ## Mean-centralization: substract gene average across all conditions
  average.gene <- rowMeans(log2.gene.tpm.norm)
  log2.gene.tpm.norm.all.conditions <- sweep(log2.gene.tpm.norm, 1, average.gene, FUN = "-")
  
  write.table(log2.gene.tpm.norm.all.conditions, 
              file = paste0("/Users/Oscar/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/data_processed_TIDE/log2mas1TPM_norm_", Cancer,".txt"), sep = "\t")
  
  # Compute TIDE
  if (Cancer == "SKCM") {
    system(paste0("/anaconda3/bin/tidepy ", 
                  "/Users/Oscar/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/data_processed_TIDE/log2mas1TPM_norm_", 
                  Cancer,".txt", " -o /Users/Oscar/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/data_processed_TIDE/output_TIDE_norm_",Cancer,".txt -c Melanoma"))
  }else{
    system(paste0("/anaconda3/bin/tidepy ", 
                  "/Users/Oscar/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/data_processed_TIDE/log2mas1TPM_norm_", 
                  Cancer,".txt", " -o /Users/Oscar/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/data_processed_TIDE/output_TIDE_norm_",Cancer,".txt -c Other"))
  }
  
  TIDE.table.norm <- read.table(file = paste0("/Users/Oscar/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/data_processed_TIDE/output_TIDE_norm_",Cancer,".txt"), sep = "\t", header = T, row.names = 1)
  TIDE.score <- data.frame(TIDE = TIDE.table.norm[,"TIDE", drop = FALSE], check.names = F)
  
  cat("TIDE score computed","\n")
  
  return(TIDE.score)
}

# ***************