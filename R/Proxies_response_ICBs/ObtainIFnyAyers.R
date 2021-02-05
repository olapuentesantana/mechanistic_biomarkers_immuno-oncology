###########################################################################################
# Script to create ObtainIFnyAyers function:
# After performance of quantile normalization, a log10 transformation was applied, and 
# signature scores were calculated by averaging of the included genes for the IFN-γ
###########################################################################################

ObtainIFnyAyers <- function(gene.tpm){
  
  # ************
  # Literature info: cytolytic markers, IFNy 6-genes signature (Ayers):
  IFNy_Ayers_read <- c("IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA") 
  
  # ***************
  # Gene expression data just for IFNy signature (Ayers)
  match_genes <- match(IFNy_Ayers_read, rownames(gene.tpm))

  if (anyNA(match_genes)){
    warning(paste0("differenty named or missing signature genes : \n", IFNy_Ayers_read[!IFNy_Ayers_read %in% rownames(gene.tpm)]))
    match_genes <- na.omit(match_genes)
  }
  
  IFNy.score <- vector("numeric", length = ncol(gene.tpm)) ; names(IFNy.score) <- colnames(gene.tpm)

  # Log2 transformation:
  log2.gene.tpm <- log2(gene.tpm + 1)
  
  # Subset Ayers IFNy signature genes
  sub_log2.gene.tpm  <- log2.gene.tpm[match_genes, ]
  
  # Average of the included genes for the IFN-y signature
  IFNy.score <- data.frame(IFny = apply(sub_log2.gene.tpm, 2, mean), check.names = FALSE)
    
  cat("IFNy (Ayers) score computed","\n")
  return(IFNy.score)
}

# ***************




