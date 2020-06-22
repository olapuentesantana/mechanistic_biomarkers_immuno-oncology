###########################################################################################
# Script to create ObtainRohISS function:
# Calculated a PCA score based on PC1 that approximates anaverage Z-score based on the 
# log2 intensities for each of the 12 chemokinegenes in the gene expression signature.
###########################################################################################

Obtain12chemokine <- function(gene.tpm){
  
  # ************
  # Literature info: 12 chemokine genes
  Chemokine_read <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21",
                      "CXCL9", "CXCL10", "CXCL11", "CXCL13")
  
  # ***************
  # Gene expression data just for 12-chemokine signature
  match_genes <- match(Chemokine_read, rownames(gene.tpm))
  
  if (anyNA(match_genes)){
    warning(paste0("differenty named or missing signature genes : \n", Chemokine_read[!Chemokine_read %in% rownames(gene.tpm)]))
    match_genes <- na.omit(match_genes)
  }
  
  Chemokine.score <- vector("numeric", length = ncol(gene.tpm)) ; names(Chemokine.score) <- colnames(gene.tpm)

  # log2(TPM +1)
  gene_expr <- log2(gene.tpm + 1)
  chemokine.pca <- prcomp(t(gene_expr[match_genes,]), center = T,scale = T) # Z-score calculated whitin prcomp
  Chemokine.score <- data.frame(chemokine = chemokine.pca$x[,1], check.names = FALSE)
  
  cat("12-chemokine (Messina) score computed","\n")
  return(Chemokine.score)
}
