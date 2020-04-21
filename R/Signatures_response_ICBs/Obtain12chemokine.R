###########################################################################################
# Script to create ObtainRohISS function:
# Calculated a PCA score based on PC1 that approximates anaverage Z-score based on the 
# log2 intensities for each of the 12 chemokinegenes in the gene expression signature.
###########################################################################################

Obtain12chemokine <- function(gene_expr){
  
  # ************
  # Literature info: 12 chemokine genes
  Chemokine_read <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21",
                      "CXCL9", "CXCL10", "CXCL11", "CXCL13")
  
  # ***************
  # Gene expression data just for 12-chemokine signature
  match_genes <- na.omit(match(Chemokine_read, rownames(gene_expr)))
  
  if (anyNA(match_genes)){cat("MISSING GENES", "\n") ; warning("Be careful, some genes are missing"); cat(!which(Chemokine_read %in% rownames(gene_expr)))}
  
  Chemokine.score <- vector("numeric", length = ncol(gene_expr)) ; names(Chemokine.score) <- colnames(gene_expr)

  chemokine.pca <- prcomp(t(gene_expr[match_genes,]), center = T,scale = T)
  Chemokine.score <- data.frame(chemokine = chemokine.pca$x[,1], check.names = FALSE)
  
  cat("12-chemokine (Messina) score computed","\n")
  return(Chemokine.score)
}
