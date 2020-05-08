#---------------------------------------------------------------------------------------------------------------#
# Compute Immune Checkpoint genes expression
#---------------------------------------------------------------------------------------------------------------#

computation.ICB_genes <- function(RNA.tpm, ....){
  
  # Extract position genes for GZMA and PRF1
  tmp = match(c("CD274","CTLA4","PDCD1"), colnames(RNA.tpm))
  
  # PDL-1 calculation
  PDL1_expr = RNA.tpm[,tmp[1]] ; names(PDL1_expr) <- rownames(RNA.tpm)

  # CTLA-4 calculation
  CTLA4_expr = RNA.tpm[,tmp[2]] ; names(CTLA4_expr) <- rownames(RNA.tpm)
  
  # PD-1 calculation
  PD1_expr = RNA.tpm[,tmp[3]] ; names(PD1_expr) <- rownames(RNA.tpm)

  
  ICB_genes_expr <- list(PDL1 = PDL1_expr , CTLA4 = CTLA4_expr , PD1 =PD1_expr )
  cat("Immune Checkpoint Genes expression computed","\n")
  return(ICB_genes_expr)
}