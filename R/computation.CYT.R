#---------------------------------------------------------------------------------------------------------------#
# Compute Cytolytic Activity
#---------------------------------------------------------------------------------------------------------------#

computation.CYT <- function(RNA.tpm, ....){

  # Extract position genes for GZMA and PRF1
  tmp = match(c("GZMA", "PRF1"), colnames(RNA.tpm))
  
  # Geometric mean of GZMA and PRF1(TPM,0.01 offset)
  RNA.tpm_0.01_offset <- RNA.tpm + 0.01
  GZMA_expr <- RNA.tpm_0.01_offset[,tmp[1]] ; names(GZMA_expr) <- rownames(RNA.tpm)
  PRF1_expr <- RNA.tpm_0.01_offset[,tmp[2]] ; names(PRF1_expr) <- rownames(RNA.tpm)
  CYT <- sqrt(GZMA_expr * PRF1_expr) 

  cat("Cytolytic Activity score computed","\n")
  return(CYT)
}