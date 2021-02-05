###########################################################################################
# Script to create ObtainCYT function:
# Computing the geometric mean of literature genes provided by the authors:
# Geometric-mean of [TPM + 0.01] cytolytic genes
###########################################################################################

ObtainCYT <- function(gene.tpm){

  # ************
  # Literature genes
  CYT_read <- c("GZMA", "PRF1")
  match_genes <- match(CYT_read, rownames(gene.tpm))

  if (anyNA(match_genes)){
    warning(paste0("differenty named or missing signature genes : \n", CYT_read[!CYT_read %in% rownames(gene.tpm)]))
    match_genes <- na.omit(match_genes)
  }

  CYT.score <- vector("numeric", length = ncol(gene.tpm)) ; names(CYT.score) <- colnames(gene.tpm)

  # Subset gene expression data
  sub_gene.tpm <- gene.tpm[match_genes,]

  # Calculated as geometric mean (so-called log-average) [TPM, 0.01 offset]
  geom_mean <- as.matrix(apply(sub_gene.tpm + 0.01, 2, function(X) exp(mean(log(X)))))

  CYT.score <- data.frame(CYT = geom_mean, check.names = FALSE)

  message("Cytolytic Activity score computed","\n")
  return(CYT.score)
}
