###########################################################################################
# Script to create ObtainExpandedIS function:
# After performance of quantile normalization, a log10 transformation was applied, and
# signature scores were calculated by averaging of the included genes for the expanded immune signature
###########################################################################################

ObtainExpandedImmuneAyers <- function(gene.tpm){

  # ************
  # Literature info: cytolytic markers, Expanded Immune 18-genes signature (Ayers):
  ExpandedImmune_Ayers_read <- c("GZMB", "GZMK", "CXCR6", "CXCL10", "CXCL13", "CCL5", "STAT1","CD3D", "CD3E",
                                 "CD2", "IL2RG" , "NKG7", "HLA-E", "CIITA","HLA-DRA", "LAG3", "IDO1", "TAGAP")

  # ***************
  # Gene expression data just for IFNy signature (Ayers)
  match_genes <- match(ExpandedImmune_Ayers_read, rownames(gene.tpm))

  if (anyNA(match_genes)){
    warning(paste0("differenty named or missing signature genes : \n", ExpandedImmune_Ayers_read[!ExpandedImmune_Ayers_read %in% rownames(gene.tpm)]))
    match_genes <- na.omit(match_genes)
  }

  ExpandedImmune.score <- vector("numeric", length = ncol(gene.tpm)) ; names(ExpandedImmune.score) <- colnames(gene.tpm)

  # Log2 transformation:
  log2.gene.tpm <- log2(gene.tpm + 1)

  # Subset Ayers IFNy signature genes
  sub_log2.gene.tpm  <- log2.gene.tpm[match_genes, ]

  # Average of the included genes for the Expanded Immune signature
  ExpandedImmune.score <- data.frame(ExpandedImmune = apply(sub_log2.gene.tpm, 2, mean), check.names = FALSE)

  cat("Expanded Immune (Ayers) score computed","\n")
  return(ExpandedImmune.score)
}

# ***************

