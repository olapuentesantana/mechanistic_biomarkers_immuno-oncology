###########################################################################################
# Script to create ObtainExpandedIS function:
# After performance of quantile normalization, a log10 transformation was applied, and 
# signature scores were calculated by averaging of the included genes for the expanded immune signature
###########################################################################################

ObtainExpandedImmuneAyers <- function(gene_expr){
  
  # ************
  # Literature info: cytolytic markers, Expanded Immune 18-genes signature (Ayers):
  ExpandedImmune_Ayers_read <- c("GZMB", "GZMK", "CXCR6", "CXCL10", "CXCL13", "CCL5", "STAT1","CD3D", "CD3E",
                                 "CD2", "IL2RG" , "NKG7", "HLA-E", "CIITA"," HLA-DRA", "LAG3", "IDO1", "TAGAP")  
  
  # ***************
  # Gene expression data just for IFNy signature (Ayers)
  match_genes <- na.omit(match(ExpandedImmune_Ayers_read, rownames(gene_expr)))
  tmp.gene_expr <- gene_expr[match_genes, ]
  
  if (anyNA(match_genes)){
    warning("Be careful, some genes are missing")
    cat("MISSING GENES", "\n") 
    cat(!which(ExpandedImmune_Ayers_read %in% rownames(gene_expr)))
  }
  
  ExpandedImmune.score <- vector("numeric", length = ncol(gene_expr)) ; names(ExpandedImmune.score) <- colnames(gene_expr)
  
  # log10 transformation
  log10.tmp.gene_expr <- log10(tmp.gene_expr + 1)
  
  # Average of the included genes for the IFN-y signature
  ExpandedImmune.score <- data.frame(ExpandedImmune = apply(log10.tmp.gene_expr, 2, mean), check.names = FALSE)
  
  
  cat("Expanded Immune (Ayers) score computed","\n")
  return(ExpandedImmune.score)
}

# ***************

