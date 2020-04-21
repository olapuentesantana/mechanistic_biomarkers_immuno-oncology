###########################################################################################
# Script to create ObtainIFnyAyers function:
# After performance of quantile normalization, a log10 transformation was applied, and 
# signature scores were calculated by averaging of the included genes for the IFN-γ
###########################################################################################

ObtainIFnyAyers <- function(gene_expr){
  
  # ************
  # Literature info: cytolytic markers, IFNy 6-genes signature (Ayers):
  IFNy_Ayers_read <- c("IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA") 
  
  # ***************
  # Gene expression data just for IFNy signature (Ayers)
  match_genes <- na.omit(match(IFNy_Ayers_read, rownames(gene_expr)))
  tmp.gene_expr <- gene_expr[match_genes, ]
  
  if (anyNA(match_genes)){
    warning("Be careful, some genes are missing")
    cat("MISSING GENES", "\n") 
    cat(!which(IFNy_Ayers_read %in% rownames(gene_expr)))
  }
  
  IFNy.score <- vector("numeric", length = ncol(gene_expr)) ; names(IFNy.score) <- colnames(gene_expr)

  # log10 transformation
  log10.tmp.gene_expr <- log10(tmp.gene_expr + 1)
  
  # Average of the included genes for the IFN-y signature
  IFNy.score <- data.frame(IFny = apply(log10.tmp.gene_expr, 2, mean), check.names = FALSE)
    
  
  cat("IFNy (Ayers) score computed","\n")
  return(IFNy.score)
}

# ***************




