###########################################################################################
# Script to create ObtainRohISS function:
# Computing the geometric mean of literature genes provided by the authors:
# Geometric-mean of TPM-normalized coounts of immune score genes
###########################################################################################

ObtainDavoliIS <- function(gene_expr){
  
  # ************
  # Literature info: Genes from immune infiltrate signature, genes specific for cytotoxic CD8+ T cells and NK cells (Davoli)
  IS_Davoli_read <- c("CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK")
  
  # ***************
  # Gene expression data just for Davoli immune signature score
  match_genes <- na.omit(match(IS_Davoli_read, rownames(gene_expr)))

  if (anyNA(match_genes)){
    warning("Be careful, some genes are missing")
    cat("MISSING GENES", "\n") 
    cat(!which(IS_Davoli_read %in% rownames(gene_expr)))
  }
  
  IS_Davoli.score <- vector("numeric", length = ncol(gene_expr)) ; names(IS_Davoli.score) <- colnames(gene_expr)
  
  # Keep matrix for Immune Signature Davoli genes
  tmp <- gene_expr[match_genes,]
  
  # Calculate rank position for each gene across samples
  ranks_tmp <- apply(tmp, 1, rank)
  
  # Get normalized rank by divided 
  ranks_tmp.norm <- (ranks_tmp - 1)/(nrow(ranks_tmp) - 1)
  
  # Average of the expression value of all the genes in the corresponding signature.
  average_ranks_tmp <- apply(ranks_tmp.norm, 1, mean)

  # For the pan-cancer analysis, the immune signature score was normalized across the samples of the same tumor type. 
  IS_Davoli.score <- data.frame(IS_Davoli = average_ranks_tmp, check.names = FALSE)
  
  cat("IS_Davoli score computed","\n")
  return(IS_Davoli.score)
}

# ***************

