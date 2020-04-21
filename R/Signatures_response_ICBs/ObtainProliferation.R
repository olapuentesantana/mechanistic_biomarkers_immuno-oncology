###########################################################################################
# Script to create ObtainRohISS function:
# Computing the geometric mean of literature genes provided by the authors:
# Geometric-mean of TPM-normalized coounts of immune score genes
###########################################################################################

ObtainProliferation <- function(gene_expr){
  
  # ************
  # Literature info: Genes from cell cycle signature, genes associated with proliferation not frequently amplified/deleted (Davoli)
  Proliferation_read <- c("CENPE", "CCNA2", "CCNB2", "MCM6", "CCNF", "BUB1", "CDC20", "CDC6", "CDK1", "PLK1")
  
  # ***************
  # Gene expression data just for Proliferation signature score
  match_genes <- na.omit(match(Proliferation_read, rownames(gene_expr)))
  
  if (anyNA(match_genes)){
    warning("Be careful, some genes are missing")
    cat("MISSING GENES", "\n") 
    cat(!which(Proliferation_read %in% rownames(gene_expr)))
  }
  
  Proliferation.score <- vector("numeric", length = ncol(gene_expr)) ; names(Proliferation.score) <- colnames(gene_expr)
  
  # Keep matrix for Immune Signature Davoli genes
  tmp <- gene_expr[match_genes,]
  
  # Calculate rank position for each gene across samples
  ranks_tmp <- apply(tmp, 1, rank)
  
  # Get normalized rank by divided 
  ranks_tmp.norm <- (ranks_tmp - 1)/(nrow(ranks_tmp) - 1)
  
  # Average of the expression value of all the genes in the corresponding signature.
  average_ranks_tmp <- apply(ranks_tmp.norm, 1, mean)
  
  # For the pan-cancer analysis, the immune signature score was normalized across the samples of the same tumor type. 
  Proliferation.score <- data.frame(Proliferation = average_ranks_tmp, check.names = FALSE)
  
  cat("Proliferation score computed","\n")
  return(Proliferation.score)
}

# ***************
