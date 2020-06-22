###########################################################################################
# Script to create ObtainRohISS function:
# Computing the geometric mean of literature genes provided by the authors:
# Geometric-mean of TPM-normalized coounts of immune score genes
###########################################################################################

ObtainDavoliIS <- function(gene.tpm){
  
  # ************
  # Literature info: Genes from immune infiltrate signature, genes specific for cytotoxic CD8+ T cells and NK cells (Davoli)
  IS_Davoli_read <- c("CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK")
  
  # ***************
  # Gene expression data just for Davoli immune signature score
  match_genes <- match(IS_Davoli_read, rownames(gene.tpm))
  
  if (anyNA(match_genes)){
    warning(paste0("differenty named or missing signature genes : \n", IS_Davoli_read[!IS_Davoli_read %in% rownames(gene.tpm)]))
    match_genes <- na.omit(match_genes)
  }
  
  IS_Davoli.score <- vector("numeric", length = ncol(gene.tpm)) ; names(IS_Davoli.score) <- colnames(gene.tpm)
  
  # Log2 transformation:
  log2.gene.tpm <- log2(gene.tpm + 1)
  
  # Keep matrix for Immune Signature Davoli genes
  sub_log2.gene.tpm <- log2.gene.tpm[match_genes,]

  # Calculate rank position for each gene across samples
  ranks_sub_log2.gene.tpm <- apply(sub_log2.gene.tpm, 1, rank)
  
  # Get normalized rank by divided 
  ranks_sub_log2.gene.tpm.norm <- (ranks_sub_log2.gene.tpm - 1)/(nrow(ranks_sub_log2.gene.tpm) - 1)
  
  # Average of the expression value of all the genes within-sample
  average_ranks_sub_log2.gene.tpm.norm<- apply(ranks_sub_log2.gene.tpm.norm, 1, mean)

  # For the pan-cancer analysis, the immune signature score was normalized across the samples of the same tumor type. 
  IS_Davoli.score <- data.frame(IS_Davoli = average_ranks_sub_log2.gene.tpm.norm, check.names = FALSE)
  
  cat("IS_Davoli score computed","\n")
  return(IS_Davoli.score)
}

# ***************

