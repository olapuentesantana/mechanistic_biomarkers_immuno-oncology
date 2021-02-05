###########################################################################################
# Script to create ObtainRohISS function:
# Computing the geometric mean of literature genes provided by the authors:
# Geometric-mean of TPM-normalized coounts of immune score genes
###########################################################################################

ObtainProliferation <- function(gene.tpm){
  
  # ************
  # Literature info: Genes from cell cycle signature, genes associated with proliferation not frequently amplified/deleted (Davoli)
  Proliferation_read <- c("CENPE", "CCNA2", "CCNB2", "MCM6", "CCNF", "BUB1", "CDC20", "CDC6", "CDK1", "PLK1")
  
  # ***************
  # Gene expression data just for Proliferation signature score
  match_genes <- match(Proliferation_read, rownames(gene.tpm))
  
  if (anyNA(match_genes)){
    warning(paste0("differenty named or missing signature genes : \n", Proliferation_read[!Proliferation_read %in% rownames(gene.tpm)]))
    match_genes <- na.omit(match_genes)
  }
  
  Proliferation.score <- vector("numeric", length = ncol(gene.tpm)) ; names(Proliferation.score) <- colnames(gene.tpm)
  
  # Log2 transformation:
  log2.gene.tpm <- log2(gene.tpm + 1)
  
  # Keep matrix for cell cycle signature genes
  sub_log2.gene.tpm <- log2.gene.tpm[match_genes,]
  
  # Calculate rank position for each gene across samples
  ranks_sub_log2.gene.tpm <- apply(sub_log2.gene.tpm, 1, rank)
  
  # Get normalized rank by divided 
  ranks_sub_log2.gene.tpm.norm <- (ranks_sub_log2.gene.tpm - 1)/(nrow(ranks_sub_log2.gene.tpm) - 1)
  
  # Average of the expression value of all the genes within-sample
  average_ranks_sub_log2.gene.tpm.norm<- apply(ranks_sub_log2.gene.tpm.norm, 1, mean)
  
  # For the pan-cancer analysis, the immune signature score was normalized across the samples of the same tumor type. 
  Proliferation.score <- data.frame(Proliferation = average_ranks_sub_log2.gene.tpm.norm, check.names = FALSE)
  
  cat("Proliferation score computed","\n")
  return(Proliferation.score)
}

# ***************
