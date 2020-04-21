###########################################################################################
# Script to create ObtainTIDE function:
# Compute TIDE using published model (tidepy)
# Input: log2(TPM + 1) + quantile normalization + mean centralization (author's suggestions)
###########################################################################################

ObtainTIDE <- function(gene_expr, Cancer){
  
  # ************
  # Pre-process 
  ## quantile normalization
  # Genes as rows
  # Sort each column of X (Genes as rows)
  sort.log2.mas1.TPM <- indices.sort <- matrix(0,nrow = nrow(gene_expr), ncol = ncol(gene_expr))
  
  for (X in 1:ncol(gene_expr)){
    tmp = sort(gene_expr[,X], decreasing = F, index.return = TRUE)
    sort.log2.mas1.TPM[,X] = tmp$x
    indices.sort[,X] = tmp$ix
  }
  # Take the means across rows of X sort
  # Assign this mean to each element in the row to get X' short 
  for (X in 1:nrow(gene_expr)){
    sort.log2.mas1.TPM[X,] <- rowMeans(sort.log2.mas1.TPM)[X]
  }
  # Get X normalized by rearranging each column of X' sort to have th esame ordering as original X
  for (X in 1:ncol(gene_expr)){
    tmp = sort(indices.sort[,X], decreasing = F, index.return = TRUE)
    sort.log2.mas1.TPM[,X] <- sort.log2.mas1.TPM[tmp$ix,X]
  }
  rownames(sort.log2.mas1.TPM) <- rownames(gene_expr)
  colnames(sort.log2.mas1.TPM) <- colnames(gene_expr)
  #boxplot.matrix(sort.TCGA.log2TPM1.pre)
  
  ## Mean-centralization: substract gene average across all conditions
  average.gene <- rowMeans(sort.log2.mas1.TPM)
  sort.log2.mas1.TPM.norm <- sweep(sort.log2.mas1.TPM,1,average.gene, FUN = "-")
  
  write.table(sort.log2.mas1.TPM.norm, file = paste0("./data/data_processed_TIDE/log2mas1TPM_", Cancer,".txt"), sep = "\t")
  
  # Compute TIDE
  if (Cancer == "SKCM") {
    system(paste0("/home/olapuent/anaconda3/bin/tidepy ", "./data/data_processed_TIDE/log2mas1TPM_", Cancer,".txt", " -o output_TIDE_",Cancer,".txt -c Melanoma"))
  }else{
    system(paste0("/home/olapuent/anaconda3/bin/tidepy ", "./data/data_processed_TIDE/log2mas1TPM_", Cancer,".txt", " -o output_TIDE_",Cancer,".txt -c Other"))
  }
  
  TIDE.table <- read.table(file = paste0("output_TIDE_",Cancer,".txt"), sep = "\t", header = T, row.names = 1)
  TIDE.score <- data.frame(TIDE = TIDE.table[,"TIDE", drop = FALSE], check.names = F)
  
  cat("TIDE score computed","\n")
  
  return(TIDE.score)
}

# ***************