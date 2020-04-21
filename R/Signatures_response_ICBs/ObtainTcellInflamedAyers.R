###########################################################################################
# Script to create ObtainTcellInflamedAyers function:
# The GEP score was computed by taking a weighted sum of the housekeeping normalized values of 
# the 18 genes on the GEP18 signature. Cristecu et al. 2018 provided the weights and also explained 
# how it is done. GEP scores listed in Table S2A were computed by first normalizing the raw counts
# by subtracting the average of the log10 counts of the house-keeping genes from the log10 count of 
# each of the predictor genes, and then a weighted sum of the normalized predictor gene values was 
# calculated using the weights for each of the 18 genes:
###########################################################################################

ObtainTcellInflamedAyers <- function(gene_expr){
  
  # ************
  # Literature info: genes related to antigen presentation, chemokine expression, cytolytic activity, and adaptive immune resistance
  # 18-gene T cell-inflamed signature (Ayers)
  T_cell_inflamed_Ayers_read <- c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1", 
                                  "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")
  
  Housekeeping.genes <- c("STK11IP", "ZBTB34", "TBC1D10B", "OAZ1", "POLR2A", "G6PD", "ABCF1", "C14orf102", "UBB", "TBP", "SDHA")
  
  weights <- data.frame(CCL5=0.008346, CD27=0.072293, CD274=0.042853, CD276=-0.0239, CD8A=0.031021 ,CMKLR1=0.151253, CXCL9=0.074135,
                        CXCR6=0.004313, `HLA-DQA1`=0.020091, `HLA-DRB1`=0.058806, `HLA-E`=0.07175, IDO1=0.060679, LAG3=0.123895, NKG7=0.075524, PDCD1LG2=0.003734,
                        PSMB10=0.032999, STAT1=0.250229, TIGIT=0.084767, check.names = FALSE)
  
  TcellInflamed.score <- vector("numeric", length = ncol(gene_expr)) ; names(TcellInflamed.score) <- colnames(gene_expr)
  
  # ***************
  # Gene expression data just for T cell inflamed signature (Ayers)
  match_genes.housekeeping <- na.omit(match(Housekeeping.genes, rownames(gene_expr)))
  match_genes.predictors <- na.omit(match(T_cell_inflamed_Ayers_read, rownames(gene_expr)))
  
  if (anyNA(match_genes.housekeeping) & anyNA(match_genes.predictors)){
    warning("Be careful, some genes are missing")
    cat("MISSING GENES", "\n") 
    cat(!which(Housekeeping.genes %in% rownames(gene_expr)))
    cat(!which(T_cell_inflamed_Ayers_read %in% rownames(gene_expr)))
  }
  
  # Get log10 counts
  log10.gene_expr <- log10(gene_expr + 1)
  
  log10.gene_expr.housekeeping <- log10.gene_expr[match_genes.housekeeping, ]
  log10.gene_expr.predictors <- log10.gene_expr[match_genes.predictors, ]
  
  # Housekeeping normalization
  average.log10.housekeeping.genes <- apply(log10.gene_expr.housekeeping,2, mean)
  norm.predictor.genes <- sweep(log10.gene_expr.predictors, 2, average.log10.housekeeping.genes, FUN = "-")
  
  # Weighted sum of the normalized predictor gene values
  tidy <- match(rownames(norm.predictor.genes), colnames(as.vector(weights)))
  TcellInflamed.score <- data.frame(T_cell_inflamed = t(norm.predictor.genes[tidy,]) %*% t(as.vector(weights)), check.names = FALSE)

  
  cat("T cell inflamed (Ayers) score computed","\n")
  return(TcellInflamed.score)
}

# ***************


