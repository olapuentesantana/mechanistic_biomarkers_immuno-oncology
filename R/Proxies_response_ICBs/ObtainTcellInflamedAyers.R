###########################################################################################
# Script to create ObtainTcellInflamedAyers function:
# The GEP score was computed by taking a weighted sum of the housekeeping normalized values of 
# the 18 genes on the GEP18 signature. Cristecu et al. 2018 provided the weights and also explained 
# how it is done. GEP scores listed in Table S2A were computed by first normalizing the raw counts
# by subtracting the average of the log10 counts of the house-keeping genes from the log10 count of 
# each of the predictor genes, and then a weighted sum of the normalized predictor gene values was 
# calculated using the weights for each of the 18 genes:
###########################################################################################

ObtainTcellInflamedAyers <- function(gene.tpm){
  
  # ************
  # Literature info: genes related to antigen presentation, chemokine expression, cytolytic activity, and adaptive immune resistance
  # 18-gene T cell-inflamed signature (Ayers)
  T_cell_inflamed_Ayers_read <- c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1", 
                                  "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")
  
  Housekeeping.genes <- c("STK11IP", "ZBTB34", "TBC1D10B", "OAZ1", "POLR2A", "G6PD", "ABCF1", "NRDE2", "UBB", "TBP", "SDHA") # C14orf102 = NRDE2
  
  # * Pending
  # EQUIVALENT : "NRDE2" = "C14orf102"
  # Some genes might have other name: case for "C14orf102", it's called "NRDE2", be carefull
  if (any(rownames(gene.tpm) %in% "C14orf102")){
    cat("Gene name changed: NRDE2 is approved symbol, not C14orf102","\n")
    rownames(gene.tpm)[rownames(gene.tpm) %in% "C14orf102"] <- "NRDE2"
  }
  
  weights <- data.frame(CCL5=0.008346, CD27=0.072293, CD274=0.042853, CD276=-0.0239, CD8A=0.031021 ,CMKLR1=0.151253, CXCL9=0.074135,
                        CXCR6=0.004313, `HLA-DQA1`=0.020091, `HLA-DRB1`=0.058806, `HLA-E`=0.07175, IDO1=0.060679, LAG3=0.123895, NKG7=0.075524, PDCD1LG2=0.003734,
                        PSMB10=0.032999, STAT1=0.250229, TIGIT=0.084767, check.names = FALSE)
  
  TcellInflamed.score <- vector("numeric", length = ncol(gene.tpm)) ; names(TcellInflamed.score) <- colnames(gene.tpm)
  
  # ***************
  # Gene expression data just for T cell inflamed signature (Ayers)
  match_genes.housekeeping <- match(Housekeeping.genes, rownames(gene.tpm))
  match_genes.predictors <- match(T_cell_inflamed_Ayers_read, rownames(gene.tpm))

  if (anyNA(c(match_genes.housekeeping, match_genes.predictors))){
    tmp <- c(T_cell_inflamed_Ayers_read, Housekeeping.genes)
    warning(paste0("differenty named or missing signature genes : \n", tmp[!tmp %in% rownames(gene.tpm)]))
    match_genes.housekeeping <- na.omit(match_genes.housekeeping)
    match_genes.predictors <- na.omit(match_genes.predictors)
  }
  
  # Log2 transformation:
  log2.gene.tpm <- log2(gene.tpm + 1)
  
  log2.gene.tpm.housekeeping <- log2.gene.tpm[match_genes.housekeeping, ]
  log2.gene.tpm.predictors <- log2.gene.tpm[match_genes.predictors, ]
  
  # Housekeeping normalization
  average.log2.gene.tpm.housekeeping <- apply(log2.gene.tpm.housekeeping,2, mean)
  log2.gene.tpm.predictors.norm <- sweep(log2.gene.tpm.predictors, 2, average.log2.gene.tpm.housekeeping, FUN = "-")
  
  # Weighted sum of the normalized predictor gene values
  tidy <- match(rownames(log2.gene.tpm.predictors.norm), colnames(as.vector(weights)))
  TcellInflamed.score <- data.frame(T_cell_inflamed = t(log2.gene.tpm.predictors.norm[tidy,]) %*% t(as.vector(weights)), check.names = FALSE)

  cat("T cell inflamed (Ayers) score computed","\n")
  return(TcellInflamed.score)
}

# ***************


