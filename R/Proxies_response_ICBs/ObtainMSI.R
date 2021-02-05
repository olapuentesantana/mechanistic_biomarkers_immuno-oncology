###########################################################################################
# Script to create ObtainMSI function:
# Computing of relation tab: comparisons of the F pairs.
# Input: TPM and log2-transformed
###########################################################################################

ObtainMSI <- function(gene.tpm){
  
  
  # ************
  # Literature info --> * (CCRN4L in tcga, NOCT approved symbol)
  MSI.pairs <- data.frame(Gene_1 = c("HNRNPL","MTA2","CALR","RASL11A","LYG1", "STRN3", "HPSE", 
                                                   "PRPF39","NOCT","AMFR"), 
                                        Gene_2 = c("CDC16","VGF","SEC22B","CAB39L","DHRS12", "TMEM192", "BCAS3", 
                                                   "ATF6","GRM8","DUSP18"))
  
  MSI_read <- unique(as.vector(as.matrix(MSI.pairs))) # 20 genes
  
  # Some genes might have other name: case for "CCRN4L", it's called "NOCT", be carefull
  if (any(rownames(gene.tpm) %in% "CCRN4L")){
    cat("Gene name changed: NOCT is approved symbol, not CCRN4L","\n")
    rownames(gene.tpm)[rownames(gene.tpm) %in% "CCRN4L"] <- "NOCT"
  }
  
  # ***************
  # Gene expression data just for the MSI genes
  match_F_1 <- match(as.character(MSI.pairs[,1]), rownames(gene.tpm))
  match_F_2 <- match(as.character(MSI.pairs[,2]), rownames(gene.tpm))

  if (anyNA(c(match_F_1,match_F_2))){
    warning(paste0("differenty named or missing signature genes : \n", MSI_read[!MSI_read %in% rownames(gene.tpm)]))
    match_F_1 <- na.omit(match_F_1)
    match_F_2 <- na.omit(match_F_2)
  }
  
  F_pair_expr <- matrix(0, nrow(MSI.pairs),2) ; colnames(F_pair_expr) <- c("A","B")
  MSI.matrix <- matrix(0, nrow(MSI.pairs), ncol(gene.tpm)) ; colnames(MSI.matrix) <- colnames(gene.tpm)
  remove_pairs <- vector("list", length = ncol(gene.tpm)) ; names(remove_pairs) <- colnames(gene.tpm)
  MSI.score <- vector("numeric", length = ncol(gene.tpm)) ; names(MSI.score) <- colnames(gene.tpm)
  
  # log2(TPM+1)
  gene_expr <- as.data.frame(log2(gene.tpm + 1))
  
  for (ii in 1:ncol(gene.tpm)){
    F_pair_expr[,1] <- gene_expr[match_F_1,ii]
    F_pair_expr[,2] <- gene_expr[match_F_2,ii]
    
    
    if(anyNA(as.vector(F_pair_expr))) remove_pairs[[ii]] = which(is.na(rowSums(F_pair_expr) == TRUE))
    
    for (x in 1:nrow(MSI.pairs)){
      if(is.element(x,remove_pairs[[ii]]) == TRUE){
        MSI.matrix[x,ii] <- NA
      }else{
        if(F_pair_expr[x,1] > F_pair_expr[x,2]){
          MSI.matrix[x,ii] <- 1
        }else{
          MSI.matrix[x,ii] <- 0
        }
      }
    }
    if(anyNA(MSI.matrix[,ii])){
      MSI.score[ii] <- sum(na.omit(MSI.matrix[,ii]))
      MSI.score[ii] <- (MSI.score[ii] *  nrow(MSI.pairs)) / (nrow(MSI.pairs) - length(remove_pairs[[ii]]))
    }else{
      MSI.score[ii] <- sum(MSI.matrix[,ii])
    }
  }
  
  # MSI.score <- colSums(MSI.matrix, na.rm = T)
  cat("MSI score computed","\n")
  
  # gene_expr_msi <- gene_expr[List.MSI.genes,]
  # answer <- apply(gene_expr_msi, 1, function(X) which(X == 0))
  # genes_samples <- unlist(answer)
  # genes <- sapply(strsplit(names(genes_samples),".", fixed = TRUE),function(X) return(X[1]))
  # samples <- sapply(strsplit(names(genes_samples),".", fixed = TRUE),function(X) return(X[2]))
  # quantify <- sapply(unique(genes), function(X) length(which(genes == X)))
  # 
  # if (is.null(quantify)) {warning("Some samples do not express MSI gene pairs") ; print(quantify)}

  return(data.frame(MSI = MSI.score, check.names = FALSE))
}

# ***************