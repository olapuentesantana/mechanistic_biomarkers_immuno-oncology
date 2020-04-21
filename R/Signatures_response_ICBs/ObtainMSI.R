###########################################################################################
# Script to create ObtainMSI function:
# Computing of relation tab: comparisons of the F pairs.
# Input: RSEM-normalized format and log2-transformed
###########################################################################################

ObtainMSI <- function(gene_expr, Cancer){
  
  
  # ************
  # Literature info
  MSI.pairs <- data.frame(Gene_1 = c("HNRNPL","MTA2","CALR","RASL11A","LYG1", "STRN3", "HPSE", 
                                                   "PRPF39","CCRN4L","AMFR"),
                                        Gene_2 = c("CDC16","VGF","SEC22B","CAB39L","DHRS12", "TMEM192", "BCAS3", 
                                                   "ATF6","GRM8","DUSP18"))
  
  List.MSI.genes <- unique(as.vector(as.matrix(MSI.pairs))) # 20 genes
  
  # * Pending
  # EQUIVALENT : "VISTA" = "C10orf54", "PDL-1" = "CD274", "TIM-3" = "HAVCR2",
  # "PD-1" = "PDCD1", "HVEM" = "TNFRSF14", "OX40L" = "TNFSF4", "CD137L" = "TNFSF9"
  
  # Some genes might have other name: case for "C10orf54", it's called "VSIR", be carefull
  
  # if (any(rownames(gene_expr) == "VSIR")){
  #   cat("Gene name changed: C10orf54 instead of VSIR","\n")
  #   rownames(gene_expr)[which(rownames(gene_expr) == "VSIR")] = "C10orf54"
  # }
  # * Pending
  
  # ***************
  # Gene expression data just for the MSI genes
  match_F_1 = match(as.character(MSI.pairs[,1]), rownames(gene_expr))
  match_F_2 = match(as.character(MSI.pairs[,2]), rownames(gene_expr))
  
  F_pair_expr <- matrix(0, nrow(MSI.pairs),2) ; colnames(F_pair_expr) <- c("A","B")
  MSI.matrix <- matrix(0, nrow(MSI.pairs), ncol(gene_expr)) ; colnames(MSI.matrix) <- colnames(gene_expr)
  remove_pairs <- vector("list", length = ncol(gene_expr)) ; names(remove_pairs) <- colnames(gene_expr)
  MSI.score <- vector("numeric", length = ncol(gene_expr)) ; names(MSI.score) <- colnames(gene_expr)
  
  for (ii in 1:ncol(gene_expr)){
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
    
    MSI.score[ii] <- sum(MSI.matrix[,ii])
    
    # if(anyNA(MSI.matrix[,ii])){
    #   MSI.score[ii] <- sum(na.omit(MSI.matrix[,ii]))
    #   MSI.score[ii] <- (MSI.score[ii] *  nrow(MSI.pairs)) / (nrow(MSI.pairs) - length(remove_pairs[[ii]]))
    # }else{
    #   MSI.score[ii] <- sum(MSI.matrix[,ii])
    # }
  }
  cat("MSI score computed","\n")
  
  gene_expr_msi <- gene_expr[List.MSI.genes,]
  answer <- apply(gene_expr_msi, 1, function(X) which(X == 0))
  genes_samples <- unlist(answer)
  genes <- sapply(strsplit(names(genes_samples),".", fixed = TRUE),function(X) return(X[1]))
  samples <- sapply(strsplit(names(genes_samples),".", fixed = TRUE),function(X) return(X[2]))
  quantify <- sapply(unique(genes), function(X) length(which(genes == X)))
  
  if (is.null(quantify)) {warning("Some samples do not express MSI gene pairs") ; print(quantify)}

  return(data.frame(MSI = MSI.score, check.names = FALSE))
}

# ***************