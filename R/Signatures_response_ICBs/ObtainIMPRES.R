###########################################################################################
# Script to create ObtainIMPRES function:
# Computing of relation tab: comparisons of the F pairs.
###########################################################################################

ObtainIMPRES <- function(gene_expr){
  
  # ************
  # Literature info
  IMPRES.checkpoint.pairs <- data.frame(Gene_1 = c("PDCD1","CD27","CTLA4","CD40","CD86", "CD28", "CD80", 
                                                   "CD274","CD86","CD40","CD86","CD40","CD28","CD40","TNFRSF14"),
                                        Gene_2 = c("TNFSF4","PDCD1","TNFSF4","CD28","TNFSF4", "CD86", "TNFSF9", 
                                                   "C10orf54","HAVCR2","PDCD1","CD200","CD80","CD276","CD274","CD86"))
  
  List.checkpoint.genes <- unique(as.vector(as.matrix(IMPRES.checkpoint.pairs))) # 15 genes
  
  # EQUIVALENT : "VISTA" = "C10orf54", "PDL-1" = "CD274", "TIM-3" = "HAVCR2",
  # "PD-1" = "PDCD1", "HVEM" = "TNFRSF14", "OX40L" = "TNFSF4", "CD137L" = "TNFSF9"
  
  # Some genes might have other name: case for "C10orf54", it's called "VSIR", be carefull
  
  if (any(rownames(gene_expr) == "VSIR")){
    cat("Gene name changed: C10orf54 instead of VSIR","\n")
    rownames(gene_expr)[which(rownames(gene_expr) == "VSIR")] = "C10orf54"
  }
    
  # ***************
  # Gene expression data just for the checkpoint genes
  match_F_1 = match(as.character(IMPRES.checkpoint.pairs[,1]), rownames(gene_expr))
  match_F_2 = match(as.character(IMPRES.checkpoint.pairs[,2]), rownames(gene_expr))
  
  F_pair_expr <- matrix(0, nrow(IMPRES.checkpoint.pairs),2) ; colnames(F_pair_expr) <- c("A","B")
  IMPRES.matrix <- matrix(0, nrow(IMPRES.checkpoint.pairs), ncol(gene_expr)) ; colnames(IMPRES.matrix) <- colnames(gene_expr)
  remove_pairs <- vector("list", length = ncol(gene_expr)) ; names(remove_pairs) <- colnames(gene_expr)
  IMPRES.score <- vector("numeric", length = ncol(gene_expr)) ; names(IMPRES.score) <- colnames(gene_expr)

  for (ii in 1:ncol(gene_expr)){
      F_pair_expr[,1] <- gene_expr[match_F_1,ii]
      F_pair_expr[,2] <- gene_expr[match_F_2,ii]
      
      if(anyNA(as.vector(F_pair_expr))) remove_pairs[[ii]] = which(is.na(rowSums(F_pair_expr) == TRUE))
      
      for (x in 1:length(List.checkpoint.genes)){
        if(is.element(x,remove_pairs[[ii]]) == TRUE){
            IMPRES.matrix[x,ii] <- NA
        }else{
          if(F_pair_expr[x,1] > F_pair_expr[x,2]){
            IMPRES.matrix[x,ii] <- 1
          }else{
            IMPRES.matrix[x,ii] <- 0
          }
        }
      }
      
    if(anyNA(IMPRES.matrix[,ii])){
      IMPRES.score[ii] <- sum(na.omit(IMPRES.matrix[,ii]))
      IMPRES.score[ii] <- (IMPRES.score[ii] *  nrow(IMPRES.checkpoint.pairs)) / (nrow(IMPRES.checkpoint.pairs) - length(remove_pairs[[ii]]))
    }else{
      IMPRES.score[ii] <- sum(IMPRES.matrix[,ii])
    }
  }
  cat("IMPRES score computed","\n")
  
  gene_expr_impres <- gene_expr[List.checkpoint.genes,]
  answer <- apply(gene_expr_impres, 1, function(X) which(X == 0))
  genes_samples <- unlist(answer)
  genes <- sapply(strsplit(names(genes_samples),".", fixed = TRUE),function(X) return(X[1]))
  samples <- sapply(strsplit(names(genes_samples),".", fixed = TRUE),function(X) return(X[2]))
  quantify <- sapply(unique(genes), function(X) length(which(genes == X)))
  
  if (is.null(quantify)) {warning("Some samples do not express IMPRES gene pairs") ; print(quantify)}
  
  return(data.frame(IMPRES = IMPRES.score, check.names = FALSE))
}

# ***************


