#' compute.LR.pairs
#'
#' This function compute LR pairs from tpm as done in Maisa's externship
#'
#' @importFrom 
#'
#' @export
#'
#' @param RNA.tpm numeric matrix of tpm values with rows=genes and columns=samples
#' @param remove.genes.ICB_proxies character vector of all those genes involved in the computation of ICB proxy's of response
#'
#' @return Ligand-Receptors pairs weights matrix in log2 tpm + 1 with rows=samples and columns=LR.pairs
#'
#--------------------------------------------------------------------
# Compute Ligand-Receptor pairs weights from transcriptomics data.
#--------------------------------------------------------------------

# This function pre-process transcriptomics data as required by DoRothEA. It matches transcript names with DoRothEA genes (regulons for the
# transcription factors) and provides information regarding how many transcripts are lost/kept in the computation. 
# TF activities are computed using DoRothEA package.

compute.LR.pairs <- function(RNA.tpm, remove.genes.ICB_proxies=TRUE,....){

  # ****************
  # packages
  

  # ****************
  # scripts

  
  # ****************
  # data
  load("../data/all_genes_ICB_proxies.RData")
  load("../data/LR.pairs.Ramilowski.RData") # Load the ligand receptor network from Ramilowski (only literature supported interactions)
 
  gene_expr <- log2(RNA.tpm + 1)
  genes <- rownames(gene_expr)
  
  # Rows with non-valid HGNC symbols were removed.
  if (any(grep("\\?",genes))) {
    gene_expr <- gene_expr[-grep("?",rownames(gene_expr), fixed = T),]
    genes <- rownames(gene_expr)
  }
  if (any(grep("\\|",genes))) {
    genes <- sapply(strsplit(rownames(gene_expr),"\\|"),function(X) return(X[1]))
  }
  
  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(genes) != 0){
    idx <- which(duplicated(genes) == TRUE)
    dup_genes <- genes[idx]
    for (ii in dup_genes){
      gene_expr[which(genes %in% ii)[1],] <- colMeans(gene_expr[which(genes %in% ii),])
      gene_expr <- gene_expr[-which(genes %in% ii)[2],]
      genes <- genes[-which(genes %in% ii)[2]]
    }
    rownames(gene_expr) <- genes 
  }
  
  # Genes to remove according to all ICB proxy's 
  if (remove.genes.ICB_proxies) {
    cat("Removing signatures genes for proxy's of ICB response:  \n")
    idy <- na.exclude(match(all_genes_to_remove, rownames(gene_expr)))
    gene_expr <- gene_expr[-idy,]
  }
  
  gene_expr <- as.data.frame(gene_expr)
  # compute the LR.pairs myself
  LR.pairs.computed <- do.call(rbind, lapply(1:nrow(LR.pairs.Ramilowski), function(x){
    ligand <- LR.pairs.Ramilowski[x, "Ligand.ApprovedSymbol"]
    receptor <- LR.pairs.Ramilowski[x, "Receptor.ApprovedSymbol"]
    bypatient <- do.call(cbind, lapply(colnames(gene_expr), function(y){
      min(gene_expr[ligand, y], gene_expr[receptor,y])
    }))
  }))
  colnames(LR.pairs.computed) <- colnames(gene_expr)
  LR.pairs.computed <- t(LR.pairs.computed)
  colnames(LR.pairs.computed) <- LR.pairs.Ramilowski$Pair.Name

  return(LR.pairs.computed)
}
