#' compute.LR.pairs
#'
#' This function compute LR and Cytokine pairs from tpm as done in Maisa's externship
#'
#' @importFrom 
#'
#' @export
#'
#' @param RNA.tpm numeric matrix of tpm values with rows=genes and columns=samples
#' @param remove.genes.ICB_proxies boolean variable to reomove all those genes involved in the computation of ICB proxy's of response
#' @param compute.cytokines.pairs boolean variable to compute cytokine pairs as well
#'
#' @return Ligand-Receptors pairs weights matrix in log2 tpm + 1 with rows=samples and columns=LR.pairs
#'
#--------------------------------------------------------------------
# Compute Ligand-Receptor and Cytokine pairs weights from transcriptomics data.
#--------------------------------------------------------------------


compute.LR.pairs <- function(RNA.tpm, remove.genes.ICB_proxies=TRUE, compute.cytokines.pairs=TRUE,....){

  
  # ****************
  # data
  load("../data/all_genes_ICB_proxies.RData")
  load("../data/Ligand_Receptors_Rdata/LR_pairs_network_only_TME.RData") # Load the ligand receptor network (only TME)
  load("../data/Ligand_Receptors_Rdata/CYTOKINE_pairs_subnetwork_only_TME.RData") # Load cytokine subnetwork (only TME)
  
  gene_expr <- log2(RNA.tpm + 1)
  genes <- rownames(gene_expr)
  
  # Genes to remove according to all ICB proxy's 
  if (remove.genes.ICB_proxies) {
    cat("Removing signatures genes for proxy's of ICB response:  \n")
    idy <- na.exclude(match(all_genes_to_remove, rownames(gene_expr)))
    gene_expr <- gene_expr[-idy,]
  }
  
  gene_expr <- as.data.frame(gene_expr)
  
  # Compute L-R pairs 
  LR.pairs.computed <- do.call(rbind, lapply(1:length(LR.pairs_network), function(x){
    ligand <- sapply(strsplit(LR.pairs_network[x], split = "_", fixed = TRUE), head, 1)
    receptor <- sapply(strsplit(LR.pairs_network[x], split = "_", fixed = TRUE), tail, 1)
    by_patient <-  apply(gene_expr[c(ligand,receptor), ], 2, min)
  }))
  colnames(LR.pairs.computed) <- colnames(gene_expr)
  LR.pairs.computed <- t(LR.pairs.computed)
  colnames(LR.pairs.computed) <- LR.pairs_network

  # Compute cytokine pairs 
  if (compute.cytokines.pairs) {
    idy <- na.exclude(match(CYTOKINE.pairs_subnetwork, colnames(LR.pairs.computed)))
    CYTOKINE.pairs.computed <- LR.pairs.computed[,idy]
  }
  
  # Remove missing pairs in LR pairs data frame
  missing.LRpairs <- colnames(LR.pairs.computed)[as.vector(apply(LR.pairs.computed, 2, anyNA))]
  LR.pairs.computed <- LR.pairs.computed[,-match(missing.LRpairs, colnames(LR.pairs.computed))]
  # Remove missing pairs in CYTOKINE pairs data frame
  missing.CYTOKINEpairs <- colnames(CYTOKINE.pairs.computed)[as.vector(apply(CYTOKINE.pairs.computed, 2, anyNA))]
  CYTOKINE.pairs.computed <- CYTOKINE.pairs.computed[,-match(missing.CYTOKINEpairs, colnames(CYTOKINE.pairs.computed))]
  
  # Output list:  
  LR.data <- list(LRpairs = as.data.frame(LR.pairs.computed), CYTOKINEpairs = as.data.frame(CYTOKINE.pairs.computed))
  
  message("L-R pairs computed \n")
  message("Cytokine pairs computed \n")
  
  return(LR.data)
}
