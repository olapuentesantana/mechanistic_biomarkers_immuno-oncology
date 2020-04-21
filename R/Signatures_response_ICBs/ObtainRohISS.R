###########################################################################################
# Script to create ObtainRohISS function:
# Computing the geometric mean of literature genes provided by the authors:
# Geometric-mean of TPM-normalized coounts of immune score genes
###########################################################################################

ObtainRohISS <- function(gene_expr){
  
  # ************
  # Literature info: cytolytic markers, HLA molecules, IFN-γ pathway genes, chemokines and adhesion molecules
  RohISS_read <- c("GZMA", "GZMB", "PRF1", "GNLY", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", 
                    "HLA-G", "HLA-H", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", 
                    "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", 
                    "IFNG", "IFNGR1", "IFNGR2", "IRF1", "STAT1", "PSMB9", "CCR5", "CCL3", "CCL4", 
                    "CCL5", "CXCL9", "CXCL10", "CXCL11", "ICAM1", "ICAM2", "ICAM3", "ICAM4", "ICAM5", "VCAM1")

  # ***************
  # Gene expression data just for roh immune signature score
  match_genes <- na.omit(match(RohISS_read, rownames(gene_expr)))
  
  if (anyNA(match_genes)){
    warning("Be careful, some genes are missing")
    cat("MISSING GENES", "\n") 
    cat(!which(RohISS_read %in% rownames(gene_expr)))
  }

  RohIS.score <- vector("numeric", length = ncol(gene_expr)) ; names(RohIS.score) <- colnames(gene_expr)
  
  RohIS.score <- sapply(colnames(gene_expr), function(id){
    tmp <- gene_expr[match_genes, id]
    # Pseudocount of 0.01 for all genes
    tmp <- tmp + 0.01
    # Pseudocount of 1 for genes with 0 expr
    if(any(tmp == 0)) tmp[which(tmp == 0)] <- tmp[which(tmp == 0)] +1

    aux <- 0
    # geometric mean
    for(X in 1:length(tmp)-1){
      if (X == 1){
        aux <- tmp[X+1] * tmp[X]
      }else{
        aux <- aux * tmp[X+1]
      }
    }
    geom_mean <- aux^(1/length(tmp))
    return(geom_mean)
  })
  RohIS.score <- data.frame(RohIS = RohIS.score, check.names = FALSE)
    
  cat("Roh IS score computed","\n")
  return(RohIS.score)
}

# ***************
# Checked with Van Allen tpm data



