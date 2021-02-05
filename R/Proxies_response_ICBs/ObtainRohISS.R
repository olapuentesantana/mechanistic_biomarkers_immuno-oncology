###########################################################################################
# Script to create ObtainRohISS function:
# Computing the geometric mean of literature genes provided by the authors:
# Geometric-mean of TPM-normalized counts of immune score genes
###########################################################################################

ObtainRohISS <- function(gene.tpm){

  # ************
  # Literature info: cytolytic markers, HLA molecules, IFN-γ pathway genes, chemokines and adhesion molecules
  RohISS_read <- c("GZMA", "GZMB", "PRF1", "GNLY", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
                    "HLA-G", "HLA-H", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1",
                    "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1",
                    "IFNG", "IFNGR1", "IFNGR2", "IRF1", "STAT1", "PSMB9", "CCR5", "CCL3", "CCL4",
                    "CCL5", "CXCL9", "CXCL10", "CXCL11", "ICAM1", "ICAM2", "ICAM3", "ICAM4", "ICAM5", "VCAM1")

  # ***************
  # Gene expression data just for roh immune signature score
  match_genes <- match(RohISS_read, rownames(gene.tpm))

  if (anyNA(match_genes)){
    warning(paste0("differenty named or missing signature genes : \n", RohISS_read[!RohISS_read %in% rownames(gene.tpm)]))
    match_genes <- na.omit(match_genes)
  }

  RohIS.score <- vector("numeric", length = ncol(gene.tpm)) ; names(RohIS.score) <- colnames(gene.tpm)

  # Subset gene expression data
  sub_gene.tpm <- gene.tpm[match_genes,]

  # Pseudocount of 0.01 for all genes
  sub_gene.tpm <- sub_gene.tpm + 0.01
  # Pseudocount of 1 for genes with 0 expr
  if(any(sub_gene.tpm == 0)) sub_gene.tpm[sub_gene.tpm == 0] <- sub_gene.tpm[sub_gene.tpm == 0] +1
  # Calculated as geometric mean (so-called log-average) [TPM, 0.01 offset]
  geom_mean <- apply(sub_gene.tpm, 2, function(X) exp(mean(log(X))))

  RohIS.score <- data.frame(RohIS = geom_mean, check.names = FALSE)

  cat("Roh IS score computed","\n")
  return(RohIS.score)
}

# ***************
# Checked with Van Allen tpm data



