#' compute.pathways.scores
#'
#' This function computes pathway activity scores from raw counts RNAseq data.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions getVarianceStabilizedData
#' @importFrom progeny progeny
#'
#' @export
#'
#' @param RNA.raw_counts numeric matrix with genes as rows and samples as columns
#'
#' @return List: scores = Pathway activity matrix with samples as rows and pathways as columns; 
#' transcripts_kept = transcripts available for computation; 
#' transcripts_left = transcripts missing for computation.
#' 
#'
#--------------------------------------------------------------------
# Compute pathways activity from transcriptomics data.
#--------------------------------------------------------------------

# This function pre-process transcriptomics data as required by PROGENy. It matches transcript names with PROGENy perturbation genes
# and provides information regarding how many transcripts are lost/kept in the computation. 
# Pathways scores are computed using PROGENy package.

compute.pathways.scores <- function(RNA.raw_counts){

  # ****************
  # packages
  if(!("BiocManager" %in% installed.packages()[,"Package"])) install.packages("BiocManager", quiet = TRUE)
  list.of.packages <- c("progeny", "DESeq2")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) BiocManager::install(new.packages, ask = FALSE)

  suppressMessages(library(DESeq2))
  suppressMessages(library(progeny))

  # ****************
  # data
  load("pathway_responsive_genes.Rdata")
  
  raw_counts <- RNA.raw_counts
  genes <- rownames(raw_counts)

  # Integers are required for "DESeq2"
  if (any(sapply(raw_counts,is.integer) == FALSE)) {
    raw_counts <- sapply(raw_counts,as.integer)
    rownames(raw_counts) <- genes
    }
  # Rows with non-valid HGNC symbols were removed.
  if (any(grep("?",genes, fixed = T))) {raw_counts <- raw_counts[-grep("?",rownames(raw_counts), fixed = T),]}
  if (any(grep("\\|",genes, fixed = T))) {rownames(raw_counts) <- sapply(strsplit(rownames(raw_counts),"\\|"),function(X) return(X[1]))}

  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(rownames(raw_counts)) != 0){
    idx <- which(duplicated(rownames(raw_counts)) == TRUE)
    dup_genes <- rownames(raw_counts)[idx]
    for (ii in dup_genes){
      raw_counts[which(rownames(raw_counts) %in% ii)[1],] <- colMeans(raw_counts[which(rownames(raw_counts) %in% ii),])
      raw_counts <- raw_counts[-which( rownames(raw_counts) %in% ii)[2],]
    }
  }

  # Variance stabilizing transformation (DESeq2 package)
  # Integer count matrix, a data frame with the sample info,  design =~1 to consider all samples as part of the same group.

  # Column data:
  colData <- data.frame(id = colnames(raw_counts))

  # Construction a DESeqDataSet: (Forced all to be data.frames($ operator))
  dset <- DESeqDataSetFromMatrix(countData = raw_counts,
                                 colData = colData,
                                 design = ~ 1)
  # Variance stabilization transformation
  dset <- estimateSizeFactors(dset)
  dset <- estimateDispersions(dset)
  gene_expr <- getVarianceStabilizedData(dset)

  # Matching transcript names with pathway responsive genes
  genes_kept <- intersect(rownames(gene_expr), pathway_responsive_genes)
  genes_left <- setdiff(pathway_responsive_genes, rownames(gene_expr))

  cat("To compute pathway scores:  \n")
  cat("Transcripts available: ",length(genes_kept),"\n")
  cat("Transcripts missing: ", length(genes_left),"\n")

  # Pathways activity (Progeny package)
  Pathway_scores <- progeny(gene_expr, scale = FALSE)

  # Output list:
  Pathways <- list(scores = Pathway_scores, transcripts_kept = length(genes_kept), transcripts_left = length(genes_left))

  return(Pathways)
}
