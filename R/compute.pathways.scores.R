#' compute.pathways.scores
#'
#' This function computes pathway activity scores from raw counts RNAseq data.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions getVarianceStabilizedData
#' @importFrom progeny progeny
#'
#' @export
#'
#' @param RNA.raw_counts numeric matrix with data
#'
#' @return Pathway activity matrix
#'
#--------------------------------------------------------------------
# Compute pathways activity from transcriptomics data.
#--------------------------------------------------------------------

# This function pre-process transcriptomics data as required by PROGENy. It matches transcript names with PROGENy perturbation genes
# and provides information regarding how many transcripts are lost/kept. Pathways scores are computed using PROGENy package.

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
  load("data/Validation/pathway_responsive_genes.Rdata")
  
  # Genes to remove according to IS and CYT
  load("data/list_genes_IS_CYT.Rdata")
  ISCYT_read <- unique(list_genes_IS_CAT)
  # Genes from IPS
  IPSG_read <- read.table("data/raw_data_tcga/IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
  # Genes from IMPRES
  IMPRES.checkpoint.pairs <- data.frame(Gene_1 = c("PDCD1","CD27","CTLA4","CD40","CD86", "CD28", "CD80", 
                                                   "CD274","CD86","CD40","CD86","CD40","CD28","CD40","TNFRSF14"),
                                        Gene_2 = c("TNFSF4","PDCD1","TNFSF4","CD28","TNFSF4", "CD86", "TNFSF9", 
                                                   "C10orf54","HAVCR2","PDCD1","CD200","CD80","CD276","CD274","CD86"))
  IMPRES_read <- unique(as.vector(as.matrix(IMPRES.checkpoint.pairs))) # 15 genes
  # Unify all genes from signatures
  list_genes_to_remove <- unique(c(ISCYT_read,as.character(IPSG_read$GENE),IMPRES_read))
  
  raw_counts <- RNA.raw_counts
  genes <- rownames(raw_counts)

  # Rows with non-valid HGNC symbols were removed.
  if (any(grep("\\?",genes))) {
    raw_counts <- raw_counts[-grep("?",rownames(raw_counts), fixed = T),]
    genes <- rownames(raw_counts)
  }
  if (any(grep("\\|",genes))) {
    genes <- sapply(strsplit(rownames(raw_counts),"\\|"),function(X) return(X[1]))
  }

  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(genes) != 0){
    idx <- which(duplicated(genes) == TRUE)
    dup_genes <- genes[idx]
    for (ii in dup_genes){
      raw_counts[which(genes %in% ii)[1],] <- colMeans(raw_counts[which(genes %in% ii),])
      raw_counts <- raw_counts[-which(genes %in% ii)[2],]
      genes <- genes[-which(genes %in% ii)[2]]
    }
    rownames(raw_counts) <- genes 
  }
  
  # Remove list of genes used to build IS and CYT
  idy <- na.exclude(match(list_genes_to_remove, rownames(raw_counts)))
  raw_counts <- raw_counts[-idy,]
  
  # Integers are required for "DESeq2"
  if (any(sapply(raw_counts,is.integer) == FALSE)) {
    raw_counts.integer <- sapply(raw_counts,as.integer)
    rownames(raw_counts.integer) <- rownames(raw_counts)
    
  }

  # Variance stabilizing transformation (DESeq2 package)
  # Integer count matrix, a data frame with the sample info,  design =~1 to consider all samples as part of the same group.

  # Column data:
  colData <- data.frame(id = colnames(raw_counts.integer))

  # Construction a DESeqDataSet: (Forced all to be data.frames($ operator))
  dset <- DESeqDataSetFromMatrix(countData = raw_counts.integer,
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
