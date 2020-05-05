#' compute.pathways.scores
#'
#' This function computes pathway activity scores from raw counts RNAseq data.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions getVarianceStabilizedData
#' @importFrom progeny progeny
#'
#' @export
#'
#' @param RNA.raw_counts numeric matrix of read counts with rows=genes and columns=samples
#' @param remove.genes.ICB_proxies character vector of all those genes involved in the computation of ICB proxy's of response
#'
#' @return Pathway activity matrix
#'
#--------------------------------------------------------------------
# Compute pathways activity from transcriptomics data.
#--------------------------------------------------------------------

# This function pre-process transcriptomics data as required by PROGENy. It matches transcript names with PROGENy perturbation genes
# and provides information regarding how many transcripts are lost/kept. Pathways scores are computed using PROGENy package.

compute.pathways.scores <- function(RNA.raw_counts, remove.genes.ICB_proxies=TRUE,....){

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
  # Genes to remove according to ICB proxies genes
  load("../data/all_genes_ICB_proxies.RData")
  # The human  matrix consists of 14 pathways that accounts for the importance of each gene on each pathway upon perturbation. 
  # load("../data/model_human_full.rda")
  # pathway_responsive_genes <- as.character(unique(model_human_full$gene))
  # MAPK_pathway <- subset(model_human_full, pathway == "MAPK")
  # MAPK_pathway <- MAPK_pathway[which(MAPK_pathway$p.value %in% sort(MAPK_pathway$p.value, decreasing = FALSE)[1:100]),]
  
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
  
  # Remove list of genes used to build proxy's of ICB response
  if (remove.genes.ICB_proxies) {
    cat("Removing signatures genes for proxy's of ICB response:  \n")
    idy <- na.exclude(match(all_genes_to_remove, rownames(raw_counts)))
    raw_counts <- raw_counts[-idy,]
  }
  
  # Integers are required for "DESeq2"
  if (any(sapply(raw_counts,is.integer) == FALSE)) {
    raw_counts.integer <- raw_counts
    mode(raw_counts.integer) <- "integer"
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
  
  # HGNC symbols are required
  try(if (any(grep("ENSG00000", rownames(gene_expr)))) stop("hgnc gene symbols are required", call. = FALSE))

  # # Matching transcript names with pathway responsive genes
  # genes_kept <- intersect(rownames(gene_expr), MAPK_pathway$gene)
  # genes_left <- setdiff(pathway_responsive_genes, rownames(gene_expr))
  # 
  cat("To compute pathway scores:  \n")
  # cat("Transcripts available: ",length(genes_kept),"\n")
  # cat("Transcripts missing: ", length(genes_left),"\n")

  # Pathways activity (Progeny package)
  Pathway_scores <- progeny(gene_expr, scale = FALSE, organism = "Human", verbose = TRUE)

  # Output list:
  # Pathways <- list(scores = Pathway_scores, transcripts_kept = length(genes_kept), transcripts_left = length(genes_left))

  return(Pathway_scores)
}
