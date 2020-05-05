#' compute.TF.activity
#'
#' This function transcription factor activity from raw counts RNAseq data.
#'
#' @importFrom dorothea dorothea
#'
#' @export
#'
#' @param RNA.tpm numeric matrix of tpm values with rows=genes and columns=samples
#' @param remove.genes.ICB_proxies character vector of all those genes involved in the computation of ICB proxy's of response
#'
#' @return TF activity matrix: matrix of normalized enrichment scores with rows=samples and columns=TFs
#'
#--------------------------------------------------------------------
# Compute TFs activity from transcriptomics data.
#--------------------------------------------------------------------

# This function pre-process transcriptomics data as required by DoRothEA It matches transcript names with DoRothEA genes (regulons for the
# transcription factors) and provides information regarding how many transcripts are lost/kept. TF activities are computed using DoRothEA package.

compute.TF.activity <- function(RNA.tpm, remove.genes.ICB_proxies=TRUE,....){

  # ****************
  # packages
  if(!("BiocManager" %in% installed.packages()[,"Package"])) install.packages("BiocManager", quiet = TRUE)
  list.of.packages <- c("dorothea", "dplyr", "viper")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) BiocManager::install(new.packages, ask = FALSE)
  
  suppressMessages(library(dorothea))
  suppressMessages(library(dplyr))
  suppressMessages(library(viper))
  
  # ****************
  # data
  load("../data/all_genes_ICB_proxies.RData")
  # acessing (human) dorothea regulons
  data(dorothea_hs, package = "dorothea")
  
  tpm <- RNA.tpm
  genes <- rownames(tpm)
  
  # Rows with non-valid HGNC symbols were removed.
  if (any(grep("\\?",genes))) {
    tpm <- tpm[-grep("?",rownames(tpm), fixed = T),]
    genes <- rownames(tpm)
  }
  if (any(grep("\\|",genes))) {
    genes <- sapply(strsplit(rownames(tpm),"\\|"),function(X) return(X[1]))
  }
  
  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(genes) != 0){
    idx <- which(duplicated(genes) == TRUE)
    dup_genes <- genes[idx]
    for (ii in dup_genes){
      tpm[which(genes %in% ii)[1],] <- colMeans(tpm[which(genes %in% ii),])
      tpm <- tpm[-which(genes %in% ii)[2],]
      genes <- genes[-which(genes %in% ii)[2]]
    }
    rownames(tpm) <- genes 
  }
  
  # Genes to remove according to all ICB proxy's 
  if (remove.genes.ICB_proxies) {
    cat("Removing signatures genes for proxy's of ICB response:  \n")
    idy <- na.exclude(match(all_genes_to_remove, rownames(tpm)))
    tpm <- tpm[-idy,]
  }
  
  # Log transformed expression matrix (log2[tpm+1]):
  start <- Sys.time()
  gene_expr <- standarization(log2(t(tpm) + 1))

  # HGNC symbols are required
  try(if (any(grep("ENSG00000", rownames(gene_expr)))) stop("hgnc gene symbols are required", call. = FALSE))
  
  # redefine gene names to match transcripts for viper
  E <- t(gene_expr)
  newNames <- sapply(rownames(E), function(x){
    # strsplit(x, "\\.")[[1]][1]
    zz_tmp <- strsplit(x, "\\.")[[1]]
    paste0(zz_tmp[1:(length(zz_tmp)-1)], collapse = "-")
  })
  rownames(E) <- newNames
  
  regulons <- dplyr::filter(dorothea_hs, confidence %in% c("A", "B"))
  all_regulated_transcripts <- unique(regulons$target)
  # all_TF <- unique(regulons$tf)
  
  # check what is the percentage of regulated transcripts and TF that we have in our data
  cat("\nTo compute TF activities: \n")
  #cat(" percentage of regulated transcripts = ", sum(all_regulated_transcripts %in% newNames)*100/length(all_regulated_transcripts), "\n")
  #cat(" percentage of TF = ", sum(all_TF %in% newNames)*100/length(all_TF), "\n")

  # Expression matrix (rows=genes; columns=samples) scaled and recentered. 
  # Note that genes need to be in a comparable scale (e.g. z-transformed). 
  # Z score expression matrix (this is calculated within the package:
  
  # TF activity: run viper
  TF_activities <- dorothea::run_viper(input = E, regulons = regulons, 
                                       options = list(method = "none", minsize = 4, eset.filter = F, cores = 1, verbose=FALSE))
  end <- Sys.time()
  
  # Samples as rows, TFs as columns
  TF_activities <- t(TF_activities)

  # Matching transcript names with TFs regulons
  genes_kept <- intersect(rownames(E), all_regulated_transcripts)
  genes_left <- setdiff(all_regulated_transcripts, rownames(E))

  # Output list:
  TFs <- list(scores = TF_activities, transcripts_kept = length(genes_kept), transcripts_left = length(genes_left))

  return(TFs)
}
