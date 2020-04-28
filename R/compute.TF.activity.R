#' compute.TF.activity
#'
#' This function transcription factor activity from raw counts RNAseq data.
#'
#' @importFrom viper viper
#'
#' @export
#'
#' @param RNA.raw_counts numeric matrix with data
#'
#' @return Transcription Factor activity matrix
#'
#--------------------------------------------------------------------
# Compute TFs activity from transcriptomics data.
#--------------------------------------------------------------------

# This function pre-process transcriptomics data as required by DoRothEA It matches transcript names with DoRothEA genes (regulons for the
# transcription factors) and provides information regarding how many transcripts are lost/kept. TF activities are computed using DoRothEA package.

compute.TF.activity <- function(RNA.tpm){

  # ****************
  # packages
  if(!("viper" %in% installed.packages()[,"Package"])) BiocManager::install("viper", ask = FALSE)
  suppressMessages(require(viper))

  # ****************
  # scripts
  source("R/scaling_function.R")
  
  # ****************
  # data
  load("data/Validation/TFs_viperRegulon.rdata") # Load TF regulon genesets in VIPER format
  load("data/all_genes_ICB_proxies.RData")
  
  # Genes to remove according to all ICB proxy's 
  idx <- na.exclude(match(all_genes_to_remove, rownames(RNA.tpm)))
  RNA.tpm <- RNA.tpm[-idx,]

  # Expression matrix (rows=genes; columns=samples) scaled and recentered. Note that genes need to be in a comparable scale (e.g. z-transformed). 
  # Z score expression matrix (log2[tpm+1]):
  gene_expr <- standarization(log2(t(RNA.tpm)+1))
  
  # redefine gene names to match transcripts for viper
  E <- t(gene_expr)
  newNames <- sapply(rownames(E), function(x){
    # strsplit(x, "\\.")[[1]][1]
    zz_tmp <- strsplit(x, "\\.")[[1]]
    paste0(zz_tmp[1:(length(zz_tmp)-1)], collapse = "-")
  })
  rownames(E) <- newNames
  
  all_regulated_transcripts <- do.call(c, lapply(viper_regulon, function(x){
    names(x$tfmode)
  }))
  all_regulated_transcripts <- unique(all_regulated_transcripts)
  all_TF <- names(viper_regulon)
  
  # check what is the percentage of regulated transcripts and TF that we have in our data
  cat("\nTo compute TF activities: \n")
  cat(" percentage of regulated transcripts = ", sum(all_regulated_transcripts %in%  newNames)*100/length(all_regulated_transcripts), "\n")
  cat(" percentage of TF = ", sum(all_TF %in% newNames)*100/length(all_TF), "\n")

  # TF activity (viper package)
  TF_activities <- viper(eset = E, regulon = viper_regulon, nes = T, 
                         method = 'none', minsize = 4, eset.filter = F)
  
  # Samples as rows, TFs as columns
  TF_activities <- t(TF_activities)

  # Matching transcript names with TFs regulons
  regulons <- unique(do.call(c,lapply(names(viper_regulon), function(X){return(names(viper_regulon[[X]]$tfmode))})))
  genes_kept <- intersect(rownames(gene_expr), regulons)
  genes_left <- setdiff(regulons, rownames(gene_expr))

  # Output list:
  TFs <- list(scores = TF_activities, transcripts_kept = length(genes_kept), transcripts_left = length(genes_left))

  return(TFs)
}
