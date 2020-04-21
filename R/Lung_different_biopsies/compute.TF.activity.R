#' compute.TF.activity
#'
#' This function transcription factor activity from raw counts RNAseq data.
#'
#' @importFrom viper viper
#'
#' @export
#'
#' @param RNA.tpm numeric matrix with genes as rows and samples as columns
#'
#' @return List: scores = TFs activity matrix samples as rows and TFs as columns; 
#' transcripts_kept = transcripts available for computation; 
#' transcripts_left = transcripts missing for computation.
#'
#'
#--------------------------------------------------------------------
# Compute TFs activity from transcriptomics data.
#--------------------------------------------------------------------

# This function pre-process transcriptomics data as required by DoRothEA. It matches transcript names with DoRothEA genes (regulons for the
# transcription factors) and provides information regarding how many transcripts are lost/kept in the computation. 
# TF activities are computed using DoRothEA package.

compute.TF.activity <- function(RNA.tpm){

  # ****************
  # packages
  if(!("viper" %in% installed.packages()[,"Package"])) BiocManager::install("viper", ask = FALSE)
  suppressMessages(require(viper))

  # ****************
  # scripts
  source("scaling_function.R")
  
  # ****************
  # data
  load("TFs_viperRegulon.rdata") # Load TF regulon genesets in VIPER format
 
  # Z score expression matrix of lo2(tpm + 1) (gene wisely as recommended by the authors)
  gene_expr <- standarization(log2(RNA.tpm + 1))
  
  # redefine gene names to match transcripts for viper
  E <- gene_expr
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
