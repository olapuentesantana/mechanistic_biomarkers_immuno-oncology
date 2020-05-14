#---------------------------------------------------------------------------------------------------------------#
# Input: Transcriptomics data in two different formats: Raw counts and TPM. Different scores required different formats.
# A list with the names of the scores to be computed has to be provided
# Output: Gold standards scores

#---------------------------------------------------------------------------------------------------------------#

computation.gold.standards <- function(RNA.raw_counts, RNA.tpm, list_gold_standards, ....){
  
  # ****************
  # scripts
  source("../R/computation.CYT.R")
  source("../R/computation.ICB_genes.R")
  
  # calculate Immune Checkpoint genes expression #
  ICB_genes <- computation.ICB_genes(RNA.tpm)
  
  gold.standards <- sapply(list_gold_standards, function(X){
    if ("impres" == X){
        # calculate Impres #
        IMPRES <- computation.IMPRES(RNA.tpm)
        return(list(IMPRES))
    }else if("IPS" == X){
        # calculate Immunophenoscore #
        IPS <- computation.IPS(RNA.tpm)
        return(list(IPS))
    }else if("CYT" == X){
        # calculate Cytolytic activity #
        CYT <- computation.CYT(RNA.tpm)
        return(list(CYT))
    }else if("CTLA4" == X){
        CTLA4 <- ICB_genes$CTLA4
        return(list(CTLA4))
    }else if("PD1" == X){
        PD1 <- ICB_genes$PD1
        return(list(PD1))
    }else if("PDL1" == X){
        PDL1 <- ICB_genes$PDL1
        return(list(PDL1))
    }
  })
  
  # calculate Immune Signature #
  
  # calculate Consensus (mean IS and CYT) #
  

  return(gold.standards)
  
}

