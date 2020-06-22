#########################################################################################################
# Script to generate ImmuneResponse scores

# Output data -->
## ImmuneResponse (13 proxy's of the response): 1

# PanCancer analysis: 

# TCGA samples available for quanTIseq and IS 
## 18 Cancer types: BLCA BRCA CESC  CRC GBM HNSC KIRC KIRP LIHC LUAD LUSC OV PAAD PRAD SKCM STAD THCA UCEC
### TCGA_samples_available_screening_with_quanTIseq_IS.RData

# TCGA samples available for quanTIseq, IS and spatial TILs  (less samples)
## 11 Cancer types: BLCA BRCA CESC  CRC LUAD LUSC PAAD PRAD SKCM STAD UCEC
### TCGA_samples_available_screening_with_quanTIseq_IS_Spat.RData

# TCGA samples available for quanTIseq, IS and proteins (less samples)
## 18 Cancer types: BLCA BRCA CESC  CRC GBM HNSC KIRC KIRP LIHC LUAD LUSC OV PAAD PRAD SKCM STAD THCA UCEC
### TCGA_samples_available_screening_with_quanTIseq_IS_Spat.RData

# * CRC = COAD + READ
#########################################################################################################

# ****************
# working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ****************
# packages
library(preprocessCore)

# *****************
# functions
source("scaling_function.R")
# Derive CYT from Rooney
source("Signatures_response_ICBs/ObtainCYT.R")
# Derive IMPRES from Auslander
source("Signatures_response_ICBs/ObtainIMPRES.R")
# Derive IPS from Charaentong
source("Signatures_response_ICBs/ObtainIPS.R")
# Derive Roh Immune Signature Score
source("Signatures_response_ICBs/ObtainRohISS.R")
# Derive 12-Chemokine Signature Score
source("Signatures_response_ICBs/Obtain12chemokine.R")
# Derive Immune Signature Davoli
source("Signatures_response_ICBs/ObtainDavoliIS.R")
# Derive Proliferation Signature 
source("Signatures_response_ICBs/ObtainProliferation.R")
# Derive IFNy Signature Ayers
source("Signatures_response_ICBs/ObtainIFnyAyers.R")
# Derive Expanded Immune Signature Ayers
source("Signatures_response_ICBs/ObtainExpandedImmuneAyers.R")
# Derive T cell Inflamd GEP Ayers
source("Signatures_response_ICBs/ObtainTcellInflamedAyers.R")
# Derive TIDE score Jiang
source("Signatures_response_ICBs/ObtainTIDE.R")
# Derive MSI score Fu
source("Signatures_response_ICBs/ObtainMSI.R")

# *****************
# Cancer types
load("../analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
load("../analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_prot.RData")
load("../analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_Spat.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

#  ******************************************************************************************************* #
# ImmuneResponse: IS, CYT, IPS, IMPRES, (check other ones available)
# we need first time
ImmuneResponse.proxy_measures <- c("CYT", "IS", "IPS", "IMPRES", "RohIS", "chemokine", "IS_Davoli", "Proliferation", "IFny",
                                   "ExpandedImmune", "T_cell_inflamed","TIDE","MSI")
ImmuneResponse.no_filter <- vector("list", length = 13) ; names(ImmuneResponse.no_filter) <- ImmuneResponse.proxy_measures
ImmuneResponse.filter_spat <- vector("list", length = 13) ; names(ImmuneResponse.filter_spat) <- ImmuneResponse.proxy_measures
ImmuneResponse.filter_prot <- vector("list", length = 13) ; names(ImmuneResponse.filter_prot) <- ImmuneResponse.proxy_measures

# ****************
# IS #
# ****************
ImmuneSignature <- read.csv("../data/raw_data_tcga/immuneSignature.csv", row.names = 1, header = T, skip = 1)
ImmuneSignature <- ImmuneSignature[-c(9082:9083),]

sapply(PanCancer.names, function(Cancer){
  
  # Sample screening:
  # load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15), rownames(ImmuneSignature))
  tmp_IS <- ImmuneSignature[keep,]
  rownames(tmp_IS) <- substr(rownames(tmp_IS), 1, 12)
  colnames(tmp_IS)[2] <- "IS"
  
  ImmuneResponse.no_filter$IS <- data.frame(IS = as.matrix(tmp_IS[,"IS", drop = FALSE]), check.names = FALSE)
  
  ### Save 
  save(ImmuneResponse.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    # load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))
    
    keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15), rownames(ImmuneSignature))
    tmp_IS <- ImmuneSignature[keep,]
    rownames(tmp_IS) <- substr(rownames(tmp_IS), 1, 12)
    colnames(tmp_IS)[2] <- "IS"
    
    ImmuneResponse.filter_spat$IS <- data.frame(IS = tmp_IS[,"IS", drop = FALSE], check.names = FALSE)
    
    ### Save 
    save(ImmuneResponse.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))
    
  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    # load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
    keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15), rownames(ImmuneSignature))
    tmp_IS <- ImmuneSignature[keep,]
    rownames(tmp_IS) <- substr(rownames(tmp_IS), 1, 12)
    colnames(tmp_IS)[2] <- "IS"
    
    ImmuneResponse.filter_prot$IS <- data.frame(IS = tmp_IS[,"IS", drop = FALSE], check.names = FALSE)
    
    ### Save 
    save(ImmuneResponse.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }
  
  
})

# ********************************
# CYT --> Input: TPM
# ********************************
# ********************************
# IPS --> Input: TPM
# ********************************
# ********************************
# IMPRES --> Input: TPM
# ********************************
# ********************************
# RohIS --> Input: TPM
# ********************************
# ********************************
# Chemokine --> Input: TPM
# ********************************
# ********************************
# Davoli IS and proliferation signature --> Input: TPM
# ********************************
# ********************************
# Ayers IFny and Expanded Immune Signature --> Input: TPM
# ********************************
# ********************************
# Ayers T cell Inflamed GEP --> Input: TPM
# ********************************
# ********************************
# TIDE --> Input: TPM (pyhton tidepy required)
# ********************************
# ********************************
# MSI score --> Input: TPM
# ********************************

# Literature
CYT.paper <- read.csv("../data/raw_data_tcga/cytolyticActivity.csv", row.names = 1, header = T)
CYT.paper <- CYT.paper[,"Cytolytic.Activity", drop = FALSE]
IPS.paper <- read.csv("../data/raw_data_tcga/patientsAll_IPS_literature.csv", row.names = 1, header = TRUE, sep = "\t")
IPS.paper <- IPS.paper[,"ips_ctla4_neg_pd1_neg", drop = FALSE]
IMPRES.paper <- read.csv("../data/raw_data_tcga/IMPRES_TCGA.csv", row.names = 1, header = T)
Davoli.paper <- read.csv("../data/raw_data_tcga/IS_Proliferation_Davoli_TCGA.csv", header = T, skip = 5, row.names = 1)
IS_Davoli.paper <- Davoli.paper[,"Immune.Signature.Score", drop = FALSE]
IS_Davoli.paper <- IS_Davoli.paper[!is.na(IS_Davoli.paper[,1]), , drop = FALSE]
Proliferation_Davoli.paper <- Davoli.paper[,"CellCycle.Signature.Score", drop = FALSE]
Proliferation_Davoli.paper <- Proliferation_Davoli.paper[!is.na(Proliferation_Davoli.paper[,1]), , drop = FALSE]

sapply(PanCancer.names, function(Cancer){
  
  cat("\n",Cancer,"\n")
  
  file <- dir(paste0("../data/raw_data_tcga/RNAseq/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
                     "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0"),
              recursive = TRUE, pattern = "data.txt", full.names = TRUE)
  
  # Extract the raw counts from the text files for each gene
  data.transcripts <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
  colnames(data.transcripts) <- substr(colnames(data.transcripts), 1, 15)
  
  # TPM/scaled_estimate values (transcriptomics data)
  get_estimates <- which(data.transcripts[1,] == "scaled_estimate")
  estimates.transcripts <- data.frame(data.transcripts[-1,get_estimates])
  colnames(estimates.transcripts) <- gsub(".","-", colnames(estimates.transcripts), fixed = TRUE)
  genes <- rownames(estimates.transcripts)
  estimates.transcripts <- sapply(estimates.transcripts, as.numeric) # numeric
  
  # Obtaining TPM (transcriptomics data)
  TPM.transcripts <- estimates.transcripts * 1e6
  rownames(TPM.transcripts) <- genes
  
  # Rows with non-valid HGNC symbols were removed.
  HGNC_symbol <- sapply(strsplit(rownames(TPM.transcripts),"\\|"),function(X) return(X[1]))
  TPM.transcripts <-  TPM.transcripts[-which(HGNC_symbol == "?"),]
  HGNC_symbol <-HGNC_symbol[-which(HGNC_symbol == "?")]
  HGNC_id <- sapply(strsplit(rownames(TPM.transcripts),"\\|"),function(X) return(X[2]))
  
  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(HGNC_symbol) != 0){
    idx <- which(duplicated(HGNC_symbol) == TRUE)
    dup_genes <- HGNC_symbol[idx]
    for (ii in dup_genes){
      TPM.transcripts[which(HGNC_symbol %in% ii)[1],] <- colMeans(TPM.transcripts[which(HGNC_symbol %in% ii),])
      TPM.transcripts <- TPM.transcripts[-which(HGNC_symbol %in% ii)[2],]
      HGNC_symbol <- HGNC_symbol[-which(HGNC_symbol %in% ii)[2]]
    }
  }
  rownames(TPM.transcripts) <- HGNC_symbol
  
  # Compute CYT
  CYT <- ObtainCYT(TPM.transcripts)
  # Compute IPS
  IPS <- ObtainIPS(TPM.transcripts)
  # Compute IMPRES
  IMPRES <- ObtainIMPRES(TPM.transcripts)
  # Compute MSI
  MSI <- ObtainMSI(TPM.transcripts)
  # Compute Roh ISS score
  Roh_ISS <- ObtainRohISS(TPM.transcripts)
  # Compute 12-Chemokine Score
  Chemokine.score <- Obtain12chemokine(TPM.transcripts)
  # Compute Immune Signature Davoli
  IS_Davoli <- ObtainDavoliIS(TPM.transcripts)
  # Compute proliferation signature
  Proliferation <- ObtainProliferation(TPM.transcripts)
  # Compute IFN Ayers
  IFny_Ayers <- ObtainIFnyAyers(TPM.transcripts)
  # Compute Expanded Immune Signature Ayers
  ExpIS_Ayers <- ObtainExpandedImmuneAyers(TPM.transcripts)
  # Compute T-cell Inflamed GEP (Ayers)
  T_cell_inflamed_Ayers <- ObtainTcellInflamedAyers(TPM.transcripts)
  # Compute TIDE score
  TIDE.score <- ObtainTIDE(TPM.transcripts, Cancer)

  # ****************
  # CYT #
  # ****************
  
  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter
  load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  ImmuneResponse.CYT.no_filter <- CYT[keep.samples.no_filter, , drop = FALSE]
  rownames(ImmuneResponse.CYT.no_filter) <- substr(rownames(ImmuneResponse.CYT.no_filter), 1, 12)
  tmp <- intersect(rownames(ImmuneResponse.CYT.no_filter), rownames(CYT.paper[!is.na(CYT.paper$Cytolytic.Activity),1, drop = FALSE]))
  
  ### Check with author's score 
  #pairs(ImmuneResponse.CYT.no_filter[tmp,] ~ CYT.paper[tmp,], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8) 
  cat("\n","Correlation_no_filter =", cor(as.vector(ImmuneResponse.CYT.no_filter[tmp,]), as.vector(CYT.paper[tmp,])), "\n")
  
  ### Save 
  ImmuneResponse.no_filter$CYT <- ImmuneResponse.CYT.no_filter
  save(ImmuneResponse.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))
    
    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    ImmuneResponse.CYT.filter_Spat <- CYT[keep.samples.filter_Spat, , drop = FALSE]
    rownames(ImmuneResponse.CYT.filter_Spat) <- substr(rownames(ImmuneResponse.CYT.filter_Spat), 1, 12)
    tmp <- intersect(rownames(ImmuneResponse.CYT.filter_Spat), rownames(CYT.paper[!is.na(CYT.paper$Cytolytic.Activity),1, drop = FALSE]))
    
    ### Check with author's score 
    #pairs(ImmuneResponse.CYT.filter_Spat[tmp,] ~ CYT.paper[tmp,], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8) 
    cat("\n","Correlation_filter_Spat =", cor(as.vector(ImmuneResponse.CYT.filter_Spat[tmp,]), as.vector(CYT.paper[tmp,])), "\n")
    
    ### Save 
    ImmuneResponse.filter_spat$CYT <- ImmuneResponse.CYT.filter_Spat
    save(ImmuneResponse.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))
    
  }
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    ImmuneResponse.CYT.filter_prot <- CYT[keep.samples.filter_prot, , drop = FALSE]
    rownames(ImmuneResponse.CYT.filter_prot) <- substr(rownames(ImmuneResponse.CYT.filter_prot), 1, 12)
    tmp <- intersect(rownames(ImmuneResponse.CYT.filter_prot), rownames(CYT.paper[!is.na(CYT.paper$Cytolytic.Activity),1, drop = FALSE]))
    
    ### Check with author's score 
    #pairs(ImmuneResponse.CYT.filter_prot[tmp,] ~ CYT.paper[tmp,], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8) 
    cat("\n","Correlation_filter_prot =", cor(as.vector(ImmuneResponse.CYT.filter_prot[tmp,]), as.vector(CYT.paper[tmp,])), "\n")
    
    ### Save 
    ImmuneResponse.filter_prot$CYT <- ImmuneResponse.CYT.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }
  
  # ****************
  # IPS #
  # ****************

  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter

  load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]],1 ,15)
  IPS.no_filter <- IPS[keep.samples.no_filter, 1, drop = FALSE]
  rownames(IPS.no_filter) <- substr(rownames(IPS.no_filter), 1, 12)
  tmp <- intersect(rownames(IPS.no_filter), substr(rownames(IPS.paper), 1, 12))
  #pairs(IPS.no_filter[tmp,] ~ IPS.paper[tmp,], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
  cat("\n","Correlation_no_filter =", cor(as.vector(IPS.no_filter[tmp,]), as.vector(IPS.paper[tmp,]), method = "spearman"), "\n")

  ImmuneResponse.no_filter$IPS <- IPS.no_filter
  save(ImmuneResponse.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))

  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat

  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    IPS.filter_Spat <- IPS[keep.samples.filter_Spat, , drop = FALSE]
    rownames(IPS.filter_Spat) <- substr(rownames(IPS.filter_Spat), 1, 12)
    tmp <- intersect(rownames(IPS.filter_Spat), substr(rownames(IPS.paper), 1, 12))
    #pairs(IPS.filter_Spat[tmp,] ~ IPS.paper[tmp,], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    cat("\n","Correlation_filter_Spat =", cor(as.vector(IPS.filter_Spat[tmp,]), as.vector(IPS.paper[tmp,]), method = "spearman"), "\n")

    ImmuneResponse.filter_spat$IPS <- IPS.filter_Spat
    save(ImmuneResponse.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

  }
  ## Filter Prot
  # ------------------------------------------------------------------------------------------------------- #

  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    IPS.filter_prot <- IPS[keep.samples.filter_prot, , drop = FALSE]
    rownames(IPS.filter_prot) <- substr(rownames(IPS.filter_prot), 1, 12)
    tmp <- intersect( rownames(IPS.filter_prot), substr(rownames(IPS.paper), 1, 12))

    #pairs(IPS.filter_prot[tmp,] ~ IPS.paper[tmp,], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    cat("\n","Correlation_filter_prot =", cor(as.vector(IPS.filter_prot[tmp,]), as.vector(IPS.paper[tmp,]), method = "spearman"), "\n")

    ImmuneResponse.filter_prot$IPS <- IPS.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

  }

  # ****************
  # IMPRES
  # ****************

  # Sample screening:

  # ------------------------------------------------------------------------------------------------------- #
  ## No filter

  load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  IMPRES.no_filter <- IMPRES[keep.samples.no_filter, , drop = FALSE]
  rownames(IMPRES.no_filter) <- substr(rownames(IMPRES.no_filter), 1, 12)
  tmp <- intersect( rownames(IMPRES.no_filter), substr(rownames(IMPRES.paper), 1, 12))
  #pairs(IMPRES.no_filter[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
  cat("\n","Correlation_no_filter =", cor(as.vector(IMPRES.no_filter[tmp,]), as.vector(IMPRES.paper[tmp,2]), method = "spearman"), "\n")

  ImmuneResponse.no_filter$IMPRES <- IMPRES.no_filter
  save(ImmuneResponse.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))

  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    IMPRES.filter_Spat <- IMPRES[keep.samples.filter_Spat, , drop = FALSE]
    rownames(IMPRES.filter_Spat) <- substr(rownames(IMPRES.filter_Spat), 1, 12)
    tmp <- intersect( rownames(IMPRES.filter_Spat), substr(rownames(IMPRES.paper), 1, 12))
    #pairs(IMPRES.filter_Spat[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    cat("\n","Correlation_no_filter =", cor(as.vector(IMPRES.filter_Spat[tmp,]), as.vector(IMPRES.paper[tmp,2]), method = "spearman"), "\n")

    ImmuneResponse.filter_spat$IMPRES <- IMPRES.filter_Spat
    save(ImmuneResponse.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

  }

  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    IMPRES.filter_prot <- IMPRES[keep.samples.filter_prot, , drop = FALSE]
    rownames(IMPRES.filter_prot) <- substr(rownames(IMPRES.filter_prot), 1, 12)
    tmp <- intersect(rownames(IMPRES.filter_prot) , substr(rownames(IMPRES.paper), 1, 12))
    #pairs(IMPRES.filter_prot[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    cat("\n","Correlation_no_filter =", cor(as.vector(IMPRES.filter_prot[tmp,]), as.vector(IMPRES.paper[tmp,2]), method = "spearman"), "\n")

    ImmuneResponse.filter_prot$IMPRES <- IMPRES.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

  }

  # ****************
  # RohIS #
  # ****************

  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter

  load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  Roh_ISS.no_filter <- Roh_ISS[keep.samples.no_filter, 1, drop = FALSE]
  rownames(Roh_ISS.no_filter) <- substr(rownames(Roh_ISS.no_filter), 1, 12)

  ImmuneResponse.no_filter$RohIS <- Roh_ISS.no_filter
  save(ImmuneResponse.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))

  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    Roh_ISS.filter_Spat <- Roh_ISS[keep.samples.filter_Spat, , drop = FALSE]
    rownames(Roh_ISS.filter_Spat) <- substr(rownames(Roh_ISS.filter_Spat), 1, 12)

    ImmuneResponse.filter_spat$RohIS <- Roh_ISS.filter_Spat
    save(ImmuneResponse.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    Roh_ISS.filter_prot <- Roh_ISS[keep.samples.filter_prot, , drop = FALSE]
    rownames(Roh_ISS.filter_prot) <- substr(rownames(Roh_ISS.filter_prot), 1, 12)

    ImmuneResponse.filter_prot$RohIS <- Roh_ISS.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

  }
  # ****************
  # Chemokine  #
  # ****************

  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter

  load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  Messina.chemokine.no_filter <- Chemokine.score[keep.samples.no_filter, , drop = FALSE]
  rownames(Messina.chemokine.no_filter) <- substr(rownames(Messina.chemokine.no_filter), 1, 12)

  ImmuneResponse.no_filter$chemokine <- Messina.chemokine.no_filter
  save(ImmuneResponse.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))

  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    Messina.chemokine.filter_Spat <- Chemokine.score[keep.samples.filter_Spat, , drop = FALSE]
    rownames(Messina.chemokine.filter_Spat) <- substr(rownames(Messina.chemokine.filter_Spat), 1, 12)

    ImmuneResponse.filter_spat$chemokine <- Messina.chemokine.filter_Spat
    save(ImmuneResponse.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    Messina.chemokine.filter_prot <- Chemokine.score[keep.samples.filter_prot, , drop = FALSE]
    rownames(Messina.chemokine.filter_prot) <- substr(rownames(Messina.chemokine.filter_prot), 1, 12)

    ImmuneResponse.filter_prot$chemokine <- Messina.chemokine.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))
  }
  # *******************************************************
  # Davoli Immune Signature and proliferation signature #
  # *******************************************************

  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter

  load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  IS_Davoli.no_filter <- IS_Davoli[keep.samples.no_filter, , drop = FALSE]
  rownames(IS_Davoli.no_filter) <- substr(rownames(IS_Davoli.no_filter), 1, 12)
  tmp <- intersect(rownames(IS_Davoli.no_filter), substr(rownames(IS_Davoli.paper), 1, 12))

  #pairs(IS_Davoli.no_filter.norm[tmp,] ~ IS_Davoli.paper[tmp,1], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
  cat("\n","Correlation_no_filter IS =", cor(as.vector(IS_Davoli.no_filter[tmp,]), as.vector(IS_Davoli.paper[tmp,1]), method = "spearman"), "\n")

  Proliferation.no_filter <- Proliferation[keep.samples.no_filter, , drop = FALSE]
  rownames(Proliferation.no_filter) <- substr(rownames(Proliferation.no_filter), 1, 12)
  tmp <- intersect(rownames(Proliferation.no_filter) , substr(rownames(Proliferation_Davoli.paper), 1, 12))

  #pairs(IS_Davoli.no_filter.norm[tmp,] ~ IS_Davoli.paper[tmp,1], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
  cat("\n","Correlation_no_filter Proliferation  =", cor(as.vector(Proliferation.no_filter[tmp,]), as.vector(Proliferation_Davoli.paper[tmp,1]), method = "spearman"), "\n")

  ImmuneResponse.no_filter$IS_Davoli <- IS_Davoli.no_filter
  ImmuneResponse.no_filter$Proliferation <- Proliferation.no_filter

  save(ImmuneResponse.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))
    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)

    IS_Davoli.filter_Spat <- IS_Davoli[keep.samples.filter_Spat, , drop = FALSE]
    rownames(IS_Davoli.filter_Spat) <- substr(rownames(IS_Davoli.filter_Spat), 1, 12)
    tmp <- intersect( rownames(IS_Davoli.filter_Spat), substr(rownames(IS_Davoli.paper), 1, 12))

    #pairs(IMPRES.filter_Spat[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    cat("\n","Correlation.filter_Spat IS =", cor(as.vector(IS_Davoli.filter_Spat[tmp,]), as.vector(IS_Davoli.paper[tmp,1]), method = "spearman"), "\n")

    Proliferation.filter_Spat <- Proliferation[keep.samples.filter_Spat, , drop = FALSE]
    rownames(Proliferation.filter_Spat) <- substr(rownames(Proliferation.filter_Spat), 1, 12)
    tmp <- intersect(rownames(Proliferation.filter_Spat) , substr(rownames(Proliferation_Davoli.paper), 1, 12))

    #pairs(IS_Davoli.no_filter.norm[tmp,] ~ IS_Davoli.paper[tmp,1], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    cat("\n","Correlation.filter_Spat Proliferation  =", cor(as.vector(Proliferation.filter_Spat[tmp,]), as.vector(Proliferation_Davoli.paper[tmp,1]), method = "spearman"), "\n")

    ImmuneResponse.filter_spat$IS_Davoli <- IS_Davoli.filter_Spat
    ImmuneResponse.filter_spat$Proliferation <- Proliferation.filter_Spat
    save(ImmuneResponse.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    IS_Davoli.filter_prot <- IS_Davoli[keep.samples.filter_prot, , drop = FALSE]
    rownames(IS_Davoli.filter_prot) <- substr(rownames(IS_Davoli.filter_prot), 1, 12)
    tmp <- intersect(rownames(IS_Davoli.filter_prot), substr(rownames(IS_Davoli.paper), 1, 12))

    #pairs(IMPRES.filter_prot[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    cat("\n","Correlation.filter_prot IS =", cor(as.vector(IS_Davoli.filter_prot[tmp,]), as.vector(IS_Davoli.paper[tmp,1]), method = "spearman"), "\n")

    Proliferation.filter_prot <- Proliferation[keep.samples.filter_prot, , drop = FALSE]
    rownames(Proliferation.filter_prot) <- substr(rownames(Proliferation.filter_prot), 1, 12)
    tmp <- intersect(rownames(Proliferation.filter_prot), substr(rownames(Proliferation_Davoli.paper), 1, 12))

    #pairs(IS_Davoli.no_filter.norm[tmp,] ~ IS_Davoli.paper[tmp,1], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    cat("\n","Correlation.filter_prot Proliferation  =", cor(as.vector(Proliferation.filter_prot[tmp,]), as.vector(Proliferation_Davoli.paper[tmp,1]), method = "spearman"), "\n")

    ImmuneResponse.filter_prot$IS_Davoli <- IS_Davoli.filter_prot
    ImmuneResponse.filter_prot$Proliferation <- Proliferation.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

  }

  # ********************************************
  # Ayers IFny and Expanded Immune Signature  #
  # ********************************************

  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter

  load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  IFny_Ayers.no_filter <- IFny_Ayers[keep.samples.no_filter, , drop = FALSE]
  rownames(IFny_Ayers.no_filter) <- substr(rownames(IFny_Ayers.no_filter), 1, 12)
  #tmp <- intersect(rownames(IFny_Ayers.no_filter), substr(rownames(IS_Davoli.paper), 1, 12))

  #pairs(IFny_Ayers.no_filter[tmp,] ~ IS_Davoli.paper[tmp,1], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
  #cat("\n","Correlation_no_filter IFny =", cor(as.vector(IFny_Ayers.no_filter[tmp,]), as.vector(IS_Davoli.paper[tmp,1]), method = "spearman"), "\n")

  ExpIS_Ayers.no_filter <- ExpIS_Ayers[keep.samples.no_filter, , drop = FALSE]
  rownames(ExpIS_Ayers.no_filter) <- substr(rownames(ExpIS_Ayers.no_filter), 1, 12)
  #tmp <- intersect(rownames(ExpIS_Ayers.no_filter) , substr(rownames(Proliferation_Davoli.paper), 1, 12))

  #pairs(ExpIS_Ayers.no_filter[tmp,] ~ IS_Davoli.paper[tmp,1], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
  #cat("\n","Correlation_no_filter Expanded IS  =", cor(as.vector(ExpIS_Ayers.no_filter[tmp,]), as.vector(Proliferation_Davoli.paper[tmp,1]), method = "spearman"), "\n")

  ImmuneResponse.no_filter$IFny <- IFny_Ayers.no_filter
  ImmuneResponse.no_filter$ExpandedImmune <- ExpIS_Ayers.no_filter

  save(ImmuneResponse.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))
    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)

    IFny_Ayers.filter_Spat <- IFny_Ayers[keep.samples.filter_Spat, , drop = FALSE]
    rownames(IFny_Ayers.filter_Spat) <- substr(rownames(IFny_Ayers.filter_Spat), 1, 12)
    #tmp <- intersect( rownames(IFny_Ayers.filter_Spat), substr(rownames(IS_Davoli.paper), 1, 12))

    #pairs(IMPRES.filter_Spat[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    #cat("\n","Correlation.filter_Spat IS =", cor(as.vector(IFny_Ayers.filter_Spat[tmp,]), as.vector(IS_Davoli.paper[tmp,1]), method = "spearman"), "\n")

    ExpIS_Ayers.filter_Spat <- ExpIS_Ayers[keep.samples.filter_Spat, , drop = FALSE]
    rownames(ExpIS_Ayers.filter_Spat) <- substr(rownames(ExpIS_Ayers.filter_Spat), 1, 12)
    #tmp <- intersect(rownames(ExpIS_Ayers.filter_Spat) , substr(rownames(Proliferation_Davoli.paper), 1, 12))

    #pairs(IS_Davoli.no_filter.norm[tmp,] ~ IS_Davoli.paper[tmp,1], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    #cat("\n","Correlation.filter_Spat Proliferation  =", cor(as.vector(ExpIS_Ayers.filter_Spat[tmp,]), as.vector(Proliferation_Davoli.paper[tmp,1]), method = "spearman"), "\n")

    ImmuneResponse.filter_spat$IFny <- IFny_Ayers.filter_Spat
    ImmuneResponse.filter_spat$ExpandedImmune <- ExpIS_Ayers.filter_Spat
    save(ImmuneResponse.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    IFny_Ayers.filter_prot <- IFny_Ayers[keep.samples.filter_prot, , drop = FALSE]
    rownames(IFny_Ayers.filter_prot) <- substr(rownames(IFny_Ayers.filter_prot), 1, 12)
    #tmp <- intersect(rownames(IFny_Ayers.filter_prot), substr(rownames(IS_Davoli.paper), 1, 12))

    #pairs(IMPRES.filter_prot[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    #cat("\n","Correlation.filter_prot IS =", cor(as.vector(IS_Davoli.filter_prot[tmp,]), as.vector(IS_Davoli.paper[tmp,1]), method = "spearman"), "\n")

    ExpIS_Ayers.filter_prot <- ExpIS_Ayers[keep.samples.filter_prot, , drop = FALSE]
    rownames(ExpIS_Ayers.filter_prot) <- substr(rownames(ExpIS_Ayers.filter_prot), 1, 12)
    #tmp <- intersect(rownames(Proliferation.filter_prot), substr(rownames(Proliferation_Davoli.paper), 1, 12))

    #pairs(IS_Davoli.no_filter.norm[tmp,] ~ IS_Davoli.paper[tmp,1], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    #cat("\n","Correlation.filter_prot Proliferation  =", cor(as.vector(Proliferation.filter_prot[tmp,]), as.vector(Proliferation_Davoli.paper[tmp,1]), method = "spearman"), "\n")

    ImmuneResponse.filter_prot$IFny <- IFny_Ayers.filter_prot
    ImmuneResponse.filter_prot$ExpandedImmune <- ExpIS_Ayers.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

  }
  # ******************************
  # Ayers T cell Inflamed GEP #
  # ******************************

  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter

  load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  T_cell_inflamed_Ayers.no_filter <- T_cell_inflamed_Ayers[keep.samples.no_filter, , drop = FALSE]
  rownames(T_cell_inflamed_Ayers.no_filter) <- substr(rownames(T_cell_inflamed_Ayers.no_filter), 1, 12)

  ImmuneResponse.no_filter$T_cell_inflamed <- T_cell_inflamed_Ayers.no_filter
  save(ImmuneResponse.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))

  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    T_cell_inflamed_Ayers.filter_Spat <- T_cell_inflamed_Ayers[keep.samples.filter_Spat, , drop = FALSE]
    rownames(T_cell_inflamed_Ayers.filter_Spat) <- substr(rownames(T_cell_inflamed_Ayers.filter_Spat), 1, 12)

    ImmuneResponse.filter_spat$T_cell_inflamed <- T_cell_inflamed_Ayers.filter_Spat
    save(ImmuneResponse.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    T_cell_inflamed_Ayers.filter_prot <- T_cell_inflamed_Ayers[keep.samples.filter_prot, , drop = FALSE]
    rownames(T_cell_inflamed_Ayers.filter_prot) <- substr(rownames(T_cell_inflamed_Ayers.filter_prot), 1, 12)

    ImmuneResponse.filter_prot$T_cell_inflamed <- T_cell_inflamed_Ayers.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))
  }

  # ****************
  # TIDE
  # ****************
  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter

  load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  TIDE.no_filter <- TIDE.score[keep.samples.no_filter, , drop = FALSE]

  ImmuneResponse.no_filter$TIDE <- TIDE.no_filter
  save(ImmuneResponse.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))


  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))
    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    TIDE.filter_Spat <- TIDE.score[keep.samples.filter_Spat, , drop = FALSE]

    ImmuneResponse.filter_spat$TIDE <- TIDE.filter_Spat
    save(ImmuneResponse.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    TIDE.filter_prot <- TIDE.score[keep.samples.filter_prot, , drop = FALSE]

    ImmuneResponse.filter_prot$TIDE <- TIDE.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

  }

  # ****************
  # MSI score
  # ****************
  # Sample screening:

  # ------------------------------------------------------------------------------------------------------- #
  ## No filter

  load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  MSI.no_filter <- MSI[keep.samples.no_filter, , drop = FALSE]
  rownames(MSI.no_filter) <- substr(rownames(MSI.no_filter), 1, 12)
  #tmp <- intersect( rownames(MSI.no_filter), substr(rownames(IMPRES.paper), 1, 12))
  #pairs(MSI.no_filter[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
  #cat("\n","Correlation_no_filter =", cor(as.vector(MSI.no_filter[tmp,]), as.vector(IMPRES.paper[tmp,2]), method = "spearman"), "\n")

  ImmuneResponse.no_filter$MSI <- MSI.no_filter
  save(ImmuneResponse.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))

  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    MSI.filter_Spat <- MSI[keep.samples.filter_Spat, , drop = FALSE]
    rownames(MSI.filter_Spat) <- substr(rownames(MSI.filter_Spat), 1, 12)
    #tmp <- intersect( rownames(MSI.filter_Spat), substr(rownames(IMPRES.paper), 1, 12))
    #pairs(MSI.filter_Spat[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    #cat("\n","Correlation_no_filter =", cor(as.vector(MSI.filter_Spat[tmp,]), as.vector(IMPRES.paper[tmp,2]), method = "spearman"), "\n")

    ImmuneResponse.filter_spat$MSI <- MSI.filter_Spat
    save(ImmuneResponse.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))

  }

  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){

    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    MSI.filter_prot <- MSI[keep.samples.filter_prot, , drop = FALSE]
    rownames(MSI.filter_prot) <- substr(rownames(MSI.filter_prot), 1, 12)
    #tmp <- intersect(rownames(MSI.filter_prot) , substr(rownames(IMPRES.paper), 1, 12))
    #pairs(MSI.filter_prot[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    #cat("\n","Correlation_no_filter =", cor(as.vector(MSI.filter_prot[tmp,]), as.vector(IMPRES.paper[tmp,2]), method = "spearman"), "\n")

    ImmuneResponse.filter_prot$MSI <- MSI.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))

  }

})


# ***************************************
# Transformation from list into matrix
# ***************************************

# sapply(PanCancer.names, function(Cancer){
# 
#   cat("\n",Cancer,"\n")
# 
#   load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
#   load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))
#   load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))
# 
#   # ------------------------------------------------------------------------------------------------------- #
#   ## No filter
#   ImmuneResponse.no_filter <- do.call(cbind, lapply(names(ImmuneResponse.no_filter), function(task){
# 
#     return(ImmuneResponse.no_filter[[task]])
# 
#   }))
#   ImmuneResponse.no_filter <- as.matrix(ImmuneResponse.no_filter)
# 
#   # ------------------------------------------------------------------------------------------------------- #
#   ## Filter Spat
#   if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
# 
#     ImmuneResponse.filter_spat <- do.call(cbind, sapply(names(ImmuneResponse.filter_spat), function(task){
# 
#       return(ImmuneResponse.filter_spat[[task]])
# 
#     }))
#     ImmuneResponse.filter_spat <- as.matrix(ImmuneResponse.filter_spat)
#   }
#   # ------------------------------------------------------------------------------------------------------- #
#   ## Filter Prot
#   if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
# 
#     ImmuneResponse.filter_prot <- do.call(cbind, sapply(names(ImmuneResponse.filter_prot), function(task){
# 
#       return(ImmuneResponse.filter_prot[[task]])
# 
#     }))
#     ImmuneResponse.filter_prot <- as.matrix(ImmuneResponse.filter_prot)
#   }
# 
#   save(ImmuneResponse.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
#   save(ImmuneResponse.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))
#   save(ImmuneResponse.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))
# 
# })

# sapply(PanCancer.names, function(Cancer){
# 
#   cat("\n",Cancer,"\n")
# 
#   # ------------------------------------------------------------------------------------------------------- #
#   ## No filter
# 
#   load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
#   load(paste0("../data/PanCancer_draft_v1/", Cancer,"/DataViews_no_filter_", Cancer,".RData"))
# 
#   #rownames(ImmuneResponse.no_filter) <- rownames(DataViews.no_filter$Pathways)
# 
#   # ------------------------------------------------------------------------------------------------------- #
#   ## Filter Spat
#   if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
# 
#     load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))
#     load(paste0("../data/PanCancer_draft_v1/", Cancer,"/DataViews_filter_spat_", Cancer,".RData"))
# 
#     rownames(ImmuneResponse.filter_spat) <- rownames(DataViews.filter_spat$ImmuneCells)
#     colnames(ImmuneResponse.filter_spat) <- sapply(strsplit(colnames(ImmuneResponse.filter_spat),
#                                                             split =".", fixed = T), head,1)
#     
#   }
#   # ------------------------------------------------------------------------------------------------------- #
#   ## Filter Prot
#   if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
# 
#     load(paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))
#     load(paste0("../data/PanCancer_draft_v1/", Cancer,"/DataViews_filter_prot_", Cancer,".RData"))
# 
#     rownames(ImmuneResponse.filter_prot) <- rownames(DataViews.filter_prot$Pathways)
#     colnames(ImmuneResponse.filter_prot) <- sapply(strsplit(colnames(ImmuneResponse.filter_prot),
#                                                             split =".", fixed = T), head,1)
#   }
#   save(ImmuneResponse.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_no_filter_",Cancer, ".RData"))
#   save(ImmuneResponse.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_",Cancer, ".RData"))
#   save(ImmuneResponse.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_prot_",Cancer, ".RData"))
# 
# })



