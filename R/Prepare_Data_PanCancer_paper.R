#########################################################################################################
# Script to generate DataViews {with (11 cancer types) and without sTIL (18 cancer types)}:

# Input data -->
  ## Mechanistic DataViews:
    ### Pathways (PROGENy)
    ### Immunecells (quanTIseq)
    ### TFs (DoRothEAv1)
    ### sTIL (SpatialTIL)
    ### LRpairs (Ligand-Receptor pairs)
    ### CYTOKINEpairs (Cytokine pairs)
  
  ## RNA_PROT DataViews:
    ### transcript (transcriptomics)
    ### Protall (proteomics)

# Output data -->
  ## ImmuneResponse (proxy's of the response)
    ### IS, CYT, IPS, IMPRES, 

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
library(gdata)
# BiocManager::install("DESeq2")
#library(DESeq2)
# BiocManager::install("progeny")
#library(progeny)
# BiocManager::install("dorothea")
#library(dorothea)


# *****************
# functions
source("scaling_function.R")
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
# Compute pathway activity
source("compute.pathways.scores.R")
# Compute TFs activity
source("compute.TF.activity.R")
# Compute LR pairs
source("compute.LR.pairs.R")

# ****************
# functions from Federica
# Plot correlations including p-value, correlation value and correlation line
panel.cor <- function(x, y, digits=2, font.cor = 1, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor.test(x,y)$estimate
  p <- cor.test(x,y)$p.value
  txt_r <- format(r, digits=digits)
  txt_p <- format(p, scientific = TRUE, digits=digits)
  txt <- paste("cor=", txt_r, "\np=", txt_p, sep="")

  if(txt_r >= 0.7 & txt_p >= 0.05) font.cor <- 2

  text(0.5, 0.5, txt, cex = 1, font = font.cor)
}

panel.lm <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 0.8, col.smooth = "#A1A1A1", ...) {

  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  # abline(stats::lm(y ~ x),  col = col.smooth, ...)
  abline(a=0, b=1,  col = col.smooth, ...)
}

# *****************
# Cancer types
load("../analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
load("../analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_prot.RData")
load("../analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_Spat.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)
PanCancer.names.some <- PanCancer.names[-c(1:6)]

# *****************
# Initialize DataViews
DataViews.no_filter <- vector("list", length = 7)
names(DataViews.no_filter) <- c("pathways", "immunecells", "TFs", "sTIL", "LRpairs", "CYTOKINEpairs","transcript")
DataViews.filter_prot <- vector("list", length = 4)
names(DataViews.filter_prot) <- c("transcript", "Protall", "pathways", "TFs")

# ***************
# Remove transcripts used to build ImmuneResponse (IS,CYT,IPS,IMPRES,RohISS,Chemokine,Proliferation,IS_Davoli,IFNy,ExpandedImmune,
# T_cell_inflamed,TIDE,MSI)

# Genes to remove according to all ICB proxy's 
load("../data/all_genes_ICB_proxies.RData")

# ***************
# RNAseq data --> Pahways, TFs

sapply(PanCancer.names, function(Cancer){
  
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
  # Log2 transformed
  log2mas1.TPM.transcripts <- log2(TPM.transcripts + 1)
  
  # Rows with non-valid HGNC symbols were removed.
  HGNC_symbol <- sapply(strsplit(rownames(log2mas1.TPM.transcripts),"\\|"),function(X) return(X[1]))
  log2mas1.TPM.transcripts <-  log2mas1.TPM.transcripts[-which(HGNC_symbol == "?"),]
  HGNC_symbol <-HGNC_symbol[-which(HGNC_symbol == "?")]
  HGNC_id <- sapply(strsplit(rownames(log2mas1.TPM.transcripts),"\\|"),function(X) return(X[2]))
  
  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(HGNC_symbol) != 0){
    idx <- which(duplicated(HGNC_symbol) == TRUE)
    dup_genes <- HGNC_symbol[idx]
    for (ii in dup_genes){
      log2mas1.TPM.transcripts[which(HGNC_symbol %in% ii)[1],] <- colMeans(log2mas1.TPM.transcripts[which(HGNC_symbol %in% ii),])
      log2mas1.TPM.transcripts <- log2mas1.TPM.transcripts[-which(HGNC_symbol %in% ii)[2],]
      HGNC_symbol <- HGNC_symbol[-which(HGNC_symbol %in% ii)[2]]
    }
  }
  
  # Remove ImmuneResponse genes
  rownames(log2mas1.TPM.transcripts) <- HGNC_symbol
  remove.genes <- na.exclude(match(all_genes_to_remove, rownames(log2mas1.TPM.transcripts)))
  log2mas1.TPM.transcripts <- log2mas1.TPM.transcripts[-remove.genes,]
  cat("Removing signatures genes for proxy's of ICB response:  \n")
  
  # Sample screening:
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
  
  transcript.no_filter <- log2mas1.TPM.transcripts[, keep.samples.no_filter]
  transcript.filter_prot <- log2mas1.TPM.transcripts[, keep.samples.filter_prot]
  
  load(paste0("../data/PanCancer/",Cancer,"/new_remove_all_genes/DataViews_no_filter_",Cancer, ".RData"))
  load(paste0("../data/PanCancer/",Cancer,"/new_remove_all_genes/DataViews_filter_prot_",Cancer, ".RData"))
  
  DataViews.no_filter$transcript <- as.data.frame(t(transcript.no_filter))
  DataViews.filter_prot$transcript <- as.data.frame(t(transcript.filter_prot))
  
  # ---------------------------------------------------------------------------------- #
  # Obtain DoRothEA scores (the function does log2(tpm+1) itself)
  # plot(rowMeans(rawcounts.transcripts), rowSds(as.matrix(rawcounts.transcripts)))
  # plot(rowMeans(log2mas1.TPM.transcripts), rowSds(as.matrix(log2mas1.TPM.transcripts)))
  # plot(rowMeans(TF_activity), rowSds(as.matrix(TF_activity)))
  
  # Remove ImmuneResponse genes (the function should take care of it)
  
  # Sample screening:
  tpm.no_filter <- TPM.transcripts[, keep.samples.no_filter]
  tpm.filter_prot <- TPM.transcripts[, keep.samples.filter_prot]

  ## TF activity computation
  
  # Computation of TF activity (input matrix [genes, samples], ouput matrix [sample, TFs])
  TF_activity.no_filter <- compute.TF.activity(RNA.tpm = tpm.no_filter, remove.genes.ICB_proxies=TRUE)
  TF_activity.filter_prot <- compute.TF.activity(RNA.tpm = tpm.filter_prot, remove.genes.ICB_proxies=TRUE)

  # Insert into DataViews
  DataViews.no_filter$TFs <- as.data.frame(TF_activity.no_filter$scores)
  DataViews.filter_prot$TFs <- as.data.frame(TF_activity.filter_prot$scores)
  
  # ----------------------------------------------------------- #
  # Obtaining raw counts (Pathways data)
  get_rawcounts <- which(data.transcripts[1,] == "raw_count")
  rawcounts.transcripts <- data.frame(data.transcripts[-1,get_rawcounts])
  colnames(rawcounts.transcripts) <- gsub(".","-", colnames(rawcounts.transcripts), fixed = TRUE)
  genes <- rownames(rawcounts.transcripts)
  rawcounts.transcripts <- sapply(rawcounts.transcripts, as.numeric) # numeric
  sapply(rawcounts.transcripts, class) # numeric
  rownames(rawcounts.transcripts) <- genes
  
  # Remove ImmuneResponse genes (the function should take care of it)
  
  # Sample screening:
  rawcounts.no_filter <- rawcounts.transcripts[, keep.samples.no_filter]
  rawcounts.filter_prot <- rawcounts.transcripts[, keep.samples.filter_prot]

  ## Pathways activity computation
  
  # Computation of TF activity (input matrix [genes, samples], ouput matrix [sample, TFs])
  Pathways.activity.no_filter <- compute.pathways.scores(RNA.raw_counts=rawcounts.no_filter, remove.genes.ICB_proxies=TRUE)
  Pathways.activity.filter_prot <- compute.pathways.scores(RNA.raw_counts=rawcounts.filter_prot, remove.genes.ICB_proxies=TRUE)
  
  # Insert into DataViews
  DataViews.no_filter$pathways <- as.data.frame(Pathways.activity.no_filter)
  DataViews.filter_prot$pathways <- as.data.frame(Pathways.activity.filter_prot)
  
  save(DataViews.no_filter, file = paste0("../data/PanCancer/",Cancer,"/new_remove_all_genes/DataViews_no_filter_",Cancer, ".RData"))
  save(DataViews.filter_prot, file = paste0("../data/PanCancer/",Cancer,"/new_remove_all_genes/DataViews_filter_prot_",Cancer, ".RData"))
  
})

# ***************
# Inter-cellular networks data --> L.R pairs and Cytokine pairs

sapply(PanCancer.names, function(Cancer){
  
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
  
  # Sample screening:
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
  
  load(paste0("../data/PanCancer/",Cancer,"/new_remove_all_genes/DataViews_no_filter_",Cancer, ".RData"))
  load(paste0("../data/PanCancer/",Cancer,"/new_remove_all_genes/DataViews_filter_prot_",Cancer, ".RData"))
  
  # Remove ImmuneResponse genes (the function should take care of it)
  
  # Sample screening:
  tpm.no_filter <- TPM.transcripts[, keep.samples.no_filter]
  tpm.filter_prot <- TPM.transcripts[, keep.samples.filter_prot]
  
  ## TF activity computation
  
  # Computation of LR pairs activity (input matrix [genes, samples], ouput matrix [sample, TFs])
  LRpairs.no_filter <- compute.LR.pairs(RNA.tpm = tpm.no_filter, remove.genes.ICB_proxies=TRUE)
  LRpairs.filter_prot <- compute.LR.pairs(RNA.tpm = tpm.filter_prot, remove.genes.ICB_proxies=TRUE)
  
  # Insert into DataViews
  DataViews.no_filter$LRpairs <- as.data.frame(LRpairs.no_filter)
  DataViews.filter_prot$LRpairs <- as.data.frame(LRpairs.filter_prot)
  
  save(DataViews.no_filter, file = paste0("../data/PanCancer/",Cancer,"/new_remove_all_genes/DataViews_no_filter_",Cancer, ".RData"))
  save(DataViews.filter_prot, file = paste0("../data/PanCancer/",Cancer,"/new_remove_all_genes/DataViews_filter_prot_",Cancer, ".RData"))
  
})

sapply(PanCancer.names, function(Cancer){
  
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
  
  # Sample screening:
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
  
  load(paste0("../data/PanCancer/",Cancer,"/new_keep_all_genes/DataViews_no_filter_",Cancer, ".RData"))
  load(paste0("../data/PanCancer/",Cancer,"/new_keep_all_genes/DataViews_filter_prot_",Cancer, ".RData"))
  
  # Remove ImmuneResponse genes (the function should take care of it)
  
  # Sample screening:
  tpm.no_filter <- TPM.transcripts[, keep.samples.no_filter]
  tpm.filter_prot <- TPM.transcripts[, keep.samples.filter_prot]
  
  ## Ligand-Receptor pairs computation
  
  # Computation of LR pairs activity (input matrix [genes, samples], ouput matrix [sample, TFs])
  LRpairs.no_filter <- compute.LR.pairs(RNA.tpm = tpm.no_filter, remove.genes.ICB_proxies=FALSE)
  LRpairs.filter_prot <- compute.LR.pairs(RNA.tpm = tpm.filter_prot, remove.genes.ICB_proxies=FALSE)
  
  # Insert into DataViews
  DataViews.no_filter$LRpairs <- as.data.frame(LRpairs.no_filter)
  DataViews.filter_prot$LRpairs <- as.data.frame(LRpairs.filter_prot)
  
  save(DataViews.no_filter, file = paste0("../data/PanCancer/",Cancer,"/new_keep_all_genes/DataViews_no_filter_",Cancer, ".RData"))
  save(DataViews.filter_prot, file = paste0("../data/PanCancer/",Cancer,"/new_keep_all_genes/DataViews_filter_prot_",Cancer, ".RData"))
  
})




# ***************
# Protein abundance (RPPA) data --> Proteomics

all.prot.tcpa <- read.csv("./data/raw_data/TCGA-PANCAN32-L4.csv", row.names = 1 , header = T, stringsAsFactors = F, check.names = T)

# Cancer types
load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_prot.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)

sapply(PanCancer.names, function(ii){
  
  load(paste0("./data/PanCancer/",ii,"/new/DataViews_filter_prot_",ii, ".RData"))
  
  keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[ii]], 1, 15), substr(rownames(all.prot.tcpa), 1, 15))
  tmp_prot <- all.prot.tcpa[keep,]
  
  DataViews.filter_prot$Protall <- as.data.frame(tmp_prot)

  #save(DataViews.filter_prot, file = paste0("./data/PanCancer/",ii,"/new/DataViews_filter_prot_",ii, ".RData"))
  
})

# ***************
# Immune cells data (quanTIseq) and Spatial TILs data (Saltz et al.)

# data
all_cell.fractions <- read.csv("./data/raw_data_tcga/quanTIseq_estimated.csv", header = TRUE, row.names = 1)
all.spatial.TILs <- read.csv("./data/raw_data_tcga/spatial_TILs_saltz.csv", header = TRUE, row.names = 1)


# ------------------------------------------------------------------------------------------------------- #
# Match samples in cell fractions data for general analysis: screening with quanTIseq_IS
# ------------------------------------------------------------------------------------------------------- #

# Cancer types
load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

for (ii in PanCancer.names){
  
  load(paste0("./data/PanCancer/",ii,"/new/DataViews_no_filter_",ii, ".RData"))
  
  keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[ii]], 1, 15), 
                rownames(all_cell.fractions))
  tmp_immunecells <- all_cell.fractions[keep,]
  
  if (nrow(tmp_immunecells) == nrow(DataViews.no_filter$pathways)){
    
    tmp_immunecells$T_cells_CD4 <- tmp_immunecells$T_cells_CD4 + tmp_immunecells$T_cells_regulatory_Tregs
    DataViews.no_filter$immunecells <- tmp_immunecells
  }else{
    
    stop("pathways sample size  != immune cells sample size")
    
    break
  }
  
  #save(DataViews.no_filter, file = paste0("./data/PanCancer/",ii,"/new/DataViews_no_filter_",ii, ".RData"))
}

# ------------------------------------------------------------------------------------------------------- #
# Match samples in cell fractions data for spatial information analysis: screening with quanTIseq_IS_Spat
# ------------------------------------------------------------------------------------------------------- #

# Cancer types
load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_Spat.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)


for (ii in PanCancer.names){
  
  DataViews.filter_Spat <- vector("list", length = 2); names(DataViews.filter_Spat) <- c("immunecells", "sTIL")
  
  keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[ii]], 1, 15), 
                substr(rownames(all_cell.fractions), 1, 15))
  
  tmp_immunecells <- all_cell.fractions[keep,]
  rownames(tmp_immunecells) <- substr(rownames(tmp_immunecells),1,12)
  
  keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[ii]], 1, 12), 
                substr(rownames(all.spatial.TILs), 1, 12))
  
  tmp_STILs <- all.spatial.TILs[keep,]
  features_to_keep <- c("til_percentage", "WCD_mean", "Ball_Hall", "Banfeld_Raftery", "C_index", "Det_Ratio")
  tmp_STILs <- tmp_STILs[,features_to_keep]
  
  tmp_immunecells$T_cells_CD4 <- tmp_immunecells$T_cells_CD4 + tmp_immunecells$T_cells_regulatory_Tregs
  DataViews.filter_Spat$immunecells <- tmp_immunecells
  DataViews.filter_Spat$sTIL <- tmp_STILs
  
  #save(DataViews.filter_Spat, file = paste0("./data/PanCancer/",ii,"/new/DataViews_filter_Spat_",ii, ".RData"))
}

#  ******************************************************************************************************* #
# ImmuneResponse: IS, CYT, IPS, IMPRES, (check other ones available)
ImmuneResponse.no_filter <- vector("list", length = 4) ; names(ImmuneResponse.no_filter) <- c("CYT", "IS", "IPS", "IMPRES")
ImmuneResponse.filter_Spat <- vector("list", length = 4) ; names(ImmuneResponse.filter_Spat) <- c("CYT", "IS", "IPS", "IMPRES")
ImmuneResponse.filter_prot <- vector("list", length = 4) ; names(ImmuneResponse.filter_prot) <- c("CYT", "IS", "IPS", "IMPRES")

# ****************
# Cytolytic Activity #

# ------------------------------------------------------------------------------------------------------- #
# Match samples in RNAseq for immuneresponse: screening with quanTIseq_IS: 
# ------------------------------------------------------------------------------------------------------- #
load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)
CYT.paper <- read.csv("data/raw_data/cytolyticActivity.csv", row.names = 1, header = T)
CYT.paper <- CYT.paper[,"Cytolytic.Activity", drop = FALSE]

sapply(PanCancer.names, function(Cancer){
  
  cat("\n",Cancer,"\n")
  ImmuneResponse.no_filter <- vector("list", length = 4) ; names(ImmuneResponse.no_filter) <- c("CYT", "IS", "IPS", "IMPRES")
  ImmuneResponse.filter_Spat <- vector("list", length = 4) ; names(ImmuneResponse.filter_Spat) <- c("CYT", "IS", "IPS", "IMPRES")
  ImmuneResponse.filter_prot <- vector("list", length = 4) ; names(ImmuneResponse.filter_prot) <- c("CYT", "IS", "IPS", "IMPRES")
  
  file <- dir(paste0("~/Desktop/PhD_TUE/Github_model/desktop/data/raw_data/RNAseq/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
                     "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0"),
              recursive = TRUE, pattern = "data.txt", full.names = TRUE)
  
  # Extract the raw counts from the text files for each gene
  data.transcripts <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
  colnames(data.transcripts) <- substr(colnames(data.transcripts), 1, 15)
  
  # TPM/scaled_estimate values (transcriptomics data)
  get_estimates <- which(data.transcripts[1,] == "scaled_estimate")
  estimates.transcripts <- data.frame(data.matrix(data.transcripts[-1,get_estimates]))
  colnames(estimates.transcripts) <- gsub(".","-", colnames(estimates.transcripts), fixed = TRUE)
  sapply(estimates.transcripts, class) # numeric
  
  # Obtaining TPM (transcriptomics data)
  TPM.transcripts <- estimates.transcripts * 1e6
  # Obtaining 0.01 offset
  TPM.transcripts <- TPM.transcripts + 0.01

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
  # Remove ImmuneResponse genes
  rownames(TPM.transcripts) <- HGNC_symbol
  
  # GZMA and PRF1 expression:
  GZMA_expr <- TPM.transcripts["GZMA",]
  PRF1_expr <- TPM.transcripts["PRF1",]
  
  # Geometric mean calculation:
  geometric_mean <- sqrt(GZMA_expr * PRF1_expr)
  geometric_mean <- t(geometric_mean)
  
  colnames(geometric_mean) <- "CYT"
  ImmuneResponse.CYT <- data.frame(CYT = geometric_mean, check.names = FALSE)

  # Sample screening:
  
  ## No filter
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  ImmuneResponse.CYT.no_filter <- ImmuneResponse.CYT[keep.samples.no_filter, , drop = FALSE]
  rownames(ImmuneResponse.CYT.no_filter) <- substr(rownames(ImmuneResponse.CYT.no_filter), 1, 12)
  tmp <- intersect(rownames(ImmuneResponse.CYT.no_filter), rownames(CYT.paper[!is.na(CYT.paper$Cytolytic.Activity),1, drop = FALSE]))
  
  ### Check with author's score 
  #pairs(ImmuneResponse.CYT.no_filter[tmp,] ~ CYT.paper[tmp,], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8) 
  cat("\n","Correlation_no_filter =", cor(as.vector(ImmuneResponse.CYT.no_filter[tmp,]), as.vector(CYT.paper[tmp,])), "\n")
  
  ### Save 
  ImmuneResponse.no_filter$CYT <- ImmuneResponse.CYT.no_filter
  save(ImmuneResponse.no_filter, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    ImmuneResponse.CYT.filter_Spat <- ImmuneResponse.CYT[keep.samples.filter_Spat, , drop = FALSE]
    rownames(ImmuneResponse.CYT.filter_Spat) <- substr(rownames(ImmuneResponse.CYT.filter_Spat), 1, 12)
    tmp <- intersect(rownames(ImmuneResponse.CYT.filter_Spat), rownames(CYT.paper[!is.na(CYT.paper$Cytolytic.Activity),1, drop = FALSE]))
    
    ### Check with author's score 
    #pairs(ImmuneResponse.CYT.filter_Spat[tmp,] ~ CYT.paper[tmp,], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8) 
    cat("\n","Correlation_filter_Spat =", cor(as.vector(ImmuneResponse.CYT.filter_Spat[tmp,]), as.vector(CYT.paper[tmp,])), "\n")
    
    ### Save 
    ImmuneResponse.filter_Spat$CYT <- ImmuneResponse.CYT.filter_Spat
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }else{
    
    ImmuneResponse.filter_Spat$CYT <- NULL
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    ImmuneResponse.CYT.filter_prot <- ImmuneResponse.CYT[keep.samples.filter_prot, , drop = FALSE]
    rownames(ImmuneResponse.CYT.filter_prot) <- substr(rownames(ImmuneResponse.CYT.filter_prot), 1, 12)
    tmp <- intersect(rownames(ImmuneResponse.CYT.filter_prot), rownames(CYT.paper[!is.na(CYT.paper$Cytolytic.Activity),1, drop = FALSE]))
    
    ### Check with author's score 
    #pairs(ImmuneResponse.CYT.filter_prot[tmp,] ~ CYT.paper[tmp,], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8) 
    cat("\n","Correlation_filter_prot =", cor(as.vector(ImmuneResponse.CYT.filter_prot[tmp,]), as.vector(CYT.paper[tmp,])), "\n")
    
    ### Save 
    ImmuneResponse.filter_prot$CYT <- ImmuneResponse.CYT.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }else{
    
    ImmuneResponse.filter_prot$CYT <- NULL
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }
  
})

# ****************
# Immune Signature #
ImmuneSignature <- read.csv("data/raw_data/immuneSignature.csv", row.names = 1, header = T, skip = 1)
ImmuneSignature <- ImmuneSignature[-c(9082:9083),]

# Cancer types
load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

sapply(PanCancer.names, function(Cancer){
  
  # Sample screening:
  load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15), rownames(ImmuneSignature))
  tmp_IS <- ImmuneSignature[keep,]
  rownames(tmp_IS) <- substr(rownames(tmp_IS), 1, 12)
  colnames(tmp_IS)[2] <- "IS"

  ImmuneResponse.no_filter$IS <- data.frame(IS = as.matrix(tmp_IS[,"IS", drop = FALSE]), check.names = FALSE)
  
  ### Save 
  save(ImmuneResponse.no_filter, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
    keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15), rownames(ImmuneSignature))
    tmp_IS <- ImmuneSignature[keep,]
    rownames(tmp_IS) <- substr(rownames(tmp_IS), 1, 12)
    colnames(tmp_IS)[2] <- "IS"
    
    ImmuneResponse.filter_Spat$IS <- data.frame(IS = tmp_IS[,"IS", drop = FALSE], check.names = FALSE)
    
    ### Save 
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }else{
    
    ImmuneResponse.filter_Spat$IS <- NULL
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){

    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
    keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15), rownames(ImmuneSignature))
    tmp_IS <- ImmuneSignature[keep,]
    rownames(tmp_IS) <- substr(rownames(tmp_IS), 1, 12)
    colnames(tmp_IS)[2] <- "IS"
    
    ImmuneResponse.filter_prot$IS <- data.frame(IS = tmp_IS[,"IS", drop = FALSE], check.names = FALSE)
    
    ### Save 
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }else{
    
    ImmuneResponse.filter_prot$IS <- NULL
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }
  

})


# ****************
# IMPRES 
# * Normalized (quantile-normalizatin) counts should be used to derive IMPRES 
# * When computing IMPRES for different cancer types, we find that the IMPRES genes 
# were not measured in a unified manner across all those, where for a few cancer types 
# they were 0 for a large percentage of the samples, which could have been a technical 
# issue for IMPRES. To maintain a uniform measure, we remove such genes for those cancer
# types and then normalize the score to 15 as in other datasets described in the main text. 

load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)
IMPRES.paper <- read.csv("./data/raw_data/IMPRES_TCGA.csv", row.names = 1, header = T)

sapply(PanCancer.names, function(Cancer){

  cat("\n",Cancer,"\n")

  file <- dir(paste0("~/Desktop/PhD_TUE/Github_model/desktop/data/raw_data/RNAseq/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
                     "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0"),
              recursive = TRUE, pattern = "data.txt", full.names = TRUE)

  # Extract the raw counts from the text files for each gene
  data.transcripts <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
  colnames(data.transcripts) <- substr(colnames(data.transcripts), 1, 15)

  get_normalized_counts <- which(data.transcripts[1,] == "normalized_count")
  normalized.transcripts <- data.frame(data.matrix(data.transcripts[-1,get_normalized_counts]))
  colnames(normalized.transcripts) <- gsub(".","-", colnames(normalized.transcripts), fixed = TRUE)
  sapply(normalized.transcripts, class) # numeric
  
  # Rows with non-valid HGNC symbols were removed.
  HGNC_symbol <- sapply(strsplit(rownames(normalized.transcripts),"\\|"),function(X) return(X[1]))
  normalized.transcripts <-  normalized.transcripts[-which(HGNC_symbol == "?"),]
  HGNC_symbol <-HGNC_symbol[-which(HGNC_symbol == "?")]
  HGNC_id <- sapply(strsplit(rownames(normalized.transcripts),"\\|"),function(X) return(X[2]))

  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(HGNC_symbol) != 0){
    idx <- which(duplicated(HGNC_symbol) == TRUE)
    dup_genes <- HGNC_symbol[idx]
    for (ii in dup_genes){
      normalized.transcripts[which(HGNC_symbol %in% ii)[1],] <- colMeans(normalized.transcripts[which(HGNC_symbol %in% ii),])
      normalized.transcripts <- normalized.transcripts[-which(HGNC_symbol %in% ii)[2],]
      HGNC_symbol <- HGNC_symbol[-which(HGNC_symbol %in% ii)[2]]
    }
  }

  rownames(normalized.transcripts) <- HGNC_symbol

  # Compute IMPRES
  IMPRES <- ObtainIMPRES(normalized.transcripts)

  # Sample screening:
  
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter

  load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  IMPRES.no_filter <- IMPRES[keep.samples.no_filter, , drop = FALSE]
  rownames(IMPRES.no_filter) <- substr(rownames(IMPRES.no_filter), 1, 12)
  tmp <- intersect( rownames(IMPRES.no_filter), substr(rownames(IMPRES.paper), 1, 12))
  #pairs(IMPRES.no_filter[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
  cat("\n","Correlation_no_filter =", cor(as.vector(IMPRES.no_filter[tmp,]), as.vector(IMPRES.paper[tmp,2]), method = "spearman"), "\n")

  ImmuneResponse.no_filter$IMPRES <- IMPRES.no_filter
  save(ImmuneResponse.no_filter, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){

    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))

    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    IMPRES.filter_Spat <- IMPRES[keep.samples.filter_Spat, , drop = FALSE]
    rownames(IMPRES.filter_Spat) <- substr(rownames(IMPRES.filter_Spat), 1, 12)
    tmp <- intersect( rownames(IMPRES.filter_Spat), substr(rownames(IMPRES.paper), 1, 12))
    #pairs(IMPRES.filter_Spat[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    cat("\n","Correlation_no_filter =", cor(as.vector(IMPRES.filter_Spat[tmp,]), as.vector(IMPRES.paper[tmp,2]), method = "spearman"), "\n")

    ImmuneResponse.filter_Spat$IMPRES <- IMPRES.filter_Spat
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }else{
    
    ImmuneResponse.filter_Spat$IMPRES <- NULL
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }
  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){

    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))

    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    IMPRES.filter_prot <- IMPRES[keep.samples.filter_prot, , drop = FALSE]
    rownames(IMPRES.filter_prot) <- substr(rownames(IMPRES.filter_prot), 1, 12)
    tmp <- intersect(rownames(IMPRES.filter_prot) , substr(rownames(IMPRES.paper), 1, 12))
    #pairs(IMPRES.filter_prot[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    cat("\n","Correlation_no_filter =", cor(as.vector(IMPRES.filter_prot[tmp,]), as.vector(IMPRES.paper[tmp,2]), method = "spearman"), "\n")

    ImmuneResponse.filter_prot$IMPRES <- IMPRES.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }else{
    
    ImmuneResponse.filter_prot$IMPRES <- NULL
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))

  }
})

# ****************
# IPS #

# data from francesca, to get same samples and order of genes
IPS.paper <- read.csv("./data/raw_data/patientsAll_IPS_literature.csv", row.names = 1, header = TRUE, sep = "\t")
IPS.paper <- IPS.paper[,"ips_ctla4_neg_pd1_neg", drop = FALSE]

sapply(PanCancer.names, function(Cancer){
  
  cat("\n",Cancer,"\n")
  
  file <- dir(paste0("~/Desktop/PhD_TUE/Github_model/desktop/data/raw_data/RNAseq/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
                     "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0"),
              recursive = TRUE, pattern = "data.txt", full.names = TRUE)
  
  # Extract the raw counts from the text files for each gene
  data.transcripts <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
  colnames(data.transcripts) <- substr(colnames(data.transcripts), 1, 15)
  
  # TPM/scaled_estimate values (transcriptomics data)
  get_estimates <- which(data.transcripts[1,] == "scaled_estimate")
  estimates.transcripts <- data.frame(data.matrix(data.transcripts[-1,get_estimates]))
  colnames(estimates.transcripts) <- gsub(".","-", colnames(estimates.transcripts), fixed = TRUE)
  sapply(estimates.transcripts, class) # numeric
  
  # Obtaining TPM (transcriptomics data)
  TPM.transcripts <- estimates.transcripts * 1e6
  # Obtaining log2(1 offset)
  log2mas1.TPM.transcripts <- log2(TPM.transcripts + 1)
  
  # Rows with non-valid HGNC symbols were removed.
  HGNC_symbol <- sapply(strsplit(rownames(log2mas1.TPM.transcripts),"\\|"),function(X) return(X[1]))
  log2mas1.TPM.transcripts <-  log2mas1.TPM.transcripts[-which(HGNC_symbol == "?"),]
  HGNC_symbol <-HGNC_symbol[-which(HGNC_symbol == "?")]
  HGNC_id <- sapply(strsplit(rownames(log2mas1.TPM.transcripts),"\\|"),function(X) return(X[2]))
  
  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(HGNC_symbol) != 0){
    idx <- which(duplicated(HGNC_symbol) == TRUE)
    dup_genes <- HGNC_symbol[idx]
    for (ii in dup_genes){
      log2mas1.TPM.transcripts[which(HGNC_symbol %in% ii)[1],] <- colMeans(log2mas1.TPM.transcripts[which(HGNC_symbol %in% ii),])
      log2mas1.TPM.transcripts <- log2mas1.TPM.transcripts[-which(HGNC_symbol %in% ii)[2],]
      HGNC_symbol <- HGNC_symbol[-which(HGNC_symbol %in% ii)[2]]
    }
  }
  
  rownames(log2mas1.TPM.transcripts) <- HGNC_symbol
  
  # Compute IPS
  IPS <- ObtainIPS(log2mas1.TPM.transcripts)

  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter
  
  load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]],1 ,15)
  IPS.no_filter <- IPS[keep.samples.no_filter, 1, drop = FALSE]
  rownames(IPS.no_filter) <- substr(rownames(IPS.no_filter), 1, 12)
  tmp <- intersect(rownames(IPS.no_filter), substr(rownames(IPS.paper), 1, 12))
  #pairs(IPS.no_filter[tmp,] ~ IPS.paper[tmp,], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8) 
  cat("\n","Correlation_no_filter =", cor(as.vector(IPS.no_filter[tmp,]), as.vector(IPS.paper[tmp,]), method = "spearman"), "\n")
  
  ImmuneResponse.no_filter$IPS <- IPS.no_filter
  save(ImmuneResponse.no_filter, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    IPS.filter_Spat <- IPS[keep.samples.filter_Spat, , drop = FALSE]
    rownames(IPS.filter_Spat) <- substr(rownames(IPS.filter_Spat), 1, 12)
    tmp <- intersect(rownames(IPS.filter_Spat), substr(rownames(IPS.paper), 1, 12))
    #pairs(IPS.filter_Spat[tmp,] ~ IPS.paper[tmp,], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8) 
    cat("\n","Correlation_filter_Spat =", cor(as.vector(IPS.filter_Spat[tmp,]), as.vector(IPS.paper[tmp,]), method = "spearman"), "\n")
    
    ImmuneResponse.filter_Spat$IPS <- IPS.filter_Spat
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }else{
    
    ImmuneResponse.filter_Spat$IPS <- NULL
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }
  ## Filter Prot
  # ------------------------------------------------------------------------------------------------------- #
  
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    IPS.filter_prot <- IPS[keep.samples.filter_prot, , drop = FALSE]
    rownames(IPS.filter_prot) <- substr(rownames(IPS.filter_prot), 1, 12)
    tmp <- intersect( rownames(IPS.filter_prot), substr(rownames(IPS.paper), 1, 12))
    
    #pairs(IPS.filter_prot[tmp,] ~ IPS.paper[tmp,], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8) 
    cat("\n","Correlation_filter_prot =", cor(as.vector(IPS.filter_prot[tmp,]), as.vector(IPS.paper[tmp,]), method = "spearman"), "\n")
    
    ImmuneResponse.filter_prot$IPS <- IPS.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }else{
    
    ImmuneResponse.filter_prot$IPS <- NULL
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }
})

# ****************
# Roh Immune Signature Score #
# data from francesca, to get same samples and order of genes

sapply(PanCancer.names, function(Cancer){
  
  cat("\n",Cancer,"\n")
  
  file <- dir(paste0("~/Desktop/PhD_TUE/Github_model/desktop/data/raw_data/RNAseq/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
                     "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0"),
              recursive = TRUE, pattern = "data.txt", full.names = TRUE)
  
  # Extract the raw counts from the text files for each gene
  data.transcripts <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
  colnames(data.transcripts) <- substr(colnames(data.transcripts), 1, 15)
  
  # TPM/scaled_estimate values (transcriptomics data)
  get_estimates <- which(data.transcripts[1,] == "scaled_estimate")
  estimates.transcripts <- data.frame(data.matrix(data.transcripts[-1,get_estimates]))
  colnames(estimates.transcripts) <- gsub(".","-", colnames(estimates.transcripts), fixed = TRUE)
  sapply(estimates.transcripts, class) # numeric
  
  # Obtaining TPM (transcriptomics data)
  TPM.transcripts <- estimates.transcripts * 1e6
  # Obtaining log2(1 offset)
  #log2mas1.TPM.transcripts <- log2(TPM.transcripts + 0.01)
  
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
  
  # Compute Roh ISS score
  Roh_ISS <- ObtainRohISS(TPM.transcripts)

  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter
  
  load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  Roh_ISS.no_filter <- Roh_ISS[keep.samples.no_filter, 1, drop = FALSE]
  rownames(Roh_ISS.no_filter) <- substr(rownames(Roh_ISS.no_filter), 1, 12)
  
  ImmuneResponse.no_filter$Roh_ISS <- Roh_ISS.no_filter
  save(ImmuneResponse.no_filter, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    Roh_ISS.filter_Spat <- Roh_ISS[keep.samples.filter_Spat, , drop = FALSE]
    rownames(Roh_ISS.filter_Spat) <- substr(rownames(Roh_ISS.filter_Spat), 1, 12)

    ImmuneResponse.filter_Spat$Roh_ISS <- Roh_ISS.filter_Spat
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))

  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    Roh_ISS.filter_prot <- Roh_ISS[keep.samples.filter_prot, , drop = FALSE]
    rownames(Roh_ISS.filter_prot) <- substr(rownames(Roh_ISS.filter_prot), 1, 12)
    
    ImmuneResponse.filter_prot$Roh_ISS <- Roh_ISS.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }
})

# ****************
# 12-Chemokine Signature  #

sapply(PanCancer.names, function(Cancer){
  
  cat("\n",Cancer,"\n")
  
  file <- dir(paste0("~/Desktop/PhD_TUE/Github_model/desktop/data/raw_data/RNAseq/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
                     "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0"),
              recursive = TRUE, pattern = "data.txt", full.names = TRUE)
  
  # Extract the raw counts from the text files for each gene
  data.transcripts <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
  colnames(data.transcripts) <- substr(colnames(data.transcripts), 1, 15)
  
  # TPM/scaled_estimate values (transcriptomics data)
  get_estimates <- which(data.transcripts[1,] == "scaled_estimate")
  estimates.transcripts <- data.frame(data.matrix(data.transcripts[-1,get_estimates]))
  colnames(estimates.transcripts) <- gsub(".","-", colnames(estimates.transcripts), fixed = TRUE)
  sapply(estimates.transcripts, class) # numeric
  
  # Obtaining TPM (transcriptomics data)
  TPM.transcripts <- estimates.transcripts * 1e6
  # Obtaining log2(1 offset)
  log2mas1.TPM.transcripts <- log2(TPM.transcripts + 1)
  
  # Rows with non-valid HGNC symbols were removed.
  HGNC_symbol <- sapply(strsplit(rownames(log2mas1.TPM.transcripts),"\\|"),function(X) return(X[1]))
  log2mas1.TPM.transcripts <-  log2mas1.TPM.transcripts[-which(HGNC_symbol == "?"),]
  HGNC_symbol <-HGNC_symbol[-which(HGNC_symbol == "?")]
  HGNC_id <- sapply(strsplit(rownames(log2mas1.TPM.transcripts),"\\|"),function(X) return(X[2]))
  
  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(HGNC_symbol) != 0){
    idx <- which(duplicated(HGNC_symbol) == TRUE)
    dup_genes <- HGNC_symbol[idx]
    for (ii in dup_genes){
      log2mas1.TPM.transcripts[which(HGNC_symbol %in% ii)[1],] <- colMeans(log2mas1.TPM.transcripts[which(HGNC_symbol %in% ii),])
      log2mas1.TPM.transcripts <- log2mas1.TPM.transcripts[-which(HGNC_symbol %in% ii)[2],]
      HGNC_symbol <- HGNC_symbol[-which(HGNC_symbol %in% ii)[2]]
    }
  }
  
  rownames(log2mas1.TPM.transcripts) <- HGNC_symbol
  
  # Compute 12-Chemokine Score
  Chemokine.score <- Obtain12chemokine(log2mas1.TPM.transcripts)

  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter
  
  load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  Messina.chemokine.no_filter <- Chemokine.score[keep.samples.no_filter, , drop = FALSE]
  rownames(Messina.chemokine.no_filter) <- substr(rownames(Messina.chemokine.no_filter), 1, 12)
  
  ImmuneResponse.no_filter$Messina.chemokine <- Messina.chemokine.no_filter
  save(ImmuneResponse.no_filter, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    Messina.chemokine.filter_Spat <- Chemokine.score[keep.samples.filter_Spat, , drop = FALSE]
    rownames(Messina.chemokine.filter_Spat) <- substr(rownames(Messina.chemokine.filter_Spat), 1, 12)
    
    ImmuneResponse.filter_Spat$Messina.chemokine <- Messina.chemokine.filter_Spat
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    Messina.chemokine.filter_prot <- Chemokine.score[keep.samples.filter_prot, , drop = FALSE]
    rownames(Messina.chemokine.filter_prot) <- substr(rownames(Messina.chemokine.filter_prot), 1, 12)
    
    ImmuneResponse.filter_prot$Messina.chemokine <- Messina.chemokine.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
  }
})

# ****************
# Immune Signature (Davoli) and proliferation signature #

load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)
Davoli.paper <- read.csv("./data/raw_data/IS_Proliferation_Davoli_TCGA.csv", header = T, skip = 5, row.names = 1)
IS_Davoli.paper <- Davoli.paper[,"Immune.Signature.Score", drop = FALSE]
IS_Davoli.paper <- IS_Davoli.paper[!is.na(IS_Davoli.paper[,1]), , drop = FALSE]
Proliferation_Davoli.paper <- Davoli.paper[,"CellCycle.Signature.Score", drop = FALSE]
Proliferation_Davoli.paper <- Proliferation_Davoli.paper[!is.na(Proliferation_Davoli.paper[,1]), , drop = FALSE]

sapply(PanCancer.names, function(Cancer){
  
  cat("\n",Cancer,"\n")
  
  file <- dir(paste0("~/Desktop/PhD_TUE/Github_model/desktop/data/raw_data/RNAseq/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
                     "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0"),
              recursive = TRUE, pattern = "data.txt", full.names = TRUE)
  
  # Extract the raw counts from the text files for each gene
  data.transcripts <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
  colnames(data.transcripts) <- substr(colnames(data.transcripts), 1, 15)
  
  # TPM/scaled_estimate values (transcriptomics data)
  get_estimates <- which(data.transcripts[1,] == "scaled_estimate")
  estimates.transcripts <- data.frame(data.matrix(data.transcripts[-1,get_estimates]))
  colnames(estimates.transcripts) <- gsub(".","-", colnames(estimates.transcripts), fixed = TRUE)
  sapply(estimates.transcripts, class) # numeric
  
  # Obtaining TPM (transcriptomics data)
  TPM.transcripts <- estimates.transcripts * 1e6
  # Obtaining log2(1 offset)
  log2mas1.TPM.transcripts <- log2(TPM.transcripts + 1)
  
  # Rows with non-valid HGNC symbols were removed.
  HGNC_symbol <- sapply(strsplit(rownames(log2mas1.TPM.transcripts),"\\|"),function(X) return(X[1]))
  log2mas1.TPM.transcripts <-  log2mas1.TPM.transcripts[-which(HGNC_symbol == "?"),]
  HGNC_symbol <-HGNC_symbol[-which(HGNC_symbol == "?")]
  HGNC_id <- sapply(strsplit(rownames(log2mas1.TPM.transcripts),"\\|"),function(X) return(X[2]))
  
  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(HGNC_symbol) != 0){
    idx <- which(duplicated(HGNC_symbol) == TRUE)
    dup_genes <- HGNC_symbol[idx]
    for (ii in dup_genes){
      log2mas1.TPM.transcripts[which(HGNC_symbol %in% ii)[1],] <- colMeans(log2mas1.TPM.transcripts[which(HGNC_symbol %in% ii),])
      log2mas1.TPM.transcripts <- log2mas1.TPM.transcripts[-which(HGNC_symbol %in% ii)[2],]
      HGNC_symbol <- HGNC_symbol[-which(HGNC_symbol %in% ii)[2]]
    }
  }
  
  rownames(log2mas1.TPM.transcripts) <- HGNC_symbol
  
  # Compute Immune Signature Davoli and proliferation signature
  IS_Davoli <- ObtainDavoliIS(log2mas1.TPM.transcripts)
  Proliferation <- ObtainProliferation(log2mas1.TPM.transcripts)

  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter
  
  load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
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
  
  save(ImmuneResponse.no_filter, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
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
    
    ImmuneResponse.filter_Spat$IS_Davoli <- IS_Davoli.filter_Spat
    ImmuneResponse.filter_Spat$Proliferation <- Proliferation.filter_Spat
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
  
  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
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
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }
  
})

# ****************
# IFny and Expanded Immune Signature (Ayers) #
# Input: normalized counts

load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

sapply(PanCancer.names, function(Cancer){
  
  cat("\n",Cancer,"\n")
  
  file <- dir(paste0("~/Desktop/PhD_TUE/Github_model/desktop/data/raw_data/RNAseq/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
                     "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0"),
              recursive = TRUE, pattern = "data.txt", full.names = TRUE)
  
  # Extract the raw counts from the text files for each gene
  data.transcripts <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
  colnames(data.transcripts) <- substr(colnames(data.transcripts), 1, 15)
  
  get_normalized_counts <- which(data.transcripts[1,] == "normalized_count")
  normalized.transcripts <- data.frame(data.matrix(data.transcripts[-1,get_normalized_counts]))
  colnames(normalized.transcripts) <- gsub(".","-", colnames(normalized.transcripts), fixed = TRUE)
  sapply(normalized.transcripts, class) # numeric
  
  # Rows with non-valid HGNC symbols were removed.
  HGNC_symbol <- sapply(strsplit(rownames(normalized.transcripts),"\\|"),function(X) return(X[1]))
  normalized.transcripts <-  normalized.transcripts[-which(HGNC_symbol == "?"),]
  HGNC_symbol <-HGNC_symbol[-which(HGNC_symbol == "?")]
  HGNC_id <- sapply(strsplit(rownames(normalized.transcripts),"\\|"),function(X) return(X[2]))
  
  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(HGNC_symbol) != 0){
    idx <- which(duplicated(HGNC_symbol) == TRUE)
    dup_genes <- HGNC_symbol[idx]
    for (ii in dup_genes){
      normalized.transcripts[which(HGNC_symbol %in% ii)[1],] <- colMeans(normalized.transcripts[which(HGNC_symbol %in% ii),])
      normalized.transcripts <- normalized.transcripts[-which(HGNC_symbol %in% ii)[2],]
      HGNC_symbol <- HGNC_symbol[-which(HGNC_symbol %in% ii)[2]]
    }
  }
  
  rownames(normalized.transcripts) <- HGNC_symbol
  
  # Compute IFN Ayers
  IFny_Ayers <- ObtainIFnyAyers(normalized.transcripts)
  # Compute Expanded Immune Signature Ayers
  ExpIS_Ayers <- ObtainExpandedImmuneAyers(normalized.transcripts)
  
  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter
  
  load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
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
  
  ImmuneResponse.no_filter$IFny_Ayers <- IFny_Ayers.no_filter
  ImmuneResponse.no_filter$ExpIS_Ayers <- ExpIS_Ayers.no_filter
  
  save(ImmuneResponse.no_filter, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
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
    
    ImmuneResponse.filter_Spat$IFny_Ayers <- IFny_Ayers.filter_Spat
    ImmuneResponse.filter_Spat$ExpIS_Ayers <- ExpIS_Ayers.filter_Spat
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
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
    
    ImmuneResponse.filter_prot$IFny_Ayers <- IFny_Ayers.filter_prot
    ImmuneResponse.filter_prot$ExpIS_Ayers <- ExpIS_Ayers.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }
  
})

# ****************
# T cell Inflamed GEP (Ayers) #
# Input: raw counts
load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

sapply(PanCancer.names, function(Cancer){
  
  cat("\n",Cancer,"\n")
  
  file <- dir(paste0("~/Desktop/PhD_TUE/Github_model/desktop/data/raw_data/RNAseq/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
                     "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0"),
              recursive = TRUE, pattern = "data.txt", full.names = TRUE)
  
  # Extract the raw counts from the text files for each gene
  data.transcripts <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
  colnames(data.transcripts) <- substr(colnames(data.transcripts), 1, 15)
  
  # raw values (transcriptomics data)
  get_rawcounts <- which(data.transcripts[1,] == "raw_count")
  rawcounts.transcripts <- data.frame(data.matrix(data.transcripts[-1,get_rawcounts]))
  colnames(rawcounts.transcripts) <- gsub(".","-", colnames(rawcounts.transcripts), fixed = TRUE)
  sapply(rawcounts.transcripts, class) # numeric
  
  # Rows with non-valid HGNC symbols were removed.
  HGNC_symbol <- sapply(strsplit(rownames(rawcounts.transcripts),"\\|"),function(X) return(X[1]))
  rawcounts.transcripts <-  rawcounts.transcripts[-which(HGNC_symbol == "?"),]
  HGNC_symbol <-HGNC_symbol[-which(HGNC_symbol == "?")]
  HGNC_id <- sapply(strsplit(rownames(rawcounts.transcripts),"\\|"),function(X) return(X[2]))
  
  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(HGNC_symbol) != 0){
    idx <- which(duplicated(HGNC_symbol) == TRUE)
    dup_genes <- HGNC_symbol[idx]
    for (ii in dup_genes){
      rawcounts.transcripts[which(HGNC_symbol %in% ii)[1],] <- colMeans(rawcounts.transcripts[which(HGNC_symbol %in% ii),])
      rawcounts.transcripts <- rawcounts.transcripts[-which(HGNC_symbol %in% ii)[2],]
      HGNC_symbol <- HGNC_symbol[-which(HGNC_symbol %in% ii)[2]]
    }
  }
  
  # Remove ImmuneResponse genes
  rownames(rawcounts.transcripts) <- HGNC_symbol
  
  # Compute T-cell Inflamed GEP (Ayers)
  T_cell_inflamed_Ayers <- ObtainTcellInflamedAyers(rawcounts.transcripts)

  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter
  
  load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  T_cell_inflamed_Ayers.no_filter <- T_cell_inflamed_Ayers[keep.samples.no_filter, , drop = FALSE]
  rownames(T_cell_inflamed_Ayers.no_filter) <- substr(rownames(T_cell_inflamed_Ayers.no_filter), 1, 12)
  
  ImmuneResponse.no_filter$T_cell_inflamed_Ayers <- T_cell_inflamed_Ayers.no_filter
  save(ImmuneResponse.no_filter, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    T_cell_inflamed_Ayers.filter_Spat <- T_cell_inflamed_Ayers[keep.samples.filter_Spat, , drop = FALSE]
    rownames(T_cell_inflamed_Ayers.filter_Spat) <- substr(rownames(T_cell_inflamed_Ayers.filter_Spat), 1, 12)
    
    ImmuneResponse.filter_Spat$T_cell_inflamed_Ayers <- T_cell_inflamed_Ayers.filter_Spat
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    T_cell_inflamed_Ayers.filter_prot <- T_cell_inflamed_Ayers[keep.samples.filter_prot, , drop = FALSE]
    rownames(T_cell_inflamed_Ayers.filter_prot) <- substr(rownames(T_cell_inflamed_Ayers.filter_prot), 1, 12)
    
    ImmuneResponse.filter_prot$T_cell_inflamed_Ayers <- T_cell_inflamed_Ayers.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
  }
})

# ****************
# TIDE (we need python package or linux function) #
# Input: log2(RPKM +1)

load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

sapply(PanCancer.names, function(Cancer){
  
  cat("\n",Cancer,"\n")
  
  file <- dir(paste0("~/Desktop/PhD_TUE/Github_model/desktop/data/raw_data/RNAseq/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
                     "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0"),
              recursive = TRUE, pattern = "data.txt", full.names = TRUE)
  
  # Extract the raw counts from the text files for each gene
  data.transcripts <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
  colnames(data.transcripts) <- substr(colnames(data.transcripts), 1, 15)
  
  # TPM/scaled_estimate values (transcriptomics data)
  get_estimates <- which(data.transcripts[1,] == "scaled_estimate")
  estimates.transcripts <- data.frame(data.matrix(data.transcripts[-1,get_estimates]))
  colnames(estimates.transcripts) <- gsub(".","-", colnames(estimates.transcripts), fixed = TRUE)
  sapply(estimates.transcripts, class) # numeric
  
  # Obtaining TPM (transcriptomics data)
  TPM.transcripts <- estimates.transcripts * 1e6
  # Obtaining log2(1 offset)
  log2mas1.TPM.transcripts <- log2(TPM.transcripts + 1)
  
  # Rows with non-valid HGNC symbols were removed.
  HGNC_symbol <- sapply(strsplit(rownames(log2mas1.TPM.transcripts),"\\|"),function(X) return(X[1]))
  log2mas1.TPM.transcripts <-  log2mas1.TPM.transcripts[-which(HGNC_symbol == "?"),]
  HGNC_symbol <-HGNC_symbol[-which(HGNC_symbol == "?")]
  HGNC_id <- sapply(strsplit(rownames(log2mas1.TPM.transcripts),"\\|"),function(X) return(X[2]))
  
  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(HGNC_symbol) != 0){
    idx <- which(duplicated(HGNC_symbol) == TRUE)
    dup_genes <- HGNC_symbol[idx]
    for (ii in dup_genes){
      log2mas1.TPM.transcripts[which(HGNC_symbol %in% ii)[1],] <- colMeans(log2mas1.TPM.transcripts[which(HGNC_symbol %in% ii),])
      log2mas1.TPM.transcripts <- log2mas1.TPM.transcripts[-which(HGNC_symbol %in% ii)[2],]
      HGNC_symbol <- HGNC_symbol[-which(HGNC_symbol %in% ii)[2]]
    }
  }
  
  rownames(log2mas1.TPM.transcripts) <- HGNC_symbol
  
  # 1. log2(TPM + 1): genes as rows
  log2.mas1.TPM <- log2mas1.TPM.transcripts
  #boxplot.matrix(TCGA.log2TPM1.pre)
  
  # 2. quantile normalization
  # Sort each column of X (values)
  sort.log2.mas1.TPM <- indices.sort <- matrix(0,nrow = nrow(log2.mas1.TPM), ncol = ncol(log2.mas1.TPM))
  
  for (X in 1:ncol(log2.mas1.TPM)){
    tmp = sort(log2.mas1.TPM[,X], decreasing = F, index.return = TRUE)
    sort.log2.mas1.TPM[,X] = tmp$x
    indices.sort[,X] = tmp$ix
  }
  # Take the means across rows of X sort
  # Assign this mean to each element in the row to get X' short 
  for (X in 1:nrow(log2.mas1.TPM)){
    sort.log2.mas1.TPM[X,] <- rowMeans(sort.log2.mas1.TPM)[X]
  }
  # Get X normalized by rearranging each column of X' sort to have th esame ordering as original X
  for (X in 1:ncol(log2.mas1.TPM)){
    tmp = sort(indices.sort[,X], decreasing = F, index.return = TRUE)
    sort.log2.mas1.TPM[,X] <- sort.log2.mas1.TPM[tmp$ix,X]
  }
  rownames(sort.log2.mas1.TPM) <- rownames(log2.mas1.TPM)
  colnames(sort.log2.mas1.TPM) <- colnames(log2.mas1.TPM)
  #boxplot.matrix(sort.TCGA.log2TPM1.pre)
  
  # 3. Substract gene average across all conditions
  average.gene <- rowMeans(sort.log2.mas1.TPM)
  sort.log2.mas1.TPM.norm <- sweep(sort.log2.mas1.TPM,1,average.gene, FUN = "-")
  
  write.table(sort.log2.mas1.TPM.norm, file = paste0("./data/data_processed_TIDE/log2mas1TPM_", Cancer,".txt"), sep = "\t")
  
  # Compute TIDE
  if (Cancer == "SKCM") {
    system(paste0("/home/olapuent/anaconda3/bin/tidepy ", "./data/data_processed_TIDE/log2mas1TPM_", Cancer,".txt", " -o output_TIDE_",Cancer,".txt -c Melanoma"))
  }else{
    system(paste0("/home/olapuent/anaconda3/bin/tidepy ", "./data/data_processed_TIDE/log2mas1TPM_", Cancer,".txt", " -o output_TIDE_",Cancer,".txt -c Other"))
  }
  
  TIDE.table <- read.table(file = paste0("output_TIDE_",Cancer,".txt"), sep = "\t", header = T, row.names = 1)
  TIDE.score <- data.frame(TIDE = TIDE.table[,"TIDE", drop = FALSE], check.names = F)
  
  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter
  
  load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  TIDE.no_filter <- TIDE.score[keep.samples.no_filter, , drop = FALSE]
  
  ImmuneResponse.no_filter$TIDE <- TIDE.no_filter
  save(ImmuneResponse.no_filter, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))

  # IFNG.score <- data.frame(IFNG = TIDE.table[,"IFNG", drop = FALSE], check.names = F)
  # IFNG.no_filter <- IFNG.score[keep.samples.no_filter, , drop = FALSE]
  # rownames(IFNG.no_filter) <- substr(rownames(IFNG.no_filter), 1, 12)
  # tmp <- intersect(rownames(IFNG.me), rownames(IFNG.no_filter))
  # 
  # 
  # pairs(IFNG.me[tmp,] ~ IFNG.no_filter[tmp,], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
  # cat("\n","Correlation =", cor(as.vector(IFNG.me[tmp,]), as.vector(IFNG.no_filter[tmp,]), method = "spearman"), "\n")
  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    TIDE.filter_Spat <- TIDE.score[keep.samples.filter_Spat, , drop = FALSE]
    
    ImmuneResponse.filter_Spat$TIDE <- TIDE.filter_Spat
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    TIDE.filter_prot <- TIDE.score[keep.samples.filter_prot, , drop = FALSE]
    
    ImmuneResponse.filter_prot$TIDE <- TIDE.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }
  
})

# ****************
# MSI score 
# Input: RSEM-normalized format and log2-transformed

load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

sapply(PanCancer.names, function(Cancer){
  
  cat("\n",Cancer,"\n")
  
  file <- dir(paste0("~/Desktop/PhD_TUE/Github_model/desktop/data/raw_data/RNAseq/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
                     "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0"),
              recursive = TRUE, pattern = "data.txt", full.names = TRUE)
  
  # Extract the raw counts from the text files for each gene
  data.transcripts <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
  colnames(data.transcripts) <- substr(colnames(data.transcripts), 1, 15)
  
  get_normalized_counts <- which(data.transcripts[1,] == "normalized_count")
  normalized.transcripts <- data.frame(data.matrix(data.transcripts[-1,get_normalized_counts]))
  colnames(normalized.transcripts) <- gsub(".","-", colnames(normalized.transcripts), fixed = TRUE)
  sapply(normalized.transcripts, class) # numeric
  
  # Rows with non-valid HGNC symbols were removed.
  HGNC_symbol <- sapply(strsplit(rownames(normalized.transcripts),"\\|"),function(X) return(X[1]))
  normalized.transcripts <-  normalized.transcripts[-which(HGNC_symbol == "?"),]
  HGNC_symbol <-HGNC_symbol[-which(HGNC_symbol == "?")]
  HGNC_id <- sapply(strsplit(rownames(normalized.transcripts),"\\|"),function(X) return(X[2]))
  
  # Rows corresponding to the same HGNC symbol were averaged.
  if(anyDuplicated(HGNC_symbol) != 0){
    idx <- which(duplicated(HGNC_symbol) == TRUE)
    dup_genes <- HGNC_symbol[idx]
    for (ii in dup_genes){
      normalized.transcripts[which(HGNC_symbol %in% ii)[1],] <- colMeans(normalized.transcripts[which(HGNC_symbol %in% ii),])
      normalized.transcripts <- normalized.transcripts[-which(HGNC_symbol %in% ii)[2],]
      HGNC_symbol <- HGNC_symbol[-which(HGNC_symbol %in% ii)[2]]
    }
  }
  
  rownames(normalized.transcripts) <- HGNC_symbol
  log2mas1.normalized <- log2(normalized.transcripts + 1)
  # Compute MSI
  MSI <- ObtainMSI(log2mas1.normalized)
  
  # Sample screening:
  
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter
  
  load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  MSI.no_filter <- MSI[keep.samples.no_filter, , drop = FALSE]
  rownames(MSI.no_filter) <- substr(rownames(MSI.no_filter), 1, 12)
  #tmp <- intersect( rownames(MSI.no_filter), substr(rownames(IMPRES.paper), 1, 12))
  #pairs(MSI.no_filter[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
  #cat("\n","Correlation_no_filter =", cor(as.vector(MSI.no_filter[tmp,]), as.vector(IMPRES.paper[tmp,2]), method = "spearman"), "\n")
  
  ImmuneResponse.no_filter$MSI <- MSI.no_filter
  save(ImmuneResponse.no_filter, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
    MSI.filter_Spat <- MSI[keep.samples.filter_Spat, , drop = FALSE]
    rownames(MSI.filter_Spat) <- substr(rownames(MSI.filter_Spat), 1, 12)
    #tmp <- intersect( rownames(MSI.filter_Spat), substr(rownames(IMPRES.paper), 1, 12))
    #pairs(MSI.filter_Spat[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    #cat("\n","Correlation_no_filter =", cor(as.vector(MSI.filter_Spat[tmp,]), as.vector(IMPRES.paper[tmp,2]), method = "spearman"), "\n")
    
    ImmuneResponse.filter_Spat$MSI <- MSI.filter_Spat
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }
  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
    MSI.filter_prot <- MSI[keep.samples.filter_prot, , drop = FALSE]
    rownames(MSI.filter_prot) <- substr(rownames(MSI.filter_prot), 1, 12)
    #tmp <- intersect(rownames(MSI.filter_prot) , substr(rownames(IMPRES.paper), 1, 12))
    #pairs(MSI.filter_prot[tmp,] ~ IMPRES.paper[tmp,2], upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8)
    #cat("\n","Correlation_no_filter =", cor(as.vector(MSI.filter_prot[tmp,]), as.vector(IMPRES.paper[tmp,2]), method = "spearman"), "\n")
    
    ImmuneResponse.filter_prot$MSI <- MSI.filter_prot
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))

  }
})

# ****************

# ***************
# Re-editing some files

# Change TCGA id of TIDE in ImmuneResponse. All samples are primary tumors, so we keep XXXX-XX-XXXX #

load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

sapply(PanCancer.names, function(Cancer){
  
  cat("\n",Cancer,"\n")
  
  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter
  
  load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  rownames(ImmuneResponse.no_filter$TIDE) <- substr(rownames(ImmuneResponse.no_filter$TIDE), 1, 12)
  
  save(ImmuneResponse.no_filter, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))
  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    rownames(ImmuneResponse.filter_Spat$TIDE) <- substr(rownames(ImmuneResponse.filter_Spat$TIDE), 1, 12)

    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    rownames(ImmuneResponse.filter_prot$TIDE) <- substr(rownames(ImmuneResponse.filter_prot$TIDE), 1, 12)

    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
  }
  
})

# Remove sTIL attribute from quanTIseq_IS and quanTIseq_IS_prot #
# Make pathways a data.frame #
# Remove pathways attribute (twice) from DataViews_filter_prot, as well as immune cells #

# ------------------------------------------------------------------------------------------------------- #
# Match samples in cell fractions data for general analysis: screening with quanTIseq_IS
# ------------------------------------------------------------------------------------------------------- #


# Cancer types
load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

sapply(PanCancer.names, function(ii){

  load(paste0("./data/PanCancer/",ii,"/new/DataViews_no_filter_",ii, ".RData"))
  
  DataViews.no_filter$sTIL <- NULL
  DataViews.no_filter$pathways <- as.data.frame(DataViews.no_filter$pathways)
  
  #save(DataViews.no_filter, file = paste0("./data/PanCancer/",ii,"/new/DataViews_no_filter_",ii, ".RData"))
  
  
})

# ------------------------------------------------------------------------------------------------------- #
# Match samples in cell fractions data for proteomics analysis: screening with quanTIseq_IS_prot
# ------------------------------------------------------------------------------------------------------- #


# Cancer types
load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_prot.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)

sapply(PanCancer.names, function(ii){
  
  load(paste0("./data/PanCancer/",ii,"/new/DataViews_filter_prot_",ii, ".RData"))

  DataViews.filter_prot$immunecells <- NULL
  DataViews.filter_prot$Pathways <- NULL
  DataViews.filter_prot$pathways <- as.data.frame(DataViews.filter_prot$pathways)
  
  #save(DataViews.filter_prot, file = paste0("./data/PanCancer/",ii,"/new/DataViews_filter_prot_",ii, ".RData"))
  
})


# ImmuneResponse should be a matrix, not a list. 

load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

sapply(PanCancer.names, function(Cancer){
  
  cat("\n",Cancer,"\n")
  
  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter
  
  load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, ".RData"))

  ImmuneResponse.no_filter <- do.call(cbind, lapply(names(ImmuneResponse.no_filter), function(proxy){
    
    tmp.proxy <- ImmuneResponse.no_filter[[proxy]]
    
    return(tmp.proxy)
    
  }))
  
  save(ImmuneResponse.no_filter, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_no_filter_",Cancer, "_matrix_format.RData"))
  
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, ".RData"))
    
    ImmuneResponse.filter_Spat <- do.call(cbind, lapply(names(ImmuneResponse.filter_Spat), function(proxy){
      
      tmp.proxy <- ImmuneResponse.filter_Spat[[proxy]]
      
      return(tmp.proxy)
      
    }))
    
    save(ImmuneResponse.filter_Spat, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_",Cancer, "_matrix_format.RData"))
    
  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, ".RData"))
    
    ImmuneResponse.filter_prot <- do.call(cbind, lapply(names(ImmuneResponse.filter_prot), function(proxy){
      
      tmp.proxy <- ImmuneResponse.filter_prot[[proxy]]
      
      return(tmp.proxy)
      
    }))
    
    save(ImmuneResponse.filter_prot, file = paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_prot_",Cancer, "_matrix_format.RData"))
    
  }
  
})










#######################################
# # Overview samples size #
# cancer.types <- names(Proteomics.pan.cancer)
# 
# patient.we.lose = do.call(rbind,lapply(cancer.types, function(X){
#  sum(is.na(rowSums(Proteomics.pan.cancer[[X]])))/nrow(Proteomics.pan.cancer[[X]])
# }))
# rownames(patient.we.lose) <- cancer.types
# 
# # % patients we will lose:
# # BLCA 0.62948207 YES
# # CESC 0.03968254 NO
# # CRC  0.16417910 YES
# # LUAD 0.40438871 YES
# # LUSC 0.43495935 YES
# # PRAD 0.40789474 YES 
# # STAD 0.39271255 YES 
# # UCEC 0.00000000 NO
# # SKCM 0.36000000 YES 
# # BRCA 0.64689655 YES
# 
# 
# # Missing values imputation #
# # We will be using pcaMethods library in r
# # if (!requireNamespace("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# # BiocManager::install("pcaMethods", version = "3.8")
# 
# # Sparse RPPA measurements were excluded, requiring that each protein is measured in at least 100 cell lines, 
# # which resulted in the inclusion of 263 protein measurements. 
# 
# #library(pcaMethods)
# par(mfrow = c(2, 5))
# #Proteomics.pan.cancer.imputated <- do.call(list,lapply(cancer.types[1:9], function(X){
# Proteomics.pan.cancer <- lapply(cancer.types, function(X){
#   matrix = Proteomics.pan.cancer[[X]]
#   # 1.Each protein is measured in at least 12.5 % of the patients
#   threshold = round(0.125*nrow(matrix),2)
#   keep_proteins = which(apply(is.na(matrix),2,sum) <= threshold)
#   remove_proteins = which(apply(is.na(matrix),2,sum) > threshold)
#   
#   # 2.Check which patients have a significant number of NA values:
#   remove_patients = apply(is.na(matrix_a),1,sum) 
#   # Some patients have 6 missing values. Why would you delete cases with missing values? They are likely to be different 
#   # to cases with complete values, so you are building bias into your sample by so doing. 
#   
#   # 3. Check how the matrix looks like
#   matrix_a = matrix[,keep_proteins]
#   
#   # 4. Plot intensity distributions with and without missing values
#   # To check whether missing values are biased to lower intense proteins, the densities  are plotted for proteins with and without missing values.
#   aux = which(apply(is.na(matrix_a),2,sum) != 0)
#   if (length(aux) > 0){
#     a = matrix_a[,-aux]
#     write.table(names(aux), file = paste0("data/Pipeline/PanCancer/",X,"/Proteins_with_NA_values_in_",X,".txt"))
#     #heatmap(as.matrix(matrix_a), Rowv = NA, Colv = NA, scale = "none", col = terrain.colors(256))
#     #heatmap(as.matrix(matrix_a), scale = "none", col = terrain.colors(256))
#   }else{
#     a = matrix_a
#     write.table("No proteins with NA values", file = paste0("data/Pipeline/PanCancer/",X,"/Proteins_with_NA_values_in_",X,".txt"))
#   }
#   b = as.vector(matrix_a)
#   b[is.na(b) == T] = 0
#   plot(density(unlist(as.vector(a))), col = "red", main = paste0("Density function:",X), lty = 2, lwd = 3)
#   lines(density(unlist(b)), col = "blue", lwd = 3)
#   text(max(matrix, na.rm = T)-4,1,labels = paste0("Proteins with NA values = ", length(aux)))
#   
#   return(matrix_a)
# })
# names(Proteomics.pan.cancer) <- cancer.types

# In our data, there is no clear difference between the two distributions.
# Thanks to the heatmap, we are able to see that certain proteins might have missing values at random (MAR).

# # matrix – Pre-processed (centered, scaled) data with variables in columns and observations in rows.
# # The data may contain missing values, denoted as NA.
# tmp <- standarization(Proteomics.pan.cancer[[X]])
# ## Perform svdImpute using the 3 largest components
# result <- pca(t(tmp), method="svdImpute", nPcs=3)
# ## Get the estimated complete observations
# cObs <- as.data.frame(completeObs(result))

###########################################################################
# Adding inter-cellular networks to DataViews
###########################################################################

load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

sapply(PanCancer.names, function(ii){
  
  cat("\n",ii,"\n")
  
  ## Ligand-receptors pairs ##
  L_R_pairs_log10 <- read.table(paste0("../data/Inter-cellular_networks_Maisa/weighted_networks_log10//",
                          ii,".txt"), header = T, sep = "\t", row.names = 2, check.names = F)
  
  L_R_pairs_log2 <- read.table(paste0("../data/Inter-cellular_networks_Maisa/weighted_networks_log2/",
                                 ii,".txt"), header = T, sep = "\t", row.names = 2, check.names = F)
  
  # Remove two string columns
  L_R_pairs_log10 <- as.data.frame(t(L_R_pairs_log10[,-c(1,2,3)]))
  L_R_pairs_log2 <- as.data.frame(t(L_R_pairs_log2[,-c(1,2,3)]))
  
  L_R_pairs_log10[1:3,1:3]
  L_R_pairs_log2[1:3,1:3]
  
 
  L_R_pairs## Cytokine pairs ##
  Cytokine_pairs <- read.table(paste0("~/Desktop/PhD_TUE/Github_model/desktop/data/Inter-cellular_networks_Maisa/weighted_networks/",
                               ii,"_cytokine.txt"), header = T, sep = "\t", row.names = 2, check.names = F)
  # Remove two string columns
  Cytokine_pairs <- as.data.frame(t(Cytokine_pairs[,-c(1,2,3)]))

  # Sample screening:
  # ------------------------------------------------------------------------------------------------------- #
  ## No filter
  
  load(paste0("./data/PanCancer/",ii,"/new/DataViews_no_filter_",ii, ".RData"))
  
  keep.samples_no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[ii]], 1, 15)
  LRpairs.no_filter <- L_R_pairs[keep.samples_no_filter, , drop = FALSE]
  Cytokinepairs.no_filter <- Cytokine_pairs[keep.samples_no_filter, , drop = FALSE]
  
  # Integrate within DataViews
  DataViews.no_filter$LRpairs <- LRpairs.no_filter
  DataViews.no_filter$CYTOKINEpairs <- Cytokinepairs.no_filter
  
  save(DataViews.no_filter, file = paste0("./data/PanCancer/",ii,"/new/DataViews_no_filter_",ii, ".RData"))

  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Spat
  if(ii %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("./data/PanCancer/",ii,"/new/DataViews_filter_Spat_",ii, ".RData"))
    
    keep.samples.filter_Spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[ii]], 1, 15)
    LRpairs.filter_Spat <- L_R_pairs[keep.samples.filter_Spat, , drop = FALSE]
    Cytokinepairs.filter_Spat <- Cytokine_pairs[keep.samples.filter_Spat, , drop = FALSE]
    
    # Integrate within DataViews
    DataViews.filter_Spat$LRpairs <- LRpairs.filter_Spat
    DataViews.filter_Spat$CYTOKINEpairs <- Cytokinepairs.filter_Spat
    
    save(DataViews.filter_Spat, file = paste0("./data/PanCancer/",ii,"/new/DataViews_filter_Spat_",ii, ".RData"))

  }
  # ------------------------------------------------------------------------------------------------------- #
  ## Filter Prot
  if(ii %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)){
    
    load(paste0("./data/PanCancer/",ii,"/new/DataViews_filter_prot_",ii, ".RData"))
    
    keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[ii]], 1, 15)
    LRpairs.filter_prot <- L_R_pairs[keep.samples.filter_prot, , drop = FALSE]
    Cytokinepairs.filter_prot <- Cytokine_pairs[keep.samples.filter_prot, , drop = FALSE]
    
    # Integrate within DataViews
    DataViews.filter_prot$LRpairs <- LRpairs.filter_prot
    DataViews.filter_prot$CYTOKINEpairs <- Cytokinepairs.filter_prot
    
    save(DataViews.filter_prot, file = paste0("./data/PanCancer/",ii,"/new/DataViews_filter_prot_",ii, ".RData"))

  }
  
})
a <- as.matrix(L_R_pairs_log10)
b <- as.matrix(L_R_pairs_log2)
c <- as.matrix(DataViews.no_filter$LRpairs)

LR_sum <- apply(c,2, sum)
remove_NA_LR_pairs <- as.numeric(na.action(na.omit(LR_sum)))
c <- c[,-remove_NA_LR_pairs]

LR_sum <- apply(a,2, sum)
remove_NA_LR_pairs <- as.numeric(na.action(na.omit(LR_sum)))
a <- a[,-remove_NA_LR_pairs]

LR_sum <- apply(b,2, sum)
remove_NA_LR_pairs <- as.numeric(na.action(na.omit(LR_sum)))
b <- b[,-remove_NA_LR_pairs]

data <- data.frame(log10 = as.vector(a), log2 = as.vector(b), compute = as.vector(c))
pairs(~., data = data, lower.panel = panel.cor, upper.panal = panel.cor)



"A2M", "LRP1"

RNA.tpm
gene_expr["A2M","TCGA-02-0047-01"]
gene_expr["LRP1","TCGA-02-0047-01"]




    