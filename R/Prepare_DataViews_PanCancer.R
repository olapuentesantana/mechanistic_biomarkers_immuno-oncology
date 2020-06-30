#########################################################################################################
# Script to generate DataViews {with (11 cancer types) and without sTIL (18 cancer types)}:

# Input data -->
## Mechanistic DataViews:
### Pathways, Immunecells, TFs, sTIL, LRpairs, CYTOKINEpairs

## RNA_PROT DataViews:
### transcript, Protall

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
# BiocManager::install("DESeq2")
#library(DESeq2)
# BiocManager::install("progeny")
#library(progeny)
# BiocManager::install("dorothea")
#library(dorothea)

# *****************
# functions
source("scaling_function.R")
# Compute pathway activity
source("compute.pathways.scores.R")
# Compute TFs activity
source("compute.TF.activity.R")
# Compute LR pairs
source("compute.LR.pairs.R")
# General pre-processing tcga gene names
source("preprocess_tcga_genes_names.R")
# Pre-processing tcga gene names to match pathway genes
source("preprocess_gene_names_to_match_pathway_genes.R")
# Pre-processing tcga gene names to match pathway genes
source("preprocess_gene_names_to_match_tf_target_genes.R")
# Pre-processing tcga gene names to match LR pairs genes
source("preprocess_gene_names_to_match_LRpairs_genes.R")

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

# *****************
# Initialize DataViews
DataViews.no_filter <- vector("list", length = 6)
names(DataViews.no_filter) <- c("Pathways", "ImmuneCells", "TFs", "LRpairs", "CYTOKINEpairs","Transcript")
DataViews.filter_prot <- vector("list", length = 4)
names(DataViews.filter_prot) <- c("Pathways", "Protall", "TFs", "Transcript")
DataViews.filter_spat <- vector("list", length = 4)
names(DataViews.filter_spat) <- c("ImmuneCells", "sTIL", "LRpairs", "CYTOKINEpairs")

# ***************
# Remove transcripts used to build ImmuneResponse (IS,CYT,IPS,IMPRES,RohISS,Chemokine,Proliferation,IS_Davoli,IFNy,ExpandedImmune,
# T_cell_inflamed,TIDE,MSI)

# Genes to remove according to all ICB proxy's 
load("../data/all_genes_ICB_proxies.RData")

# ***************
# RNAseq data --> Pahways, TFs

sapply(PanCancer.names[5:length(PanCancer.names)], function(Cancer){
  
  message("Cancer type: ", Cancer, "\n")
  
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
  # # Log2 transformed
  # log2mas1.TPM.transcripts <- log2(TPM.transcripts + 1)
  # 
  # # Correct gene names in TCGA data
  # log2mas1.TPM.transcripts_process_v1 <- preprocess_tcga_genes_names(log2mas1.TPM.transcripts)
  # 
  # remove.genes <- na.exclude(match(all_genes_to_remove, rownames(log2mas1.TPM.transcripts_process_v1)))
  # log2mas1.TPM.transcripts_process_v1 <- log2mas1.TPM.transcripts_process_v1[-remove.genes,]
  # message("Removing signatures genes for proxy's of ICB response:  \n")
  # 
  # Sample screening:
  keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
  keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
  keep.samples.filter_spat <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)
  # 
  # transcript.no_filter <- log2mas1.TPM.transcripts_process_v1[, keep.samples.no_filter]
  # transcript.filter_prot <- log2mas1.TPM.transcripts_process_v1[, keep.samples.filter_prot]
  # 
  # DataViews.no_filter$Transcript <- as.data.frame(t(transcript.no_filter))
  # DataViews.filter_prot$Transcript <- as.data.frame(t(transcript.filter_prot))
  
  # ---------------------------------------------------------------------------------- #
  ## TF activity computation
  # ---------------------------------------------------------------------------------- #
  # 
  # # Correct gene names in TCGA data
  # TPM.transcripts_process_v1 <- preprocess_tcga_genes_names(TPM.transcripts)
  # 
  # # Correct gene names in order to match correctly tf target genes
  # TPM.transcripts_process_v2 <- preprocess_gene_names_to_match_tf_target_genes(TPM.transcripts_process_v1)
  # 
  # # Sample screening:
  # tpm.no_filter <- TPM.transcripts_process_v2[, keep.samples.no_filter]
  # tpm.filter_prot <- TPM.transcripts_process_v2[, keep.samples.filter_prot]
  # tpm.filter_spat <- TPM.transcripts_process_v2[, keep.samples.filter_spat]
  # 
  # # Computation of TF activity (input matrix [genes, samples], ouput matrix [sample, TFs])
  # TF_activity.no_filter <- compute.TF.activity(RNA.tpm = tpm.no_filter, remove.genes.ICB_proxies=TRUE)
  # TF_activity.filter_prot <- compute.TF.activity(RNA.tpm = tpm.filter_prot, remove.genes.ICB_proxies=TRUE)
  # 
  # # Insert into DataViews
  # DataViews.no_filter$TFs <- as.data.frame(TF_activity.no_filter$scores)
  # DataViews.filter_prot$TFs <- as.data.frame(TF_activity.filter_prot$scores)
  
  # ---------------------------------------------------------------------------------- #
  ## Ligand-Receptor and cytokine pairs computation
  # ---------------------------------------------------------------------------------- #
  
  load(paste0("../data/PanCancer_draft_v1/",Cancer,"/DataViews_no_filter_",Cancer, ".RData"))

  # Correct gene names in TCGA data
  TPM.transcripts_process_v1 <- preprocess_tcga_genes_names(TPM.transcripts)
  
  # Correct gene names in order to match correctly pathway responsive genes
  TPM.transcripts_process_v2 <- preprocess_gene_names_to_match_LRpairs_genes(TPM.transcripts_process_v1)
  
  # Sample screening:
  tpm.no_filter <- TPM.transcripts_process_v2[, keep.samples.no_filter]

  # Computation of LR pairs (input matrix [genes, samples], ouput matrix [sample, LR])
  LRpairs.no_filter <- compute.LR.pairs(RNA.tpm = tpm.no_filter, remove.genes.ICB_proxies=TRUE, 
                                        compute.cytokines.pairs=TRUE)
  
  # Insert into DataViews
  DataViews.no_filter$LRpairs <- as.data.frame(LRpairs.no_filter$LRpairs)
  DataViews.no_filter$CYTOKINEpairs <- as.data.frame(LRpairs.no_filter$CYTOKINEpairs)
  
  ## Filter Spat
  if(Cancer %in% names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)){
    
    load(paste0("../data/PanCancer_draft_v1/",Cancer,"/DataViews_filter_spat_",Cancer, ".RData"))
    
    # Sample screening:
    tpm.filter_spat <- TPM.transcripts_process_v2[, keep.samples.filter_spat]
    
    # Computation of LR pairs (input matrix [genes, samples], ouput matrix [sample, LR])
    LRpairs.filter_spat <- compute.LR.pairs(RNA.tpm = tpm.filter_spat, remove.genes.ICB_proxies=TRUE, 
                                            compute.cytokines.pairs=TRUE)
    # Insert into DataViews
    DataViews.filter_spat$LRpairs <- as.data.frame(LRpairs.filter_spat$LRpairs)
    DataViews.filter_spat$CYTOKINEpairs <- as.data.frame(LRpairs.filter_spat$CYTOKINEpairs)
    save(DataViews.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/DataViews_filter_spat_",Cancer, ".RData"))
    
  }
  
  # ---------------------------------------------------------------------------------- #
  ## Pathways activity computation
  # ---------------------------------------------------------------------------------- #
  # Obtaining raw counts 
  # get_rawcounts <- which(data.transcripts[1,] == "raw_count")
  # rawcounts.transcripts <- data.frame(data.transcripts[-1,get_rawcounts])
  # colnames(rawcounts.transcripts) <- gsub(".","-", colnames(rawcounts.transcripts), fixed = TRUE)
  # genes <- rownames(rawcounts.transcripts)
  # rawcounts.transcripts <- sapply(rawcounts.transcripts, as.numeric) # numeric
  # sapply(rawcounts.transcripts, class) # numeric
  # rownames(rawcounts.transcripts) <- genes
  # 
  # # Correct gene names in TCGA data
  # rawcounts.transcripts_process_v1 <- preprocess_tcga_genes_names(rawcounts.transcripts)
  # 
  # # Correct gene names in order to match correctly pathway responsive genes
  # rawcounts.transcripts_processed_v2 <- preprocess_gene_names_to_match_pathways_genes(rawcounts.transcripts_process_v1)
  # 
  # # Sample screening:
  # rawcounts.no_filter <- rawcounts.transcripts_processed_v2[, keep.samples.no_filter]
  # rawcounts.filter_prot <- rawcounts.transcripts_processed_v2[, keep.samples.filter_prot]
  # 
  # # Computation of Pathway scores (input matrix [genes, samples], ouput matrix [sample, pathways])
  # Pathways.activity.no_filter <- compute.pathways.scores(RNA.raw_counts=rawcounts.no_filter, 
  #                                                        remove.genes.ICB_proxies=TRUE)
  # Pathways.activity.filter_prot <- compute.pathways.scores(RNA.raw_counts=rawcounts.filter_prot, 
  #                                                          remove.genes.ICB_proxies=TRUE)
  # # Insert into DataViews
  # DataViews.no_filter$Pathways <- as.data.frame(Pathways.activity.no_filter$scores)
  # DataViews.filter_prot$Pathways <- as.data.frame(Pathways.activity.filter_prot$scores)
  
  save(DataViews.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/DataViews_no_filter_",Cancer, ".RData"))
  # save(DataViews.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/DataViews_filter_prot_",Cancer, ".RData"))
  # save(DataViews.filter_spat, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/DataViews_filter_spat_",Cancer, ".RData"))
  
})


## remove.genes.ICB_proxies=FALSE
# sapply(PanCancer.names.some, function(Cancer){
#   
#   file <- dir(paste0("../data/raw_data_tcga/RNAseq/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
#                      "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0"),
#               recursive = TRUE, pattern = "data.txt", full.names = TRUE)
#   
#   # Extract the raw counts from the text files for each gene
#   data.transcripts <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
#   colnames(data.transcripts) <- substr(colnames(data.transcripts), 1, 15)
#   
#   # TPM/scaled_estimate values (transcriptomics data)
#   get_estimates <- which(data.transcripts[1,] == "scaled_estimate")
#   estimates.transcripts <- data.frame(data.transcripts[-1,get_estimates])
#   colnames(estimates.transcripts) <- gsub(".","-", colnames(estimates.transcripts), fixed = TRUE)
#   genes <- rownames(estimates.transcripts)
#   estimates.transcripts <- sapply(estimates.transcripts, as.numeric) # numeric
#   
#   # Obtaining TPM (transcriptomics data)
#   TPM.transcripts <- estimates.transcripts * 1e6
#   rownames(TPM.transcripts) <- genes
#   
#   # Sample screening:
#   keep.samples.no_filter <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)
#   keep.samples.filter_prot <- substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[Cancer]], 1, 15)
#   
#   load(paste0("../data/PanCancer/",Cancer,"/PanCancer_draft_v1/DataViews_no_filter_",Cancer, ".RData"))
#   load(paste0("../data/PanCancer/",Cancer,"/PanCancer_draft_v1/DataViews_filter_prot_",Cancer, ".RData"))
#   
#   # Remove ImmuneResponse genes (the function should take care of it)
#   
#   # Sample screening:
#   tpm.no_filter <- TPM.transcripts[, keep.samples.no_filter]
#   tpm.filter_prot <- TPM.transcripts[, keep.samples.filter_prot]
#   
  # ---------------------------------------------------------------------------------- #
  ## Ligand-Receptor and cytokine pairs computation
  # ---------------------------------------------------------------------------------- #
#   
#   # Computation of LR pairs activity (input matrix [genes, samples], ouput matrix [sample, LR])
#   LRpairs.no_filter <- compute.LR.pairs(RNA.tpm = tpm.no_filter, remove.genes.ICB_proxies=FALSE, compute.cytokines.pairs=TRUE)
#   LRpairs.filter_prot <- compute.LR.pairs(RNA.tpm = tpm.filter_prot, remove.genes.ICB_proxies=FALSE, compute.cytokines.pairs=TRUE)
#   
#   # Insert into DataViews
#   DataViews.no_filter$LRpairs <- as.data.frame(LRpairs.no_filter$LRpairs)
#   DataViews.no_filter$CYTOKINEpairs <- as.data.frame(LRpairs.no_filter$CYTOKINEpairs)
#   
#   DataViews.filter_prot$LRpairs <- as.data.frame(LRpairs.filter_prot$LRpairs)
#   DataViews.filter_prot$CYTOKINEpairs <- as.data.frame(LRpairs.filter_prot$CYTOKINEpairs)
#   
#   save(DataViews.no_filter, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/DataViews_no_filter_",Cancer, ".RData"))
#   save(DataViews.filter_prot, file = paste0("../data/PanCancer_draft_v1/",Cancer,"/DataViews_filter_prot_",Cancer, ".RData"))
#   
# })

# ***************
# Protein abundance (RPPA) data --> Proteomics

all.prot.tcpa <- read.csv("../data/raw_data_tcga/TCGA-PANCAN32-L4.csv", row.names = 1 , header = T, stringsAsFactors = F, check.names = T)

# Cancer types
load("../analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_prot.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)

sapply(PanCancer.names, function(ii){
  
  load(paste0("../data/PanCancer_draft_v1/",ii,"/DataViews_filter_prot_",ii, ".RData"))
  
  keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot[[ii]], 1, 15), substr(rownames(all.prot.tcpa), 1, 15))
  tmp_prot <- all.prot.tcpa[keep,]
  tmp_prot <- tmp_prot[, 3:ncol(tmp_prot)]
  rownames(tmp_prot) <- substr(rownames(tmp_prot), 1, 15)
  
  DataViews.filter_prot$Protall <- as.data.frame(tmp_prot)
  
  save(DataViews.filter_prot, file = paste0("../data/PanCancer_draft_v1/",ii,"/DataViews_filter_prot_",ii, ".RData"))
  
})

# ***************
# Immune cells data (quanTIseq) and Spatial TILs data (Saltz et al.)

# data
all_cell.fractions <- read.csv("../data/raw_data_tcga/quanTIseq_estimated.csv", header = TRUE, row.names = 1)
all.spatial.TILs <- read.csv("../data/raw_data_tcga/spatial_TILs_saltz.csv", header = TRUE, row.names = 1)


# ------------------------------------------------------------------------------------------------------- #
# Match samples in cell fractions data for general analysis: screening with quanTIseq_IS
# ------------------------------------------------------------------------------------------------------- #

# Cancer types
load("../analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

for (ii in PanCancer.names){

  load(paste0("../data/PanCancer_draft_v1/",ii,"/DataViews_no_filter_",ii, ".RData"))

  keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[ii]], 1, 15),
                rownames(all_cell.fractions))
  tmp_immunecells <- all_cell.fractions[keep,]

  if (nrow(tmp_immunecells) == nrow(DataViews.no_filter$Pathways)){

    tmp_immunecells$T_cells_CD4 <- tmp_immunecells$T_cells_CD4 + tmp_immunecells$T_cells_regulatory_Tregs
    DataViews.no_filter$ImmuneCells <- tmp_immunecells
  }else{

    stop("pathways sample size  != immune cells sample size")

    break
  }
  save(DataViews.no_filter, file = paste0("../data/PanCancer_draft_v1/",ii,"/DataViews_no_filter_",ii, ".RData"))
}

# ------------------------------------------------------------------------------------------------------- #
# Match samples in cell fractions data for spatial information analysis: screening with quanTIseq_IS_Spat
# ------------------------------------------------------------------------------------------------------- #
# ***************
# Immune cells data (quanTIseq) and Spatial TILs data (Saltz et al.)

# data
all_cell.fractions <- read.csv("../data/raw_data_tcga/quanTIseq_estimated.csv", header = TRUE, row.names = 1)
all.spatial.TILs <- read.csv("../data/raw_data_tcga/spatial_TILs_saltz.csv", header = TRUE, row.names = 1)

# Cancer types
load("../analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_Spat.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)

for (ii in PanCancer.names){

  load(paste0("../data/PanCancer_draft_v1/",ii,"/DataViews_filter_spat_", ii,".RData"))
  keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[ii]], 1, 15),
                substr(rownames(all_cell.fractions), 1, 15))

  tmp_immunecells <- all_cell.fractions[keep,]
  rownames(tmp_immunecells) <- substr(rownames(tmp_immunecells),1,12)

  keep <- match(substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[ii]], 1, 12),
                substr(rownames(all.spatial.TILs), 1, 12))

  tmp_STILs <- all.spatial.TILs[keep,]
  features_to_keep <- c("til_percentage", "WCD_mean", "Ball_Hall", "Banfeld_Raftery", "C_index", "Det_Ratio")
  tmp_STILs <- tmp_STILs[, features_to_keep]
  tmp_STILs$Det_Ratio <- as.numeric(tmp_STILs$Det_Ratio)

  tmp_immunecells$T_cells_CD4 <- tmp_immunecells$T_cells_CD4 + tmp_immunecells$T_cells_regulatory_Tregs
  DataViews.filter_spat$ImmuneCells <- tmp_immunecells
  DataViews.filter_spat$sTIL <- tmp_STILs

  save(DataViews.filter_spat, file = paste0("../data/PanCancer_draft_v1/",ii,"/DataViews_filter_spat_",ii, ".RData"))
}




