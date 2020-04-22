###########################################################################################################
# Script to match samples from transcriptomics data (TCGA) with:

# Immune Signature
# Spatial Information of TILs

# This data is the limiting factor for the amount of patients considered

# #########################################################################################################

# ****************
# working directory
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ****************
# data
PanCancer.names <- dir("~/Desktop/Master_TUE/MT_Oscar/Raw_Data/Transcriptomics/20160128_version/stddata__2016_01_28/")

## screening data
IS <- read.csv("~/Desktop/Master_TUE/MT_Oscar/Raw_Data/Patient_Response/Aprox/Output_Data/IS.csv", sep=";", header = T, skip = 1)
quanTIseq <- read.csv("~/Desktop/Master_TUE/MT_Oscar/Raw_Data/Patient_Response/Aprox/Input_Data/quanTIseq_estimated.csv", sep=",", header = T)
spatial_TIL <- read.csv("~/Desktop/PhD_TUE/Github_model/desktop/data/raw_data/spatial_TILs_saltz.csv", sep=",", header = T)

# samples tcga id for IS, quanTISeq and Spatial_TIls
IS.samples <- as.vector(as.character(IS$TCGA_ID))
quanTIseq.samples <- as.vector(as.character(quanTIseq$X))
spatial_TIL.samples <- as.vector(as.character(spatial_TIL$ParticipantBarcode)) # All primary umors I think, except for skin cancer (metastatis can be included)

## RNA-seq data
TCGA.samples.pancancer <- do.call(c, lapply(PanCancer.names, function(Cancer){

  TCGA.samples.cancer <- vector("list", length = 1)
  names(TCGA.samples.cancer) <- Cancer
  
  file <- dir(paste0("~/Desktop/Master_TUE/MT_Oscar/Raw_Data/Transcriptomics/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
                     "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0"),
                     recursive = TRUE, pattern = "data.txt", full.names = TRUE)
            
  # Extract the raw counts from the text files for each gene
  data <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
  get_rawcounts <- which(data[1,] == "raw_count")
  # Extract tcga samples per cancer type
  TCGA.samples.cancer[[Cancer]] <- colnames(data)[get_rawcounts]
  
  return(TCGA.samples.cancer)
  
}))


# ****************
# Select cancer type: (* do not considered cancer types with n<100) check this with federica


# ****************
# Sample screening:
# 1. quanTIseq data (it considered all cancer types and primary tumors)

TCGA.samples.pancancer_with_screen_quantiseg <- do.call(c, lapply(PanCancer.names, function(Cancer){
  
  TCGA.samples.cancer_with_screen_quantiseg <- vector("list", length = 1)
  names(TCGA.samples.cancer_with_screen_quantiseg) <- Cancer
  
  keep <- na.omit(match(quanTIseq.samples, substr(TCGA.samples.pancancer[[Cancer]], 1, 15)))
  TCGA.samples.cancer_with_screen_quantiseg[[Cancer]] <- TCGA.samples.pancancer[[Cancer]][keep]
  
  return(TCGA.samples.cancer_with_screen_quantiseg)
  
}))

remove_cancers_no_present_in_quantiseq <- sapply(PanCancer.names, function(X) is_empty(TCGA.samples.pancancer_with_screen_quantiseg[[X]]))
TCGA.samples.pancancer_with_screen_quantiseg$ESCA <- NULL
TCGA.samples.pancancer_with_screen_quantiseg$LAML <- NULL
TCGA.samples.pancancer_with_screen_quantiseg$PCPG <- NULL
TCGA.samples.pancancer_with_screen_quantiseg$SARC <- NULL
PanCancer.names <- PanCancer.names[-which(PanCancer.names %in% c("ESCA", "LAML", "PCPG", "SARC"))]

# 2. Two directions:
## 2a. Immune Signature

TCGA.samples.pancancer_with_screen_quantiseg_IS <- do.call(c, lapply(PanCancer.names, function(Cancer){
  
  TCGA.samples.cancer_with_screen_quantiseg_IS <- vector("list", length = 1)
  names(TCGA.samples.cancer_with_screen_quantiseg_IS) <- Cancer
  
  keep <- na.omit(match(IS.samples, substr(TCGA.samples.pancancer_with_screen_quantiseg[[Cancer]], 1, 15)))
  TCGA.samples.cancer_with_screen_quantiseg_IS[[Cancer]] <- TCGA.samples.pancancer_with_screen_quantiseg[[Cancer]][keep]
  
  return(TCGA.samples.cancer_with_screen_quantiseg_IS)
  
}))

## 2b. Immune Signature and Spatial Information on TILs

# Cases were staged according to the American Joint Committee on Cancer (AJCC). 
# Each frozen primary tumor specimen had a companion normal tissue specimen (blood or blood components, 
# including DNA extracted at the tissue source site). Adjacent tissue was submitted for some cases. 

TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs <- do.call(c, lapply(PanCancer.names, function(Cancer){
  
  TCGA.samples.cancer_with_screen_quantiseg_IS_SpatialTILs <- vector("list", length = 1)
  names(TCGA.samples.cancer_with_screen_quantiseg_IS_SpatialTILs) <- Cancer
  
  keep <- na.omit(match(spatial_TIL.samples, substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 12)))
  TCGA.samples.cancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]] <- TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]][keep]
  
  return(TCGA.samples.cancer_with_screen_quantiseg_IS_SpatialTILs)
  
}))

# # Remove COAD and READ as separate cancer types in both datasets:
# 
# if (all(c(TCGA.samples.pancancer_with_screen_quantiseg_IS$COAD, TCGA.samples.pancancer_with_screen_quantiseg_IS$READ)
#         %in% TCGA.samples.pancancer_with_screen_quantiseg_IS$COADREAD)){
#   
#   TCGA.samples.pancancer_with_screen_quantiseg_IS$COAD <- NULL
#   TCGA.samples.pancancer_with_screen_quantiseg_IS$READ <- NULL
#   #names(TCGA.samples.pancancer_with_screen_quantiseg_IS)[which(names(TCGA.samples.pancancer_with_screen_quantiseg_IS) == "COADREAD")] <- "CRC"
#   
# }
#               
# if (all(c(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs$COAD, TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs$READ)
#         %in% TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs$COADREAD)){
#   
#   TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs$COAD <- NULL
#   TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs$READ <- NULL
#   #names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)[which(names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs) == "COADREAD")] <- "CRC"
#   
# }             
# PanCancer.names <- PanCancer.names[-which(PanCancer.names %in% c("COAD", "READ"))]

# Remove empty cancer types
remove_cancers_no_present_in_quantiseq_IS_SpatialTILs <- sapply(PanCancer.names, function(X) is_empty(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[X]]))
TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs$GBM <- NULL
TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs$HNSC <- NULL
TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs$KIRC <- NULL
TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs$KIRP <- NULL
TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs$LIHC <- NULL
TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs$OV <- NULL
TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs$THCA <- NULL

# names(TCGA.samples.pancancer_with_screen_quantiseg_IS)[which(names(TCGA.samples.pancancer_with_screen_quantiseg_IS) == "COADREAD")] <- "CRC"
# names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)[which(names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs) == "COADREAD")] <- "CRC"
# 
# PanCancer.names <- PanCancer.names[-which(PanCancer.names %in% c("COAD", "READ"))]
# PanCancer.names[which(PanCancer.names %in% c("COADREAD"))] <- "CRC"


# save(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs, file = "TCGA_samples_available_screening_with_quanTIseq_IS_Spat.RData")
# save(TCGA.samples.pancancer_with_screen_quantiseg_IS, file = "TCGA_samples_available_screening_with_quanTIseq_IS.RData")

# Check if the proteomics data was the limitation for the small amount of samples
# Extract the raw counts from the text files for each gene
data.prot <- read.csv("~/Desktop/PhD_TUE/Github_Elastic_Net/data/Proteomics/TCGA-PANCAN32-L4.csv", sep=",", header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)

## Protein abundance data
TCGA.samples.pancancer.prot <- do.call(c, lapply(PanCancer.names, function(Cancer){
  
  TCGA.samples.cancer.prot <- vector("list", length = 1)
  names(TCGA.samples.cancer.prot) <- Cancer

  # Extract tcga samples per cancer type
  TCGA.samples.cancer.prot[[Cancer]] <- rownames(data.prot)[which(data.prot$Cancer_Type == Cancer)]
  
  return(TCGA.samples.cancer.prot)
  
}))

# 3. Two directions:
## 3a. Immune Signature + Protein data

TCGA.samples.pancancer_with_screen_quantiseg_IS_prot <- do.call(c, lapply(PanCancer.names, function(Cancer){
  
  TCGA.samples.cancer_with_screen_quantiseg_IS_prot <- vector("list", length = 1)
  names(TCGA.samples.cancer_with_screen_quantiseg_IS_prot) <- Cancer
  
  keep <- na.omit(match(substr(TCGA.samples.pancancer.prot[[Cancer]],1,15), substr(TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]], 1, 15)))
  TCGA.samples.cancer_with_screen_quantiseg_IS_prot[[Cancer]] <- TCGA.samples.pancancer_with_screen_quantiseg_IS[[Cancer]][keep]
  
  return(TCGA.samples.cancer_with_screen_quantiseg_IS_prot)
  
}))

TCGA.samples.pancancer_with_screen_quantiseg_IS_prot$COADREAD <- c(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot$COAD, TCGA.samples.pancancer_with_screen_quantiseg_IS_prot$READ)
TCGA.samples.pancancer_with_screen_quantiseg_IS_prot$COAD <- NULL
TCGA.samples.pancancer_with_screen_quantiseg_IS_prot$READ <- NULL

## 3b. Immune Signature + Spatial Information on TILs + Protein data

TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs_prot <- do.call(c, lapply(PanCancer.names, function(Cancer){
  
  TCGA.samples.cancer_with_screen_quantiseg_IS_SpatialTILs_prot <- vector("list", length = 1)
  names(TCGA.samples.cancer_with_screen_quantiseg_IS_SpatialTILs_prot) <- Cancer
  
  keep <- na.omit(match(substr(TCGA.samples.pancancer.prot[[Cancer]],1,15), substr(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]], 1, 15)))
  TCGA.samples.cancer_with_screen_quantiseg_IS_SpatialTILs_prot[[Cancer]] <- TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs[[Cancer]][keep]
  
  return(TCGA.samples.cancer_with_screen_quantiseg_IS_SpatialTILs_prot)
  
}))

TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs_prot$COADREAD <- c(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs_prot$COAD, TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs_prot$READ)
TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs_prot$COAD <- NULL
TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs_prot$READ <- NULL


# save(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs, file = "TCGA_samples_available_screening_with_quanTIseq_IS_Spat.RData")
# save(TCGA.samples.pancancer_with_screen_quantiseg_IS, file = "TCGA_samples_available_screening_with_quanTIseq_IS.RData")
save(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs_prot, file = "TCGA_samples_available_screening_with_quanTIseq_IS_Spat_prot.RData")
save(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot, file = "TCGA_samples_available_screening_with_quanTIseq_IS_prot.RData")


