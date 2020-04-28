###########################################################################
# Genes involved in transcriptomics signatures used as proxy's of response
###########################################################################

# ****************
# working directory
setwd("~/Desktop/PhD_TUE/Github_model/desktop/")

# Genes to remove according to IS and CYT
load("./data/list_genes_IS_CYT.Rdata")
# Genes that were not included in 'Predictive gene signature' IS:
ISCYT_read <- unique(list_genes_IS_CAT)

# Genes from IPS
IPSG_read <- read.table("data/raw_data_tcga/IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)

# Genes from IMPRES
IMPRES.checkpoint.pairs <- data.frame(Gene_1 = c("PDCD1","CD27","CTLA4","CD40","CD86", "CD28", "CD80", 
                                                 "CD274","CD86","CD40","CD86","CD40","CD28","CD40","TNFRSF14"),
                                      Gene_2 = c("TNFSF4","PDCD1","TNFSF4","CD28","TNFSF4", "CD86", "TNFSF9", 
                                                 "C10orf54","HAVCR2","PDCD1","CD200","CD80","CD276","CD274","CD86"))
IMPRES_read <- unique(as.vector(as.matrix(IMPRES.checkpoint.pairs))) # 15 genes

# Genes from Roh Immune Signature Score
RohISS_read <- c("GZMA", "GZMB", "PRF1", "GNLY", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", 
                 "HLA-G", "HLA-H", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", 
                 "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", 
                 "IFNG", "IFNGR1", "IFNGR2", "IRF1", "STAT1", "PSMB9", "CCR5", "CCL3", "CCL4", 
                 "CCL5", "CXCL9", "CXCL10", "CXCL11", "ICAM1", "ICAM2", "ICAM3", "ICAM4", "ICAM5", "VCAM1")

# Genes from 12-chemokine signature (Messina)
Chemokine_read <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21",
                    "CXCL9", "CXCL10", "CXCL11", "CXCL13")

# Genes from cell cycle signature, genes associated with proliferation not frequently amplified/deleted (Davoli)
Proliferation_read <- c("CENPE", "CCNA2", "CCNB2", "MCM6", "CCNF", "BUB1", "CDC20", "CDC6", "CDK1", "PLK1")

# Genes from immune infiltrate signature, genes specific for cytotoxic CD8+ T cells and NK cells (Davoli)
IS_Davoli_read <- c("CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK")

# # Genes from T-effector and interferon-γ gene signature (Feherenbacher): authors do not provide any computational measure
# Teff_IFny_read <- c("CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21") 

# IFNy 6-genes signature (Ayers):
IFNy_Ayers_read <- c("IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA") 
  
# Expanded immune 18-genes signature (Ayers):
ExpandedImmune_Ayers_read <- c("GZMB", "GZMK", "CXCR6", "CXCL10", "CXCL13", "CCL5", "STAT1","CD3D", "CD3E",
                               "CD2", "IL2RG" , "NKG7", "HLA-E", "CIITA"," HLA-DRA", "LAG3", "IDO1", "TAGAP") 

# 18-gene IFNγ T cell-inflamed signature (Ayers):
T_cell_inflamed_Ayers_read <- c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1", 
                          "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")

Housekeeping.genes <- c("STK11IP", "ZBTB34", "TBC1D10B", "OAZ1", "POLR2A", "G6PD", "ABCF1", "C14orf102", "UBB", "TBP", "SDHA")

weights <- data.frame(CCL5=0.008346, CD27=0.072293, CD274=0.042853, CD276=-0.0239, CD8A=0.031021 ,CMKLR1=0.151253, CXCL9=0.074135,
                      CXCR6=0.004313, HLA.DQA1=0.020091, HLA.DRB1=0.058806, HLA.E=0.07175, IDO1=0.060679, LAG3=0.123895, NKG7=0.075524, PDCD1LG2=0.003734,
                      PSMB10=0.032999, STAT1=0.250229, TIGIT=0.084767)

# TIDE
#TIDE_read <- 

# MSI score: 10 gene pairs 
MSI.pairs <- data.frame(Gene_1 = c("HNRNPL","MTA2","CALR","RASL11A","LYG1", "STRN3", "HPSE", 
                                                 "PRPF39","CCRN4L","AMFR"),
                        Gene_2 = c("CDC16","VGF","SEC22B","CAB39L","DHRS12", "TMEM192", "BCAS3", 
                                                 "ATF6","GRM8","DUSP18"))

MSI_read <- unique(as.vector(as.matrix(MSI.pairs))) # 20 genes


# VEGF-dependent vasculature (VDV) genes]
#angio_read <- c("EGFA", "KDR", "ESM1", "PECAM1", "ANGPTL4", "CD34")

# Teff 
#Teff_read <- c("CD8A", "EOMES", "PRF1", "IFNG", "CD274")

# Teff (Rosenberg)
Teff_read <- c("CD8A", "CXCL9", "CXCL10", "PRF1", "IFNG", "GZMA", "GZMB", "TBX21") 

# myeloid inflammation
#Myeloid_inflammation_read <-  c("IL-6", "CXCL1", "CXCL2", "CXCL3", "CXCL8", "PTGS2")

# Unify all genes from signatures

IS_read <- ISCYT_read[-length(ISCYT_read)]
ICB.proxies.genes <- list(CYT = c("GZMA", "PRF1"),
                                IS = IS_read,
                                IPS = as.character(IPSG_read$GENE),
                                IMPRES = IMPRES_read,
                                RohISS = RohISS_read,
                                Chemokine = Chemokine_read,
                                Proliferation = Proliferation_read,
                                IS_Davoli = IS_Davoli_read,
                                IFny = IFNy_Ayers_read,
                                ExpandedImmune = ExpandedImmune_Ayers_read,
                                T_cell_inflamed = T_cell_inflamed_Ayers_read,
                                TIDE = NULL,
                                MSI = MSI_read)
save(ICB.proxies.genes, file = "data/list_each_ICB_proxy_with_involved_genes.RData")


