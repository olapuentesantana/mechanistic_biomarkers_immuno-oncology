# ***************************************************************
# Prepara data for Maisa #
# TPM values for each cancer type in a separated .txt file #
# ***************************************************************

# ****************
# working directory
setwd("~/Desktop/PhD_TUE/Github_model/desktop/")

###########################################################################
# TCGA data
###########################################################################

# *****************
# Cancer types
load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

sapply(PanCancer.names, function(Cancer){
  
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
  # Log2 transformed
  #log2mas1.TPM.transcripts <- log2(TPM.transcripts + 1)
  
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
  
  write.table(TPM.transcripts, file = paste0("./data/Maisa/TCGA_tpm_",Cancer,".txt"), sep = "\t")
  
})

###########################################################################
# Validation data
###########################################################################

# ***************************************************************
# Hugo: melanoma (PD-1)
# ***************************************************************
sampleinfo <- read.table("data/Validation/Francesca/Hugo/Hugo_patients.txt", 
                         sep = "\t", header = TRUE, row.names = 1)
tpm_hugo <- as.data.frame(t(read.table("data/Validation/Francesca/Hugo/Hugo_gene_tpm.txt", 
                            sep = "\t", header = TRUE, row.names = 1)))
## Match samples
tpm_hugo <- tpm_hugo[match(rownames(sampleinfo), rownames(tpm_hugo)), ]

write.table(tpm_hugo, file = paste0("./data/Maisa/hugo_tpm.txt"), sep = "\t")

# ***************************************************************
# Riaz: melanoma (PD-1)
# ***************************************************************
load("data/Validation/Francesca/Riaz/Riaz_sampleinfo.rdata")
tpm_riaz <- as.data.frame(t(read.table("data/Validation/Francesca/Riaz/Riaz_gene_tpm.txt", 
                         sep = "\t", header = TRUE, row.names = 1)))
## Match samples
tpm_riaz <- tpm_riaz[match(sampleinfo$SRA, rownames(tpm_riaz)), ]

write.table(tpm_riaz, file = paste0("./data/Maisa/riaz_tpm.txt"), sep = "\t")

# ***************************************************************
# Auslander: melanoma (PD-1, CTLA-4)
# ***************************************************************
load("data/Validation/Francesca/Auslander/Auslander_sampleinfo.rdata")
tpm_auslander <- as.data.frame(t(read.table("data/Validation/Francesca/Auslander/Auslander_gene_tpm.txt", 
                              sep = "\t", header = TRUE, row.names = 1)))
## Match samples: Only 36 patients have RNAseq data
tpm_auslander <- tpm_auslander[na.omit(match(sampleinfo$SRA, rownames(tpm_auslander))), ]

write.table(tpm_auslander, file = paste0("./data/Maisa/auslander_tpm.txt"), sep = "\t")

# ***************************************************************
# Liu: melanoma (PD-1)
# ***************************************************************
sampleinfo <- read.csv(file = "data/Validation/Internet/LIU/LIU_clinical_data.csv",
                       sep = ",", header = TRUE, row.names = 2, nrows = 144)

tpm_liu <- as.data.frame(read.csv(file = "data/Validation/Internet/LIU/LIU_tpm_RNAseq.csv",
                        sep = ",",  header = TRUE, row.names = 1))

## Match samples: Only 121 patients have RNAseq data
sampleinfo <- sampleinfo[match(rownames(tpm_liu), rownames(sampleinfo)),]

write.table(tpm_liu, file = paste0("./data/Maisa/liu_tpm.txt"), sep = "\t")
# ***************************************************************
# Gide: melanoma (PD-1, PD-L1, CTLA-4)
# ***************************************************************

sampleinfo <- read.csv(file = "data/Validation/Francesca/Gide/PRJEB23709/PRJEB23709_SraRunTable.txt",
                       sep = ",", header = TRUE, row.names = 1)

tpm_gide <- as.data.frame(t(read.csv(file = "data/Validation/Francesca/Gide/PRJEB23709/PRJEB23709_gene_tpm.txt",
                                  sep = "\t",  header = TRUE, row.names = 1)))

## Match samples: 
tpm_gide <- tpm_gide[match(rownames(sampleinfo), rownames(tpm_gide)), ]

write.table(tpm_gide, file = paste0("./data/Maisa/gide_tpm.txt"), sep = "\t")


# ***************************************************************
# Zhao: GBM (PD-1)
# ***************************************************************

sampleinfo <- read.csv(file = "data/Validation/Francesca/Zhao/PRJNA482620/PRJNA482620_SraRunTable.txt",
                       sep = ",", header = TRUE, row.names = 1)

tpm_zhao <- as.data.frame(t(read.csv(file = "data/Validation/Francesca/Zhao/PRJNA482620/PRJNA482620_gene_tpm.txt",
                                     sep = "\t",  header = TRUE, row.names = 1)))

## Match samples: Only 17  patients have RNAseq data (34 samples)
tpm_zhao <- tpm_zhao[na.omit(match(rownames(sampleinfo), rownames(tpm_zhao))), ]

write.table(tpm_zhao, file = paste0("./data/Maisa/zhao_tpm.txt"), sep = "\t")

# ***************************************************************
# Kim: STAD (PD-1)
# ***************************************************************

sampleinfo <- read.csv(file = "data/Validation/Francesca/Kim/PRJEB25780/PRJEB25780_SraRunTable.txt",
                       sep = ",", header = TRUE, row.names = 1)

tpm_kim <- as.data.frame(t(read.csv(file = "data/Validation/Francesca/Kim/PRJEB25780/PRJEB25780_gene_tpm.txt",
                                     sep = "\t",  header = TRUE, row.names = 1)))

## Match samples: check why info data says 78 samples with RNA-seq data, while we find 188 samples in tpm matrix
tpm_kim <- tpm_kim[na.omit(match(rownames(sampleinfo), rownames(tpm_kim))), ]

write.table(tpm_kim, file = paste0("./data/Maisa/kim_tpm.txt"), sep = "\t")




