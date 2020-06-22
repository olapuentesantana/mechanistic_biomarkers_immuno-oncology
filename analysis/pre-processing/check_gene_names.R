# #########################################################################################################
# Script to check and update gene hugo symbols using multi-symbol checker
# #########################################################################################################

# *****************
# working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# *****************
# Write gene names into a .csv file
Cancer = "GBM"
cat("\n",Cancer,"\n")

file <- dir(paste0("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/raw_data_tcga/RNAseq/20160128_version/stddata__2016_01_28/", Cancer,"/20160128/",
                   "gdac.broadinstitute.org_", Cancer,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0"),
            recursive = TRUE, pattern = "data.txt", full.names = TRUE)

# Extract the raw counts from the text files for each gene
data.transcripts <- read.csv(file, sep="\t",header=T,  row.names=1, stringsAsFactors=F, check.names=FALSE)
colnames(data.transcripts) <- substr(colnames(data.transcripts), 1, 15)

# TPM/scaled_estimate values (transcriptomics data)
get_estimates <- which(data.transcripts[1,] == "raw_count") #scaled_estimate
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
HGNC_symbol <- HGNC_symbol[-which(HGNC_symbol == "?")]
entrez_id <- sapply(strsplit(rownames(TPM.transcripts),"\\|"),function(X) return(X[2]))
anyDuplicated(entrez_id)
rownames(TPM.transcripts) <- entrez_id

library(org.Hs.eg.db)
library(AnnotationDbi)

# Entrez id hugo symbols
Symbols_from_id <- AnnotationDbi::select(org.Hs.eg.db, keys=entrez_id, columns='SYMBOL', keytype='ENTREZID')
id_from_symbols <- AnnotationDbi::select(org.Hs.eg.db, keys=all_regulated_transcripts, columns='ENTREZID', keytype='SYMBOL')

# TCGA hugo symbols
pos_1 <- match(Symbols_from_id[!is.na(Symbols_from_id$SYMBOL),"ENTREZID"], entrez_id)
rownames(TPM.transcripts)[pos_1] <- Symbols_from_id[!is.na(Symbols_from_id$SYMBOL),"SYMBOL"]
pos_2 <- match(Symbols_from_id[is.na(Symbols_from_id$SYMBOL),"ENTREZID"], entrez_id)
rownames(TPM.transcripts)[pos_2] <- HGNC_symbol[pos_2]

# Rows corresponding to the same HGNC symbol were averaged.
gene.symbol <- rownames(TPM.transcripts)
if(anyDuplicated(gene.symbol) != 0){
  idx <- which(duplicated(gene.symbol) == TRUE)
  dup_genes <- gene.symbol[idx]
  for (ii in dup_genes){
    TPM.transcripts[which(gene.symbol %in% ii)[1],] <- colMeans(TPM.transcripts[which(gene.symbol %in% ii),])
    TPM.transcripts <- TPM.transcripts[-which(gene.symbol %in% ii)[2],]
    gene.symbol <- gene.symbol[-which(gene.symbol %in% ii)[2]]
  }
}
rownames(TPM.transcripts) <- gene.symbol

# write.csv(rownames(TPM.transcripts), file = paste0("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/multi_symbol_checker/gene_names_to_be_checked_TCGA.csv"),
#           quote = FALSE, row.names = FALSE, col.names = FALSE)
  

# *****************
# Read correct gene names from multi-symbol checker

gene.names_corrected <- read.csv(file = paste0("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/multi_symbol_checker/hgnc-symbol-check_TCGA.csv"),
                                 header = TRUE, row.names = NULL, skip = 1)

tmp <- gene.names_corrected[duplicated(gene.names_corrected$Input),]


trial.data <- log2(TPM.transcripts + 1)

## Quantile normalization
correct_symbols_trial.norm <- normalize.quantiles(as.matrix(trial.data))
dimnames(correct_symbols_trial.norm) <- dimnames(trial.data)

## Mean-centralization: substract gene average across all conditions
average.gene <- rowMeans(correct_symbols_trial.norm)
correct_symbols_trial.norm.all <- sweep(correct_symbols_trial.norm, 1, average.gene, FUN = "-")

write.table(correct_symbols_trial.norm.all, file = paste0("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/data_processed_TIDE/correct_symbols_trial.txt"), sep = "\t")
system(paste0("/anaconda3/bin/tidepy ", "~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/data_processed_TIDE/correct_symbols_trial.txt", " -o ~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/data_processed_TIDE/output_trial.txt -c Other"))

  
TIDE.table.trial <- read.table(file = paste0("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/data_processed_TIDE/output_trial.txt"), sep = "\t", header = T, row.names = 1)
TIDE.table.norm <- read.table(file = paste0("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/data_processed_TIDE/output_tide_norm_GBM.txt"), sep = "\t", header = T, row.names = 1)

TIDE.score.gbm <- data.frame(TIDE = TIDE.table.norm[,"TIDE", drop = FALSE], check.names = F)
TIDE.score.trial <- data.frame(TIDE = TIDE.table.trial[,"TIDE", drop = FALSE], check.names = F)

cor(TIDE.score.gbm, TIDE.score.trial)



