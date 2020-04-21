# #########################################################################################################
# Script to test algorithm on validation datasets
# #########################################################################################################

# ****************
# working directory
if(!("rstudioapi" %in% installed.packages()[,"Package"])) install.packages("rstudioapi",  quiet = TRUE)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ****************
# scripts
source("compute.pathways.scores.R")
source("compute.TF.activity.R")

# ****************
# RNA-seq data 

# Raw counts values
gene.count <- read.table(file = paste0("_gene_count.txt"),  
                         sep = "\t", header = TRUE, row.names = 1)
# Tpm values
gene.tpm <- read.table(file = paste0("_gene_tpm.txt"),  
                       sep = "\t", header = TRUE, row.names = 1)

# ****************
# Computation of pathways scores (input matrix [genes, samples], ouput matrix [sample, pathways])
Pathway_activities <- compute.pathways.scores(RNA.raw_counts = gene.count)

# ****************
# Computation of TF activity (input matrix [genes, samples], ouput matrix [sample, TFs])
TF_activities <- compute.TF.activity(RNA.tpm = gene.tpm)


