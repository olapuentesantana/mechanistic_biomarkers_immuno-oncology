#
# Run optimization for model training 
# BEMKL, RMLTR
#
#
#
#
#
#
#########################################################################################################

# working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# packages
library(devtools)
library(matrixStats)
library(easier) # install easier package from github

# source files
source("./tcga_training/randomized_crossvalidation.R")

# model prior parameters
load("../data/model_prior_parameters.RData")

# cancer type 
cancer_type <- "SKCM"

# ----------------------------------------- #
# RNA-seq data: tpm and count data
# ----------------------------------------- #
tpm <- read.table()
counts <- read.table()
# ----------------------------------------- #
# compute system-based derived signatures
# ----------------------------------------- #

# Computation of cell fractions
cell_fractions <- compute_cell_fractions(RNA.tpm=tpm)
# Computation of pathway activity
pathways_activity <- compute_pathways_scores(RNA.countss=counts, remove.genes.ICB_proxies=TRUE)
# Computation of TF activity
TF_activity <- compute_TF_activity(RNA.tpm=tpm, remove.genes.ICB_proxies=FALSE)
# Computation of LR pairs weights
lrpairs_weights <- compute_LR_pairs(RNA.tpm=tpm, remove.genes.ICB_proxies=FALSE, compute.cytokines.pairs=FALSE, cancertype="pancan")
# Computation of Cell-Cell scores
ccpairsgrouped_scores <- compute_CC_pairs_grouped(lrpairs=lrpairs_weights$LRpairs, cancertype="pancan")

# ----------------------------------------- #
# compute scores of immune response
# ----------------------------------------- #

tasks <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS") # Ock_IS was derived from publication
tmp_file_path <- c("")
tasks_values <- compute_gold_standards(RNA.tpm=tpm, list_gold_standards=tasks, cancertype=cancer_type, output_file_path=tmp_file_path)

# returns a list, convert to matrix
# Into matrix
immune_response <- do.call(cbind, lapply(tasks, function(X){
  immune_response <- t(tasks_values[[X]])
  return(immune_response)
}))
# Assess correlation between chemokines and the other correlated tasks
tasks_cormat <- cor(immune_response)
cor_sign <- sign(tasks_cormat[,"chemokines"])
cor_sign <- cor_sign[names(cor_sign) != "chemokines"]
if (all(cor_sign == -1)){
  immune_response[,"chemokines"] <- -immune_response[,"chemokines"]
}

# select view
view_name <- list(Pathways = "gaussian",
                  ImmuneCells = "gaussian")
# view data
DataViews <- list(Pathways = Pathway_activities,
                  ImmuneCells = cell_fractions)

# select algorithm 
algorithm <- "RMTLR" # BEMKL

# randomized cross-validation
rand_cv_res <- randomized_crossvalidation(drug_source = immune_response,
                                          views_source = DataViews,
                                          view_combination = view_name,
                                          algorithm = algorithm,
                                          standardize_any = T,
                                          standardize_response = T, 
                                          parameters = model_prior_param[[algorithm]],
                                          k_fold=5,
                                          random=100)












