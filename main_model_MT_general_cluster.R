# #########################################################################################################
# Script to apply multi-task learning
# #########################################################################################################

# ****************
# working directory
setwd("/home/feduati/Oscar/PhD_Oscar/new/")

# ****************
# packages
library(devtools, lib.loc= Sys.getenv("R_LIBS"))
library(matrixStats,lib.loc= Sys.getenv("R_LIBS"))

# ****************
# Input: select view
argsJob <- commandArgs(trailingOnly=TRUE)
View <- print(argsJob[1])

View_input <- vector("character", length = length(strsplit(View,"_", fixed = TRUE)[[1]]))
View_input <- do.call("c",lapply(1:length(View_input), function(X){
    tmp <- "gaussian"
    names(tmp) <- sapply(strsplit(View,"_", fixed = TRUE), function(Y) return(Y[X]))
    return(tmp)
}))

# ****************
# Input: select cancer type
CancerType <- print(argsJob[2])

# ****************
# scripts
source("R/PanCancer_draft_v1/cross_validation_MT_cluster.R")

# ****************
# data
load("data/parameters_4_all.RData")

# ****************
# algorithms
input_algorithm <- c("Multi_Task_EN", "BEMKL")

# ****************
# Uncorrelated tasks
filter_tasks <- c("IPS","IMPRES", "Proliferation", "TIDE", "MSI")

load(paste0("data/PanCancer_draft_v1/", CancerType,"/DataViews_filter_spat_", CancerType,".RData"))
load(paste0("data/PanCancer_draft_v1/", CancerType,"/ImmuneResponse_filter_spat_", CancerType,".RData"))
  
# load(paste0("data/PanCancer_draft_v1/DataViews_no_filter_", Cancer,".RData"))
# load(paste0("data/PanCancer_draft_v1/ImmuneResponse_no_filter_", Cancer,".RData"))

# Multi-gaussian regression on validation dataset
ImmuneResponse.filter_spat <- ImmuneResponse.filter_spat[, !colnames(ImmuneResponse.filter_spat) %in% filter_tasks]

# Remove NA values in DataViews sTIL (sample-wise)
if(anyNA(DataViews.filter_spat$sTIL$Det_Ratio)){
  DataViews.filter_spat$sTIL <- DataViews.filter_spat$sTIL[!is.na(DataViews.filter_spat$sTIL$Det_Ratio),]
  DataViews.filter_spat$ImmuneCells <- DataViews.filter_spat$ImmuneCells[match(rownames(DataViews.filter_spat$sTIL),
                                                                               rownames(DataViews.filter_spat$ImmuneCells)),]
  ImmuneResponse.filter_spat <- ImmuneResponse.filter_spat[match(rownames(DataViews.filter_spat$sTIL),
                                                                 rownames(ImmuneResponse.filter_spat)),]
}

# Remove NA values in DataViews LRpairs and CYTOKINE pairs (feature-wise)
if (anyNA(DataViews.filter_spat$LRpairs)){
  LR_sum <- apply(DataViews.filter_spat$LRpairs, 2, sum)
  DataViews.filter_spat$LRpairs <- DataViews.filter_spat$LRpairs[,!is.na(LR_sum)]
  Cyt_sum <- apply(DataViews.filter_spat$CYTOKINEpairs, 2, sum)
  DataViews.filter_spat$CYTOKINEpairs <- DataViews.filter_spat$CYTOKINEpairs[,!is.na(Cyt_sum)]
}

cat("View:", names(View_input), "\n")
cat("Cancer Type:", CancerType, "\n")

# Select k fold:
K = 5

all_cv_res <- lapply(input_algorithm, function(Y){
  cv_res <- cross_validation.MT(drug_source = ImmuneResponse.filter_spat,
                                views_source = DataViews.filter_spat,
                                view_combination = View_input,
                                algorithm = Y,
                                standardize_any = T,
                                standardize_response = T,
                                parameters = parameters,
                                k_fold=K,
                                random=100)
  
  return(cv_res)
})
names(all_cv_res) <- input_algorithm
# Save results: name file changes together with input file
if (length(View_input) == 1){
  save(all_cv_res, file = paste0("output/PanCancer_draft_v1/",CancerType,"/all_cv_res_",CancerType,
                                 "_train_rand100_with_all_tasks_", names(View_input)[1],".RData"))
}else if(length(View_input) == 2){
  save(all_cv_res, file = paste0("output/PanCancer_draft_v1/",CancerType,"/all_cv_res_",CancerType,
                                 "_train_rand100_with_all_tasks_", names(View_input)[1],"_",names(View_input)[2],".RData"))
}
 
    
    
    
    