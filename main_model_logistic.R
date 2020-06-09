# #########################################################################################################
# Script to logistic regression on validation data
# #########################################################################################################

# ****************
# working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ****************
# packages
library(devtools)
library(matrixStats)

# ****************
# scripts
source("./R/cross_validation_multi-task.R")

# ****************
# Select validation dataset 
load(paste0("./data/Validation/All_DataViews_test_pre.RData"))
load(paste0("./data/Validation/All_Labels_test_pre.RData"))
Datasets.names <- names(All.DataViews.test)
names(Datasets.names) <- c("SKCM", "SKCM", "SKCM", "STAD", "SKCM", "SKCM", "GBM", "SKCM", "SKCM")

# ****************
# views
views <- c(pathways = 'binomial', #1
           Protall = 'binomial', #2
           immunecells = 'binomial', #3
           TFs = 'binomial', #4
           transcript = 'binomial', #5
           sTIL = 'binomial', #6
           LRpairs = 'binomial', #7
           CYTOKINEpairs = 'binomial')  #8)

# ****************
view_combinations <- list(views[1], views[3], views[c(1,3)], views[4], views[5], views[7], views[8])

# ****************
# data
load("./data/parameters_4_all.RData")

# ****************
#input_algorithm = names(parameters)
input_algorithm <- c("Logistic_EN")

for (dataset in Datasets.names){
  
  DataViews.no_filter <- All.DataViews.test[["comb.Gide_Auslander"]]
  ImmuneResponse.no_filter <- All.Labels.test[["comb.Gide_Auslander"]]
  rownames(ImmuneResponse.no_filter) <- ImmuneResponse.no_filter$Sample
  ImmuneResponse.no_filter <- ImmuneResponse.no_filter[,-1, drop = FALSE]

  ImmuneResponse.no_filter$label <- gsub(" ","", ImmuneResponse.no_filter$label)
  ImmuneResponse.no_filter$label <- gsub("CR|PR|MR|PRCR","R", ImmuneResponse.no_filter$label)
  ImmuneResponse.no_filter$label <- gsub("SD|PD|PD ","NR", ImmuneResponse.no_filter$label)
  ImmuneResponse.no_filter$label <- factor(ImmuneResponse.no_filter$label, levels = c("NR", "R"))

  # Remove NA values in Dataviews LRpairs and CYTOKINE pairs
  if (anyNA(DataViews.no_filter$LRpairs)){
    
    LR_sum <- apply(DataViews.no_filter$LRpairs,2, sum)
    remove_NA_LR_pairs <- as.numeric(na.action(na.omit(LR_sum)))
    DataViews.no_filter$LRpairs <- DataViews.no_filter$LRpairs[,-remove_NA_LR_pairs]
    Cyt_sum <- apply(DataViews.no_filter$CYTOKINEpairs,2, sum)
    remove_NA_Cyt_pairs <- as.numeric(na.action(na.omit(Cyt_sum)))
    DataViews.no_filter$CYTOKINEpairs <- DataViews.no_filter$CYTOKINEpairs[,-remove_NA_Cyt_pairs]
    
  }
  
  sapply(1:length(view_combinations), function(X) {
    
    # Input data:
    i <- view_combinations[[X]]
    cat("View:", names(i), "\n")
    
    # Select k fold:
    K = 5
    
    all_cv_res <- lapply(input_algorithm, function(Y){
      
      # Coefficients estimation: optimized by cross validation
      cv_res <- cross_validation(drug_source = ImmuneResponse.no_filter,
                                 views_source = DataViews.no_filter,
                                 view_combination = i,
                                 algorithm = Y,
                                 standardize_any = T,
                                 standardize_response = T, # False in logistic regression
                                 parameters = parameters,
                                 k_fold=K,
                                 random=100)
      
      return(cv_res)
    })
    names(all_cv_res) <- input_algorithm
    # Save results: name file changes together with input file
    if (length(i) == 1){
      # save(all_cv_res, file = paste0("./output/BIM_cluster_presentation/","all_cv_res_",Cancer, "_train_rand100_with_combo_gide_auslander_",names(i)[1],".RData"))
      save(all_cv_res, file = paste0("./output/BIM_cluster_presentation/","all_cv_res_",Cancer, "_train_rand100_with_cor_tasks_",names(i)[1],".RData"))
      
    }else if(length(i) == 2){
      # save(all_cv_res, file = paste0("./output/BIM_cluster_presentation/","all_cv_res_",Cancer, "_train_rand100_with_combo_gide_auslander_",names(i)[1],"_",names(i)[2],".RData"))
      save(all_cv_res, file = paste0("./output/BIM_cluster_presentation/","all_cv_res_",Cancer, "_train_rand100_with_cor_tasks_",names(i)[1],"_",names(i)[2],".RData"))
      
    }
    
  })
  
}



