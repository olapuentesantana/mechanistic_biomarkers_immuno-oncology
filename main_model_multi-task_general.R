# #########################################################################################################
# Script to apply multi-task learning on TCGA data
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
# Select cancer type
load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
# load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_prot.RData")
# load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_Spat.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)
#PanCancer.names.some <- PanCancer.names[3:length(PanCancer.names)]

# ****************
# views
views <- c(pathways = 'gaussian', #1
           Protall = 'gaussian', #2
           immunecells = 'gaussian', #3
           TFs = 'gaussian', #4
           transcript = 'gaussian', #5
           sTIL = 'gaussian', #6
           LRpairs = 'gaussian', #7
           CYTOKINEpairs = 'gaussian')  #8)

# ****************
view_combinations <- list(views[1], views[3], views[c(1,3)], views[4], views[5], views[7], views[8])

# ****************
# data
load("./data/parameters_4_all.RData")

# ****************
#input_algorithm = names(parameters)
input_algorithm <- c("Multi_Task_EN")
filter_tasks <- c("IPS","IMPRES", "Proliferation", "TIDE", "MSI","chemokine")

for (Cancer in PanCancer.names){
  Cancer = "SKCM"
  
  load(paste0("./data/PanCancer/",Cancer,"/new/DataViews_filter_Spat_", Cancer,".RData"))
  load(paste0("./data/PanCancer/",Cancer,"/new/ImmuneResponse_filter_Spat_", Cancer,"_matrix_format.RData"))
  
  load(paste0("./data/BIM_cluster_presentation/","DataViews_no_filter_", Cancer,".RData"))
  load(paste0("./data/BIM_cluster_presentation/","ImmuneResponse_no_filter_", Cancer,"_matrix_format.RData"))
  
  # Multi-gaussian regression on validation dataset
  ImmuneResponse.no_filter <- ImmuneResponse.no_filter[, -which(colnames(ImmuneResponse.no_filter) %in% filter_tasks)]
  
  # sTIL_sum <- apply(DataViews.filter_Spat$sTIL,1, sum)
  # 
  # if(anyNA(sTIL_sum)){
  #   remove_NA_sTIL <- as.numeric(na.action(na.omit(sTIL_sum)))
  #   DataViews.filter_Spat$sTIL <- DataViews.filter_Spat$sTIL[-remove_NA_sTIL,]
  #   DataViews.filter_Spat$immunecells <- DataViews.filter_Spat$immunecells[-remove_NA_sTIL,]
  #   
  # }

  # Remove NA values in Dataviews LRpairs and CYTOKINE pairs
  if (anyNA(DataViews.no_filter$LRpairs)){
    
    LR_sum <- apply(DataViews.no_filter$LRpairs,2, sum)
    remove_NA_LR_pairs <- as.numeric(na.action(na.omit(LR_sum)))
    DataViews.no_filter$LRpairs <- DataViews.no_filter$LRpairs[,-remove_NA_LR_pairs]
    Cyt_sum <- apply(DataViews.no_filter$CYTOKINEpairs,2, sum)
    remove_NA_Cyt_pairs <- as.numeric(na.action(na.omit(Cyt_sum)))
    DataViews.no_filter$CYTOKINEpairs <- DataViews.no_filter$CYTOKINEpairs[,-remove_NA_Cyt_pairs]
    
  }
  # Remove 0 variance features 
  
  # Consensus define as mean betweeen both responses. Used just for testing
  # ImmuneResponse$consensus <- apply(scale(ImmuneResponse[,1:2]),1,mean)
  # Remove IPS and IMPRES
  #ImmuneResponse$IPS <- NULL
  #ImmuneResponse$impres <- NULL
  
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
 

   
