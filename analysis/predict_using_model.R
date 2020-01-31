# #########################################################################################################
# Script to obtain predictions based on the models obtained during training..
## Validation datasets: Riaz, Auslander, Hugo, LIU.
# #########################################################################################################

# ****************
# working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ****************
# packages
library(reshape2)
library(pdist)

# ****************
# scripts
source("./R/BEMKL/predict_BEMKL.R")
source("./R/GLMs/predict_GLMs.R")

# ****************
# Select cancer type
PanCancer.names = c("BLCA",  "BRCA","CESC", "CRC", "LUAD", "LUSC", "PRAD", "SKCM", "STAD", "UCEC")

Cancer = "SKCM"

# ****************
# output data
files <- dir(paste0("output/model/",Cancer), pattern = "all_cv_res", full.names = T, recursive = T)

# ****************
# Initialize variables
views <- c(PROGENy = 'gaussian', #1
           Protall = 'gaussian', #2
           quanTIseq = 'gaussian', #3
           DoRothEAv1 = 'gaussian', #4
           transcript = 'gaussian', #5
           SpatialTIL = 'gaussian')  #6

view_combinations <- list(views[3],
                          views[c(3,6)])

label_IS = list()
predictions_IS.min = predictions_IS.1se = list()
label_CYT = list()
predictions_CYT.min = predictions_CYT.1se = list()
label_IPS = list()
predictions_IPS.min = predictions_IPS.1se = list()
label_impres = list()
predictions_impres.min = predictions_impres.1se = list()
label_consensus = list()
predictions_consensus.min = predictions_consensus.1se = list()
label_common = list()
predictions_common.min = predictions_common.1se = list()

standardize_any = T
K = 100
input_algorithm =  c("Lasso", "BEMKL", "Elastic_Net", "CV_linear_reg_L1&L2", "Ridge")
k_fold <- c(5,10)
summary_test_views_alg_kfold <- vector("list", length = length(view_combinations))
summary_alg_kfold <- vector("list", length = length(k_fold))

# ****************
# data
# Training sets Dataviews data: TCGA
load(paste0("./data/PanCancer/",Cancer,"/DataViews_", Cancer,"_all_train_with_sTIL.RData"))
# Validation sets DataViews data: pre-treatment
load("data/Validation/DataViews_pre_treatment_test_Riaz_Auslander_Hugo.RData")
# Validation sets ImmuneResponse data: pre-treatment
load("data/Validation/ImmuneResponse_pre_treatment_test_Riaz_Auslander_Hugo.RData")

DataViews.test.pre$Auslander$transcript = log2(DataViews.test.pre$Auslander$transcript + 1)
DataViews.test.pre$Riaz$transcript = log2(DataViews.test.pre$Riaz$transcript + 1)
DataViews.test.pre$Hugo$transcript = log2(DataViews.test.pre$Hugo$transcript + 1)

# Normalization should be done taking into account the train set. #
for (j in 1:length(view_combinations)){
  
  input_name = paste(names(view_combinations[[j]]), collapse ="_")
  names(summary_test_views_alg_kfold)[j] <- input_name
  
  cat("View:",input_name, "\n")
  
  summary_test_alg_kfold <- lapply(names(ImmuneResponse.test.pre), function(X){
    
    # check data type is available for validation data
    check = all(names(view_combinations[[j]]) %in% names(DataViews.test.pre[[X]]))
    
    if(check != TRUE) stop("View: ",input_name, " not available for validation dataset")
    #stopifnot(check == TRUE)
    
    cat("Validation data set:",X, "\n")
    
    summary_alg_kfold <- lapply(input_algorithm, function(alg){
      
      cat("Model -->", alg, "\n")
      
      view_combination <- view_combinations[[j]]
      
      if (alg %in% c("Lasso", "Elastic_Net", "CV_linear_reg_L1&L2","Ridge")){
        
        summary_kfold <- lapply(k_fold, function(fold){
        
          cat("kfold =", fold, "\n")
        
          pos_file = grep(paste0("_",input_name,"_kfold_",fold), files, fixed = TRUE)
          load(files[pos_file])

          # Generalized linear models:
          summary_kfold <- predict_GLMs(DataViews, DataViews.test.pre, ImmuneResponse.test.pre, j, K,
                                        all_cv_res, alg, view_combination, standardize_any, X) 
          return(summary_kfold)
        })
        
        names(summary_kfold) <- k_fold
        summary_alg_kfold <- summary_kfold
        return(summary_alg_kfold)
        
      }else{
        
        pos_file = grep(paste0("_",input_name,"_kfold_",5), files, fixed = TRUE)
        load(files[pos_file])
        
        # BEMKL:
        summary_BEMKL <- predict_BEMKL(DataViews, DataViews.test.pre, ImmuneResponse.test.pre, j, K,
                                         all_cv_res, alg, view_combination, standardize_any, X) 
      }
      summary_alg_kfold$BEMKL <- summary_BEMKL
    })
    names(summary_alg_kfold) <- input_algorithm
    summary_test_alg_kfold <- summary_alg_kfold
    return(summary_test_alg_kfold)
  })
  names(summary_test_alg_kfold) <- names(ImmuneResponse.test.pre)
  summary_test_views_alg_kfold[[input_name]] <- summary_test_alg_kfold
}
save(summary_test_views_alg_kfold, file = "output/Predictions_validation/SKCM/Labels_predictions_test_pre_treatment_Riaz_Auslander_Hugo_Pan_Cancer.RData")

