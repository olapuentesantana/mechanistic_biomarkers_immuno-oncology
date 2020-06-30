# #########################################################################################################
# Script to predict response on external datasets
# #########################################################################################################

# ****************
# working directory
setwd("/home/feduati/Oscar/PhD_Oscar/new/")

# ****************
# packages
library(devtools, lib.loc= Sys.getenv("R_LIBS"))
library(matrixStats,lib.loc= Sys.getenv("R_LIBS"))
library(pdist,lib.loc= Sys.getenv("R_LIBS"))

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
source("R/PanCancer_draft_v1/MT_GLM_cluster/predict_MT_GLM_cluster.R")
source("R/PanCancer_draft_v1/BEMKL_cluster/predict_BEMKL_cluster.R")

# ****************
# algorithms
input_algorithm <- c("Multi_Task_EN", "BEMKL")
analysis <- c("cor")

# ****************
# External datasets

# Validation sets DataViews data: pre-treatment
load("data/PanCancer_draft_v1/Validation/All_DataViews_test_pre.RData")

# Validation sets ImmuneResponse data: pre-treatment
load("data/PanCancer_draft_v1/Validation/All_Labels_test_pre.RData")

Datasets.names <- names(All.DataViews.test)
names(Datasets.names) <- c("SKCM", "SKCM", "SKCM", "STAD", "SKCM", "SKCM", "GBM", "SKCM", "SKCM")

# ****************
# Training sets Dataviews data: TCGA 
load(paste0("data/PanCancer_draft_v1/",CancerType,"/DataViews_no_filter_", CancerType,".RData"))

# ****************
# Output model from training
file <- dir(path = paste0("output/PanCancer_draft_v1/", CancerType), pattern = "all_cv_res_", full.names = T, recursive = F)

# check data type is available for validation data
filter.cancer <- which(names(Datasets.names) == CancerType)

predictions.test <- lapply(Datasets.names[filter.cancer], function(dataset){

  summary_analysis <- lapply(analysis, function(anal){
    
      # Load file
      which_file <- grep(pattern = paste0("_with_", anal,"_tasks_", View,".RData"), file, fixed = T)
      load(file[which_file])
      
      # combo info
      view_combination <- View_input
      
      summary_alg <- lapply(input_algorithm, function(alg){
        
        if (alg %in% c("BEMKL")){
          
          pred_alg <- predict_BEMKL(DataViews.train = DataViews.no_filter, DataViews.test = All.DataViews.test[[dataset]],
                                    Label.test = All.Labels.test[[dataset]], View = View, Trained.model = all_cv_res[[alg]],
                                    View.info = view_combination)
          
        }else if (alg %in% c("Multi_Task_EN")){
          
          pred_alg <- predict_MT_GLM(DataViews.train = DataViews.no_filter, DataViews.test = All.DataViews.test[[dataset]],
                                     Label.test = All.Labels.test[[dataset]], View = View, Trained.model = all_cv_res[[alg]],
                                     View.info = view_combination)
        }
        return(pred_alg)
      })
      names(summary_alg) <- input_algorithm
      return(summary_alg)
    })
  names(summary_analysis) <- analysis
  return(summary_analysis)
})
names(predictions.test) <- Datasets.names[filter.cancer]

save(predictions.test, file = paste0("output/PanCancer_draft_v1/Validation/", CancerType, "/Labels_predictions_test_pre_treatment_", View,
                                     ".RData"))





