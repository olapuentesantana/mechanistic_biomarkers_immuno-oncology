# #########################################################################################################
# Script to use lasso to predict immune response: valid for laptop and cluster
# #########################################################################################################

# ****************
# working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ****************
# packages
library(devtools)

# ****************
# scripts
source("./R/cross_validation.R")

# ****************
# Select cancer type
PanCancer.names = c( "BLCA", "BRCA","CESC", "CRC", "LUAD", "LUSC", "PRAD", "SKCM", "STAD", "UCEC")

# ****************
# views
views <- c(PROGENy = 'gaussian', #1
           Protall = 'gaussian', #2
           quanTIseq = 'gaussian', #3
           DoRothEAv1 = 'gaussian', #4
           transcript = 'gaussian', #5
           SpatialTIL = 'gaussian')  #6

# ****************
view_combinations <- list(views[3],
                         views[c(3,6)])
# ****************
# data
load("./data/parameters_4_all.RData")

# ****************
input_algorithm = names(parameters)

for (Cancer in PanCancer.names){

  load(paste0("./data/PanCancer/",Cancer,"/DataViews_", Cancer,"_all_train_with_sTIL.RData"))
  load(paste0("./data/PanCancer/",Cancer,"/ImmuneResponse_", Cancer,"_all_train_with_sTIL_included_IPS_IMPRES.RData"))
  
  # Dataviews SpatialTIL
  dont_keep <- c("Number_of_TIL_Patches", "number_of_clusters", "WCD_mean","CE_mean")
  DataViews$SpatialTIL <- DataViews$SpatialTIL[,-which(colnames(DataViews$SpatialTIL) %in% dont_keep)]
  
  # Consensus define as mean betweeen both responses. Used just for testing
  ImmuneResponse$consensus <- apply(scale(ImmuneResponse[,1:2]),1,mean)

  sapply(1:length(view_combinations), function(X) {
    
    # Input data:
    i <- view_combinations[[X]]
    cat("View:", names(i), "\n")
    
    # Select k fold:
    K = 5
    
    all_cv_res <- lapply(input_algorithm, function(Y){
      
      # Coefficients estimation: optimized by cross validation
      cv_res <- cross_validation(drug_source = ImmuneResponse,
                                 views_source = DataViews,
                                 view_combination = i,
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
    if (length(i) == 1){
      save(all_cv_res, file = paste0("./output/model/",Cancer,"/","all_cv_res_",Cancer, "_train_rand100_withConsensus_",names(i)[1],"_kfold_10",".RData"))
    }else if(length(i) == 2){
      save(all_cv_res, file = paste0("./output/model/",Cancer,"/","all_cv_res_",Cancer, "_train_rand100_withConsensus_",names(i)[1],"_",names(i)[2],"_kfold_10",".RData"))
    }else if(length(i) == 3){
      save(all_cv_res, file = paste0("./output/model/",Cancer,"/","all_cv_res_",Cancer, "_train_rand100_withConsensus_",names(i)[1],"_",names(i)[2],"_",names(i)[3],"_kfold_10",".RData"))
    }else if(length(i) == 4){
      save(all_cv_res, file = paste0("./output/model/",Cancer,"/","all_cv_res_",Cancer, "_train_rand100_withConsensus_",
                                 names(i)[1],"_",names(i)[2],"_",names(i)[3], names(i)[4],"_kfold_10", ".RData"))
    }
      
  })
  
}
  
