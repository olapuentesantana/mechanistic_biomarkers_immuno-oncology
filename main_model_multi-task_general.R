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
#load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_prot.RData")
#load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_Spat.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)

# ****************
# views
views <- c(Pathways = 'gaussian', #1
           Protall = 'gaussian', #2
           ImmuneCells = 'gaussian', #3
           TFs = 'gaussian', #4
           Transcript = 'gaussian', #5
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
  
  load(paste0("./data/PanCancer_draft_v1/",Cancer,"/DataViews_filter_spat_", Cancer,".RData"))
  load(paste0("./data/PanCancer_draft_v1/",Cancer,"/ImmuneResponse_filter_spat_", Cancer,".RData"))
  
  # Multi-gaussian regression on validation dataset
  ImmuneResponse.filter_spat <- ImmuneResponse.filter_spat[, !colnames(ImmuneResponse.filter_spat) %in% filter_tasks]
  
  # Remove NA values in DataViews sTIL (sample-wise)
  if(anyNA(DataViews.filter_spat$sTIL$Det_Ratio)){
    DataViews.filter_spat$sTIL <- DataViews.filter_spat$sTIL[!is.na(DataViews.filter_spat$sTIL$Det_Ratio),]
    DataViews.filter_spat$ImmuneCells <- DataViews.filter_spat$ImmuneCells[match(rownames(DataViews.filter_spat$sTIL),
                                                                                 rownames(DataViews.filter_spat$ImmuneCells)),]
  }

  # Remove NA values in DataViews LRpairs and CYTOKINE pairs (feature-wise)
  if (anyNA(DataViews.filter_spat$LRpairs)){
    LR_sum <- apply(DataViews.filter_spat$LRpairs, 2, sum)
    DataViews.filter_spat$LRpairs <- DataViews.filter_spat$LRpairs[,!is.na(LR_sum)]
    Cyt_sum <- apply(DataViews.filter_spat$CYTOKINEpairs, 2, sum)
    DataViews.filter_spat$CYTOKINEpairs <- DataViews.filter_spat$CYTOKINEpairs[,!is.na(Cyt_sum)]
  }
  
  sapply(1:length(view_combinations), function(X) {
    
    # Input data:
    i <- view_combinations[[X]]
    cat("View:", names(i), "\n")
    
    # Select k fold:
    K = 5
    
    all_cv_res <- lapply(input_algorithm, function(Y){
      
      # Coefficients estimation: optimized by cross validation
      cv_res <- cross_validation(drug_source = ImmuneResponse.filter_spat,
                                 views_source = DataViews.filter_spat,
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
      save(all_cv_res, file = paste0("./output/PanCancer_draft_v1/","all_cv_res_",Cancer, "_train_rand100_with_cor_tasks_",names(i)[1],".RData"))
    }else if(length(i) == 2){
      save(all_cv_res, file = paste0("./output/PanCancer_draft_v1/","all_cv_res_",Cancer, "_train_rand100_with_cor_tasks_",names(i)[1],"_",names(i)[2],".RData"))
    }
    
  })
  
}
 

   
