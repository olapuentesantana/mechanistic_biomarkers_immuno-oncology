# #########################################################################################################
# Script to validate the model on external datasets obtained from literature -->
# Data types: pathways; immune cells; 
#             proteins; transcript;
#             L-R pairs; TFs; 
#             sTIL;cytokine pairs;
#
# PanCancer: BLCA, BRCA, CESC, CRC, GBM, HNSC, KIRC, KIRP,
#            LIHC, LUAD, LUSC, OV, PAAD, PRAD, SKCM, STAD, THCA, UCEC
#
# Algorithms: L21, Elastic_Net
# #########################################################################################################

# ****************
# working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ****************
# packages
library(ggplot2)
library(reshape)
library(dplyr)
library(grid) 
library(gridExtra)
library("ggsignif")
library(ggpubr)
library(pdist)

# ****************
# scripts
source("../R/BEMKL/predict_BEMKL.R")
source("../R/GLMs/predict_GLMs.R")
source("../R/L21_regularization_EN/predict_L21.R")
source("../R/predict_Multi_Task_EN.R")

# ****************
# Select cancer type
## no filter
load("./pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)
## filter spat
#load("./pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_Spat.RData")
#PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)
## filter prot
#load("./pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_prot.RData")
#PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)

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

# all
view_combinations <- list(views[c(1)], views[c(3)], views[c(1,3)], views[8], views[5])

# ****************
# Initialize variables
input_algorithm <- c("Multi_Task_EN")
analysis <- c("all", "cor")

# ****************
# Cancer types
# PanCancer.names <- c("GBM", "SKCM", "STAD")
PanCancer.names <- c("SKCM")
all.predictions.test <- list()

# ****************
# External datasets

# Validation sets DataViews data: pre-treatment
load("../data/Validation/All_DataViews_test_pre.RData")

# Validation sets ImmuneResponse data: pre-treatment
load("../data/Validation/All_Labels_test_pre.RData")

Datasets.names <- names(All.DataViews.test)
names(Datasets.names) <- c("SKCM", "SKCM", "SKCM", "STAD", "SKCM", "SKCM", "GBM", "SKCM", "SKCM")

#for (Cancer in PanCancer.names){
Cancer = "SKCM"
    cat("Cancer: ", Cancer, "\n")
  
    # Training sets Dataviews data: TCGA 
    # load(paste0("../data/PanCancer/",Cancer,"/new/DataViews_no_filter_", Cancer,".RData"))
    # load(paste0("../data/Federica_presentation_colab/DataViews_no_filter_",Cancer,".RData"))
    load(paste0("../data/BIM_cluster_presentation/","DataViews_no_filter_", Cancer,".RData"))
    
    # Output model from training
    # file <- dir(paste0("../output/Federica_presentation_colab/"), full.names = T, recursive = F)
    file <- dir(path = paste0("../output/BIM_cluster_presentation"), pattern = "all_cv_res_", full.names = T, recursive = F)
    
    # check data type is available for validation data
    filter.cancer <- which(names(Datasets.names) == Cancer)
    filter.cancer <- 6 
    predictions.test <- lapply(Datasets.names[filter.cancer], function(dataset){
      
      cat("Validation data set:",dataset, "\n")
      
      summary_analysis <- lapply(analysis, function(anal){
        
        cat("Analysis:",anal, " tasks \n")
        
        summary_view <- lapply(1:length(view_combinations), function(view){
          
          input_name <- paste(names(view_combinations[[view]]), collapse ="_")
          cat("View:",input_name, "\n")
          
          # Load file
          which_file <- grep(pattern = paste0("_with_", anal,"_tasks_", input_name,".RData"), file, fixed = T)
          load(file[which_file])
          
          # combo info
          view_combination <- view_combinations[[view]]
          
          summary_alg <- lapply(input_algorithm, function(alg){
            
            cat("Model -->", alg, "\n")
          
            if (alg %in% c("Elastic_Net")){
          
              # Generalized linear models:
              pred_alg <- predict_GLMs(DataViews.train = DataViews.no_filter, DataViews.test = All.DataViews.test[[dataset]],
                                       Label.test = All.Labels.test[[dataset]], View = input_name, Trained.model = all_cv_res, 
                                       Algorithm = alg, View.info = view_combination) 
              
            }else if (alg %in% c("L21")){
              
              # L21 regularization:
              pred_alg <- predict_L21(DataViews.train = DataViews.no_filter, DataViews.test = All.DataViews.test[[dataset]],
                                      Label.test = All.Labels.test[[dataset]], View = input_name, Trained.model = all_cv_res, 
                                      Algorithm = alg, View.info = view_combination)
              
            }else if (alg %in% c("BEMKL")){
            
              # BEMKL:
              # pred_alg <- predict_BEMKL(DataViews.train = DataViews.no_filter, DataViews.test = All.DataViews.test[[dataset]],
              #                           Label.test = All.Labels.test[[dataset]], View = input_name, Trained.model = all_cv_res, 
              #                           Algorithm = alg, View.info = view_combination)
            }else if (alg %in% c("Multi_Task_EN")){
            
            pred_alg <- predict_Multi_Task_EN(DataViews.train = DataViews.no_filter, DataViews.test = All.DataViews.test[[dataset]],
                                              Label.test = All.Labels.test[[dataset]], View = input_name, Trained.model = all_cv_res,
                                              Algorithm = alg, View.info = view_combination)
          }
          return(pred_alg)
            
        })
          
        names(summary_alg) <- input_algorithm
        return(summary_alg)
          
      })
        names(summary_view) <- do.call(c, lapply(1:length(view_combinations), function(X){
                              name_view <- paste(names(view_combinations[[X]]), collapse = "_"); return(name_view)}))
      return(summary_view)
        
    })
    names(summary_analysis) <- analysis
    return(summary_analysis)
    
  })
  names(predictions.test) <- Datasets.names[filter.cancer]
  all.predictions.test[[Cancer]] <- predictions.test
    
#}
    
save(all.predictions.test, file = "../output/BIM_cluster_presentation/Labels_predictions_test_pre_treatment_all_datasets.RData")


#--------------------------------------------------------------------
# Plot ROC curve and barplot with the area under the ROC curve.
# No combination of test datasets. Each is evaluated separately.
#--------------------------------------------------------------------

# ****************
# load predictions
load("../output/BIM_cluster_presentation/Labels_predictions_test_pre_treatment_all_datasets.RData")
predictions_immune_response <- all.predictions.test

# ****************
# initialize variables
# Datasets.names <-  dir("../data/Validation/Francesca", full.names = F, recursive = F)
# names(Datasets.names) <- c("SKCM", "SKCM", "SKCM", "STAD", "SKCM", "SKCM", "GBM")

# ****************
# Colors for visualization
# colors.DataType <- toupper(c("#8c9f3e","#b65cbf","#4eac7c","#c95574","#747fca","#ca743e"))
# names(colors.DataType) <- names(predictions_immune_response$GBM$Zhao$all)


colors.DataType <- toupper(c("#c85878","#59aa54","#b65ebd","#999a3e","#6f7dcb"))# "#cc5136", "#45b0a4","#c88645"))
names(colors.DataType) <- names(predictions_immune_response$SKCM$comb.Gide_Auslander$all)

colors.tasks <- toupper(c("#b47645","#ad58c5","#6cb643","#d24787","#52ad7b","#cf4740", "#4bafd0",
                          "#dc7b31","#6776cb","#c1ad46","#b975b1","#6c7b33","#c26671"))
names(colors.tasks) <- names(predictions_immune_response$SKCM$comb.Gide_Auslander$all$pathways$Multi_Task_EN$pred)

colors.algorithm <- toupper(c("#ff7433"))#,"#853760"))
names(colors.algorithm) <- names(predictions_immune_response$SKCM$comb.Gide_Auslander$all$pathways)

alpha.algorithm <- c(1)
names(alpha.algorithm) <- names(predictions_immune_response$SKCM$comb.Gide_Auslander$all$pathways)

alpha.analysis_alg <- c(0.93,0.4)
names(alpha.analysis_alg) <- names(predictions_immune_response$SKCM$comb.Gide_Auslander)

# Check labels, we need to be consistent --> label.ordering = c("NR", "R") #
predictions_immune_response <- lapply(PanCancer.names, function(Cancer){
  
  # check data type is available for validation data
  filter.cancer <- which(names(Datasets.names) == Cancer) 
  filter.cancer <- 9
  predictions_immune_response[[Cancer]] <- lapply(Datasets.names[filter.cancer], function(dataset){
    predictions_immune_response[[Cancer]][[dataset]] <- lapply(names(alpha.analysis_alg), function(anal){
     predictions_immune_response[[Cancer]][[dataset]][[anal]] <- lapply(names(colors.DataType), function(view){
       predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]] <- lapply(names(alpha.algorithm), function(algorithm){ 
          predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]][[algorithm]]$lab <- lapply(names(colors.tasks), function(task){
            
              df.label <- predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]][[algorithm]]$lab[[task]][[1]]
              
              if (dataset %in% c("Kim","Gide","Riaz","Liu", "comb.Gide_Auslander")){
                df.label <- gsub(" ","", df.label)
                df.label <- gsub("CR|PR|MR|PRCR","R", df.label)
                df.label <- gsub("SD|PD|PD ","NR", df.label)
                
              } else if (dataset == "Hugo"){
                df.label <- gsub("Complete Response|Partial Response", "R", df.label)
                df.label <- gsub("Progressive Disease","NR", df.label)
                
              } else if (dataset == "Zhao"){
                df.label <- gsub("Yes|Yes, > 6 months stable","R", df.label)
                df.label <- gsub("No","NR", df.label)
                
              } else if (dataset == "comb.Hugo_Liu_Riaz"){
                df.label <- gsub(" ","", df.label)
                df.label <- gsub("CR|PR|MR|PRCR|CompleteResponse|PartialResponse","R", df.label)
                df.label <- gsub("SD|PD|PD|ProgressiveDisease","NR", df.label)
                
              } 
              predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]][[algorithm]]$lab[[task]][[1]] <- df.label
              return(predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]][[algorithm]]$lab[[task]])
            })
          names(predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]][[algorithm]]$lab) <- names(colors.tasks)
          return(predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]][[algorithm]])
          })
       names(predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]]) <- names(alpha.algorithm)
       return(predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]])
      })
     names(predictions_immune_response[[Cancer]][[dataset]][[anal]]) <- names(colors.DataType)
     return(predictions_immune_response[[Cancer]][[dataset]][[anal]])
    })
    names(predictions_immune_response[[Cancer]][[dataset]]) <- names(alpha.analysis_alg)
    return(predictions_immune_response[[Cancer]][[dataset]])
  })
  names(predictions_immune_response[[Cancer]]) <- Datasets.names[filter.cancer]
  return(predictions_immune_response[[Cancer]])
})
names(predictions_immune_response) <- PanCancer.names

ROC_info <- list()
# filter.analysis.alg <- list(all = c("L21", "Elastic_Net"), all_top = "L21")
filter.analysis.alg <- list(all = c("Multi_Task_EN"), cor = c("Multi_Task_EN"))
filter.analysis.task <- list(all = names(colors.tasks), cor = c("CYT", "IS","RohIS","IS_Davoli","IFny","ExpandedImmune", "T_cell_inflamed")) #, all_top = names(colors.tasks)[-c(3,4,8,12,13)])
cv_model <- names(predictions_immune_response$SKCM$comb.Gide_Auslander$all$pathways_immunecells$Multi_Task_EN$pred$CYT)
# Collect predictions #
ROC_info <-  lapply(PanCancer.names, function(Cancer){
  
  cat("Cancer:",Cancer, "\n")
  
  # check data type is available for validation data
  filter.cancer <- which(names(Datasets.names) == Cancer) 
  filter.cancer <- 9
  ROC_info <-  lapply(Datasets.names[filter.cancer], function(dataset){
    
    cat("Dataset:",dataset, "\n")
    
    ROC_info <-  lapply(names(alpha.analysis_alg), function(anal){
      
      cat("Analysis:",anal, "\n")
      
      ROC_info <-  lapply(names(colors.DataType), function(view){
        
        cat("View:",view, "\n")
        
        ROC_info <-  lapply(filter.analysis.alg[[anal]], function(algorithm){ 
          
          cat("Algorithm:",algorithm, "\n")
          
          ROC_info <- lapply(c(filter.analysis.task[[anal]],"common_mean","common_median"), function(task){
            
            cat("Task:",task, "\n")
            
            df <- predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]][[algorithm]]
            
            # if (task == "common_all"){
              
             df$pred[["common_mean"]][["min.mse"]][[1]] <- matrix(0,nrow(df$pred$CYT[["min.mse"]][[1]]), ncol = 100)
             df$pred[["common_mean"]][["1se.mse"]][[1]] <- matrix(0,nrow(df$pred$CYT[["1se.mse"]][[1]]), ncol = 100)
             
             df$lab[["common_mean"]][[1]] <- df$lab$CYT[[1]]
             
             df$pred[["common_median"]][["min.mse"]][[1]] <- matrix(0,nrow(df$pred$CYT[["min.mse"]][[1]]), ncol = 100)
             df$pred[["common_median"]][["1se.mse"]][[1]] <- matrix(0,nrow(df$pred$CYT[["1se.mse"]][[1]]), ncol = 100)
             
             df$lab[["common_median"]][["min.mse"]][[1]] <- df$lab$CYT[[1]]

            if (anal == "cor"){
              for (j in 1:100) {
                for(i in 1:nrow(df$pred$CYT$min.mse[[1]])){
                  
                  df$pred[["common_mean"]]$min.mse[[1]][i,j] <- mean(c(df$pred$CYT$min.mse[[1]][i,j], df$pred$IS$min.mse[[1]][i,j],df$pred$RohIS$min.mse[[1]][i,j], 
                                                         df$pred$IS_Davoli$min.mse[[1]][i,j], df$pred$IFny$min.mse[[1]][i,j],df$pred$ExpandedImmune$min.mse[[1]][i,j], 
                                                         df$pred$T_cell_inflamed$min.mse[[1]][i,j]))
                  
                  df$pred[["common_mean"]]$`1se.mse`[[1]][i,j] <- mean(c(df$pred$CYT$`1se.mse`[[1]][i,j], df$pred$IS$`1se.mse`[[1]][i,j],df$pred$RohIS$`1se.mse`[[1]][i,j], 
                                                                       df$pred$IS_Davoli$`1se.mse`[[1]][i,j], df$pred$IFny$`1se.mse`[[1]][i,j],df$pred$ExpandedImmune$`1se.mse`[[1]][i,j], 
                                                                       df$pred$T_cell_inflamed$`1se.mse`[[1]][i,j]))
                  
                  df$pred[["common_median"]]$min.mse[[1]][i,j] <- median(c(df$pred$CYT$min.mse[[1]][i,j], df$pred$IS$min.mse[[1]][i,j],df$pred$RohIS$min.mse[[1]][i,j], 
                                                                         df$pred$IS_Davoli$min.mse[[1]][i,j], df$pred$IFny$min.mse[[1]][i,j],df$pred$ExpandedImmune$min.mse[[1]][i,j], 
                                                                         df$pred$T_cell_inflamed$min.mse[[1]][i,j]))
                  
                  df$pred[["common_median"]]$`1se.mse`[[1]][i,j] <-  median(c(df$pred$CYT$`1se.mse`[[1]][i,j], df$pred$IS$`1se.mse`[[1]][i,j],df$pred$RohIS$`1se.mse`[[1]][i,j], 
                                                                            df$pred$IS_Davoli$`1se.mse`[[1]][i,j], df$pred$IFny$`1se.mse`[[1]][i,j],df$pred$ExpandedImmune$`1se.mse`[[1]][i,j], 
                                                                            df$pred$T_cell_inflamed$`1se.mse`[[1]][i,j]))
                  
                  
                }
              }
             rownames(df$pred[["common_mean"]]$min.mse[[1]]) <- rownames(df$lab$CYT[[1]])
             rownames(df$pred[["common_mean"]]$`1se.mse`[[1]]) <- rownames(df$lab$CYT[[1]])
             rownames(df$pred[["common_median"]]$min.mse[[1]]) <- rownames(df$lab$CYT[[1]])
             rownames(df$pred[["common_median"]]$`1se.mse`[[1]]) <- rownames(df$lab$CYT[[1]])
             
            }else{
              
              for (j in 1:100) {
                for(i in 1:nrow(df$pred$CYT$min.mse[[1]])){
                  
                  df$pred[["common_mean"]]$min.mse[[1]][i,j] <- mean(c(df$pred$CYT$min.mse[[1]][i,j],df$pred$IS$min.mse[[1]][i,j], df$pred$IPS$min.mse[[1]][i,j],
                                                                       df$pred$IMPRES$min.mse[[1]][i,j], df$pred$RohIS$min.mse[[1]][i,j], df$pred$chemokine$min.mse[[1]][i,j],
                                                                       df$pred$IS_Davoli$min.mse[[1]][i,j], df$pred$Proliferation$min.mse[[1]][i,j], df$pred$IFny$min.mse[[1]][i,j],
                                                                       df$pred$ExpandedImmune$min.mse[[1]][i,j], 
                                                                       df$pred$T_cell_inflamed$min.mse[[1]][i,j], df$pred$TIDE$min.mse[[1]][i,j], df$pred$MSI$min.mse[[1]][i,j]))
                  
                  df$pred[["common_mean"]]$`1se.mse`[[1]][i,j] <- mean(c(df$pred$CYT$`1se.mse`[[1]][i,j],df$pred$IS$`1se.mse`[[1]][i,j], df$pred$IPS$`1se.mse`[[1]][i,j],
                                                                         df$pred$IMPRES$`1se.mse`[[1]][i,j], df$pred$RohIS$`1se.mse`[[1]][i,j], df$pred$chemokine$`1se.mse`[[1]][i,j],
                                                                         df$pred$IS_Davoli$`1se.mse`[[1]][i,j], df$pred$Proliferation$`1se.mse`[[1]][i,j], df$pred$IFny$`1se.mse`[[1]][i,j],
                                                                         df$pred$ExpandedImmune$`1se.mse`[[1]][i,j], 
                                                                         df$pred$T_cell_inflamed$`1se.mse`[[1]][i,j], df$pred$TIDE$`1se.mse`[[1]][i,j], df$pred$MSI$`1se.mse`[[1]][i,j]))
                  
                  df$pred[["common_median"]]$min.mse[[1]][i,j] <- median(c(df$pred$CYT$min.mse[[1]][i,j],df$pred$IS$min.mse[[1]][i,j], df$pred$IPS$min.mse[[1]][i,j],
                                                                           df$pred$IMPRES$min.mse[[1]][i,j], df$pred$RohIS$min.mse[[1]][i,j], df$pred$chemokine$min.mse[[1]][i,j],
                                                                           df$pred$IS_Davoli$min.mse[[1]][i,j], df$pred$Proliferation$min.mse[[1]][i,j], df$pred$IFny$min.mse[[1]][i,j],
                                                                           df$pred$ExpandedImmune$min.mse[[1]][i,j], 
                                                                           df$pred$T_cell_inflamed$min.mse[[1]][i,j], df$pred$TIDE$min.mse[[1]][i,j], df$pred$MSI$min.mse[[1]][i,j]))
                  
                  df$pred[["common_median"]]$`1se.mse`[[1]][i,j] <-  median(c(df$pred$CYT$`1se.mse`[[1]][i,j],df$pred$IS$`1se.mse`[[1]][i,j], df$pred$IPS$`1se.mse`[[1]][i,j],
                                                                             df$pred$IMPRES$`1se.mse`[[1]][i,j], df$pred$RohIS$`1se.mse`[[1]][i,j], df$pred$chemokine$`1se.mse`[[1]][i,j],
                                                                             df$pred$IS_Davoli$`1se.mse`[[1]][i,j], df$pred$Proliferation$`1se.mse`[[1]][i,j], df$pred$IFny$`1se.mse`[[1]][i,j],
                                                                             df$pred$ExpandedImmune$`1se.mse`[[1]][i,j], 
                                                                             df$pred$T_cell_inflamed$`1se.mse`[[1]][i,j], df$pred$TIDE$`1se.mse`[[1]][i,j], df$pred$MSI$`1se.mse`[[1]][i,j]))
                  
                  
                }
              }
              rownames(df$pred[["common_mean"]]$min.mse[[1]]) <- rownames(df$lab$CYT[[1]])
              rownames(df$pred[["common_mean"]]$`1se.mse`[[1]]) <- rownames(df$lab$CYT[[1]])
              rownames(df$pred[["common_median"]]$min.mse[[1]]) <- rownames(df$lab$CYT[[1]])
              rownames(df$pred[["common_median"]]$`1se.mse`[[1]]) <- rownames(df$lab$CYT[[1]])
              
              
              
              
              
            }
            # }else if (task == "common_top"){
            #   
            #   df$pred[["common_top"]][[1]] <- matrix(0,nrow(df$pred$CYT[[1]]), ncol = 100)
            #   df$lab[["common_top"]][[1]] <- df$lab$CYT[[1]]
            #   
            #   for (j in 1:100) {
            #     for(i in 1:nrow(df$pred$CYT[[1]])){
            #       
            #       df$pred[["common_top"]][[1]][i,j]<- median(c(df$pred$CYT[[1]][i,j], df$pred$IS[[1]][i,j], 
            #                                             df$pred$RohIS[[1]][i,j], df$pred$chemokine[[1]][i,j], 
            #                                             df$pred$IS_Davoli[[1]][i,j], df$pred$IFny[[1]][i,j],
            #                                             df$pred$ExpandedImmune[[1]][i,j], df$pred$T_cell_inflamed[[1]][i,j]))
            #     }
            #   }
            #   rownames(df$pred[["common_top"]][[1]]) <- rownames(df$lab$CYT[[1]])
            #   
            # }
            
            pred.min.mse <- ROCR::prediction(df$pred[[task]]$min.mse[[1]],df$lab[["CYT"]][[1]], label.ordering = c("NR", "R"))
            pred.1se.mse <- ROCR::prediction(df$pred[[task]]$`1se.mse`[[1]],df$lab[["CYT"]][[1]], label.ordering = c("NR", "R"))
            
            perf.min.mse <- ROCR::performance(pred.min.mse,"tpr","fpr")
            perf.1se.mse <- ROCR::performance(pred.1se.mse,"tpr","fpr")
            
            AUC.min.mse <- unlist(ROCR::performance(pred.min.mse, "auc")@y.values)
            AUC.1se.mse <- unlist(ROCR::performance(pred.1se.mse, "auc")@y.values)
            
            data_ROC <- list(perf.min.mse = perf.min.mse, perf.1se.mse = perf.1se.mse)
            Barplot <- list(AUC.min.mse = AUC.min.mse, AUC.1se.mse = AUC.1se.mse)
            return(list(Curve = data_ROC, Barplot = Barplot))
          })
          names(ROC_info) <- c(filter.analysis.task[[anal]],"common_mean","common_median")
          return(ROC_info)
        })
        names(ROC_info) <- filter.analysis.alg[[anal]]
        return(ROC_info)
      })
      names(ROC_info) <- names(colors.DataType)
      return(ROC_info)
    })
    names(ROC_info) <- names(alpha.analysis_alg)
    return(ROC_info)
  })
  names(ROC_info) <- Datasets.names[filter.cancer]
  return(ROC_info)
})
names(ROC_info) <- PanCancer.names


#--------------------------------------------------------------------
# Barplot with the area under the ROC curve.
# No combination of test datasets. Each is evaluated separately.
#--------------------------------------------------------------------

AUC.mean.sd <- do.call(rbind, lapply(PanCancer.names, function(Cancer){
  cat("Cancer:",Cancer, "\n")
  # check data type is available for validation data
  filter.cancer <- which(names(Datasets.names) == Cancer) 
  filter.cancer <- 9
  AUC.mean.sd <- do.call(rbind, lapply(Datasets.names[filter.cancer], function(dataset){
    cat("Dataset:",dataset, "\n")
    AUC.mean.sd <- do.call(rbind, lapply(names(alpha.analysis_alg), function(anal){
      cat("Analysis:",anal, "\n")
      AUC.mean.sd <- do.call(rbind, lapply(names(colors.DataType), function(view){
        cat("View:",view, "\n")
        AUC.mean.sd <- do.call(rbind, lapply(filter.analysis.alg[[anal]], function(algorithm){ 
          cat("Algorithm:",algorithm, "\n")
          AUC.mean.sd <- do.call(rbind, lapply(c(filter.analysis.task[[anal]],"common_mean","common_median"), function(task){

            AUC.data <- data.frame(Cancer = Cancer,
                                   Dataset = dataset,
                                   Analysis = anal,
                                   View = view, 
                                   Model = algorithm,
                                   cv_model = rep(c("min.mse", "1se.mse"), each = 100),
                                   Task = task,
                                   iteration = seq(1,100),
                                   AUC = c(as.numeric(ROC_info[[Cancer]][[dataset]][[anal]][[view]][[algorithm]][[task]]$Barplot[["AUC.min.mse"]]),
                                           as.numeric(ROC_info[[Cancer]][[dataset]][[anal]][[view]][[algorithm]][[task]]$Barplot[["AUC.1se.mse"]])))
            
            AUC.mean.sd <- do.call(data.frame, aggregate(AUC ~ Cancer + Dataset + Analysis + View + Model + cv_model + Task, 
                                               data = AUC.data, FUN = function(x) c(median = median(x), sd = sd(x))))
          }))
          return(AUC.mean.sd)
        }))
        return(AUC.mean.sd)
      }))
      return(AUC.mean.sd)
    }))
    return(AUC.mean.sd)
  }))
  return(AUC.mean.sd)
}))

# AUC.mean.sd$Analysis_alg <- factor(paste0(AUC.mean.sd$Model,"_", AUC.mean.sd$Analysis))
# alpha.analysis_alg <- c(1)
# names(alpha.analysis_alg) <- levels(AUC.mean.sd$Analysis_alg)

# Subset # 
CancerType <- "SKCM"
tmp.barplot <- subset(AUC.mean.sd, Cancer == CancerType)
tmp.barplot$Task <- factor(tmp.barplot$Task)
tmp.barplot$cv_model <- factor(tmp.barplot$cv_model, levels = c("min.mse", "1se.mse"))

#
# CancerType <- "STAD"
# AlgorithmType <- "L21"
# Tasks.keep <- filter.analysis.task$all_top
tmp.barplot <- subset(AUC.mean.sd, cv_model == "1se.mse")
tmp.barplot$Task <- factor(tmp.barplot$Task, levels = c("CYT","IS","IPS","IMPRES","RohIS", "chemokine","IS_Davoli","Proliferation",
                           "IFny","ExpandedImmune", "T_cell_inflamed","TIDE","MSI", "common_mean","common_median"))
tmp.barplot$Analysis <- factor(tmp.barplot$Analysis)

tmp.barplot <- subset(AUC.mean.sd, Task %in% c("common_mean","common_median") & cv_model == "1se.mse")
tmp.barplot$Analysis <- factor(tmp.barplot$Analysis)
tmp.barplot$Task <- factor(tmp.barplot$Task)

# Barplot #
ggplot2::ggplot(tmp.barplot, aes(x=View, y=round(AUC.median,2), fill=Task,  
                                 colour = Task, alpha = Analysis)) +
  ggplot2::geom_bar(stat="identity", position = position_dodge()) +
  ggplot2::scale_fill_manual(name = "Task", 
                             labels = levels(tmp.barplot$Task),
                              values = c("red","green")) +
  ggplot2::scale_color_manual(name = "Task", 
                     labels = levels(tmp.barplot$Task),
                     values =  c("red","green")) + 
  ggplot2::scale_alpha_manual(name = "Analysis",
                     labels = levels(tmp.barplot$Analysis),
                     values = c(1,0.5))  +
  ggplot2::theme_bw() +
  #ggplot2::facet_grid(Dataset ~ .) +
  ggplot2::ylim(0, 1) +
  ggplot2::ylab("AUC") +
  ggplot2::geom_errorbar(aes(ymin = AUC.median - AUC.sd, ymax = AUC.median + AUC.sd), 
                         width=.5, color="gray", position = position_dodge(0.9)) +
  ggplot2::geom_text(aes(label= round(AUC.median,2)),stat = "identity", color="black", size = 3, angle = 90,
                     position = position_dodge(0.9)) +
  scale_x_discrete(labels = c("immunecells"= "ImmuneCells",
                              "LRpairs"= "L-R pairs",
                              "CYTOKINEpairs"= "Cytokine pairs",
                              "pathways"="Pathways",
                              "pathways_immunecells"="Pathways \n + \n Immunecells",
                              "TFs" = "TFs",
                              "transcript" = "Transcriptomics")) +  
  theme(axis.text.x = element_text(size=10,face="bold", angle = 0, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,face="bold"), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_line(size = 1),
        axis.title.y = element_text(size=10,face="bold",vjust = 0.9), strip.background = element_blank(), panel.spacing.y =  unit(1, "lines"),
        legend.position = "right", legend.direction = "vertical",
        legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
        panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(size=0.5, colour = "black")) +
# ggtitle(paste0("Mechanistic signatures performance across ",CancerType, " datasets (",AnalysisType,")"))
  ggtitle(paste0("Mechanistic signatures AUC values for each ",CancerType, " datasets"))

ggsave(paste0("../figures/BIM_cluster_presentation/","Barplot_with_aucs_on_SKCM_datasets.pdf"), width = 12, height = 12)


#--------------------------------------------------------------------
# Plot ROC curve
# No combination of test datasets. Each is evaluated separately.
#--------------------------------------------------------------------
library(ROCR)

source("../R/computation.gold.standards.R")
list.gold.standards <- c("CYT", "PD1", "PDL1","CTLA4")
gold.standards <- computation.gold.standards(RNA.tpm = 2^All.DataViews.test$comb.Gide_Auslander$transcript -1, list_gold_standards = list.gold.standards)

df <- predictions_immune_response[["SKCM"]][["comb.Gide_Auslander"]][["cor"]][["pathways"]][["Multi_Task_EN"]]$lab[[1]][[1]]
gold.standards.curve <- list()
gold.standards.curve<- do.call(c,lapply(list.gold.standards, function(X){

  pred <- ROCR::prediction(gold.standards[[X]],df[,1], label.ordering = c("NR", "R"))
  perf <- ROCR::performance(pred,"tpr","fpr")
  AUC <- unlist(ROCR::performance(pred, "auc")@y.values)
  
  gold.standards.curve[[X]]$Curve <- perf
  gold.standards.curve[[X]]$Barplot <- AUC
  return(gold.standards.curve)
}))

# lapply(PanCancer.names, function(Cancer){
#   cat("Cancer:",Cancer, "\n")
#   # check data type is available for validation data
#   filter.cancer <- which(names(Datasets.names) == Cancer) 
#   lapply(Datasets.names[filter.cancer], function(dataset){
#     cat("Dataset:",dataset, "\n")
#     lapply(names(alpha.analysis_alg), function(anal){
#       cat("Analysis:",anal, "\n")
#       lapply(names(colors.DataType), function(view){
#         cat("View:",view, "\n")
#         lapply(filter.analysis.alg[[anal]], function(algorithm){ 
#           cat("Algorithm:",algorithm, "\n")
#           lapply(c(filter.analysis.task[[anal]],"common_all","common_top"), function(task){
            
            # pathways_immunecells - common 
            dataset = "comb.Gide_Auslander"
            pdf("../figures/BIM_cluster_presentation/ROC_curves_predictions_on_SKCM_comb_Gide_Auslander_validation_dataset_pathways_immunecells.pdf", width = 12, height = 12)
            par(cex.axis=2.2, mar = c(5, 5, 5, 5))
            ROCR::plot(ROC_info[[Cancer]][[dataset]][["cor"]][["pathways_immunecells"]][[algorithm]][["common_median"]]$Curve$perf.1se.mse,
                       avg = "threshold", col = colors.DataType["CYTOKINEpairs"], lwd = 8, type = "S", cex.lab=2.2, ylab="True Positive Rate",
                       xlab="False Positive Rate")
            
            # # immunecells - common
            # ROCR::plot(ROC_info[[Cancer]][[dataset]][["cor"]][["immunecells"]][[algorithm]][["common_median"]]$Curve$perf.1se.mse,
            #            avg = "threshold", col = colors.DataType[[2]], lwd = 5, type = "S", add = TRUE)
            # # pathways - common
            # ROCR::plot(ROC_info[[Cancer]][[dataset]][["cor"]][["pathways"]][[algorithm]][["common_median"]]$Curve$perf.1se.mse,
            #            avg = "threshold", col = colors.DataType[[3]], lwd = 5, type = "S", add = TRUE)

            # transcript - common
            ROCR::plot(ROC_info[[Cancer]][[dataset]][["cor"]][["transcript"]][[algorithm]][["common_median"]]$Curve$perf.1se.mse,
                        avg = "threshold", col = colors.DataType["transcript"], lwd = 8, type = "s", add = TRUE)
            
            # gold standard - CYT
            ROCR::plot(gold.standards.curve$CYT$Curve, col = "gray", lwd = 2, type = "S", lty = 1, add = TRUE)
            
            # # gold standard - PD1
            ROCR::plot(gold.standards.curve$PD1$Curve,col = "gray50", lwd = 2, type = "S", lty = 1, add = TRUE)
            
            # # gold standard - PDL1
            # ROCR::plot(gold.standards.curve$PDL1$Curve,col = "chocolate1", lwd = 2, type = "S", lty = 1, add = TRUE)
            
            # # # gold standard - CTLA4
            # ROCR::plot(gold.standards.curve$CTLA4$Curve,col = "gray50", lwd = 2, type = "S", lty = 1, add = TRUE)
            
            abline(a=0, b=1, lty = 3, lwd = 1, col = "antiquewhite4")
            
            
            #title(paste0("Dataset: ", dataset, "(PD-1, PD-1 & CTLA-4)") )

            # AUC_median <- round(subset(AUC.mean.sd, Cancer == Cancer & Dataset == dataset & View == view &
            #                              Analysis == anal & Model == algorithm & Task == task, AUC.median),2)
            
            legend(x = 0.49,y = 0.42, legend = c(paste0("Pathways \n       +          ","(AUC=",
                                                       round(subset(AUC.mean.sd, Task == "common_median" & Analysis == "cor" & cv_model == "1se.mse" &
                                                                             View == "pathways_immunecells")$AUC.median,3),")"," \nImmuneCells"),
                                                 # c(paste0("Cytokine pairs", "(AUC=",
                                                 #        round(subset(AUC.mean.sd, Task == "common_median" & Analysis == "cor" & cv_model == "1se.mse" &
                                                 #                       View == "CYTOKINEpairs")$AUC.median,3),")"),
                                                 paste0("Transcriptomics", "(AUC=",
                                                 round(subset(AUC.mean.sd, Task == "common_median" & Analysis == "cor" & cv_model == "1se.mse" &
                                                                View == "transcript")$AUC.median,3),")"), 
                                                 paste0("Cytolytic Ativity", "(AUC=", round(gold.standards.curve[["CYT"]][["Barplot"]],3),")"),
                                                        
                                                 paste0("PD1", "(AUC=", round(gold.standards.curve[["PD1"]][["Barplot"]],3),")")),
                                                        
                   col = c(colors.DataType["CYTOKINEpairs"],colors.DataType["transcript"],"gray","gray50"),
                   lty = c(1,1,1,1), lwd = c(8,8,2,2), cex = 1.8, bty = "n")
            

            dev.off()
            
#           })
#         })
#       })
#     })
#   })
# })
###################
# Patient score
###################
            
All.scores <- do.call(rbind, lapply(PanCancer.names, function(Cancer){
  cat("Cancer:",Cancer, "\n")
  # check data type is available for validation data
  filter.cancer <- which(names(Datasets.names) == Cancer) 
  filter.cancer <- 9
  All.scores <- do.call(rbind, lapply(Datasets.names[filter.cancer], function(dataset){
    cat("Dataset:",dataset, "\n")
    All.scores <- do.call(rbind, lapply(names(alpha.analysis_alg), function(anal){
      cat("Analysis:",anal, "\n")
      All.scores <- do.call(rbind, lapply(names(colors.DataType), function(view){
        cat("View:",view, "\n")
        All.scores <- do.call(rbind, lapply(filter.analysis.alg[[anal]], function(algorithm){ 
          cat("Algorithm:",algorithm, "\n")
          All.scores <- do.call(rbind, lapply(filter.analysis.task[[anal]], function(task){
            
            df <- predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]][[algorithm]][["pred"]][[task]][["1se.mse"]][[view]]
            label <- predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]][[algorithm]][["lab"]][[task]][[view]]
            
            All.scores <- do.call(rbind, lapply(rownames(df), function(patient){
            
              df.patient <- df[patient,]
              label.patient <- label[patient,]
              
              All.scores <- data.frame(Cancer = Cancer,
                                   Dataset = dataset,
                                   Analysis = anal,
                                   View = view, 
                                   Model = algorithm,
                                   cv_model = "1se.mse",
                                   Task = task,
                                   iteration = seq(1,100),
                                   patient = patient,
                                   label = label.patient,
                                   pred = df.patient)
            }))
            return(All.scores)
          }))
          return(All.scores)
        }))
        return(All.scores)
      }))
      return(All.scores)
    }))
    return(All.scores)
  }))
  return(All.scores)
}))

All.scores.view <- subset(All.scores, View == "pathways_immunecells" & Analysis == "cor")

# All.scores.view$Task <- factor(All.scores.view$Task, levels = c("CYT", "IS", "IPS", "IMPRES", "RohIS", "chemokine", "IS_Davoli", "Proliferation", "IFny", "ExpandedImmune", 
#                                                                   "T_cell_inflamed", "MSI", "TIDE", "common_median"))

All.scores.view$Task <- factor(All.scores.view$Task, levels = c("CYT", "IS","RohIS", "IS_Davoli", "IFny", "ExpandedImmune", 
                                                                "T_cell_inflamed"))
median_pred <- aggregate(pred ~ patient + Task, 
                                     FUN = "median", na.rm = T, data = All.scores.view)

colors.patients <-  predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]][[algorithm]][["lab"]][[task]][[view]][,1]
colors.patients <- gsub("NR","red", colors.patients)
colors.patients <- gsub("R","blue", colors.patients)

# Choose three patients: "SRR7344575", "SRR7344574", "SRR7344567"
median_pred.three.patients <- subset(median_pred, patient %in% c("SRR7344575", "SRR7344574", "SRR7344546"))
median_pred.tasks <- aggregate(pred ~ patient, FUN = "median", na.rm = T, data = median_pred)

order.patients <- median_pred.tasks$patient[order(abs(median_pred.tasks$pred), decreasing = TRUE)]

median_pred$patient <- factor(median_pred$patient, levels = order.patients)
colors.patients <- colors.patients[order.patients]

ggplot2::ggplot(median_pred, aes(x=pred, y=patient, fill = Task, color = Task)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_fill_manual(name = "Task",
                             labels = levels(All.scores.view$Task),
                             values = as.character(colors.tasks)) +
  ggplot2::scale_color_manual(name = "Task",
                              labels = levels(All.scores.view$Task),
                              values = as.character(colors.tasks)) +
  # ggplot2::scale_alpha_manual(name = "Analysis",
  #                             labels = levels(tmp.barplot$Analysis),
  #                             values = c(1,0.5))  +
  # ggplot2::theme_linedraw() +
  #ggplot2::xlim(-1, 1) +
  theme(axis.text.x = element_text(size=18,face="bold", angle = 0, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size=10, colour = colors.patients), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), panel.background = element_rect(fill = "white"),
        legend.position = "right", legend.direction = "vertical", panel.grid.major.y = element_line(colour = "black"),
        legend.text=element_text(size=16), legend.title = element_text(size = 16, face="bold", vjust = 0.5)) +
  labs(x = "prediction", y = "patients")

ggsave(paste0("../figures/BIM_cluster_presentation/patient_score_comb_gide_auslander_pathways_immunecells.pdf"), width = 12, height = 12)

                        
                        
                        
                  
