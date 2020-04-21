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
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop")

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
source("R/BEMKL/predict_BEMKL.R")
source("R/GLMs/predict_GLMs.R")
source("R/L21_regularization_EN/predict_L21.R")

# ****************
# Select cancer type
## no filter
load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)
## filter spat
#load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_Spat.RData")
#PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)
## filter prot
#load("./analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_prot.RData")
#PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)

# ****************
# views
views <- c(pathways = 'gaussian', #1
           Protall = 'gaussian', #2
           immunecells = 'gaussian', #3
           TFs = 'gaussian', #4
           transcript = 'gaussian', #5
           LRpairs = 'gaussian', #6
           CYTOKINEpairs = 'gaussian')  #7)  

# view_combinations <- list(views[c(1,3)], views[1], views[3], views[4], views[8])
view_combinations <- list(views[c(1,3)], views[1], views[3], views[4], views[7])

# ****************
# Initialize variables
input_algorithm <- c("L21", "Elastic_Net")
analysis <- c("all","all_top")

# ****************
# Cancer types
PanCancer.names <- c("GBM", "SKCM", "STAD")
all.predictions.test <- list()

# ****************
# External datasets
Datasets.names <-  dir("data/Validation/Francesca", full.names = F, recursive = F)
names(Datasets.names) <- c("SKCM", "SKCM", "SKCM", "STAD", "SKCM", "SKCM", "GBM")

for (Cancer in PanCancer.names){
    
  cat("Cancer: ", Cancer, "\n")
  
    # Training sets Dataviews data: TCGA 
    load(paste0("data/PanCancer/",Cancer,"/new/DataViews_no_filter_", Cancer,".RData"))
  
    # Output model from training
    file <- dir(paste0("output/new/",Cancer), full.names = T, recursive = F)
    
    # check data type is available for validation data
    filter.cancer <- which(names(Datasets.names) == Cancer) 
    
    predictions.test <- lapply(Datasets.names[filter.cancer], function(dataset){
      
      cat("Validation data set:",dataset, "\n")
    
      # Validation sets DataViews data: pre-treatment
      load("data/Validation/All_DataViews_test_pre.RData")
    
      # Validation sets ImmuneResponse data: pre-treatment
      load("data/Validation/All_Labels_test_pre.RData")
      
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
    
}
    
save(all.predictions.test, file = "output/new_v2/Validation/Labels_predictions_test_pre_treatment_all_datasets.RData")


#--------------------------------------------------------------------
# Plot ROC curve and barplot with the area under the ROC curve.
# No combination of test datasets. Each is evaluated separately.
#--------------------------------------------------------------------

# ****************
# load predictions
load("output/new_v2/Validation/Labels_predictions_test_pre_treatment_all_datasets.RData")
predictions_immune_response <- all.predictions.test

# ****************
# initialize variables
Datasets.names <-  dir("data/Validation/Francesca", full.names = F, recursive = F)
names(Datasets.names) <- c("SKCM", "SKCM", "SKCM", "STAD", "SKCM", "SKCM", "GBM")

# ****************
# Colors for visualization
colors.DataType <- toupper(c("#8c9f3e","#b65cbf","#4eac7c","#c95574","#747fca","#ca743e"))
names(colors.DataType) <- names(predictions_immune_response$GBM$Zhao$all)

colors.tasks <- toupper(c("#b47645","#ad58c5","#6cb643","#d24787","#52ad7b","#cf4740", "#4bafd0",
                          "#dc7b31","#6776cb","#c1ad46","#b975b1","#6c7b33","#c26671"))
names(colors.tasks) <- names(predictions_immune_response$GBM$Zhao$all$pathways$L21$pred)

colors.algorithm <- toupper(c("#ff7433","#853760"))
names(colors.algorithm) <- names(predictions_immune_response$GBM$Zhao$all$pathways)

alpha.algorithm <- c(1,0.4)
names(alpha.algorithm) <- names(predictions_immune_response$GBM$Zhao$all$pathways)

alpha.analysis_alg <- c(1,0.4)
names(alpha.analysis_alg) <- names(predictions_immune_response$GBM$Zhao)

# Check labels, we need to be consistent --> label.ordering = c("NR", "R") #
predictions_immune_response <- lapply(PanCancer.names, function(Cancer){
  
  # check data type is available for validation data
  filter.cancer <- which(names(Datasets.names) == Cancer) 
  
  predictions_immune_response[[Cancer]] <- lapply(Datasets.names[filter.cancer], function(dataset){
    predictions_immune_response[[Cancer]][[dataset]] <- lapply(names(alpha.analysis_alg), function(anal){
     predictions_immune_response[[Cancer]][[dataset]][[anal]] <- lapply(names(colors.DataType), function(view){
       predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]] <- lapply(names(alpha.algorithm), function(algorithm){ 
          predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]][[algorithm]]$lab <- lapply(names(colors.tasks), function(task){
            
              df.label <- predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]][[algorithm]]$lab[[task]][[1]]
              
              if (dataset %in% c("Kim","Gide","Riaz","Liu")){
                df.label <- gsub(" ","", df.label)
                df.label <- gsub("CR|PR|MR|PRCR","R", df.label)
                df.label <- gsub("SD|PD|PD ","NR", df.label)
              } else if (dataset == "Hugo"){
                df.label <- gsub("Complete Response|Partial Response", "R", df.label)
                df.label <- gsub("Progressive Disease","NR", df.label)
              } else if (dataset == "Zhao"){
                df.label <- gsub("Yes|Yes, > 6 months stable","R", df.label)
                df.label <- gsub("No","NR", df.label)
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
filter.analysis.alg <- list(all = c("L21", "Elastic_Net"), all_top = "L21")
filter.analysis.task <- list(all = names(colors.tasks), all_top = names(colors.tasks)[-c(3,4,8,12,13)])
# Collect predictions #
ROC_info <-  lapply(PanCancer.names, function(Cancer){
  
  cat("Cancer:",Cancer, "\n")
  
  # check data type is available for validation data
  filter.cancer <- which(names(Datasets.names) == Cancer) 
  
  ROC_info <-  lapply(Datasets.names[filter.cancer], function(dataset){
    
    cat("Dataset:",dataset, "\n")
    
    ROC_info <-  lapply(names(alpha.analysis_alg), function(anal){
      
      cat("Analysis:",anal, "\n")
      
      ROC_info <-  lapply(names(colors.DataType), function(view){
        
        cat("View:",view, "\n")
        
        ROC_info <-  lapply(filter.analysis.alg[[anal]], function(algorithm){ 
          
          cat("Algorithm:",algorithm, "\n")
          
          ROC_info <- lapply(filter.analysis.task[[anal]], function(task){
            
            cat("Task:",task, "\n")
            
            df <- predictions_immune_response[[Cancer]][[dataset]][[anal]][[view]][[algorithm]]
            
            pred <- ROCR::prediction(df$pred[[task]][[1]],df$lab[[task]][[1]], label.ordering = c("NR", "R"))
            perf <- ROCR::performance(pred,"tpr","fpr")
            AUC <- unlist(ROCR::performance(pred, "auc")@y.values)
            
            data_ROC <- list(perf)
            Barplot <- list(AUC)
            return(list(Curve = data_ROC, Barplot = Barplot))
          })
          names(ROC_info) <- filter.analysis.task[[anal]]
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
# Plot ROC curve
# No combination of test datasets. Each is evaluated separately.
#--------------------------------------------------------------------
library(pROC)

lapply(PanCancer.names, function(Cancer){
  cat("Cancer:",Cancer, "\n")
  # check data type is available for validation data
  filter.cancer <- which(names(Datasets.names) == Cancer) 
  lapply(Datasets.names[filter.cancer], function(dataset){
    cat("Dataset:",dataset, "\n")
    lapply(names(alpha.analysis_alg), function(anal){
      cat("Analysis:",anal, "\n")
      lapply(names(colors.DataType), function(view){
        cat("View:",view, "\n")
        lapply(filter.analysis.alg[[anal]], function(algorithm){ 
          cat("Algorithm:",algorithm, "\n")
          lapply(filter.analysis.task[[anal]], function(task){

          ROCR::plot(ROC_info[[Cancer]][[dataset]][[anal]][[view]][[algorithm]][[task]]$Curve[[1]],
                     avg = "threshold", col = colors.tasks[[task]], lwd = 2)
          title(paste0("Cancer: ", Cancer, " ,Dataset: ", dataset, " ,Analysis: ", anal, " ,View: ", view, 
                       " ,Algorithm: ", algorithm," ,Task: ",task), cex.main = 0.8)
          AUC_median <- round(subset(AUC.mean.sd, Cancer == Cancer & Dataset == dataset & View == view &
                                            Analysis == anal & Model == algorithm & Task == task, AUC.median),2)
   
          legend("bottomright", legend = c(paste0(view,",AUC=", AUC_median)), col = colors.tasks[[task]],
                 lty = 1, cex = 0.8,  text.width = 0.1, bty = "n")
          })
        })
      })
    })
  })
})









#--------------------------------------------------------------------
# Barplot with the area under the ROC curve.
# No combination of test datasets. Each is evaluated separately.
#--------------------------------------------------------------------

AUC.mean.sd <- do.call(rbind, lapply(PanCancer.names, function(Cancer){
  cat("Cancer:",Cancer, "\n")
  # check data type is available for validation data
  filter.cancer <- which(names(Datasets.names) == Cancer) 
  AUC.mean.sd <- do.call(rbind, lapply(Datasets.names[filter.cancer], function(dataset){
    cat("Dataset:",dataset, "\n")
    AUC.mean.sd <- do.call(rbind, lapply(names(alpha.analysis_alg), function(anal){
      cat("Analysis:",anal, "\n")
      AUC.mean.sd <- do.call(rbind, lapply(names(colors.DataType), function(view){
        cat("View:",view, "\n")
        AUC.mean.sd <- do.call(rbind, lapply(filter.analysis.alg[[anal]], function(algorithm){ 
          cat("Algorithm:",algorithm, "\n")
          AUC.mean.sd <- do.call(rbind, lapply(filter.analysis.task[[anal]], function(task){

            AUC.data <- data.frame(Cancer = rep(Cancer, times = 100),
                                   Dataset = rep(dataset, times = 100),
                                   Analysis = rep(anal, times = 100),
                                   View = rep(view, times = 100), 
                                   Model = rep(algorithm, times = 100),
                                   Task = rep(task, times = 100), 
                                   iteration =  seq(1,100),
                                   AUC = as.numeric(unlist(ROC_info[[Cancer]][[dataset]][[anal]][[view]][[algorithm]][[task]]$Barplot)))
            
            AUC.mean.sd <- do.call(data.frame, aggregate(AUC ~ Cancer + Dataset + Analysis + View + Model + Task, 
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

AUC.mean.sd$Analysis_alg <- factor(paste0(AUC.mean.sd$Model,"_", AUC.mean.sd$Analysis))
alpha.analysis_alg <- c(1,0.6,0.15)
names(alpha.analysis_alg) <- levels(AUC.mean.sd$Analysis_alg)

# Subset # 
CancerType <- "GBM"
tmp.barplot <- subset(AUC.mean.sd, Cancer == CancerType)
#
# CancerType <- "STAD"
# AlgorithmType <- "L21"
# Tasks.keep <- filter.analysis.task$all_top
# tmp.barplot <- subset(AUC.mean.sd, Cancer == CancerType & Model == AlgorithmType & Task %in% Tasks.keep)

# Barplot #
ggplot2::ggplot(tmp.barplot, aes(x=Task, y=round(AUC.median,2), fill=View,  
                                 colour = View, alpha = Analysis_alg)) +
  ggplot2::geom_bar(stat="identity", position = position_dodge()) +
  ggplot2::scale_fill_manual(name = "Mechanistic signatures", 
                             labels = names(colors.DataType),
                              values = colors.DataType) +
  ggplot2::scale_color_manual(name = "Mechanistic signatures", 
                     labels = names(colors.DataType),
                     values = colors.DataType) + 
  ggplot2::scale_alpha_manual(name = "Algorithms", 
                     labels = names(alpha.analysis_alg),
                     values = alpha.analysis_alg)  +
  ggplot2::theme_minimal() +
  ggplot2::facet_grid(Dataset ~ .) +
  ggplot2::ylim(0, 1) +
  ggplot2::ylab("AUC") +
  ggplot2::geom_errorbar(aes(ymin = AUC.median - AUC.sd, ymax = AUC.median + AUC.sd), 
                         width=.5, color="gray", position = position_dodge(0.9)) +
  ggplot2::geom_text(aes(label= round(AUC.median,2)),stat = "identity", color="black", size = 3, angle = 90,
                     position = position_dodge(0.9)) +
  scale_x_discrete(labels = function(x) substr(x,1,8)) +
  theme(axis.text.x = element_text(size=10,face="bold", angle = 45, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,face="bold"), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_line(size = 1),
        axis.title.y = element_text(size=10,face="bold",vjust = 0.9), strip.background = element_blank(), panel.spacing.y =  unit(1, "lines"),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
        panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(size=0.5, colour = "black")) +
# ggtitle(paste0("Mechanistic signatures performance across ",CancerType, " datasets (",AnalysisType,")"))
  ggtitle(paste0("Mechanistic signatures performance across ",CancerType, " datasets (all vs all_top)"))

ggsave(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/new_v2/",CancerType,"/Validation/",
              "L21_with_all_datatypes_and_all_top_tasks.pdf"), width = 12, height = 12)
