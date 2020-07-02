# #########################################################################################################
# Script to plot predictions on external datasets:
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
library(ggpubr)
library(pdist)

# ****************
# initialize variables
PanCancer.names <- c("GBM","SKCM", "STAD") 

# ****************
# views
views <- c(Pathways = 'gaussian', ImmuneCells = 'gaussian', TFs = 'gaussian',
           LRpairs = 'gaussian', CYTOKINEpairs = 'gaussian', Transcript = 'gaussian')   

# **************** 
# Select data to examine

## Pathways and immune cells ## 
view_combinations <- list(views[1], views[2],views[c(1,2)], views[3], views[4], views[5])

#--------------------------------------------------------------------
# Manipulate labels for each dataset, we need just two levels
#--------------------------------------------------------------------

summary_cancer <- list()

for (CancerType in PanCancer.names){
  
  message(CancerType,"\n")

  # load predictions
  file <- dir(path = paste0("../output/PanCancer_draft_v1/Validation/", CancerType,"/"),
              pattern = "Labels", full.names = T, recursive = F)
  
  summary_predictions <- do.call(rbind, lapply(1:length(view_combinations), function(view){
    
    # Load file
    input_name <- paste(names(view_combinations[[view]]), collapse ="_")
    which_file <- grep(pattern = paste0("pre_treatment_", input_name,".RData"), file, fixed = T)
    load(file[which_file])
    
    predictions_immune_response <- predictions.test
    input_dataset <- names(predictions_immune_response)
    input_analysis <- names(predictions_immune_response[[1]])
    input_algorithm <- names(predictions_immune_response[[1]][[1]])
    
    predictions_immune_response <- lapply(input_dataset, function(dataset){
      predictions_immune_response[[dataset]] <- lapply(input_analysis, function(anal){
        predictions_immune_response[[dataset]][[anal]] <- lapply(input_algorithm, function(algorithm){ 
          
          # Check labels, we need to be consistent --> label.ordering = c("NR", "R") #

              df.label <- predictions_immune_response[[dataset]][[anal]][[algorithm]]$lab[[1]]
              
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
              predictions_immune_response[[dataset]][[anal]][[algorithm]]$lab[[1]] <- df.label
              return(predictions_immune_response[[dataset]][[anal]][[algorithm]])
          })
        names(predictions_immune_response[[dataset]][[anal]]) <- input_algorithm
        return(predictions_immune_response[[dataset]][[anal]])
        })
      names(predictions_immune_response[[dataset]]) <- input_analysis
      return(predictions_immune_response[[dataset]])
    })
    names(predictions_immune_response) <- input_dataset
    
    summary_predictions <- list(predictions_immune_response)
    return(summary_predictions)
  }))
  names(summary_predictions) <- sapply(1:length(view_combinations), function(X) paste(names(view_combinations[[X]]), collapse ="_"))
  
  summary_cancer[[CancerType]] <- summary_predictions
}

#--------------------------------------------------------------------
# Use of ROCR to prepare data to plot ROC curves
#--------------------------------------------------------------------

ROC_info <- list()

ROC_info <-  lapply(names(summary_cancer), function(Cancer){
  
  cat("Cancer:", Cancer, "\n")
  
  ROC_info <-  lapply(names(summary_cancer[[Cancer]]), function(view){
    
    cat("View:",view, "\n")
  
  ROC_info <-  lapply(names(summary_cancer[[Cancer]][[view]]), function(dataset){
    
    cat("Dataset:",dataset, "\n")
    
    ROC_info <-  lapply(names(summary_cancer[[Cancer]][[view]][[dataset]]), function(anal){
      
      cat("Analysis:",anal, "\n")
      
      ROC_info <-  lapply(names(summary_cancer[[Cancer]][[view]][[dataset]][[anal]]), function(algorithm){
        
        cat("Algorithm:",algorithm, "\n")
          
          ROC_info <- lapply(c(names(summary_cancer[[Cancer]][[view]][[dataset]][[anal]][[algorithm]]$pred),"common_mean","common_median"), function(task){
            
            cat("Task:",task, "\n")
            
            df <- summary_cancer[[Cancer]][[view]][[dataset]][[anal]][[algorithm]]
            
            # Build commond_mean and common_median task (according to algorithm)
            if (algorithm == "BEMKL"){
            
              df$pred[["common_mean"]][[1]] <- df$pred[["common_median"]][[1]] <- matrix(0,nrow(df$pred$CYT[[1]]), ncol = 100)

            }else if (algorithm == "Multi_Task_EN"){
            
              df$pred[["common_mean"]][["min.mse"]][[1]] <- matrix(0,nrow(df$pred$CYT[["min.mse"]][[1]]), ncol = 100)
              df$pred[["common_mean"]][["1se.mse"]][[1]] <- matrix(0,nrow(df$pred$CYT[["1se.mse"]][[1]]), ncol = 100)
              
              df$pred[["common_median"]][["min.mse"]][[1]] <- matrix(0,nrow(df$pred$CYT[["min.mse"]][[1]]), ncol = 100)
              df$pred[["common_median"]][["1se.mse"]][[1]] <- matrix(0,nrow(df$pred$CYT[["1se.mse"]][[1]]), ncol = 100)
              
            }
            
            if (anal == "cor"){
                if (algorithm == "BEMKL"){
                  
                  for (j in 1:100) {
                  
                    for(i in 1:nrow(df$pred$CYT[[1]])){
                    
                      df$pred[["common_mean"]][[1]][i,j] <- mean(c(df$pred$CYT[[1]][i,j], df$pred$IS[[1]][i,j],df$pred$RohIS[[1]][i,j], 
                                                                   df$pred$IS_Davoli[[1]][i,j], df$pred$IFny[[1]][i,j],df$pred$ExpandedImmune[[1]][i,j], 
                                                                   df$pred$chemokine[[1]][i,j], df$pred$T_cell_inflamed[[1]][i,j]))
                    
                      df$pred[["common_median"]][[1]][i,j] <- median(c(df$pred$CYT[[1]][i,j], df$pred$IS[[1]][i,j],df$pred$RohIS[[1]][i,j], 
                                                                       df$pred$IS_Davoli[[1]][i,j], df$pred$IFny[[1]][i,j],df$pred$ExpandedImmune[[1]][i,j], 
                                                                       df$pred$chemokine[[1]][i,j], df$pred$T_cell_inflamed[[1]][i,j]))
                    }
                  }
                  rownames(df$pred[["common_mean"]][[1]]) <- rownames(df$lab[[1]])
                  rownames(df$pred[["common_mean"]][[1]]) <- rownames(df$lab[[1]])
                  rownames(df$pred[["common_median"]][[1]]) <- rownames(df$lab[[1]])
                  rownames(df$pred[["common_median"]][[1]]) <- rownames(df$lab[[1]])
                  
                  
                }else if (algorithm == "Multi_Task_EN"){
                  
                  for (j in 1:100) {
                  
                    for(i in 1:nrow(df$pred$CYT$min.mse[[1]])){
                      
                      df$pred[["common_mean"]]$min.mse[[1]][i,j] <- mean(c(df$pred$CYT$min.mse[[1]][i,j], df$pred$IS$min.mse[[1]][i,j],df$pred$RohIS$min.mse[[1]][i,j], 
                                                                           df$pred$IS_Davoli$min.mse[[1]][i,j], df$pred$IFny$min.mse[[1]][i,j],df$pred$ExpandedImmune$min.mse[[1]][i,j], 
                                                                           df$pred$chemokine$min.mse[[1]][i,j], df$pred$T_cell_inflamed$min.mse[[1]][i,j]))
                      
                      df$pred[["common_mean"]]$`1se.mse`[[1]][i,j] <- mean(c(df$pred$CYT$`1se.mse`[[1]][i,j], df$pred$IS$`1se.mse`[[1]][i,j],df$pred$RohIS$`1se.mse`[[1]][i,j], 
                                                                             df$pred$IS_Davoli$`1se.mse`[[1]][i,j], df$pred$IFny$`1se.mse`[[1]][i,j],df$pred$ExpandedImmune$`1se.mse`[[1]][i,j], 
                                                                             df$pred$chemokine$`1se.mse`[[1]][i,j], df$pred$T_cell_inflamed$`1se.mse`[[1]][i,j]))
                      
                      df$pred[["common_median"]]$min.mse[[1]][i,j] <- median(c(df$pred$CYT$min.mse[[1]][i,j], df$pred$IS$min.mse[[1]][i,j],df$pred$RohIS$min.mse[[1]][i,j], 
                                                                               df$pred$IS_Davoli$min.mse[[1]][i,j], df$pred$IFny$min.mse[[1]][i,j],df$pred$ExpandedImmune$min.mse[[1]][i,j], 
                                                                               df$pred$chemokine$min.mse[[1]][i,j], df$pred$T_cell_inflamed$min.mse[[1]][i,j]))
                      
                      df$pred[["common_median"]]$`1se.mse`[[1]][i,j] <-  median(c(df$pred$CYT$`1se.mse`[[1]][i,j], df$pred$IS$`1se.mse`[[1]][i,j],df$pred$RohIS$`1se.mse`[[1]][i,j], 
                                                                                  df$pred$IS_Davoli$`1se.mse`[[1]][i,j], df$pred$IFny$`1se.mse`[[1]][i,j],df$pred$ExpandedImmune$`1se.mse`[[1]][i,j], 
                                                                                  df$pred$chemokine$`1se.mse`[[1]][i,j], df$pred$T_cell_inflamed$`1se.mse`[[1]][i,j]))
                      
                   }
                }
                rownames(df$pred[["common_mean"]]$min.mse[[1]]) <- rownames(df$lab[[1]])
                rownames(df$pred[["common_mean"]]$`1se.mse`[[1]]) <- rownames(df$lab[[1]])
                rownames(df$pred[["common_median"]]$min.mse[[1]]) <- rownames(df$lab[[1]])
                rownames(df$pred[["common_median"]]$`1se.mse`[[1]]) <- rownames(df$lab[[1]])
              }
              
          }
            #else{
            #   
            #   for (j in 1:100) {
            #     for(i in 1:nrow(df$pred$CYT$min.mse[[1]])){
            #       
            #       df$pred[["common_mean"]]$min.mse[[1]][i,j] <- mean(c(df$pred$CYT$min.mse[[1]][i,j],df$pred$IS$min.mse[[1]][i,j], df$pred$IPS$min.mse[[1]][i,j],
            #                                                            df$pred$IMPRES$min.mse[[1]][i,j], df$pred$RohIS$min.mse[[1]][i,j], df$pred$chemokine$min.mse[[1]][i,j],
            #                                                            df$pred$IS_Davoli$min.mse[[1]][i,j], df$pred$Proliferation$min.mse[[1]][i,j], df$pred$IFny$min.mse[[1]][i,j],
            #                                                            df$pred$ExpandedImmune$min.mse[[1]][i,j], 
            #                                                            df$pred$T_cell_inflamed$min.mse[[1]][i,j], df$pred$TIDE$min.mse[[1]][i,j], df$pred$MSI$min.mse[[1]][i,j]))
            #       
            #       df$pred[["common_mean"]]$`1se.mse`[[1]][i,j] <- mean(c(df$pred$CYT$`1se.mse`[[1]][i,j],df$pred$IS$`1se.mse`[[1]][i,j], df$pred$IPS$`1se.mse`[[1]][i,j],
            #                                                              df$pred$IMPRES$`1se.mse`[[1]][i,j], df$pred$RohIS$`1se.mse`[[1]][i,j], df$pred$chemokine$`1se.mse`[[1]][i,j],
            #                                                              df$pred$IS_Davoli$`1se.mse`[[1]][i,j], df$pred$Proliferation$`1se.mse`[[1]][i,j], df$pred$IFny$`1se.mse`[[1]][i,j],
            #                                                              df$pred$ExpandedImmune$`1se.mse`[[1]][i,j], 
            #                                                              df$pred$T_cell_inflamed$`1se.mse`[[1]][i,j], df$pred$TIDE$`1se.mse`[[1]][i,j], df$pred$MSI$`1se.mse`[[1]][i,j]))
            #       
            #       df$pred[["common_median"]]$min.mse[[1]][i,j] <- median(c(df$pred$CYT$min.mse[[1]][i,j],df$pred$IS$min.mse[[1]][i,j], df$pred$IPS$min.mse[[1]][i,j],
            #                                                                df$pred$IMPRES$min.mse[[1]][i,j], df$pred$RohIS$min.mse[[1]][i,j], df$pred$chemokine$min.mse[[1]][i,j],
            #                                                                df$pred$IS_Davoli$min.mse[[1]][i,j], df$pred$Proliferation$min.mse[[1]][i,j], df$pred$IFny$min.mse[[1]][i,j],
            #                                                                df$pred$ExpandedImmune$min.mse[[1]][i,j], 
            #                                                                df$pred$T_cell_inflamed$min.mse[[1]][i,j], df$pred$TIDE$min.mse[[1]][i,j], df$pred$MSI$min.mse[[1]][i,j]))
            #       
            #       df$pred[["common_median"]]$`1se.mse`[[1]][i,j] <-  median(c(df$pred$CYT$`1se.mse`[[1]][i,j],df$pred$IS$`1se.mse`[[1]][i,j], df$pred$IPS$`1se.mse`[[1]][i,j],
            #                                                                   df$pred$IMPRES$`1se.mse`[[1]][i,j], df$pred$RohIS$`1se.mse`[[1]][i,j], df$pred$chemokine$`1se.mse`[[1]][i,j],
            #                                                                   df$pred$IS_Davoli$`1se.mse`[[1]][i,j], df$pred$Proliferation$`1se.mse`[[1]][i,j], df$pred$IFny$`1se.mse`[[1]][i,j],
            #                                                                   df$pred$ExpandedImmune$`1se.mse`[[1]][i,j], 
            #                                                                   df$pred$T_cell_inflamed$`1se.mse`[[1]][i,j], df$pred$TIDE$`1se.mse`[[1]][i,j], df$pred$MSI$`1se.mse`[[1]][i,j]))
            #       
            #       
            #     }
            #   }
            #   rownames(df$pred[["common_mean"]]$min.mse[[1]]) <- rownames(df$lab$CYT[[1]])
            #   rownames(df$pred[["common_mean"]]$`1se.mse`[[1]]) <- rownames(df$lab$CYT[[1]])
            #   rownames(df$pred[["common_median"]]$min.mse[[1]]) <- rownames(df$lab$CYT[[1]])
            #   rownames(df$pred[["common_median"]]$`1se.mse`[[1]]) <- rownames(df$lab$CYT[[1]])
            #   
            # }
       
            # Obtain data for ROCR package (according to algorithm)
            
            if (algorithm == "BEMKL"){
              
              pred.min.mse <- ROCR::prediction(df$pred[[task]][[1]],df$lab[[1]], label.ordering = c("NR", "R"))
              perf.min.mse <- ROCR::performance(pred.min.mse,"tpr","fpr")
              AUC.min.mse <- unlist(ROCR::performance(pred.min.mse, "auc")@y.values)
              
              perf.1se.mse <- NULL
              AUC.1se.mse <- NULL
              
            }else if (algorithm == "Multi_Task_EN"){
            
              pred.min.mse <- ROCR::prediction(df$pred[[task]]$min.mse[[1]],df$lab[[1]], label.ordering = c("NR", "R"))
              pred.1se.mse <- ROCR::prediction(df$pred[[task]]$`1se.mse`[[1]],df$lab[[1]], label.ordering = c("NR", "R"))
            
              perf.min.mse <- ROCR::performance(pred.min.mse,"tpr","fpr")
              perf.1se.mse <- ROCR::performance(pred.1se.mse,"tpr","fpr")
            
              AUC.min.mse <- unlist(ROCR::performance(pred.min.mse, "auc")@y.values)
              AUC.1se.mse <- unlist(ROCR::performance(pred.1se.mse, "auc")@y.values)
            }
            
            data_ROC <- list(perf.min.mse = perf.min.mse, perf.1se.mse = perf.1se.mse)
            Barplot <- list(AUC.min.mse = AUC.min.mse, AUC.1se.mse = AUC.1se.mse)
            
            return(list(Curve = data_ROC, Barplot = Barplot))
          })
          names(ROC_info) <- c(names(summary_cancer[[Cancer]][[view]][[dataset]][[anal]][[algorithm]]$pred),
                               "common_mean","common_median")
          return(ROC_info)
        })
        names(ROC_info) <- names(summary_cancer[[Cancer]][[view]][[dataset]][[anal]])
        return(ROC_info)
      })
      names(ROC_info) <- names(summary_cancer[[Cancer]][[view]][[dataset]])
      return(ROC_info)
    })
    names(ROC_info) <- names(summary_cancer[[Cancer]][[view]])
    return(ROC_info)
  })
  names(ROC_info) <- names(summary_cancer[[Cancer]])
  return(ROC_info)
})
names(ROC_info) <- names(summary_cancer)


#--------------------------------------------------------------------
# Calculate mean and sd from AUC to plot barplots.
#--------------------------------------------------------------------

AUC.mean.sd <- do.call(rbind, lapply(names(summary_cancer), function(Cancer){
  
  cat("Cancer:", Cancer, "\n")
  
  AUC.mean.sd <- do.call(rbind, lapply(names(summary_cancer[[Cancer]]), function(view){
    
    cat("View:",view, "\n")
    
    AUC.mean.sd <- do.call(rbind, lapply(names(summary_cancer[[Cancer]][[view]]), function(dataset){
      
      cat("Dataset:",dataset, "\n")
      
      AUC.mean.sd <- do.call(rbind, lapply(names(summary_cancer[[Cancer]][[view]][[dataset]]), function(anal){
        
        cat("Analysis:",anal, "\n")
        
        AUC.mean.sd <- do.call(rbind, lapply(names(summary_cancer[[Cancer]][[view]][[dataset]][[anal]]), function(algorithm){
          
          cat("Algorithm:",algorithm, "\n")
          
          AUC.mean.sd <- do.call(rbind, lapply(c(names(summary_cancer[[Cancer]][[view]][[dataset]][[anal]][[algorithm]]$pred),"common_mean","common_median"), function(task){
            
            cat("Task:",task, "\n")
            
            if (algorithm == "BEMKL"){
              
              AUC.data <- data.frame(Cancer = Cancer,
                                     Dataset = dataset,
                                     Analysis = anal,
                                     View = view, 
                                     Model = algorithm,
                                     cv_model = rep(c("min.mse"), each = 100),
                                     Task = task,
                                     iteration = seq(1,100),
                                     AUC = c(as.numeric(ROC_info[[Cancer]][[view]][[dataset]][[anal]][[algorithm]][[task]]$Barplot[["AUC.min.mse"]])))
              
            }else if (algorithm == "Multi_Task_EN"){
              
              AUC.data <- data.frame(Cancer = Cancer,
                                     Dataset = dataset,
                                     Analysis = anal,
                                     View = view, 
                                     Model = algorithm,
                                     cv_model = rep(c("min.mse", "1se.mse"), each = 100),
                                     Task = task,
                                     iteration = seq(1,100),
                                     AUC = c(as.numeric(ROC_info[[Cancer]][[view]][[dataset]][[anal]][[algorithm]][[task]]$Barplot[["AUC.min.mse"]]),
                                             as.numeric(ROC_info[[Cancer]][[view]][[dataset]][[anal]][[algorithm]][[task]]$Barplot[["AUC.1se.mse"]])))
              
            }
            
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

# Subset # 
CancerType <- "STAD" # "SKCM", "STAD" "GBM
tmp.barplot <- subset(AUC.mean.sd, Cancer == CancerType)
tmp.barplot$Task <- factor(tmp.barplot$Task, levels = c("CYT","IS","RohIS", "chemokine","IS_Davoli",
                                                        "IFny","ExpandedImmune", "T_cell_inflamed",
                                                        "common_mean","common_median"))
tmp.barplot$Analysis <- factor(tmp.barplot$Analysis)

tmp.barplot.EN <- subset(tmp.barplot, Model == "Multi_Task_EN" & cv_model == "1se.mse")
tmp.barplot.BE <- subset(tmp.barplot, Model == "BEMKL")
  
tmp.barplot.both <- rbind(tmp.barplot.EN, tmp.barplot.BE)
tmp.barplot.both$Model <- factor(tmp.barplot.both$Model)

# ****************
# Colors for visualization
# colors.DataType <- toupper(c("#8c9f3e","#b65cbf","#4eac7c","#c95574","#747fca","#ca743e"))
# names(colors.DataType) <- names(predictions_immune_response$GBM$Zhao$all)

colors.DataType <- toupper(c("#99007a","#01851e","#ffa5e3",
                              "#d5a318","#0078d1","#c21a35"))
names(colors.DataType) <- names(ROC_info$SKCM)

colors.tasks <- toupper(c("#b47645","#ad58c5","#6cb643","#d24787","#52ad7b","#cf4740", "#4bafd0",
                          "#dc7b31","#6776cb","#c1ad46"))
names(colors.tasks) <- names(ROC_info$SKCM$Pathways$Auslander$cor$Multi_Task_EN)

alpha.model <- toupper(c(0.93,0.5))
names(colors.algorithm) <- names(ROC_info$SKCM$Pathways$Auslander$cor)

dataset <- "Kim" # "comb.Hugo_Liu_Riaz"
tmp.barplot.both.tmp <- subset(tmp.barplot.both, Dataset == dataset)
# Barplot #
ggplot2::ggplot(tmp.barplot.both.tmp, aes(x=View, y=round(AUC.median,2), fill=Task,  
                                 colour = Task, alpha = Model)) +
  ggplot2::geom_bar(stat="identity", position = position_dodge()) +
  ggplot2::scale_fill_manual(name = "Task", 
                             labels = levels(tmp.barplot.both$Task),
                             values = as.vector(colors.tasks)) +
  ggplot2::scale_color_manual(name = "Task", 
                              labels = levels(tmp.barplot.both$Task),
                              values =  as.vector(colors.tasks)) + 
  ggplot2::scale_alpha_manual(name = "Model",
                              labels = levels(tmp.barplot.both$Model),
                              values = c(0.93,0.4))  +
  ggplot2::theme_classic() +
  ggplot2::facet_grid( . ~ Dataset) +
  ggplot2::ylim(0, 1) +
  ggplot2::ylab("AUC") +
  ggplot2::geom_errorbar(aes(ymin = AUC.median - AUC.sd, ymax = AUC.median + AUC.sd), 
                         width=.5, color="gray", position = position_dodge(0.9)) +
  ggplot2::geom_text(aes(label= round(AUC.median,2)),stat = "identity", color="black", size = 3, angle = 45,
                     position = position_dodge(0.9)) +
  scale_x_discrete(labels = c("ImmuneCells"= "IC",
                              "LRpairs"= "L-R",
                              "CYTOKINEpairs"= "CK-CK",
                              "Pathways"="PA",
                              "Pathways_ImmuneCells"="PA + IC",
                              "TFs" = "TFs",
                              "Transcript" = "Transcriptomics")) +  
  ggplot2::geom_hline(yintercept = 0.5, linetype="dashed") + 
  theme(axis.text.x = element_text(size=10,face="bold", angle = 45, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,face="bold"), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_line(size = 1),
        axis.title.y = element_text(size=10,face="bold",vjust = 0.9), panel.spacing.y =  unit(1, "lines"),
        legend.box.background = element_rect(color="black", size=0.3),
        legend.position = c(0.74,0.88), legend.direction = "horizontal",
        legend.box.margin = margin(1, 1, 1, 1),
        legend.text=element_text(size=8), 
        legend.title = element_text(size = 8, vjust = 0.5),
        strip.background = element_rect(fill = "lightgray", colour = "white"),
        panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(size=0.5, colour = "black")) 
  #ggtitle(paste0("Mechanistic signatures AUC values  ",CancerType, " datasets"))

ggsave(paste0("../figures/PanCancer_draft_v1/External_validation/","Barplot_with_aucs_on_", CancerType, "_", dataset, "_dataset.pdf"),
       width = 12, height = 12)


#--------------------------------------------------------------------
# Plot ROC curve
#--------------------------------------------------------------------
dataset = "Liu"

source("../R/computation.gold.standards.R")
list.gold.standards <- c("CYT", "PD1", "PDL1","CTLA4")
load("../data/PanCancer_draft_v1/Validation/All_DataViews_test_pre.RData")

gold.standards <- computation.gold.standards(RNA.tpm = 2^All.DataViews.test$Liu$Transcript -1,
                                             list_gold_standards = list.gold.standards)

df <- summary_cancer$SKCM$Pathways_ImmuneCells$Liu$cor$Multi_Task_EN$lab[[1]]
gold.standards.curve <- list()
gold.standards.curve<- do.call(c,lapply(list.gold.standards, function(X){
  
  pred <- ROCR::prediction(gold.standards[[X]],df[,1], label.ordering = c("NR", "R"))
  perf <- ROCR::performance(pred,"tpr","fpr")
  AUC <- unlist(ROCR::performance(pred, "auc")@y.values)
  
  gold.standards.curve[[X]]$Curve <- perf
  gold.standards.curve[[X]]$Barplot <- AUC
  return(gold.standards.curve)
}))

# pathways_immunecells - common 
pdf("../figures/PanCancer_draft_v1/External_validation/ROC_curve_SKCM_comb_Liu_pathways_immunecells.pdf",
    width = 12, height = 12)
par(cex.axis=2.2, mar = c(5, 5, 5, 5))
ROCR::plot(ROC_info$SKCM$Pathways_ImmuneCells$Liu$cor$Multi_Task_EN[["common_median"]]$Curve$perf.1se.mse,
           avg = "threshold", col = colors.DataType["Pathways_ImmuneCells"], lwd = 8, type = "S", cex.lab=2.2, ylab="True Positive Rate",
           xlab="False Positive Rate")

# # immunecells - common
# ROCR::plot(ROC_info[[Cancer]][[dataset]][["cor"]][["immunecells"]][[algorithm]][["common_median"]]$Curve$perf.1se.mse,
#            avg = "threshold", col = colors.DataType[[2]], lwd = 5, type = "S", add = TRUE)
# # pathways - common
# ROCR::plot(ROC_info[[Cancer]][[dataset]][["cor"]][["pathways"]][[algorithm]][["common_median"]]$Curve$perf.1se.mse,
#            avg = "threshold", col = colors.DataType[[3]], lwd = 5, type = "S", add = TRUE)

# transcript - common
ROCR::plot(ROC_info$SKCM$Pathways$Liu$cor$Multi_Task_EN[["common_median"]]$Curve$perf.1se.mse,
           avg = "threshold", col = colors.DataType["Pathways"], lwd = 8, type = "s", add = TRUE)

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
                                            round(subset(AUC.mean.sd, Dataset == dataset & Task == "common_median" & Analysis == "cor" & cv_model == "1se.mse" &
                                                           View == "Pathways_ImmuneCells" & algorithm == "Multi_Task_EN")$AUC.median,3),")"," \nImmuneCells"),
                                     # c(paste0("Cytokine pairs", "(AUC=",
                                     #        round(subset(AUC.mean.sd, Task == "common_median" & Analysis == "cor" & cv_model == "1se.mse" &
                                     #                       View == "CYTOKINEpairs")$AUC.median,3),")"),
                                     paste0("Pathways", "(AUC=",
                                            round(subset(AUC.mean.sd, Dataset == dataset & Task == "common_median" & Analysis == "cor" & cv_model == "1se.mse" &
                                                           View == "Pathways" & algorithm == "Multi_Task_EN")$AUC.median,3),")"), 
                                     paste0("Cytolytic Ativity", "(AUC=", round(gold.standards.curve[["CYT"]][["Barplot"]],3),")"),
                                     
                                     paste0("PD1", "(AUC=", round(gold.standards.curve[["PD1"]][["Barplot"]],3),")")),
       
       col = c(colors.DataType["Pathways_ImmuneCells"],colors.DataType["Pathways"],"gray","gray50"),
       lty = c(1,1,1,1), lwd = c(8,8,2,2), cex = 1.4, bty = "n")


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





