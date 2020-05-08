# #########################################################################################################
# Script to plot model feature selection obtained from cross-validation -->
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

# ****************
# Select cancer type
## no filter
load("./pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
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
           sTIL = 'gaussian', #6
           LRpairs = 'gaussian', #7
           CYTOKINEpairs = 'gaussian')  #8)  

# **************** 
# Select data to examine

# All: comparison
view_combinations <- list(views[c(1,3)], views[1], views[3])

# ------------------------------------------------------------------------------------------------------------ #
# Collect data from L21 cross-validation results --> kfold = 5
# ------------------------------------------------------------------------------------------------------------ #

# **************** 
# Initialize variable to collect results
summary_view_features <- NULL
algorithms <- c("Elastic_Net")
analysis <- c("all")

# for (Cancer in PanCancer.names){

  # load(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/data/PanCancer/",Cancer,"/new/DataViews_no_filter_", Cancer,".RData"))
  load(paste0("../data/Federica_presentation_colab/DataViews_no_filter_SKCM.RData"))

  # file <- dir(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/output/new/",Cancer), full.names = T, recursive = F)
  file <- dir(paste0("../output/Federica_presentation_colab"), full.names = T, recursive = F)

  summary_analysis <- do.call(rbind, lapply(analysis, function(anal){
    
    summary_view <- do.call(rbind, lapply(1:length(view_combinations), function(view){
  
      # Load file
      input_name <- paste(names(view_combinations[[view]]), collapse ="_")
      which_file <- grep(pattern = paste0("_with_", anal,"_tasks_", input_name,".RData"), file, fixed = T)
      load(file[which_file])
      
      # 2 algorithms
      summary_alg <- do.call(rbind, lapply(names(all_cv_res), function(alg){
        
        tmp_alg <- all_cv_res[[alg]]
    
        # 100 iterations
        summary_iter <- do.call(rbind, lapply(1:length(tmp_alg), function(iteration){
        
          tmp_iter_model <- tmp_alg[[iteration]]$model
          
          if (alg == "L21") {tasks <- colnames(tmp_iter_model$cv.MTL.features$min.mse)}
          if (alg == "Elastic_Net") {tasks <- names(tmp_iter_model$Coef)}
          
          # All tasks
          summary_task <- do.call(rbind, lapply(tasks, function(task){
            
            if (alg == "L21") {
              tmp_iter_model_coef_task <- tmp_iter_model$cv.MTL.features$min.mse[,task]
              tmp_iter_model_hyp_task <- paste0("min-mse (",tmp_iter_model$cv.MTL.hyperparameters$min.mse$Lam1,",",
                                                round(as.numeric(tmp_iter_model$cv.MTL.hyperparameters$min.mse$Lam2),3),")")
              feature_names <- as.character(names(tmp_iter_model_coef_task)[2:length(names(tmp_iter_model_coef_task))])
            }
            
            if (alg == "Elastic_Net") {
              tmp_iter_model_coef_task <- tmp_iter_model$Coef[[task]]$min.mse
              tmp_iter_model_hyp_task <- paste0("min-mse (", tmp_iter_model$hyperparameters[[task]]$min.mse$alpha, ",",
                                                round(as.numeric(tmp_iter_model$hyperparameters[[task]]$min.mse$lambda),3), ")")
              feature_names <- as.character(rownames(tmp_iter_model_coef_task)[2:length(rownames(tmp_iter_model_coef_task))])
            }
  
            info <- data.frame(algorithm = alg,
                               iteration = rep(iteration, times = length(feature_names)),
                               model = rep(tmp_iter_model_hyp_task, times = length(feature_names)),
                               task = rep(task, times = length(feature_names)),
                               feature = feature_names,
                               estimate = as.numeric(tmp_iter_model_coef_task[feature_names]))
            
            return(info)
            
          }))
          
          return(summary_task)
          
        }))
        
        return(summary_iter)
      }))

      n_algorithms <- length(all_cv_res)
      n_iterations <- length(all_cv_res$Elastic_Net)
      n_tasks <- length(all_cv_res$Elastic_Net[[1]]$performances$MSE)
      n_features <- length(unique(summary_alg$feature))
      
      summary_alg$AnalysisType <- rep(anal, len = n_algorithms * n_iterations * n_tasks * n_features)
      summary_alg$DataType <- rep(input_name, len = n_algorithms * n_iterations * n_tasks  * n_features)
      summary_alg$CancerType <- rep(paste0(Cancer,"(n=",nrow(DataViews.no_filter$pathways),")"),
                                     len = n_algorithms * n_iterations * n_features * n_tasks)
      return(summary_alg)
    }))
    
    return(summary_view)
  }))
  
  summary_view_features <- rbind(summary_view_features, summary_analysis)
  
# }

##############################################################################
# Visualization
##############################################################################

# Making sure every variable is transformed to factor
summary_view_features$DataType <- factor(summary_view_features$DataType)
summary_view_features$CancerType <- factor(summary_view_features$CancerType)
summary_view_features$task <- factor(summary_view_features$task)
summary_view_features$AnalysisType <- factor(summary_view_features$AnalysisType)
summary_view_features$Analysis_Alg <- paste0(summary_view_features$algorithm,"_",
                                             summary_view_features$AnalysisType)
summary_view_features$Analysis_Alg <- factor(summary_view_features$Analysis_Alg)

# Colors for visualization
colors.cancer_types <- toupper(c("#d29d00","#9445cc","#8cce2f","#ef49c3","#1b8c00","#998fff",
                                 "#e78400","#0063a4","#fc4745","#00c98b","#b6006b","#006c43",
                                 "#ff8376","#604588", "#dcc666","#912f47","#7c5b00","#90350e"))

names(colors.cancer_types) <- levels(summary_view_features$CancerType)

colors.DataType <- toupper(c("#8c9f3e","#b65cbf","#4eac7c","#c95574","#747fca"))#,"#ca743e"))
names(colors.DataType) <- levels(summary_view_features$DataType)

colors.tasks <- toupper(c("#b47645","#ad58c5","#6cb643","#d24787","#52ad7b","#cf4740", "#4bafd0",
                          "#dc7b31","#6776cb","#c1ad46","#b975b1","#6c7b33","#c26671"))
names(colors.tasks) <- levels(summary_view_features$task)

colors.algorithm <- toupper(c("#ff7433","#853760"))
names(colors.algorithm) <- levels(summary_view_features$algorithm)

alpha.algorithm <- c(1,0.2)
names(alpha.algorithm) <- levels(summary_view_features$algorithm)

alpha.analysis_alg <- c(1,0.6,0.15)
names(alpha.analysis_alg) <- levels(summary_view_features$Analysis_Alg)

## ------------------------------------------------------------------------ ##
# 1. PanCancer comparison:
# One plot per task: compare different data types across types of tumors (L21 and EN)
## ----------------- ------------------------------------------------------ ##

# Sort features per median value and all tasks (~ feature + dataType + CancerType + algorithm)
## Get median features
median.features <- aggregate(estimate ~ feature + DataType  + CancerType , 
                                          FUN = "median", na.rm = TRUE, data = summary_view_features)
## Sort features by median value:
median.features.sort <- median.features[order(abs(median.features$estimate), decreasing = TRUE),]
median.features.sort$feature <- as.character(median.features.sort$feature)

## Restructuring data.frame
summary_view_features.tmp <- summary_view_features
# estimates
summary_view_features_tmp.order <- summary_view_features.tmp[order(match(summary_view_features.tmp$feature, 
                                                                         median.features.sort$feature)),]
summary_view_features_tmp.order$feature <- factor(summary_view_features_tmp.order$feature, 
                                                  levels =  unique(median.features.sort$feature))

## ------------------------------------------------------------------------ ##
# 3. Comparison of running L21 with 8 or 13 tasks:
# a) One plot per cancer type: compare different data types across tasks (L21 and EN)
# b) One plot per task: compare different data types across types of tumors (L21 and EN)
## ------------------------------------------------------------------------ ##

# Per task, DataType
# a)

# sapply(names(colors.DataType)[c(2:4)], function(input){
  
input <- "pathways_immunecells"

  summary_view_features.mech_sig <- subset(summary_view_features_tmp.order, DataType == input)

  # Per task
  sapply(names(colors.tasks), function(output){

    summary_view_features.mech_sig.task <- subset(summary_view_features.mech_sig, task == output)

      ggplot(summary_view_features.mech_sig.task, aes(x = feature, y = estimate, fill = CancerType,
                                                       colour = CancerType, alpha = Analysis_Alg)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(name = "Cancer Type",
                          labels = names(colors.cancer_types),
                          values = colors.cancer_types) +
        scale_colour_manual(name = "Cancer Type",
                            labels = names(colors.cancer_types),
                            values = colors.cancer_types)  +
        scale_alpha_manual(name = "Analysis",
                           labels = names(alpha.analysis_alg),
                           values = alpha.analysis_alg)  +
        theme_minimal() +
        scale_x_discrete(labels = function(x) format(x, width = 5)) +
        theme(axis.text.x = element_text(size=10,face="bold", angle = 45, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
              axis.text.y = element_text(size=10,face="bold"), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_line(size = 1),
              axis.title.y = element_text(size=10,face="bold",vjust = 0.9), strip.background = element_blank(), panel.spacing.y =  unit(1, "lines"),
              legend.position = c(0.75,0.85), legend.direction = "horizontal",
              legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
              panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(size=0.5, colour = "black")) +
        geom_hline(yintercept = 0, linetype = 2) +
        ggtitle(paste0(input, ": EN feature selection for ", output))

    ggsave(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/new_v2/PanCancer/Feature_selection/",
                  "PanCancer_feature_selection_EN_vs_L21_with_", input,"_for_", output,".pdf"), width = 12, height = 12)

  })

#})

# b)
# Per CancerType
sapply(names(colors.cancer_types), function(Cancer){
  
  # Per dataType
  sapply(names(colors.DataType), function(input){
    
    # Per task
      summary_view_features.mech_sig.cancer <- subset(summary_view_features_tmp.order, DataType == input &
                                                      CancerType == Cancer)
  
      ggplot(summary_view_features.mech_sig.cancer, aes(x = feature, y = estimate, fill = task,
                                               colour = task, alpha = Analysis_Alg)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(name = "Task", 
                          labels = names(colors.tasks),
                          values = colors.tasks) +
        scale_color_manual(name = "Task", 
                           labels = names(colors.tasks),
                           values = colors.tasks) + 
        scale_alpha_manual(name = "Analysis", 
                           labels = names(alpha.analysis_alg),
                           values = alpha.analysis_alg)  +
        theme_minimal() +
        #ylim(c(0,1)) +
        coord_fixed(ratio = 4) +
        scale_x_discrete(labels = function(x) substr(x,1,13)) +
        theme(axis.text.x = element_text(size=10,face="bold", angle = 45, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
              axis.text.y = element_text(size=10,face="bold"), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_line(size = 1),
              axis.title.y = element_text(size=10,face="bold",vjust = 0.9), strip.background = element_blank(), panel.spacing.y =  unit(1, "lines"),
              legend.position = "bottom", legend.direction = "horizontal",
              legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
              panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(size=0.5, colour = "black")) + 
        geom_hline(yintercept = 0, linetype = 2) +
        ggtitle(paste0(input," FE - all vs all_top tasks for ", Cancer))
    
    ggsave(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/new_v2/", 
                sapply(strsplit(Cancer, split = "(", fixed = T), function(X) {return(X[1])}), 
                "/Feature_selection/", sapply(strsplit(Cancer, split = "(", fixed = T), function(X) {return(X[1])}),
                "_EN_vs_L21_with_", input,"_all_vs_top_tasks.pdf"), width = 12, height = 12)
  
  })
})


# c)
# Stability selection
color.analysis_alg <- c("#d8005e","#a09f00","#ffad72")
names(color.analysis_alg) <- levels(summary_view_features_tmp.order$Analysis_Alg)

frequency <- aggregate(estimate ~ feature + task + Analysis_Alg + DataType + CancerType, FUN = function(X){sum(X != 0)},
                       data = summary_view_features_tmp.order)

sign.estimate <- aggregate(estimate ~ feature + task + Analysis_Alg + DataType + CancerType, FUN = function(X){sign(median(X))},
                           data = summary_view_features_tmp.order)

frequency <- cbind(frequency, sign.feature = sign.estimate$estimate)
frequency$sign.feature <- gsub(-1,"-", frequency$sign.feature)
frequency$sign.feature <- gsub(1,"+", frequency$sign.feature)
frequency$sign.feature <- factor(frequency$sign.feature)

# # Issue: Balance classes of directions due to coefficients close to zero
# table.tmp <- table(sign(try$estimate))
# sign <- ifelse(table.tmp[c(1:3)] > sum(table.tmp[c(1:3)])*0.8,1, 0)
# direction <- ifelse(any(sign == 1), names(sign)[which(sign == 1)], 0)

# Per CancerType
# sapply(names(colors.cancer_types), function(Cancer){
Cancer = "SKCM(n=467)"
  # Per dataType
  # sapply(names(colors.DataType), function(input){
input = "pathways_immunecells"
    frequency.cancer.input <- subset(frequency, DataType == input & CancerType == Cancer)
  
      ggplot(frequency.cancer.input, aes(x = feature, y = estimate)) +# fill = Analysis_Alg)) +
      geom_bar(stat="identity", color="black") +
      geom_text(aes(label= sign.feature), stat = "identity", color="black", size = 3, 
                position = position_stack(vjust = 0.5)) +
      # scale_fill_manual(name = "Analysis", 
      #                    labels = names(color.analysis_alg),
      #                    values = color.analysis_alg)  +
      theme_minimal() +
      facet_grid(task ~ .) +
      scale_x_discrete(labels = function(x) substr(x,1,20)) +
      theme(axis.text.x = element_text(size=5,face="bold", angle = 45, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
            axis.text.y = element_text(size=9,face="bold"), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_line(size = 1),
            axis.title.y = element_text(size=10,face="bold",vjust = 0.9), strip.background = element_blank(), panel.spacing.y =  unit(1, "lines"),
            legend.position = "bottom", legend.direction = "horizontal",  strip.text.y = element_text(angle = 360),
            legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
            panel.border = element_blank(), panel.background = element_blank(), 
            axis.line = element_line(size=0.5, colour = "black"), aspect.ratio = 0.05) + 
      labs(y = "Number of\n times") +
      ggtitle(paste0(input," - ", Cancer,": feature selection all vs top tasks"), 
              subtitle = "* Balance classes of directions due to coefficients close to zero are resolved with sign = 0")
    
    ggsave(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/new_v2/",
                  sapply(strsplit(Cancer, split = "(", fixed = T), function(X) {return(X[1])}), 
                  "/Feature_selection/", sapply(strsplit(Cancer, split = "(", fixed = T), function(X) {return(X[1])}),
                  "_EN_vs_L21_barplot_with_", input,"_all_vs_top_tasks.pdf"), width = 12, height = 12)
  })
})

frequency.L21 <- subset(frequency, Analysis_Alg %in% c("L21_all", "L21_all_top") & DataType == "TFs")

## ---------------------------------------------------------------------------------------------- ##
# 2. Cancer specific comparison:
# One plot per cancer type: compare features across all tasks (L21 and EN)
## ----------------- ---------------------------------------------------------------------------- ##




## ---------------------------------------------------------------------------------------------- ##
# 3. Combo improve single data
# Heatmap: showing coefficients when considering each individual feature set (e.g. pathways and
# immune cells separately) and then when considering them together. Additional heatmap can be used 
# to show the difference between looking at features individually or together.
## ----------------- -----------------------------------------------------------------------------##

color.combos <- c("#b3669e", "#98984d")
names(color.combos) <- c("separate", "together")

summary_view_features$Combo <- summary_view_features$dataType

summary_view_features$Combo <- gsub("pathways_immunecells","together", summary_view_features$Combo, fixed = TRUE)
summary_view_features$Combo <- gsub("immunecells","separate", summary_view_features$Combo, fixed = TRUE)
summary_view_features$Combo <- gsub("pathways","separate", summary_view_features$Combo, fixed = TRUE)
summary_view_features$Combo <- gsub("TFs","separate", summary_view_features$Combo, fixed = TRUE)


