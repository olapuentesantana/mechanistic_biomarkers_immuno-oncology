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
#load("./pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_Spat.RData")
#PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS_SpatialTILs)
## filter prot
#load("./pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS_prot.RData")
#PanCancer.names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS_prot)

# ****************
# views
views <- c(Pathways = 'gaussian')

# **************** 
# Select data to examine

## Pathways and immune cells ## 
view_combinations <- list(views[1])

# --------------------------------------------- #
# Collect performance info from cross-validation 
# --------------------------------------------- #

# Initialize variables
summary_view_all <- NULL
algorithms <- c("Multi_Task_EN", "BEMKL")
analysis <- c("cor") # all

for (Cancer in PanCancer.names){
  
  message(Cancer,"\n")
  
  load(paste0("../data/PanCancer_draft_v1/",Cancer,"/DataViews_no_filter_", Cancer,".RData"))
  
  file <- dir(path = paste0("../output/PanCancer_draft_v1/", Cancer,"/"), pattern = "all_cv_res_",
              full.names = T, recursive = F)
  
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
          
          tmp_iter <- tmp_alg[[iteration]]$performances
          
          # 4 metrics
          summary_measure <- do.call(rbind, lapply(names(tmp_iter), function(measure){
            
            tmp_measure <- tmp_iter[[measure]]
            
            if (alg == "BEMKL") {
              
              info <- data.frame(algorithm = alg,
                                 cv = "min.mse",
                                 iteration = iteration,
                                 metric = measure,
                                 task = c(names(tmp_measure), "common_median"),
                                 perf_min = c(tmp_measure, median(tmp_measure)))
              return(info)
              
            }else if (alg == "Multi_Task_EN") {
              
              # cv- 1se.mse and min.mse
              summary_cv <- do.call(rbind,lapply(names(tmp_measure), function(cv){
                
                tmp_cv <- tmp_measure[[cv]]
                
                info <- data.frame(algorithm = alg,
                                   cv = cv,
                                   iteration = iteration,
                                   metric = measure,
                                   task = c(names(tmp_cv),"common_median"),
                                   perf_min = c(tmp_cv, median(tmp_cv)))
                
                return(info) 
              }))
              return(summary_cv)
            }
          }))
          return(summary_measure)
        }))
        return(summary_iter) 
      }))
      len_summary <- nrow(summary_alg)
      
      summary_alg$AnalysisType <- rep(anal, len = len_summary)
      summary_alg$DataType <- rep(input_name, len = len_summary)
      summary_alg$CancerType <- rep(paste0(Cancer,"(n=",nrow(DataViews.no_filter$ImmuneCells),")"),
                                    len =  len_summary)
      
      return(summary_alg)
    }))
    
    return(summary_view)
  }))
  
  summary_view_all <- rbind(summary_view_all, summary_analysis)
}

# --------------------------------------------- #
# Visualization: DataType performances
# --------------------------------------------- #

# Theoretically, we cannot assume that the response follows a normal distribution, so it would be better 
# to use spearman correlation as metric

# Convert data.frame variables into factors
summary_view_all$DataType <- factor(summary_view_all$DataType)
summary_view_all$CancerType <- factor(summary_view_all$CancerType)
summary_view_all$cv <- factor(summary_view_all$cv, levels = c("min.mse", "1se.mse"))
summary_view_all$metric <- factor(summary_view_all$metric, levels = c("MSE", "PeCorr", "SpCorr", "CI"))
summary_view_all$AnalysisType <- factor(summary_view_all$AnalysisType)
summary_view_all$task <- factor(summary_view_all$task, levels = unique(summary_view_all$task))

# Find color palettes for optimal distincition
# Cancer type
colors.cancer_types <- toupper(c("#31efb8","#6138a9","#b59000","#0175e7","#c54c17","#be8fff",
                                 "#005211","#45005b","#00a369","#b5003b","#f4d07e","#694e8e",
                                 "#9b1e00","#ff7ea8","#6d2000","#ff7e5b","#580e08","#ca815b"))
names(colors.cancer_types) <- levels(summary_view_all$CancerType)
# View
colors.DataType <- toupper(c("#01dca1"))
names(colors.DataType) <- levels(summary_view_all$DataType)
# Task
colors.tasks <- toupper(c("#760068","#df9d1e","#121278","#616400",
                          "#0162c3","#811200","#ff7ec2","#ffb37e","#ef4959"))
names(colors.tasks) <- levels(summary_view_all$task)
# Algorithm
colors.algorithm <- toupper(c("#ff7433","#853760"))
names(colors.algorithm) <- levels(summary_view_all$algorithm)
alpha.algorithm <- c(1,0.2)
names(alpha.algorithm) <- levels(summary_view_all$algorithm)

# Elastic Net
summary_view_all.SpCorr_1se.mse.EN <- subset(summary_view_all, metric == "SpCorr" & algorithm == "Multi_Task_EN" 
                                             & cv == "1se.mse"
                                             & task == "common_median" 
                                             & AnalysisType == "cor")
summary_view_all.SpCorr_1se.mse.EN$task <- factor(summary_view_all.SpCorr_1se.mse.EN$task)
# BEMKL
summary_view_all.SpCorr_1se.mse.BA <- subset(summary_view_all, metric == "SpCorr" & algorithm == "BEMKL"
                                             & task == "common_median" 
                                             & AnalysisType == "cor")
summary_view_all.SpCorr_1se.mse.BA$task <- factor(summary_view_all.SpCorr_1se.mse.BA$task)

# Both algorithms
summary_view_all.SpCorr.both <- rbind(summary_view_all.SpCorr_1se.mse.EN, summary_view_all.SpCorr_1se.mse.BA)
summary_view_all.SpCorr.both$task <- factor(summary_view_all.SpCorr.both$task)
summary_view_all.SpCorr.both$algorithm <- factor(summary_view_all.SpCorr.both$algorithm)

# Using ggplot
ggplot(summary_view_all.SpCorr.both, aes(x = CancerType, y = perf_min, fill = algorithm,
                                         colour = algorithm)) + # , alpha = algorithm)) +
  geom_boxplot() +
  scale_fill_manual(name = "Algorithm",
                    labels = levels(summary_view_all.SpCorr.both$algorithm),
                    values = as.vector(colors.algorithm)) +
  scale_color_manual(name = "Algorithm",
                     labels = levels(summary_view_all.SpCorr.both$algorithm),
                     values = as.vector(colors.algorithm)) +
  #geom_text(aes(CancerType, 0.99, label=lab), nudge_x = 0.1, check_overlap = TRUE,
  #          color = "black", size = 4, data=summary_view_all.SpCorr.both) + 
  #geom_segment(data = dLines, aes(x = x, xend = xend, y = y, yend = y), color = "black") +
  ylim(0.25,1) +
  theme_bw() + theme(panel.grid = element_blank()) +
  #facet_grid(algorithm ~ .) +
  theme(axis.text.x = element_text(size=8,face="bold", angle = 55, vjust = 0.75), axis.title.x = element_blank(),
        axis.text.y = element_text(size=8,face="bold"), axis.title.y = element_text(size=10,face="bold"), 
        axis.ticks.x = element_blank(),
        legend.position = "top", legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
        panel.border = element_blank(), strip.background = element_rect(fill = "lightgray", colour = "white")) + #, axis.line.y = element_line(size=0.5, colour = "gray"))    
  labs(y = "Spearman Correlation")
# ggtitle("Multi-task Lasso regression (using common median)")

ggsave(paste0("../figures/PanCancer_draft_v1/Pathways/PanCancer_both_alg_Pathways_training_performance_SpCorr_1SE_MSE.pdf"),
       width = 10, height = 10)

# --------------------------------------------- #
# Collect features info from cross-validation 
# --------------------------------------------- #

# Initialize variables
summary_view_features <- NULL
algorithms <- c("Multi_Task_EN")
analysis <- c("cor") # all

for (Cancer in PanCancer.names){
  message(Cancer,"\n")
  
  load(paste0("../data/PanCancer_draft_v1/",Cancer,"/DataViews_no_filter_", Cancer,".RData"))
  
  file <- dir(path = paste0("../output/PanCancer_draft_v1/", Cancer,"/"), pattern = "all_cv_res_", 
              full.names = T, recursive = F)
  
  summary_analysis <- do.call(rbind, lapply(analysis, function(anal){
    
    summary_view <- do.call(rbind, lapply(1:length(view_combinations), function(view){
      
      # Load file
      input_name <- paste(names(view_combinations[[view]]), collapse ="_")
      which_file <- grep(pattern = paste0("_with_", anal,"_tasks_", input_name,".RData"), file, fixed = T)
      load(file[which_file])
      
      tmp_alg <- all_cv_res[["Multi_Task_EN"]]
      
      # 100 iterations
      summary_iter <- do.call(rbind, lapply(1:length(tmp_alg), function(iteration){
        
        tmp_iter_model <- tmp_alg[[iteration]]$model
        
        # cv = 1se.mse
        tmp_iter_model_cv <- tmp_iter_model$cv.glmnet.features[["1se.mse"]]
        
        # # All tasks
        summary_task <- do.call(rbind, lapply(colnames(tmp_iter_model_cv), function(task){
          
          info <- data.frame(iteration = iteration,
                             hyp_model = paste0(tmp_iter_model$cv.glmnet.hyperparameters[["1se.mse"]]$alpha,",",
                                                round(as.numeric(tmp_iter_model$cv.glmnet.hyperparameters[["1se.mse"]]$lambda),3)),
                             task = task,
                             feature = rownames(tmp_iter_model_cv)[2:nrow(tmp_iter_model_cv)],
                             estimate = tmp_iter_model_cv[2:nrow(tmp_iter_model_cv),task])
          
          return(info)
        }))
        
        return(summary_task)
      }))
      
      len_summary <- nrow(summary_iter)
      
      summary_iter$AnalysisType <- rep(anal, len = len_summary)
      summary_iter$DataType <- rep(input_name, len = len_summary)
      summary_iter$CancerType <- rep(paste0(Cancer,"(n=",nrow(DataViews.no_filter$ImmuneCells),")"),
                                     len =  len_summary)
      return(summary_iter)
    }))
    return(summary_view)
  }))
  summary_view_features <- rbind(summary_view_features, summary_analysis)
}

# --------------------------------------------- #
# Visualization: Features selected
# --------------------------------------------- #

# Convert data.frame variables into factors
summary_view_features$DataType <- factor(summary_view_features$DataType)
summary_view_features$CancerType <- factor(summary_view_features$CancerType)
summary_view_features$task <- factor(summary_view_features$task)
summary_view_features$AnalysisType <- factor(summary_view_features$AnalysisType)

# Find color palettes for optimal distincition
# Cancer type
colors.cancer_types <- toupper(c("#31efb8","#6138a9","#b59000","#0175e7","#c54c17","#be8fff",
                                 "#005211","#45005b","#00a369","#b5003b","#f4d07e","#694e8e",
                                 "#9b1e00","#ff7ea8","#6d2000","#ff7e5b","#580e08","#ca815b"))
names(colors.cancer_types) <- levels(summary_view_all$CancerType)
# View
colors.DataType <- toupper(c("#01dca1"))
names(colors.DataType) <- levels(summary_view_all$DataType)
# Task
colors.tasks <- toupper(c("#760068","#df9d1e","#121278","#616400",
                          "#0162c3","#811200","#ff7ec2","#ffb37e"))
names(colors.tasks) <- levels(summary_view_features$task)

# Using ggplot:
library(gtable)
# 1. We show the distribution of the coefficients for each feature across tasks
summary_view_features.comb <- subset(summary_view_features, DataType == "Pathways")

# # One task
# summary_view_features.task <- subset(summary_view_features, task == "CYT")
# median.features <- aggregate(estimate ~ feature + DataType + CancerType, 
#                              FUN = "median", na.rm = T, data = summary_view_features.task)
# median.features <- median.features[order(abs(median.features$estimate), decreasing = TRUE),]
# summary_view_features.task$feature <- factor(summary_view_features.task$feature, 
#                                              levels = unique(median.features$feature))

# Each cancer type with all tasks #
sapply(as.character(unique(summary_view_features$CancerType)), function(X){
  
  summary_view_features.comb.cancer <- subset(summary_view_features.comb, CancerType == X)
  # Sort features by median estimate value
  median.features <- aggregate(estimate ~ feature + DataType + CancerType,
                               FUN = "median", na.rm = T, data = summary_view_features.comb.cancer)
  median.features <- median.features[order(abs(median.features$estimate), decreasing = TRUE),]
  summary_view_features.comb.cancer$feature <- factor(summary_view_features.comb.cancer$feature,
                                                      levels = unique(median.features$feature))
  
  # keep only features that appear more than 50 times
  keep <- as.character(median.features$feature[median.features$estimate !=0])
  tmp <- summary_view_features.comb.cancer[summary_view_features.comb.cancer$feature %in% keep,]
  
  # boxplot
  boxplot <- ggplot(summary_view_features.comb.cancer, aes(x = feature, y = estimate, fill = task, color = task)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(name = "Task",
                      labels = names(colors.tasks),
                      values = colors.tasks) +
    scale_color_manual(name = "Task",
                       labels = names(colors.tasks),
                       values = colors.tasks) +
    theme_minimal() +
    #facet_grid(task ~.) + 
    theme(panel.grid = element_blank()) +
    #coord_flip() + 
    # ylim(c(-0.1,0.9)) +
    #coord_fixed(ratio = 2) +
    geom_hline(yintercept = 0, linetype="dashed") + 
    theme(axis.text.x = element_text(size=6, angle = 55, vjust = 0.75), axis.title.x = element_blank(),
          axis.text.y = element_text(size=8,face="bold"), 
          axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          legend.position = c(0.75, 0.75), legend.direction = "horizontal",
          legend.box.background = element_rect(color="black", size=0.3),
          legend.box.margin = margin(4, 4, 4, 4),
          legend.text=element_text(size=8), 
          legend.title = element_text(size = 10, face="bold", vjust = 0.5)) +
    labs(y = "estimated coefficient") 
  
  # barplot
  median.features.across_iterations <- aggregate(estimate ~ feature + DataType + CancerType, FUN = "median",
                                                 na.rm = T, data = summary_view_features.comb.cancer)
  
  frequency.across_iterations <- aggregate(estimate ~ feature + DataType + CancerType, 
                                           FUN = function(X){sum(X != 0)}, 
                                           data = summary_view_features.comb.cancer)
  
  frequency.across_iterations$estimate <- frequency.across_iterations$estimate/8
  levels(frequency.across_iterations$feature) <- levels(summary_view_features.comb.cancer$feature)
  
  tmp <- frequency.across_iterations[frequency.across_iterations$feature %in% keep,]
  
  barplot <- ggplot(frequency.across_iterations, aes(x = feature, y = estimate)) +
    geom_bar(stat="identity", fill = "red", color="red") +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    # coord_fixed(ratio = 1)+
    #facet_grid(task ~ .) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=6), 
          axis.ticks.y = element_line(size = 0.5),
          legend.position = "none", 
          panel.border = element_blank(), panel.background = element_blank(), 
          axis.line.y = element_line(size=0.5, colour = "black"), 
          aspect.ratio = 0.05) + 
    labs(y = "Frequency", x = "features")
  
  # Combine plots
  g1 <- ggplotGrob(boxplot)
  g2 <- ggplotGrob(barplot)
  g <- rbind(g1, g2, size = "first")
  g$widths <- unit.pmax(g1$widths, g2$widths)
  grid.newpage()
  pdf(paste0("../Figures/PanCancer_draft_v1/Pathways/", sapply(strsplit(X, split = "(", fixed = T), head, 1),
             "_Pathways_1se_mse_cor_model.pdf"), width = 12, height = 12)
  grid.draw(g)
  dev.off()
  
  
})

# All cancer types with all tasks #
summary_view_features.comb$feature <- factor(summary_view_features.comb$feature,
                                             levels = unique(summary_view_features.comb$feature))
## Sort features by median estimate value
median.features <- aggregate(estimate ~ feature + DataType + CancerType + task,
                             FUN = "median", na.rm = T, data = summary_view_features.comb)
median.features <- median.features[order(abs(median.features$estimate), decreasing = TRUE),]

median.features$feature <- factor(median.features$feature,
                                  levels = unique(median.features$feature))
# Heatmap
ggplot(median.features, aes(x = feature, y = CancerType, fill = estimate)) +
  geom_tile(color = "gray", size = 0.25) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "median estimate") +
  theme_minimal() +
  facet_grid(task ~.) + 
  theme(panel.grid = element_blank()) +
  # ylim(c(-0.1,0.9)) +
  #coord_fixed(ratio = 3) +
  theme(axis.text.x = element_text(size=6, angle = 55, vjust = 0.65), axis.title.x = element_blank(),
        axis.text.y = element_text(size=6), axis.title.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "right", legend.direction = "vertical",
        legend.box.background = element_rect(color="black", size=0.3),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.text=element_text(size=8), 
        legend.title = element_text(size = 8, vjust = 0.5),
        strip.background = element_rect(fill = "lightgray", colour = "white"))

ggsave(paste0("../Figures/PanCancer_draft_v1/Pathways/PanCancer_heatmap_1se_mse_model_cor_tasks.pdf"), width = 12, height = 12)


