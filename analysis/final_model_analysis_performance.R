# #########################################################################################################
# Script to plot model performances obtained from cross-validation -->
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
# setwd("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop")

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

## Immune cells and spatial info ## 
# view_combinations <- list(views[3],views[c(3,6)])

## TFs ##
#view_combinations <- list(views[4])

## DoRothEAv1, PROGENy, PROGENy_quanTIseq, quanTIseq, quanTIseq_spatialTILs##
# view_combinations <- list(views[4], views[1], views[c(1,3)], views[3],views[c(3,6)])

# Pathways, immune cells
# view_combinations <- list(views[c(1,3)], views[1], views[3])

# All: comparison
view_combinations <- list(views[c(1,3)], views[1], views[3], views[4], views[7], views[5])

# ------------------------------------------------------------------------------------------------------------ #
# Collect data from L21 cross-validation results --> kfold = 5
# ------------------------------------------------------------------------------------------------------------ #

# **************** 
# Initialize variable to collect results
summary_view_all <- NULL
algorithms <- c("Elastic_Net")
# algorithms <- c("L21")
analysis <- c("all")

#for (Cancer in PanCancer.names){
  Cancer = "SKCM"
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

            tmp_iter <- tmp_alg[[iteration]]$performances

            # 5 metrics
            summary_measure <- do.call(rbind, lapply(names(tmp_iter), function(measure){

              tmp_measure <- tmp_iter[[measure]]

              # all tasks
              summary_task <- do.call(rbind,lapply(names(tmp_measure), function(task){

                tmp_task <- tmp_measure[[task]]

                if (alg == "Elastic_Net") {
                  info <- data.frame(algorithm = alg,
                                     iteration = iteration,
                                     metric = measure,
                                     task = task,
                                     perf_min = tmp_task$min.mse)

                }else {
                  info <- data.frame(algorithm = alg,
                                     iteration = iteration,
                                     metric = measure,
                                     task = task,
                                     perf_min = tmp_task)
                }

                return(info)
              }))
              summary_subtask <- subset(summary_task, task %in% c("CYT","IS","RohIS","chemokine","IS_Davoli", "IFny",
                                                                 "ExpandedImmune","T_cell_inflamed"))
              common <- data.frame(algorithm = rep(alg, times = 2),
                                   iteration = rep(iteration, times = 2),
                                   metric = rep(measure, times = 2),
                                   task = c("common_all", "common_top"),
                                   perf_min = c(median(summary_task$perf_min), median(summary_subtask$perf_min)))
              
              summary_task <- rbind(summary_task, common)
              
              return(summary_task)
            }))
            return(summary_measure)
          }))
          return(summary_iter)
        }))

      n_views <- length(view_combinations)
      n_algorithms <- length(all_cv_res)
      n_iterations <- length(all_cv_res$Elastic_Net)
      n_measures <- length(all_cv_res$Elastic_Net[[1]]$performances)
      n_tasks <- length(all_cv_res$Elastic_Net[[1]]$performances$MSE) + 2

      summary_alg$AnalysisType <- rep(anal, len = n_algorithms * n_iterations * n_measures * n_tasks)
      summary_alg$DataType <- rep(input_name, len = n_algorithms * n_iterations * n_measures * n_tasks)
      summary_alg$CancerType <- rep(paste0(Cancer,"(n=",nrow(DataViews.no_filter$pathways),")"),
                                    len =  n_algorithms * n_iterations * n_measures * n_tasks)

      return(summary_alg)
    }))

    return(summary_view)
 }))

  summary_view_all <- rbind(summary_view_all, summary_analysis)
#}

# ------------------------------------------------------------------------------------------------------------ #
# Collect data from other algorithms cross-validation results --> kfold = 5
# ------------------------------------------------------------------------------------------------------------ #

# # **************** 
# # Initialize variable to collect results
# summary_view_all.rest <- NULL
# algorithms <- c("CV_linear_reg_L1&L2", "BEMKL","Lasso", "Elastic_Net")
# 
# for (Cancer in PanCancer.names){
#   
#   load(paste0("../data/PanCancer/", Cancer, "/DataViews_", Cancer,"_all_train_with_sTIL.RData"))
#   file <- dir(paste0("../output/model/",Cancer), pattern = paste0("_kfold"), full.names = T, recursive = F)
#   
#   summary_view <- do.call(rbind, lapply(1:length(view_combinations), function(view){
#     
#     # Load file
#     input_name <- paste(names(view_combinations[[view]]), collapse ="_")
#     which_file <- grep(pattern = paste0("_withConsensus_", input_name, "_kfold_5",".RData"), file, fixed = T)
#     load(file[which_file])
#     
#     # 2 algorithms
#     summary_alg <- do.call(rbind, lapply(algorithms, function(alg){
#       
#       tmp_alg <- all_cv_res[[alg]]
#       
#       # 100 iterations
#       summary_iter <- do.call(rbind, lapply(1:length(tmp_alg), function(iteration){
#         
#         tmp_iter <- tmp_alg[[iteration]]$performances
#         
#         # 5 metrics
#         summary_measure <- do.call(rbind, lapply(names(tmp_iter), function(measure){
#           
#           tmp_measure <- tmp_iter[[measure]]
#           
#           # 5 tasks
#           summary_task <- do.call(rbind,lapply(names(tmp_measure), function(task){
#             
#             tmp_task <- tmp_measure[[task]]
#             
#             if (alg == "BEMKL") {
#               info <- data.frame(algorithm = alg,
#                                  iteration = iteration,
#                                  metric = measure,
#                                  task = task,
#                                  perf_min = tmp_task,
#                                  perf_min_1SE = NA)
#               
#             }else {
#               info <- data.frame(algorithm = alg,
#                                  iteration = iteration,
#                                  metric = measure,
#                                  task = task,
#                                  perf_min = tmp_task$min.mse,
#                                  perf_min_1SE = tmp_task$`1se.mse`)
#             }
#             return(info)
#           }))
#           return(summary_task)
#         }))
#         return(summary_measure)
#       }))
#       return(summary_iter)
#     }))
#     
#     n_views <- length(view_combinations)
#     n_algorithms <- length(algorithms)
#     n_iterations <- length(all_cv_res$Lasso)
#     n_measures <- length(all_cv_res$Lasso[[1]]$performances)
#     n_tasks <- length(all_cv_res$Lasso[[1]]$performances$MSE)
#     
#     summary_alg$DataType <- rep(input_name, len = n_algorithms * n_iterations * n_measures * n_tasks)
#     summary_alg$CancerType <- rep(paste0(Cancer,"(n=",nrow(DataViews$PROGENy),")"),
#                                   len =  n_algorithms * n_iterations * n_measures * n_tasks)
#     
#     return(summary_alg)
#   }))
#   
#   summary_view_all.rest <- rbind(summary_view_all.rest, summary_view)
# }

##############################################################################
# Visualization
##############################################################################

# Theoretically, we cannot assume that the response follows a normal distribution, so it would be better 
# to use spearman correlation as metric

summary_view_all.SpCorr <- subset(summary_view_all, metric == "SpCorr")
summary_view_all.SpCorr$DataType <- factor(summary_view_all.SpCorr$DataType)
summary_view_all.SpCorr$CancerType <- factor(summary_view_all.SpCorr$CancerType)
summary_view_all.SpCorr$task <- factor(summary_view_all.SpCorr$task)
summary_view_all.SpCorr$AnalysisType <- factor(summary_view_all.SpCorr$AnalysisType)

summary_view_all.SpCorr$Analysis_Alg <- paste0(summary_view_all.SpCorr$algorithm,"_",
                                                    summary_view_all.SpCorr$AnalysisType)

summary_view_all.SpCorr$Analysis_Alg <- factor(summary_view_all.SpCorr$Analysis_Alg)

# Colors for visualization
colors.cancer_types <- toupper(c("#d29d00","#9445cc","#8cce2f","#ef49c3","#1b8c00","#998fff",
                                 "#e78400","#0063a4","#fc4745","#00c98b","#b6006b","#006c43",
                                 "#ff8376","#604588", "#dcc666","#912f47","#7c5b00","#90350e"))

names(colors.cancer_types) <- levels(summary_view_all.SpCorr$CancerType)

colors.DataType <- toupper(c("#8c9f3e","#b65cbf","#4eac7c","#c95574","#747fca","#ca743e"))
names(colors.DataType) <- levels(summary_view_all.SpCorr$DataType)

colors.tasks <- toupper(c("#b47645","#ad58c5","#6cb643","#d24787","#52ad7b","#cf4740", "#4bafd0",
                          "#dc7b31","#6776cb","#c1ad46","#b975b1","#6c7b33","#c26671"))
names(colors.tasks) <- levels(summary_view_all.SpCorr$task)

colors.algorithm <- toupper(c("#ff7433","#853760"))
names(colors.algorithm) <- levels(summary_view_all.SpCorr$algorithm)

alpha.algorithm <- c(1,0.2)
names(alpha.algorithm) <- levels(summary_view_all.SpCorr$algorithm)

alpha.analysis_alg <- c(1,0.4,0.05)
names(alpha.analysis_alg) <- levels(summary_view_all.SpCorr$Analysis_Alg)

## ------------------------------------------------------------------------ ##
# 1. PanCancer comparison:
# One plot per task: compare different data types across types of tumors (L21 and EN)
## ----------------- ------------------------------------------------------ ##


# Per CancerType
# sapply(names(colors.tasks), function(ImmuneResponse){
#   
summary_view_all.SpCorr.task <- subset(summary_view_all.SpCorr, task %in% c("common_top") & 
                                       DataType %in% c("pathways", "immunecells", "pathways_immunecells", "transcript"))
  
ggplot(summary_view_all.SpCorr.task, aes(x = DataType, y = perf_min, fill = DataType,
                                           colour = DataType, alpha = task)) +
    geom_boxplot() +
    scale_fill_manual(name = "Mechanistic signatures", 
                      labels = names(colors.DataType[c(1,3,4,6)]),
                      values = colors.DataType[c(1,3,4,6)]) +
    scale_color_manual(name = "Mechanistic signatures",
                       labels = names(colors.DataType[c(1,3,4,6)]),
                       values = colors.DataType[c(1,3,4,6)]) +
    scale_alpha_manual(name = "Task",
                     labels = c("common_top"),
                     values = c(0.5)) +
    theme_bw() +
    ylim(c(0.6,1)) +
    #coord_fixed(ratio = 4) +
    #facet_grid(algorithm ~ .) +
    scale_x_discrete(labels = c("immunecells"= "ImmuneCells",
                                "pathways"="Pathways",
                                "pathways_immunecells"="Pathways \n + \n Immunecells",
                                "transcript" = "Transcriptomics")) +
    theme(axis.text.x = element_text(size=14,face="bold", angle = 0, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
          axis.text.y = element_text(size=14,face="bold"), axis.ticks.x = element_line(size=1), axis.ticks.y = element_line(size=1),
          axis.title.y = element_text(size=14,face="bold",vjust = 0.9), strip.background = element_blank(), 
          legend.position = "right") +
          # legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
          # panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(size=0.5, colour = "black")) +
    labs(y = "Spearman Correlation") 
    # ggtitle(paste0("PanCancer comparison of mechanistic signatures performance using ", ImmuneResponse, "(all tasks)"))
  
  ggsave(paste0("../figures/Federica_presentation_colab/PROGENy_updated/Elastic_Net_mechanistic_signatures_perfomance_common_across_top_tasks.pdf"), width = 12, height = 12)
  
#})

## ------------------------------------------------------------------------ ##
# 2. Cancer specific comparison:
# One plot per cancer type: compare different data types across tasks (L21 and EN)
## ----------------- ------------------------------------------------------ ##

# Per CancerType
sapply(names(colors.cancer_types), function(Cancer){

  summary_view_all.SpCorr.Cancer <- subset(summary_view_all.SpCorr, CancerType == Cancer)

  ggplot(summary_view_all.SpCorr.Cancer, aes(x = task, y = perf_min, fill = DataType,
                                             colour = DataType, alpha = Analysis_Alg)) +
  geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(name = "Mechanistic signatures", 
                      labels = names(colors.DataType[1:5]),
                      values = colors.DataType[1:5]) +
    scale_color_manual(name = "Mechanistic signatures", 
                       labels = names(colors.DataType[1:5]),
                       values = colors.DataType[1:5]) + 
    scale_alpha_manual(name = "Algorithms", 
                       labels = names(alpha.analysis_alg),
                       values = alpha.analysis_alg)   +
    theme_minimal() +
    #ylim(c(0,1)) +
    coord_fixed(ratio = 4) +
    #facet_grid(algorithm ~ .) +
    scale_x_discrete(labels = function(x) substr(x,1,13)) +
    theme(axis.text.x = element_text(size=10,face="bold", angle = 45, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
          axis.text.y = element_text(size=10,face="bold"), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_line(size = 1),
          axis.title.y = element_text(size=10,face="bold",vjust = 0.9), strip.background = element_blank(), panel.spacing.y =  unit(1, "lines"),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
          panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(size=0.5, colour = "black")) +
    labs(y = "SpCorr") + 
  ggtitle(paste0("Comparison mechanistic signatures performance across all tasks for ",Cancer))
  
ggsave(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/new_v2/", 
              sapply(strsplit(Cancer, split = "(", fixed = T), function(X) {return(X[1])}), 
              "/Performance/", sapply(strsplit(Cancer, split = "(", fixed = T), function(X) {return(X[1])}),
              "_EN_vs_L21_with_all_datatypes_tasks.pdf"), width = 12, height = 12)

})

## ------------------------------------------------------------------------ ##
# 3. Comparison of running L21 with 8 or 13 tasks:
# a) One plot per cancer type: compare different data types across tasks (L21 and EN)
# b) One plot per task: compare different data types across types of tumors (L21 and EN)
## ------------------------------------------------------------------------ ##
# summary_view_all.SpCorr.anal.all_vs_top <- subset(summary_view_all.SpCorr, algorithm == "L21")
# 
# # Per task
# # a)
# sapply(names(colors.tasks), function(ImmuneResponse){
#   
#   summary_view_all.SpCorr.task <- subset(summary_view_all.SpCorr.anal.all_vs_top, task == ImmuneResponse)
#   
#   ggplot(summary_view_all.SpCorr.task, aes(x = CancerType, y = perf_min, fill = DataType,
#                                            colour = DataType, alpha = Analysis_Alg)) +
#     geom_boxplot(outlier.shape = NA) +
#     scale_fill_manual(name = "Mechanistic signatures", 
#                       labels = names(colors.DataType[1:5]),
#                       values = colors.DataType[1:5]) +
#     scale_color_manual(name = "Mechanistic signatures", 
#                        labels = names(colors.DataType[1:5]),
#                        values = colors.DataType[1:5]) + 
#     scale_alpha_manual(name = "Analysis", 
#                        labels = names(alpha.analysis_alg)[2:3],
#                        values = alpha.analysis_alg[2:3])  +
#     theme_minimal() +
#     #ylim(c(0,1)) +
#     coord_fixed(ratio = 4) +
#     #facet_grid(algorithm ~ .) +
#     scale_x_discrete(labels = function(x) substr(x,1,13)) +
#     theme(axis.text.x = element_text(size=10,face="bold", angle = 45, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
#           axis.text.y = element_text(size=10,face="bold"), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_line(size = 1),
#           axis.title.y = element_text(size=10,face="bold",vjust = 0.9), strip.background = element_blank(), panel.spacing.y =  unit(1, "lines"),
#           legend.position = "bottom", legend.direction = "horizontal",
#           legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
#           panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(size=0.5, colour = "black")) +
#     labs(y = "SpCorr") + 
#     ggtitle(paste0("PanCancer mechanistic signatures performance using ", ImmuneResponse, 
#                    " (all vs all_top tasks)"))
#   
#   ggsave(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/new_v2/PanCancer/Performance/", 
#                 "PanCancer_EN_vs_L21_with_all_datatypes_", ImmuneResponse ,"_all_vs_top_tasks.pdf"), width = 12, height = 12)
#   
# })

# b)
# Per CancerType
sapply(names(colors.cancer_types), function(Cancer){
  
  summary_view_all.SpCorr.Cancer <- subset(summary_view_all.SpCorr.anal.all_vs_top, CancerType == Cancer)
  
  ggplot(summary_view_all.SpCorr.Cancer, aes(x = task, y = perf_min, fill = DataType,
                                             colour = DataType, alpha = Analysis_Alg)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(name = "Mechanistic signatures", 
                      labels = names(colors.DataType[1:5]),
                      values = colors.DataType[1:5]) +
    scale_color_manual(name = "Mechanistic signatures", 
                       labels = names(colors.DataType[1:5]),
                       values = colors.DataType[1:5]) + 
    scale_alpha_manual(name = "Analysis", 
                       labels = names(alpha.analysis_alg)[2:3],
                       values = alpha.analysis_alg[2:3])  +
    theme_minimal() +
    #ylim(c(0,1)) +
    coord_fixed(ratio = 4) +
    #facet_grid(algorithm ~ .) +
    scale_x_discrete(labels = function(x) substr(x,1,13)) +
    theme(axis.text.x = element_text(size=10,face="bold", angle = 45, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
          axis.text.y = element_text(size=10,face="bold"), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_line(size = 1),
          axis.title.y = element_text(size=10,face="bold",vjust = 0.9), strip.background = element_blank(), panel.spacing.y =  unit(1, "lines"),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
          panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(size=0.5, colour = "black")) +
    labs(y = "SpCorr") + 
    ggtitle(paste0("Mechanistic signatures performance - all vs all_top tasks for ",Cancer))
  
  ggsave(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/new_v2/", 
                sapply(strsplit(Cancer, split = "(", fixed = T), function(X) {return(X[1])}), 
                "/Performance/", sapply(strsplit(Cancer, split = "(", fixed = T), function(X) {return(X[1])}),
                "_EN_vs_L21_with_all_datatypes_all_vs_top_tasks.pdf"), width = 12, height = 12)
  
})
