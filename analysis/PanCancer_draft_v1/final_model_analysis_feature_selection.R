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
view_combinations <- list(views[c(1,3)], views[8])
# ------------------------------------------------------------------------------------------------------------ #
# Collect data from L21 cross-validation results --> kfold = 5
# ------------------------------------------------------------------------------------------------------------ #

# **************** 
# Initialize variable to collect results
summary_view_features <- NULL
algorithms <- c("Multi_Task_EN")
analysis <- c("all","cor")

# for (Cancer in PanCancer.names){
  Cancer = "SKCM"
  # load(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/data/PanCancer/",Cancer,"/new/DataViews_no_filter_", Cancer,".RData"))
  load(paste0("../data/BIM_cluster_presentation/DataViews_no_filter_", Cancer,".RData"))

  # file <- dir(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/output/new/",Cancer), full.names = T, recursive = F)
  file <- dir(path = paste0("../output/BIM_cluster_presentation"), pattern = "all_cv_res_", full.names = T, recursive = F)
  
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
        
          # if (alg == "L21") {tasks <- colnames(tmp_iter_model$cv.MTL.features)}
          if (alg == "Multi_Task_EN")   tmp_iter_model <- tmp_alg[[iteration]]$model
          # if (alg == "Elastic_Net") {tasks <- names(tmp_iter_model$Coef)}
          
          # cv- 1se.mse and min.mse
          summary_cv <- do.call(rbind,lapply(names(tmp_iter_model$cv.glmnet.features), function(cv){
            
            if (alg == "Multi_Task_EN")   tmp_iter_model_cv <- tmp_iter_model$cv.glmnet.features[[cv]]
            
          # # All tasks
          summary_task <- do.call(rbind, lapply(colnames(tmp_iter_model_cv), function(task){
            
            # if (alg == "L21") {
            #   tmp_iter_model_coef_task <- tmp_iter_model$cv.MTL.features$min.mse[,task]
            #   tmp_iter_model_hyp_task <- paste0("min-mse (",tmp_iter_model$cv.MTL.hyperparameters$min.mse$Lam1,",",
            #                                     round(as.numeric(tmp_iter_model$cv.MTL.hyperparameters$min.mse$Lam2),3),")")
            #   feature_names <- as.character(names(tmp_iter_model_coef_task)[2:length(names(tmp_iter_model_coef_task))])
            # }
            # 
            # if (alg == "Elastic_Net") {
            #   tmp_iter_model_coef_task <- tmp_iter_model$Coef[[task]]$min.mse
            #   tmp_iter_model_hyp_task <- paste0("min-mse (", tmp_iter_model$hyperparameters[[task]]$min.mse$alpha, ",",
            #                                     round(as.numeric(tmp_iter_model$hyperparameters[[task]]$min.mse$lambda),3), ")")
            #   feature_names <- as.character(rownames(tmp_iter_model_coef_task)[2:length(rownames(tmp_iter_model_coef_task))])
            # }
            # 
            info <- data.frame(algorithm = alg,
                               iteration = iteration,
                               cv_model = cv,
                               hyp_model = paste0(tmp_iter_model$cv.glmnet.hyperparameters[[cv]]$alpha,",",
                                                     round(as.numeric(tmp_iter_model$cv.glmnet.hyperparameters[[cv]]$lambda),3)),
                               task = task,
                               feature = rownames(tmp_iter_model_cv)[2:nrow(tmp_iter_model_cv)],
                               estimate = tmp_iter_model_cv[2:nrow(tmp_iter_model_cv),task])
            
            return(info)
            
          }))
          
          return(summary_task)
          
        }))
          return(summary_cv)
        
        }))
        
        return(summary_iter)
      }))

      n_algorithms <- length(all_cv_res)
      n_iterations <- length(all_cv_res[[algorithms[1]]])
      n_measures <- length(all_cv_res[[algorithms[1]]][[1]]$performances)
      n_cv <- length(all_cv_res[[algorithms[1]]][[1]]$performances$MSE) 
      n_tasks <- length(all_cv_res[[algorithms[1]]][[1]]$performances$MSE$`1se.mse`)
      n_features <- length(rownames(all_cv_res[[algorithms[1]]][[1]]$model$cv.glmnet.features$`1se.mse`)) -1
      
      summary_alg$AnalysisType <- rep(anal, len = n_algorithms * n_iterations * n_tasks * n_features * n_cv)
      summary_alg$DataType <- rep(input_name, len = n_algorithms * n_iterations * n_tasks  * n_features * n_cv)
      summary_alg$CancerType <- rep(paste0(Cancer,"(n=",nrow(DataViews.no_filter$pathways),")"),
                                     len = n_algorithms * n_iterations * n_features * n_tasks * n_cv)
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
summary_view_features$cv_model <- factor(summary_view_features$cv_model)

# summary_view_features$Analysis_Alg <- paste0(summary_view_features$algorithm,"_",
#                                              summary_view_features$AnalysisType)
# summary_view_features$Analysis_Alg <- factor(summary_view_features$Analysis_Alg)

# Colors for visualization
colors.cancer_types <- toupper(c("#d29d00","#9445cc","#8cce2f","#ef49c3","#1b8c00","#998fff",
                                 "#e78400","#0063a4","#fc4745","#00c98b","#b6006b","#006c43",
                                 "#ff8376","#604588", "#dcc666","#912f47","#7c5b00","#90350e"))

# names(colors.cancer_types) <- levels(summary_view_features$CancerType)

colors.DataType <- toupper(c("#b65ebd","#c88645"))
names(colors.DataType) <- levels(summary_view_features$DataType)

colors.tasks <- toupper(c("#b47645","#ad58c5","#6cb643","#d24787","#52ad7b","#cf4740", "#4bafd0",
                          "#dc7b31","#6776cb","#c1ad46","#b975b1","#6c7b33","#c26671"))
names(colors.tasks) <- levels(summary_view_features$task)

# colors.algorithm <- toupper(c("#ff7433","#853760"))
# names(colors.algorithm) <- levels(summary_view_features$algorithm)

# alpha.algorithm <- c(1,0.2)
# names(alpha.algorithm) <- levels(summary_view_features$algorithm)
# 
# alpha.analysis_alg <- c(1,0.6,0.15)
# names(alpha.analysis_alg) <- levels(summary_view_features$Analysis_Alg)

## ------------------------------------------------------------------------ ##
# 1. PanCancer comparison:
# One plot per task: compare different data types across types of tumors (L21 and EN)
## ----------------- ------------------------------------------------------ ##

# Sort features per median value and all tasks (~ feature + cv_model + dataType + CancerType + algorithm)
## Get median features
median.features <- aggregate(estimate ~ feature + DataType + CancerType, 
                                          FUN = "median", na.rm = T, data = summary_view_features)
## Sort features by median value:
median.features.sort <- median.features[order(abs(median.features$estimate), decreasing = TRUE),]
median.features.sort$feature <- as.character(median.features.sort$feature)

## Restructuring data.frame
# estimates
# summary_view_features.1se.mse <- subset(summary_view_features, cv_model = "1se.mse")
# summary_view_features.sort <- summary_view_features.1se.mse[order(match(summary_view_features.1se.mse$feature, 
#                                                                          median.features.sort$feature)),]
# summary_view_features.sort$feature <- factor(summary_view_features.sort$feature, 
#                                                   levels =  unique(median.features.sort$feature))

## ------------------------------------------------------------------------ ##
# 3. Comparison of running L21 with 8 or 13 tasks:
# a) One plot per cancer type: compare different data types across tasks (L21 and EN)
# b) One plot per task: compare different data types across types of tumors (L21 and EN)
## ------------------------------------------------------------------------ ##

# Per task, DataType
# a)

# sapply(names(colors.DataType)[c(2:4)], function(input){
  
input <- "pathways_immunecells"

  summary_view_features.mech_sig <- subset(summary_view_features.sort, DataType == input)

  # Per task
  sapply(names(colors.tasks), function(output){

    # summary_view_features.mech_sig.task <- subset(summary_view_features.mech_sig, task == output)
      summary_view_features.mech_sig$task <- factor(summary_view_features.mech_sig$task)
      ggplot(summary_view_features.mech_sig, aes(x = feature, y = estimate, fill = task,
                                                       colour = task)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(name = "Task",
                          labels = levels(summary_view_features.mech_sig$task),
                          values = colors.tasks[1:7]) +
        scale_colour_manual(name = "Task",
                            labels = levels(summary_view_features.mech_sig$task),
                            values = colors.tasks[1:7])  +
        # scale_alpha_manual(name = "Task",
        #                    labels = levels(summary_view_features.mech_sig$task),
        #                    values = colors.tasks[1:7])  +
        theme_minimal() +
        scale_x_discrete(labels = function(x) format(x, width = 5)) +
        theme(axis.text.x = element_text(size=10,face="bold", angle = 45, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
              axis.text.y = element_text(size=10,face="bold"), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_line(size = 1),
              axis.title.y = element_text(size=10,face="bold",vjust = 0.9), strip.background = element_blank(), panel.spacing.y =  unit(1, "lines"),
              legend.position = c(0.75,0.85), legend.direction = "horizontal",
              legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
              panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(size=0.5, colour = "black")) +
        labs(y = "estimated coefficient") + 
        ggtitle(paste0(input, ": EN feature selection for ", output))

    ggsave(paste0("../Figures/BIM_cluster_presentation/",), width = 12, height = 12)

  })

#})

# b)
# Per CancerType
sapply(names(colors.cancer_types), function(Cancer){
  
  # Per dataType
  sapply(names(colors.DataType), function(input){
    
    # Per task
      summary_view_features.mech_sig.cancer <- subset(summary_view_features, DataType == "pathways_immunecells" &
                                                      CancerType == "SKCM(n=467)" & cv_model == "1se.mse" & AnalysisType == "cor" & task == "CYT")
      median.features <- aggregate(estimate ~ feature + iteration + DataType, FUN = "median", na.rm = T, data = summary_view_features.mech_sig.cancer)
      median.features.iterations <- aggregate(estimate ~ feature, FUN = "median", na.rm = T, data = summary_view_features.mech_sig.cancer)
      
      frequency.iterations <- aggregate(estimate ~ feature + task, FUN = function(X){sum(X != 0)}, data = summary_view_features.mech_sig.cancer)
      ## Sort features by median value:
      median.features.iterations <- median.features.iterations[order(abs(median.features.iterations$estimate), decreasing = TRUE),]
      remove_features <- median.features.iterations$feature[which(median.features.iterations$estimate == 0)]
      # frequency.iterations <- frequency.iterations[!frequency.iterations$feature %in% remove_features,]
                        
      non_zero_summary <- summary_view_features.mech_sig.cancer[!summary_view_features.mech_sig.cancer$feature %in% remove_features,]
      
      summary_view_features.mech_sig.cancer$feature <- factor(summary_view_features.mech_sig.cancer$feature,
                                                              levels = unique(median.features.iterations$feature))
      
      non_zero_summary$feature <- factor(non_zero_summary$feature, levels = unique(median.features.iterations$feature))
      # One task
      summary_view_features.mech_sig.cancer.task <- subset(summary_view_features.mech_sig.cancer, task == "CYT")
      median.features <- aggregate(estimate ~ feature + DataType + task + AnalysisType, FUN = "median", na.rm = T, data = summary_view_features.mech_sig.cancer.task)
      median.features <- median.features[order(abs(median.features$estimate), decreasing = TRUE),]
      
      summary_view_features.mech_sig.cancer.task$feature <- factor(summary_view_features.mech_sig.cancer.task$feature, 
                                                                   levels = as.character(median.features$feature))
      ggplot(summary_view_features.mech_sig.cancer.task, aes(x = feature, y = estimate, fill = task, color = task)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(name = "Task",
                          labels = "Cytolytic Activity",
                          values = "#999a3e") +
        scale_color_manual(name = "Task",
                           labels = "Cytolytic Activity",
                           values = "#999a3e") +

        theme_minimal() +
        # ylim(c(-0.1,0.9)) +
        # coord_fixed(ratio = 2) +
        geom_hline(yintercept = 0, linetype="dashed") + 
        scale_x_discrete(labels = c("T_cells_CD8" = "CD8 T cells",
                                    "Macrophages_M2"= "M2",
                                    "B_cells" = "B cells",
                                    "T_cells_regulatory_Tregs" = "Tregs",
                                    "Macrophages_M1"= "M1",
                                    "T_cells_CD4" = "CD4 T cells",
                                    "NK_cells" = "NK cells",
                                    "Dendritic_cells" = "DC cells")) +
        theme(axis.text.x = element_text(size=20,face="bold", angle = 45, vjust = 0.5, hjust=0.5), axis.title.x = element_text(size=16,face="bold",vjust = 0.9),
              axis.text.y = element_text(size=20,face="bold"), axis.ticks.x = element_blank(), axis.ticks.y = element_line(size = 1),
              axis.title.y = element_text(size=16,face="bold",vjust = 0.9), strip.background = element_blank(), panel.spacing.y =  unit(1, "lines"),
              legend.position = c(0.85,0.8), legend.direction = "vertical",
              legend.text=element_text(size=18), legend.title = element_text(size = 12, face="bold", vjust = 0.5),
              panel.border = element_blank(), panel.background = element_blank(), axis.line = element_blank()) + 
        labs(y = "estimated coefficient", x = "features") 
        #ggtitle(paste0("Pathways + ImmuneCells cross-validation coefficients (min MSE + 1SE model)"))
    
    ggsave(paste0("../Figures/BIM_cluster_presentation/Boxplot_displaying_coef_pathways_immunecells_1se_mse_model_cor_CYT.pdf"), width = 12, height = 12)
   
    frequency.iterations.task <- subset(frequency.iterations, task == "CYT")
    ggplot(frequency.iterations.task, aes(x = feature, y = estimate, fill = task, color = task)) +
      geom_bar(stat="identity", color="red") +
      # geom_text(aes(label= sign.feature), stat = "identity", color="black", size = 1, 
      #           position = position_stack(vjust = 0.5)) +
      scale_fill_manual(name = "Task", 
                        labels = "CYT",
                        values = "red",
                        guide = F) +
      scale_color_manual(name = "Task", 
                         labels = "CYT",
                         values = "red",
                         guide = F) + 
      theme_minimal() +
      # coord_fixed(ratio = 1)+
      # facet_grid(task ~ .) +
      scale_x_discrete(labels = c("T_cells_CD8" = "CD8 T cells",
                                  "Macrophages_M2"= "M2",
                                  "B_cells" = "B cells",
                                  "T_cells_regulatory_Tregs" = "Tregs",
                                  "Macrophages_M1"= "M1",
                                  "T_cells_CD4" = "CD4 T cells",
                                  "NK_cells" = "NK cells",
                                  "Dendritic_cells" = "DC cells")) +
      theme(axis.text.x = element_blank(), axis.title.x =  element_blank(),
            axis.text.y = element_text(size=12,face="bold"), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_line(size = 1),
            axis.title.y = element_text(size=12,face="bold",vjust = 0.9), strip.background = element_blank(), panel.spacing.y =  unit(1, "lines"),
            legend.position = "bottom", legend.direction = "horizontal",  strip.text.y = element_text(angle = 360),
            legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
            panel.border = element_blank(), panel.background = element_blank(), 
            #axis.line = element_line(size=0.5, colour = "black"), 
            aspect.ratio = 0.05) + 
      labs(y = "Number of\n times") +
      # ggtitle(paste0(input," - ", Cancer,": feature selection all tasks in min MSE + 1SE model"), 
      #         subtitle = "* Balance classes of directions due to coefficients close to zero are resolved with sign = 0")
    
    ggsave(paste0("../figures/BIM_cluster_presentation/Features_stability_pathways_immunecells_cor_tasks.pdf"), width = 12, height = 12)
    
    
    
    
    
  })
})


# c)
# Stability selection
# color.analysis_alg <- c("#d8005e","#a09f00","#ffad72")
# names(color.analysis_alg) <- levels(summary_view_features_tmp.order$Analysis_Alg)

summary_view_features.mech_sig.cancer <- subset(summary_view_features, DataType == "CYTOKINEpairs" &
                                                  CancerType == "SKCM(n=467)" & cv_model == "1se.mse" & AnalysisType == "cor")
median.features <- aggregate(estimate ~ feature, FUN = "median", na.rm = T, data = summary_view_features.mech_sig.cancer)
## Sort features by median value:
median.features.sort <- median.features[order(abs(median.features$estimate), decreasing = TRUE),]
median.features.sort$feature <- as.character(median.features.sort$feature)
features_order <- median.features.sort$feature
summary_view_features.mech_sig.cancer$feature <- factor(summary_view_features.mech_sig.cancer$feature, levels = unique(features_order))

frequency <- aggregate(estimate ~ feature + cv_model + task + DataType + AnalysisType + CancerType, FUN = function(X){sum(X != 0)},
                       data = summary_view_features.mech_sig.cancer)

sign.estimate <- aggregate(estimate ~ feature + cv_model + task  + DataType + AnalysisType + CancerType, FUN = function(X){sign(median(X))},
                           data = summary_view_features.mech_sig.cancer)

frequency <- cbind(frequency, sign.feature = sign.estimate$estimate)
frequency$sign.feature <- gsub(-1,"-", frequency$sign.feature)
frequency$sign.feature <- gsub(1,"+", frequency$sign.feature)
frequency$sign.feature <- factor(frequency$sign.feature)

# # Issue: Balance classes of directions due to coefficients close to zero
# table.tmp <- table(sign(try$estimate))
# sign <- ifelse(table.tmp[c(1:3)] > sum(table.tmp[c(1:3)])*0.8,1, 0)
# direction <- ifelse(any(sign == 1), names(sign)[which(sign == 1)], 0)

# # Per CancerType
# Cancer = "SKCM(n=467)"
# # Per DataType
input = "CYTOKINEpairs"
# 
# frequency.cancer.input <- subset(frequency, DataType == input & CancerType == Cancer & cv_model == "1se.mse")

ggplot(frequency, aes(x = feature, y = estimate, fill = task, color = task)) +
  geom_bar(stat="identity", color="black") +
  geom_text(aes(label= sign.feature), stat = "identity", color="black", size = 1, 
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(name = "Task", 
                    labels = names(colors.tasks),
                    values = as.character(colors.tasks),
                    guide = F) +
  scale_color_manual(name = "Task", 
                     labels = names(colors.tasks),
                     values = as.character(colors.tasks),
                     guide = F) + 
  theme_minimal() +
  coord_fixed(ratio = 2)+
  facet_grid(task ~ .) +
  scale_x_discrete(labels = c("T_cells_CD8" = "CD8 T cells",
                              "Macrophages_M2"= "M2",
                              "B_cells" = "B cells",
                              "T_cells_regulatory_Tregs" = "Tregs",
                              "Macrophages_M1"= "M1",
                              "T_cells_CD4" = "CD4 T cells",
                              "NK_cells" = "NK cells",
                              "Dendritic_cells" = "DC cells")) +
  theme(axis.text.x = element_text(size=7,face="bold", angle = 90, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size=9,face="bold"), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_line(size = 1),
        axis.title.y = element_text(size=10,face="bold",vjust = 0.9), strip.background = element_blank(), panel.spacing.y =  unit(1, "lines"),
        legend.position = "bottom", legend.direction = "horizontal",  strip.text.y = element_text(angle = 360),
        legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
        panel.border = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=0.5, colour = "black"), aspect.ratio = 0.05) + 
  labs(y = "Number of\n times") +
  ggtitle(paste0(input," - ", Cancer,": feature selection all tasks in min MSE + 1SE model"), 
          subtitle = "* Balance classes of directions due to coefficients close to zero are resolved with sign = 0")

ggsave(paste0("../figures/BIM_cluster_presentation/Features_stability_CYTOKINEpairs_all_tasks.pdf"), width = 12, height = 12)


## ---------------------------------------------------------------------------------------------- ##
# 2. Combo improve single data
# Heatmap: showing coefficients when considering each individual feature set (e.g. pathways and
# immune cells separately) and then when considering them together. Additional heatmap can be used 
# to show the difference between looking at features individually or together.
## ----------------- -----------------------------------------------------------------------------##

color.combos <- c("#ad9fff","#009841")
names(color.combos) <- c("single", "combo")

# Per CancerType
Cancer = "SKCM(n=467)"
# Per DataType combinations (Combo)
input <- levels(frequency$DataType)

frequency.cancer.input <- subset(frequency, DataType %in% input & CancerType == Cancer & cv_model == "1se.mse")
frequency.cancer.input$Combo <- frequency.cancer.input$DataType
frequency.cancer.input$Combo <- gsub("pathways_immunecells","combo", frequency.cancer.input$Combo, fixed = TRUE)
frequency.cancer.input$Combo <- gsub("immunecells","single", frequency.cancer.input$Combo, fixed = TRUE)
frequency.cancer.input$Combo <- gsub("pathways","single", frequency.cancer.input$Combo, fixed = TRUE)
# frequency.cancer.input$Combo <- gsub("TFs","separate", frequency.cancer.input$Combo, fixed = TRUE)

frequency.cancer.input$Combo <- factor(frequency.cancer.input$Combo)

ggplot(frequency.cancer.input, aes(x = feature, y = estimate, fill = Combo)) +
  geom_bar(stat="identity", color="black", position = position_dodge()) +
  geom_text(aes(label= sign.feature), stat = "identity", color="black", size = 3, 
            position = position_dodge(width = 1)) +
  scale_fill_manual(name = "DataType alone or in combination",
                     labels = names(color.combos),
                     values = color.combos)  +
  theme_minimal() +
  facet_grid(task ~ .) +
  scale_x_discrete(labels = function(x) substr(x,1,18)) +
  theme(axis.text.x = element_text(size=9,face="bold", angle = 45, vjust = 0.5, hjust=0.5), axis.title.x = element_blank(),
        axis.text.y = element_text(size=9,face="bold"), axis.ticks.x = element_line(size = 1), axis.ticks.y = element_line(size = 1),
        axis.title.y = element_text(size=10,face="bold",vjust = 0.9), strip.background = element_blank(), panel.spacing.y =  unit(1, "lines"),
        legend.position = "bottom", legend.direction = "horizontal",  strip.text.y = element_text(angle = 360),
        legend.text=element_text(size=9), legend.title = element_text(size = 10, face="bold", vjust = 0.5),
        panel.border = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=0.5, colour = "black"), aspect.ratio = 0.05) + 
  labs(y = "Number of\n times") +
  ggtitle(paste0("Pathways, ImmuneCells"," - ", Cancer,": feature selection all tasks in min MSE + 1SE model"), 
          subtitle = "* Balance classes of directions due to coefficients close to zero are resolved with sign = 0")

ggsave(paste0("../figures/BIM_cluster_presentation/Features_stability_Pathways_ImmmuneCells_all_tasks.pdf"), width = 12, height = 12)


