# #########################################################################################################
# Script to visualize the data
# #########################################################################################################

# ****************
# working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ****************
# packages
library(ggplot2)
library(reshape)
library(corrplot)
library(matrixStats)
library(grid)
library(tidyverse)
library(caret)

# ****************
# Select cancer type
#PanCancer.names = c("BLCA",  "BRCA","CESC", "CRC", "LUAD", "LUSC", "PRAD", "SKCM", "STAD")
PanCancer.names = c("CRC", "LUAD", "SKCM")
# ****************
# views
views <- c(PROGENy = 'gaussian', #1
           Protall = 'gaussian', #2
           quanTIseq = 'gaussian', #3
           DoRothEAv1 = 'gaussian', #4
           transcript = 'gaussian', #5
           SpatialTIL = 'gaussian',  #6
           LRpairs = 'gaussian', #7
           CYTOKINEpairs = 'gaussian')  #8


# ****************
# functions from Federica
# Plot correlations including p-value, correlation value and correlation line 
panel.cor <- function(x, y, digits=2, font.cor = 1, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor.test(x,y)$estimate
  p <- cor.test(x,y)$p.value
  txt_r <- format(r, digits=digits)
  txt_p <- format(p, scientific = TRUE, digits=digits)
  txt <- paste("cor=", txt_r, "\np=", txt_p, sep="")
  
  if(txt_r >= 0.7 & txt_p >= 0.05) font.cor <- 2
  
  text(0.5, 0.5, txt, cex = 1, font = font.cor)
}

panel.lm <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 0.8, col.smooth = "#A1A1A1", ...) {
  
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  # abline(stats::lm(y ~ x),  col = col.smooth, ...)
  abline(a=0, b=1,  col = col.smooth, ...)
}

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


# ****************
# views combinations

# # quanTIseq_spatialTILs
#view_combinations <- list(views[3],views[c(3,6)])

# PROGENy_quanTIseq
#view_combinations <- list(views[1], views[c(1,3)])

# PROGENy_DoRothEAv1
#view_combinations <- list(views[1], views[c(1,4)])

# PROGENy_LRpairs
#view_combinations <- list(views[1], views[c(1,7)])

# PROGENy_CYTOKINEpairs
#view_combinations <- list(views[1], views[c(1,8)])

# quanTIseq_LRpairs
#view_combinations <- list(views[3], views[c(3,7)])

# quanTIseq_CYTOKINEpairs
#view_combinations <- list(views[3], views[c(3,8)])

# LRpairs
#view_combinations <- list(views[7])

# CYTOKINEpairs
view_combinations <- list(views[8])
# ------------------------------------------------------------------------------------------------------------ #
# Analysis of spatial information features correlation
# ------------------------------------------------------------------------------------------------------------ #

# Check: spatialTILs correlation
# All cancer types together
comb <- do.call(rbind, lapply(PanCancer.names, function(Cancer){

  # Load data
  load(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/data/PanCancer/",
              Cancer, "/DataViews_", Cancer,"_all_train_with_sTIL.RData"))
  
  tmp_data <- DataViews[[names(view_combinations[[1]])]]
  # Remove NA values in Dataviews LRpairs and CYTOKINE pairs
  sum <- apply(tmp_data,2, sum)
  remove_NA <- as.numeric(na.action(na.omit(sum)))
  tmp_data <- tmp_data[,-remove_NA]
  
  return(tmp_data)
  
}))
  
  # Remove NA values in Dataviews LRpairs and CYTOKINE pairs
  # LR_sum <- apply(DataViews$LRpairs,2, sum)
  # remove_NA_LR_pairs <- as.numeric(na.action(na.omit(LR_sum)))
  # DataViews$LRpairs <- DataViews$LRpairs[,-remove_NA_LR_pairs]
  # Cyt_sum <- apply(DataViews$CYTOKINEpairs,2, sum)
  # remove_NA_Cyt_pairs <- as.numeric(na.action(na.omit(Cyt_sum)))
  # DataViews$CYTOKINEpairs <- DataViews$CYTOKINEpairs[,-remove_NA_Cyt_pairs]
  
  # input_name <- paste(names(view_combinations[[2]]), collapse ="_")
  # one <- names(view_combinations[[2]][1])
  # two <-  names(view_combinations[[2]][2])
  # 
  # # dataframe
  # data <- data.frame(x = cbind(DataViews[[one]], DataViews[[two]]))
  # colnames(data) <- substr(colnames(data), 3, 2000)
  
  data <- comb
  
  # remove zero variance variables
  zero_var <- apply(as.matrix(data),2, var)
  if (any(zero_var == 0)) data <- data[,-which(zero_var == 0)]
  hist(zero_var)
  
  # matrix of the p-value of the correlation
  p.mat <- cor.mtest(data)
  input_name = paste0("combo_cancers_", names(view_combinations[[1]]))
  # pdf(file = paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/pre-processing/DataTypes_correlation/", 
  #     input_name, "/", input_name,"_for_", Cancer,"_corrplot.pdf"),height = 18, width = 18)
  
  pdf(file = paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/pre-processing/DataTypes_correlation/", 
                    input_name, "/",input_name,"_corrplot.pdf"),height = 18, width = 18)
  
  corrplot(cor(data, method = "pearson"), method = "color", number.cex = 0.55, tl.cex=0.55,
           order = "hclust",
           p.mat = p.mat, sig.level = 0.05, insig = c("blank"),
           na.label = "square", na.label.col = "white")
  dev.off()

# })

# **********************************************************************************
# Check transcript levels of certain cytokines that look to correlate perfectly
# According to corrplot: IL10, IL10RA, IL10RB
# Cytokines weights:
IL10_IL10RA <- data$IL10_IL10RA
IL10_IL10RB <- data$IL10_IL10RB

IL4_IL2RG <- data$IL4_IL2RG
IL4_CD53 <- data$IL4_CD53
IL4_IL13RA1 <- data$IL4_IL13RA1

plot(IL4_CD53, IL4_IL13RA1)

# transcripts:
comb <- do.call(rbind, lapply(PanCancer.names, function(Cancer){
  
  # Load data
  load(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/data/PanCancer/",
              Cancer, "/DataViews_", Cancer,"_all_train_with_sTIL.RData"))
  
  tmp_data <- DataViews$transcript
  
  return(tmp_data)
  
}))

IL10RA <- comb$IL10RA
IL10RB <- comb$IL10RB
IL10 <- comb$IL10
min(IL10[1], IL10RA[1])
min(IL10[1], IL10RB[1])
IL10_IL10RA[1]
IL10_IL10RB[1]

all <- data.frame(IL10RA = IL10RA, IL10RB = IL10RB, IL10 = IL10)
pairs( ~ . , data = all, upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8) 

#IL2RG <- comb$IL2RG
CD53  <- comb$CD53 
IL4 <- comb$IL4
IL13RA1 <- comb$IL13RA1

min(IL4[1], CD53[1])
min(IL4[1], IL13RA1[1])
IL10_IL10RA[1]
IL10_IL10RB[1]

all <- data.frame(CD53 = CD53, IL4 = IL4, IL13RA1 = IL13RA1)
pairs( ~ . , data = all, upper.panel = panel.cor,lower.panel = panel.lm, cex.labels = 0.8) 


  
# ********************************************************************************** 
# Pairs
# pdf(file = "/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/pre-processing/DataTypes_correlation/quanTIseq_SpatialTIL_for_all_cancer_types_pairs.pdf",
#     height = 10, width = 10)
# pairs( ~ . , data = data, upper.panel = panel.cor,lower.panel = panel.lm, 
#        cex.labels = 0.8) 
# dev.off()



# Decision: spatialTILs correlation

for (Cancer in PanCancer.names){

  # dont_keep <- c("Number_of_TIL_Patches", "NP_mean", "number_of_clusters","Ball_Hall_Adjusted", "Banfeld_Raftery_Adjusted", "Det_Ratio_Adjusted")
  # DataViews$SpatialTIL <- DataViews$SpatialTIL[,-which(colnames(DataViews$SpatialTIL) %in% dont_keep)]
  # corrplot(cor(DataViews$quanTIseq,DataViews$SpatialTIL, method = "spearman"))
}

# 
# Collect results from models performance with 4 different set of features
# 

# **************** 
# Initialize variable to collect results
summary_view_all <- NULL
algorithms <- c("CV_linear_reg_L1&L2", "BEMKL","Lasso", "Elastic_Net")
PanCancer.names = c("BLCA", "SKCM")
# quanTIseq_spatialTILs
view_combinations <- list(views[3], views[3], views[3], views[3],views[c(3,6)])


for (Cancer in PanCancer.names){
  
  load(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/",
              "data/PanCancer/", Cancer, "/DataViews_", Cancer,"_all_train_with_sTIL.RData"))
  
  file <- dir(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/output/model/",Cancer)
              , full.names = T, recursive = F)

  summary_view <- do.call(rbind, lapply(1:length(view_combinations), function(view){
    
    input_name <- paste(names(view_combinations[[view]]), collapse ="_")
    
    if (input_name == "quanTIseq_SpatialTIL") {
      pattern <- paste0("_withConsensus_", input_name, "_just_4_TILs")
    }else{
      pattern <- paste0("_withConsensus_", input_name, "_kfold_5")
    }
  
    which_file <- grep(pattern = pattern, file, fixed = T)
    
    summary_try <- do.call(rbind, lapply(which_file, function(pos){
      
      if (input_name == "quanTIseq_SpatialTIL") {
        try <- substr(file[pos],141, nchar(file[pos])-6)
      }else{
        try <- "alone"
      }
      
      # Load file
      load(file[pos])
      
      # 2 algorithms
      summary_alg <- do.call(rbind, lapply(algorithms, function(alg){
        
        tmp_alg <- all_cv_res[[alg]]
        
        # 100 iterations
        summary_iter <- do.call(rbind, lapply(1:length(tmp_alg), function(iteration){
          
          tmp_iter <- tmp_alg[[iteration]]$performances
          
          # 5 metrics
          summary_measure <- do.call(rbind, lapply(names(tmp_iter), function(measure){
            
            tmp_measure <- tmp_iter[[measure]]
            
            # 5 tasks
            summary_task <- do.call(rbind,lapply(names(tmp_measure), function(task){
              
              tmp_task <- tmp_measure[[task]]
              
              if (alg == "BEMKL") {
                info <- data.frame(algorithm = alg,
                                   iteration = iteration,
                                   metric = measure,
                                   task = task,
                                   perf_min = tmp_task,
                                   try = try)
    
              }else {
                info <- data.frame(algorithm = alg,
                                   iteration = iteration,
                                   metric = measure,
                                   task = task,
                                   perf_min = tmp_task$min.mse,
                                   try = try)
              }
              return(info)
            }))
            return(summary_task)
          }))
          return(summary_measure)
        }))
        return(summary_iter)
      }))
      
      n_views <- length(view_combinations)
      n_algorithms <- length(algorithms)
      n_iterations <- length(all_cv_res$Lasso)
      n_measures <- length(all_cv_res$Lasso[[1]]$performances)
      n_tasks <- length(all_cv_res$Lasso[[1]]$performances$MSE)
      
      summary_alg$dataType <- rep(input_name, len = n_algorithms * n_iterations * n_measures * n_tasks)
      summary_alg$CancerType <- rep(paste0(Cancer,"(n=",nrow(DataViews$PROGENy),")"),
                                    len =  n_algorithms * n_iterations * n_measures * n_tasks)
      
      return(summary_alg)
    }))
    return(summary_try)
  }))

  summary_view_all <- rbind(summary_view_all, summary_view)
}
# manual trick for plotting #
summary_view.tmp <- subset(summary_view_all, task == "consensus")
tries <- levels(summary_view.tmp$try)[-1]
# BLCA
step <- (nrow(summary_view.tmp)/2)/8
sequence <- seq(1, (nrow(summary_view.tmp)/2)/2, by = step)
for (ii in 1:length(sequence)){
  summary_view.tmp$try[sequence[ii]:(sequence[ii]+step - 1)] <- 
    gsub("alone", tries[ii], summary_view.tmp$try[sequence[ii]:(sequence[ii]+step - 1)], fixed = T)
}

# SKCM
step <- (nrow(summary_view.tmp)/2)/8
sequence <- seq((nrow(summary_view.tmp)/2) + 1, 3*nrow(summary_view.tmp)/2/2, by = step)
for (ii in 1:length(sequence)){
  summary_view.tmp$try[sequence[ii]:(sequence[ii]+step - 1)] <- 
    gsub("alone", tries[ii], summary_view.tmp$try[sequence[ii]:(sequence[ii]+step - 1)], fixed = T)
}

colors.Cancer_type <- toupper(c("#82c8b8","#8447c0","#7fcf5c","#c95091","#cbb652","#7e88c0","#c15237","#53633a","#4e2e45","#cb9c98"))
colors.Cancer_type <- colors.Cancer_type[c(1,8)]

plot_performance.tries <- ggplot(summary_view.tmp, aes(x = algorithm, y = perf_min, fill = dataType, colour = CancerType)) +
  geom_boxplot(outlier.shape = NA) +
  scale_colour_manual(values=colors.Cancer_type) +
  theme_bw() +
  ylim(c(0,1)) +
  coord_fixed(ratio = 4) +
  facet_grid(metric ~ try) +
  theme(axis.text.x = element_text(size=7,face="bold", angle = 25, vjust = 0.75), axis.title.x = element_blank(),
        strip.text.y = element_text(size=7,face="bold"), strip.text.x = element_text(size=7,face="bold"),
        axis.title.y = element_blank(),
        strip.background = element_blank()) 

plot_performance.tries

pdf("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/pre-processing/spatial_TILs_correlation/performance_each_subset_features.pdf",
    height = 10 , width = 10)
plot_performance.tries
dev.off()


# **************** 
# Initialize variable to collect results
summary_view_features <- NULL
algorithms <- c("CV_linear_reg_L1&L2","Lasso", "Elastic_Net")

# quanTIseq_spatialTILs
view_combinations <- list(views[3], views[3], views[3], views[3],views[c(3,6)])
  
for (Cancer in PanCancer.names){
  
  load(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/",
              "data/PanCancer/", Cancer, "/DataViews_", Cancer,"_all_train_with_sTIL.RData"))
  
  file <- dir(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/output/model/",Cancer)
              , full.names = T, recursive = F)
  
  summary_view <- do.call(rbind, lapply(1:length(view_combinations), function(view){
    
    input_name <- paste(names(view_combinations[[view]]), collapse ="_")
    
    if (input_name == "quanTIseq_SpatialTIL") {
      pattern <- paste0("_withConsensus_", input_name, "_just_4_TILs")
    }else{
      pattern <- paste0("_withConsensus_", input_name, "_kfold_5")
    }
    
    which_file <- grep(pattern = pattern, file, fixed = T)
    
    summary_try <- do.call(rbind, lapply(which_file, function(pos){
      
      if (input_name == "quanTIseq_SpatialTIL") {
        try <- substr(file[pos],141, nchar(file[pos])-6)
      }else{
        try <- "alone"
      }
      # Load file
      load(file[pos])
      
      # 3 algorithms
      summary_alg <- do.call(rbind, lapply(algorithms, function(alg){
      
      tmp_alg <- all_cv_res[[alg]]
      
        # 100 iterations
        summary_iter <- do.call(rbind, lapply(1:length(tmp_alg), function(iteration){
          
          tmp_iter_model <- tmp_alg[[iteration]]$model
          
          # 5 tasks
          summary_task <- do.call(rbind, lapply(names(tmp_iter_model$Coef), function(task){
            
            tmp_iter_model_coef_task <- tmp_alg[[iteration]]$model$Coef[[task]]$min.mse
            tmp_iter_model_hyp_task <- tmp_alg[[iteration]]$model$hyperparameters[[task]]$min.mse
            tmp_iter_model_freq_task <- tmp_alg[[iteration]]$model$freq_features[[task]]
            
            feature_names <- as.character(rownames(tmp_iter_model_coef_task)[2:nrow(tmp_iter_model_coef_task)])
            
            info <- data.frame(algorithm = alg,
                               iteration = rep(iteration, times = length(feature_names)),
                               model = rep(paste0("min-mse (",tmp_iter_model_hyp_task$alpha,",",
                                                  round(as.numeric(tmp_iter_model_hyp_task$lambda),3),")"),
                                           times = length(feature_names)),
                               task = rep(task, times = length(feature_names)),
                               feature = feature_names,
                               estimate = as.numeric(tmp_iter_model_coef_task[feature_names,1]),
                               try = try)
            
            return(info)
            
          }))
          
          return(summary_task)
          
        }))
        
        n_iterations <- 100
        n_tasks <- 5
        n_features <- length(unique(summary_iter$feature))
        
        #summary_iter$algorithm <- rep(alg, len = n_algorithms * n_iterations * n_tasks * n_features)
        summary_iter$dataType <- rep(input_name, len =  n_iterations * n_tasks  * n_features)
        summary_iter$CancerType <- rep(paste0(Cancer,"(n=",nrow(DataViews$PROGENy),")"),
                                       len = n_iterations * n_tasks * n_features)
        return(summary_iter)
      
      }))
      return(summary_alg)
    }))
    return(summary_try)
  }))
  summary_view_features <- rbind(summary_view_features, summary_view)
}

# manual trick for plotting #
summary_view_features_tmp <- subset(summary_view_features, task == "consensus")
tries <- levels(summary_view_features_tmp$try)[-1]

# BLCA
end <- 3 * 100 * 11 * 4
sequence <- seq(1, end, by = end/4)
for (ii in 1:length(sequence)){
  summary_view_features_tmp$try[sequence[ii]:(sequence[ii]+end/4 - 1)] <- 
    gsub("alone", tries[ii], summary_view_features_tmp$try[sequence[ii]:(sequence[ii]+end/4 - 1)], fixed = T)
}


# SKCM
step <- 3 * 100 * 11 * 4
start <- step + 3 * 100 * 16 * 4 
sequence <- seq(start + 1, start + step, by = step/4)
for (ii in 1:length(sequence)){
  summary_view_features_tmp$try[sequence[ii]:(sequence[ii]+step/4 - 1)] <- 
    gsub("alone", tries[ii], summary_view_features_tmp$try[sequence[ii]:(sequence[ii]+step/4 - 1)], fixed = T)
}

# Per feature across cancer types -->
## Get median value for each feature
median_values_features <- aggregate(estimate ~ feature + dataType + try + algorithm + CancerType,
                                    FUN = "median", na.rm = TRUE, data = summary_view_features_tmp)

## Order features by median value:
median_values_features <- median_values_features[order(abs(median_values_features$estimate), decreasing = TRUE),]
median_values_features$feature <- as.character(median_values_features$feature)

### Restructuring data.frame
# estimates
summary_view_features_tmp.order <- summary_view_features_tmp[order(match(summary_view_features_tmp$feature, median_values_features$feature)),]
summary_view_features_tmp.order$feature <- factor(summary_view_features_tmp.order$feature, levels =  unique(median_values_features$feature))

# Colors for data type labels
labs <- as.vector(unique(summary_view_features_tmp.order$feature))
color_points <- cbind(labs, rep("royalblue2",length(labs)))
input <- strsplit("quanTIseq_SpatialTIL","_", fixed = TRUE)[[1]]
#input = strsplit("PROGENy_quanTIseq","_", fixed = TRUE)[[1]]
set.seed(1234)
if (length(input) >1) {
  color_points <- do.call(cbind,lapply(1:(length(input)-1), function(X){
    Y = input[X+1]
    colnames(DataViews[[Y]]) = sapply(strsplit(colnames(DataViews[[Y]]), ".", fixed = T), head,1)
    keep = which(color_points[,1] %in% colnames(DataViews[[Y]]))
    color_points[keep,2] =  rep(colours()[sample(1:500,1)], len = length(keep))
    return(color_points)
  }))
}
# Colors for cancer types, derived from website tool
colors.Cancer_type <- toupper(c("#82c8b8","#8447c0","#7fcf5c","#c95091","#cbb652","#7e88c0","#c15237","#53633a","#4e2e45","#cb9c98"))
colors.Cancer_type <- colors.Cancer_type[c(1,8)]
names(colors.Cancer_type) <- as.vector(unique(summary_view_features_tmp.order$CancerType))

summary_view_features_tmp.order$color <- rep(NA, times = nrow(summary_view_features_tmp.order))

color <- do.call(c,lapply(names(colors.Cancer_type) , function(X){
  tmp <- which(summary_view_features_tmp.order$CancerType == X)
  summary_view_features_tmp.order$color[tmp] <- as.vector(rep(colors.Cancer_type[X], times = length(tmp)))
}))

summary_view_features_tmp.order$color <- color
summary_view_features_tmp.order$color <- factor(summary_view_features_tmp.order$color, levels = colors.Cancer_type)


# Features coef estimates #
summary_view_features_tmp.order.BLCA <- subset(summary_view_features_tmp.order, CancerType == "BLCA(n=249)")
summary_view_features_tmp.order.SKCM <- subset(summary_view_features_tmp.order, CancerType == "SKCM(n=286)")

plot_estimates <- ggplot(summary_view_features_tmp.order.SKCM, aes(x = feature, y = estimate, fill = try, colour = try)) +
  geom_boxplot(outlier.shape = NA) + # geom_point() + 
  #scale_fill_manual(values=colors.Cancer_type) +
  #scale_colour_manual(values=colors.Cancer_type) +
  theme_bw() +
  facet_grid(algorithm ~ CancerType, scales = "free_y") + 
  #ylim(-0.5, 1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,  colour = color_points[,2])) +
  theme(axis.text.x = element_text(size=8,face="bold"), axis.text.y = element_text(size=8,face="bold"),
        axis.title.y = element_blank(), axis.title.x =  element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal", legend.background = element_rect(size=0.5),  legend.text=element_text(size=7, hjust = 0.3),
        legend.title = element_text(size = 7, vjust = 0.75)) + guides(color = FALSE) 

plot_estimates

# Features frequency #
features_frequency <- aggregate(estimate ~ dataType + try + feature + + algorithm + CancerType, FUN = function(x){length(which(x !=0 ))}, 
                                data = summary_view_features_tmp.order)
length(which(median_values_features$estimate != 0))

features_frequency.tmp <- subset(features_frequency, dataType == "quanTIseq_SpatialTIL")
features_frequency.tmp.BLCA <- subset(features_frequency.tmp, CancerType == "BLCA(n=249)")
features_frequency.tmp.SKCM <- subset(features_frequency.tmp, CancerType == "SKCM(n=286)")

features.freq <- ggplot(features_frequency.tmp.SKCM, aes(x=feature, y=algorithm, fill=estimate, group = try)) +
  geom_tile(colour="black",size=0.1) +
  scale_fill_gradient2(name = "Number of\n times", low="white", high= "red") + #midpoint = 50
  theme_bw() +
  coord_fixed(ratio = 0.5) +
  facet_grid(algorithm ~ .) + 
  theme(panel.border=element_blank(), panel.grid = element_blank(), panel.grid.major = element_blank(), 
        strip.background = element_blank(),strip.text.x = element_blank(), axis.ticks.y = element_blank()) +
  geom_text(aes(label= estimate), stat = "identity", color="black", size=3, angle = 90, fontface=1, vjust= 0.6, position=position_dodge(1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,  colour = color_points[,2])) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.title.x = element_blank(),legend.position = "top", legend.direction = "horizontal",
        legend.text=element_text(size=7), legend.title = element_text(size =7, vjust = 0.75)) +
  theme(plot.margin=unit(c(4,0,-2,0.8), "cm")) 

features.freq

# TOGETHER: BEMKL & ELASTIC NET
p1 <- plot_estimates
p2 <- features.freq

# pdf("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/pre-processing/spatial_TILs_correlation/stability_per_each_subset_features_SKCM.pdf",
#     height = 12 , width = 12)
# grid.draw(grobTree(rectGrob(gp=gpar(fill="white", lwd=0)), arrangeGrob(p2,p1, ncol=1, widths = 1)))
# dev.off()


# ------------------------------------------------------------------------------------------------------------ #
# Check multicollinearity through variance inflation factor (VIF)
# ------------------------------------------------------------------------------------------------------------ #

# ****************
# Select cancer type
PanCancer.names = c("BLCA",  "BRCA", "CRC", "LUAD", "LUSC", "PRAD", "SKCM", "STAD") # omit cesc because of lack of samples
#PanCancer.names = c("CRC", "LUAD", "SKCM")

# ****************
# views
views <- c(PROGENy = 'gaussian', #1
           Protall = 'gaussian', #2
           quanTIseq = 'gaussian', #3
           DoRothEAv1 = 'gaussian', #4
           transcript = 'gaussian', #5
           SpatialTIL = 'gaussian',  #6
           LRpairs = 'gaussian', #7
           CYTOKINEpairs = 'gaussian')  #8

# ****************
# views combinations

# # quanTIseq_spatialTILs
#view_combinations <- list(views[3],views[c(3,6)])

# PROGENy_quanTIseq
#view_combinations <- list(views[1], views[c(1,3)])

# PROGENy_DoRothEAv1
#view_combinations <- list(views[1], views[c(1,4)])

# PROGENy_LRpairs
#view_combinations <- list(views[1], views[c(1,7)])

# PROGENy_CYTOKINEpairs
#view_combinations <- list(views[1], views[c(1,8)])

# quanTIseq_LRpairs
#view_combinations <- list(views[3], views[c(3,7)])

# quanTIseq_CYTOKINEpairs
#view_combinations <- list(views[3], views[c(3,8)])

# view_combinations <- list(views[1], views[c(1,4)], views[c(1,8)], views[c(1,7)],
#                           views[3], views[c(3,6)], views[c(3,8)], views[c(3,7)],
#                           views[4], views[8], views[7])


view_combinations <- list(views[1], views[c(1,4)],
                          views[3], views[c(3,6)], 
                          views[4], views[6])
performance <- NULL

cancer_multicollinearity <- do.call(rbind, lapply(PanCancer.names, function(Cancer){
  
  cat(Cancer)
  # Load data
  load(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/data/PanCancer/",
              Cancer, "/DataViews_", Cancer,"_all_train_with_sTIL.RData"))
  
  if (Cancer == "CRC") DataViews$quanTIseq$Monocytes <- NULL
  
  load(paste0("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/data/PanCancer/",
              Cancer, "/ImmuneResponse_", Cancer,"_all_train_with_sTIL_included_IPS_IMPRES.RData"))
  
  # Consensus define as mean betweeen both responses. Used just for testing
  ImmuneResponse$consensus <- apply(scale(ImmuneResponse[,1:2]),1,mean)
  
  N <- nrow(DataViews$PROGENy)
  Ndrug <- ncol(ImmuneResponse)
  
  set.seed(123)
  prediction_indices = sample(1:N, round(N*20/100))
  learning_indices = setdiff(1:N, prediction_indices)
  
  # if (Cancer %in% c("CRC", "LUAD","SKCM")){
  # 
  #   # Remove NA values in Dataviews LRpairs and CYTOKINE pairs
  #   LR_sum <- apply(DataViews$LRpairs,2, sum)
  #   remove_NA_LR_pairs <- as.numeric(na.action(na.omit(LR_sum)))
  #   DataViews$LRpairs <- DataViews$LRpairs[,-remove_NA_LR_pairs]
  #   Cyt_sum <- apply(DataViews$CYTOKINEpairs,2, sum)
  #   remove_NA_Cyt_pairs <- as.numeric(na.action(na.omit(Cyt_sum)))
  #   DataViews$CYTOKINEpairs <- DataViews$CYTOKINEpairs[,-remove_NA_Cyt_pairs]
  # 
  # }else{
  # 
  #   view_combinations <- list(views[1], views[c(1,4)],
  #                             views[3], views[c(3,6)],
  #                             views[4],  views[6])
  # }

  
  view_multicollinearity <- do.call(rbind, lapply(1:length(view_combinations), function(view){
    
    cat(view)
    input_name <- names(view_combinations[[view]])
    view_name <- paste(names(view_combinations[[view]]), collapse = "_")
    
    drug_multicollinearity <- do.call(rbind, lapply(colnames(ImmuneResponse)[5], function(drug){
      
      # separate input data in learing and prediction
      learning.X <- as.matrix(do.call(cbind, lapply(input_name, function(x){DataViews[[x]][learning_indices,]})))
      prediction.X <- as.matrix(do.call(cbind,lapply(input_name, function(x){DataViews[[x]][prediction_indices,]})))
      
      # separate output data in learing and prediction
      learning.Y <- data.frame(Y = ImmuneResponse[learning_indices, drug])
      validation.Y <-  data.frame(Y = ImmuneResponse[prediction_indices, drug])
      
      data.learning <- as.data.frame(cbind(learning.X, learning.Y))
      data.test <- as.data.frame(cbind(prediction.X, validation.Y))
      
      # Build  a simple regression model
      model1 <- lm(Y ~., data = data.learning)
      
      # if (anyNA(model1$coefficients) == TRUE) {
      #   coef <- model1$coefficients[-1]
      #   data.learning <- data.learning[, !is.na(coef)]
      #   data.test <- data.test[,!is.na(coef)]
      #   
      # }
      # 
      # # Build  a simple regression model after na removed
      # model1 <- lm(Y ~., data = data.learning)
      
      # Make predictions
      predictions <- model1 %>% predict(data.test)
      
      # Model performance before
      performance_before <- data.frame(MSE = mean(t((as.matrix(predictions) - validation.Y)^2)),
                               R2 = R2(as.matrix(predictions), validation.Y),
                               drug = drug,
                               view = view_name,
                               CancerType = Cancer)
      
      # Detecting multicollinearity
      collinearity <- data.frame(features = names(car::vif(model1)),
                                 vif = as.numeric(car::vif(model1)),
                                 drug = drug,
                                 view = view_name,
                                 CancerType = Cancer)
      
      high_vif <- which(collinearity$vif == max(collinearity$vif))
      

      # Build a model excluding the tax variable
      data.learning <- data.learning[,-high_vif]
      data.test <- data.test[,-high_vif]
      model2 <- lm(Y ~. , data = data.learning)
      
      # Make predictions
      predictions <- model2 %>% predict(data.test)
      
      performance_after <- data.frame(MSE = mean(t((as.matrix(predictions) - validation.Y)^2)),
                                       R2 = R2(as.matrix(predictions), validation.Y),
                                       drug = drug,
                                       view = view_name,
                                       CancerType = Cancer)
      
      # 
      # summary_performance <- data.frame(error_measure = rep(c("MSE","R2"), times = 2),
      #                       error = c(performance_before$MSE, performance_before$R2, performance_after$MSE, performance_after$R2),
      #                       error_type = rep(c("before","after"), each = 2),
      #                       view = rep(view_name, times = 4),
      #                       drug = rep(drug, times = 4),
      #                       CancerType = rep(Cancer, times = 4))

      summary_collinearity <- data.frame(collinearity = collinearity$vif,
                                         features = collinearity$features,
                                         view = rep(view_name, times = length(collinearity$features)),
                                         drug = rep(drug, times = length(collinearity$features)),
                                         CancerType = rep(Cancer, times = length(collinearity$features)))

      return(summary_collinearity)
    }))
    
    return(drug_multicollinearity)
  }))

  return(view_multicollinearity)
}))
    
# 
# cancer_multicollinearity.MSE <- subset(cancer_multicollinearity, error_measure == "MSE")
# cancer_multicollinearity.MSE$error_type <- relevel(cancer_multicollinearity.MSE$error_type, ref = "before")
# 
# features.perf <- ggplot(cancer_multicollinearity.MSE, aes(x=error_type, y=CancerType, fill=error)) +
#   geom_tile(colour="black",size=0.1) +
#   scale_fill_gradient(name = "MSE", low="white", high= "red") + #midpoint = 50
#   theme_bw() +
#   #coord_fixed(ratio = 1) +
#   facet_grid(. ~ view) + 
#   theme(panel.border=element_blank(), panel.grid = element_blank(), panel.grid.major = element_blank(), 
#         strip.background = element_blank(), axis.ticks.y = element_blank()) +
#   #geom_text(aes(label= estimate), stat = "identity", color="black", size=3, angle = 90, fontface=1, vjust= 0.6, position=position_dodge(1)) +
#   theme(axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust=1)) +
#   theme(axis.title.y = element_blank(), axis.title.x = element_blank(),legend.position = "right", legend.direction = "vertical",
#         legend.text=element_text(size=9), legend.title = element_text(size =9, vjust = 0.75)) +
#   ggtitle("Performance after removing feature with highest collinearity;  lm(Y ~. , data = data.learning)")
# 
# features.perf
# 
# pdf("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/pre-processing/Multicollinearity/views_performance_after_removing_highest_vif_feature.pdf",
#     height = 12 , width = 12)
# features.perf
# dev.off()

cancer_multicollinearity.view <- subset(cancer_multicollinearity, view %in% c("SpatialTIL"))

features.coll <- ggplot(cancer_multicollinearity.view, aes(x=features, y=CancerType, fill=collinearity)) +
  geom_tile(colour="black",size=0.1) +
  scale_fill_gradient(name = "VIF", low="white", high= "red") + #midpoint = 50
  theme_bw() +
  #coord_fixed(ratio = 1) +
  #facet_grid(. ~ CancerType) +
  theme(panel.border=element_blank(), panel.grid = element_blank(), panel.grid.major = element_blank(),
        strip.background = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.text.x = element_text(size = 9, face = 2, angle = 45, vjust = 1, hjust=1)) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),legend.position = "right", legend.direction = "vertical",
        legend.text=element_text(size=9), legend.title = element_text(size =9, vjust = 0.75)) +
  ggtitle(paste0(unique(cancer_multicollinearity.view$view),": features collinearity"))

features.coll

pdf("/home/olapuent/Desktop/PhD_TUE/Github_model/desktop/figures/pre-processing/Multicollinearity/SpatialTIL_features_collinearity.pdf",
    height = 12 , width = 12)
features.coll
dev.off()








