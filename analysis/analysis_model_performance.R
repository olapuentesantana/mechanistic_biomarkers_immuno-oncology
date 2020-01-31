# #########################################################################################################
# Script to plot model performances obtained from cross-validation
# #########################################################################################################

# ****************
# working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ****************
# packages
library(ggplot2)
library(reshape)

# ****************
# Select cancer type
PanCancer.names = c("BLCA",  "BRCA","CESC", "CRC", "LUAD", "LUSC", "PRAD", "SKCM", "STAD","UCEC")

# ****************
# views
views <- c(PROGENy = 'gaussian', #1
           Protall = 'gaussian', #2
           quanTIseq = 'gaussian', #3
           DoRothEAv1 = 'gaussian', #4
           transcript = 'gaussian', #5
           SpatialTIL = 'gaussian')  #6


# ****************
# algorithms
algorithms <- c("Lasso",
                "BEMKL",
                "Elastic_Net",
                "CV_linear_reg_L1&L2",
                "Ridge")

# ****************
# quanTIseq_spatialTILs
view_combinations <- list(views[3],views[c(3,6)])

# PROGENy_quanTIseq
#view_combinations <- list(views[1],views[c(1,3)])

# ------------------------------------------------------------------------------------------------------------ #
# Collect data from cross-validation results
# ------------------------------------------------------------------------------------------------------------ #

#for (Cancer in PanCancer.names){
Cancer <- "SKCM"

  load(paste0("../data/PanCancer/", Cancer, "/DataViews_", Cancer,"_all_train_with_sTIL.RData"))
  file <- dir(paste0("../output/model/",Cancer,""), pattern = paste0("_kfold_5"), full.names = T, recursive = F)
 
  summary_view <- do.call(rbind, lapply(1:length(view_combinations), function(view){
      
    # Load file
    input_name <- paste(names(view_combinations[[view]]), collapse ="_")
    which_file <- grep(pattern = paste0(input_name, "_kfold_5",".RData"), file, fixed = T)
    load(file[which_file])
    algorithms <- names(all_cv_res)
    
    # 5 algorithms
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
                                 perf_min_1SE = NA)
            }else {
              info <- data.frame(algorithm = alg,
                                 iteration = iteration,
                                 metric = measure, 
                                 task = task, 
                                 perf_min = tmp_task$min.mse,
                                 perf_min_1SE = tmp_task$`1se.mse`)
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
    
    summary_alg$dataType = rep(input_name, len = n_algorithms * n_iterations * n_measures * n_tasks)
    summary_alg$CancerType = rep(paste0(Cancer,"(n=",nrow(DataViews$PROGENy),")"), 
                                  len =  n_algorithms * n_iterations * n_measures * n_tasks)
    return(summary_alg)
  }))
  
#}
  
# ------------------------------------------------------------------------------------------------------------ #
# Statistical analysis: 
# Does the integration of spatial information improve signficantly the performance of the model?
# Does it provide complementary information?
# ------------------------------------------------------------------------------------------------------------ #
  
analysis_view <- summary_view

Combination <- paste(names(view_combinations[[2]]), collapse ="_")
Single <- paste(names(view_combinations[[1]]), collapse ="_")

# variables
cancer_types <- as.vector(unique(analysis_view$CancerType))
metrics_types <- as.vector(unique(analysis_view$metric))
tasks_types <- as.vector(unique(analysis_view$task))
algorithms_types <- as.vector(unique(analysis_view$algorithm))

results_tasks <- do.call(list,lapply(tasks_types, function(task_type){
    
  results_metrics <- do.call(list,lapply(metrics_types, function(metric_type){
      
    results_cancers <- do.call(list,lapply(cancer_types, function(cancer_type){
        
      results_algorithms <- do.call(list,lapply(algorithms_types, function(alg){
        
          # Screening for task, metric, cancer and algorithm
          analysis_tmp <- subset(analysis_view, task == task_type & CancerType == cancer_type & metric == metric_type & algorithm == alg)
          Combination_cancer_task_metric <- subset(analysis_tmp, dataType %in% Combination)
          Single_cancer_task_metric <- subset(analysis_tmp, dataType %in% Single)
          
          # initialize variables
          P.values.t <- P.values.w <- P.values.t_adj <- P.values.w_adj <- vector("list", length(Combination))
          names(P.values.t) <- names(P.values.w) <- names(P.values.t_adj) <- names(P.values.w_adj) <- Single
          
          # mse.min
          x.min <- Combination_cancer_task_metric$perf_min
          y.min <- Single_cancer_task_metric$perf_min
          
          # t.test test and wilcoxon test
          
          # alternative = "greater" is the alternative that x has a larger mean than y.
          # t-test
          test.t.min = t.test(x.min, y.min, alternative = c("greater"))
          P.values.t[[Single]]$perf_min = test.t.min$p.value
          # wilcoxon test
          test.w.min = wilcox.test(x.min, y.min, alternative = c("greater"))
          P.values.w[[Single]]$perf_min = test.w.min$p.value
          
          if (alg != "BEMKL") {
            # mse.min + 1SE
            x.1se <- Combination_cancer_task_metric$perf_min_1SE
            y.1se <- Single_cancer_task_metric$perf_min_1SE
            
            test.t.1se <- t.test(x.1se, y.1se, alternative = c("greater"))
            test.w.1se <- wilcox.test(x.1se, y.1se, alternative = c("greater"))
            P.values.t[[Single]]$perf_min_1SE <- test.t.1se$p.value
            P.values.w[[Single]]$perf_min_1SE <- test.w.1se$p.value
            
          }else{
            # mse.min + 1SE
           P.values.t[[Single]]$perf_min_1SE <- NA
           P.values.w[[Single]]$perf_min_1SE <- NA
          }
          
          P.values <- list(t_test = P.values.t, wilcox_test = P.values.w)
          return(P.values)
        }))
        names(results_algorithms) <- algorithms_types
        return(results_algorithms)
      }))
      names(results_cancers) <- cancer_types
      return(results_cancers)
    }))
    names(results_metrics) <- metrics_types
    return(results_metrics)
  }))
names(results_tasks) <- tasks_types
  
# Next: organize the results in a dataframe
  
results_heatmap <- do.call(rbind,lapply(tasks_types, function(task_type){
    
  results_heatmap <- do.call(rbind,lapply(metrics_types, function(metric_type){
      
    results_heatmap <- do.call(rbind,lapply(cancer_types, function(cancer_type){
      
      results_heatmap <- do.call(rbind,lapply(algorithms_types, function(alg){
        
        tmp <- results_tasks[[task_type]][[metric_type]][[cancer_type]][[alg]]
        P.values.t.test <- unlist(tmp$t_test)
        P.values.w.test <- unlist(tmp$wilcox_test)
        
        Ready.to.plot <- data.frame(algorithm = alg,
                                    task = task_type,
                                    metric = metric_type,
                                    CancerType = cancer_type,
                                    Combination = Combination,
                                    Single = Single,
                                    t.test.min = as.vector(P.values.t.test)[1],
                                    t.test.1se = as.vector(P.values.t.test)[2],
                                    w.test.min = as.vector(P.values.w.test)[1],
                                    w.test.1se = as.vector(P.values.w.test)[2])
        
        # Adding significance level:
        Ready.to.plot <- mutate(Ready.to.plot, sig.min=ifelse(Ready.to.plot$w.test.min < 0.05  , "p<0.05", "Not Sig"))
        Ready.to.plot <- mutate(Ready.to.plot, sig.1se=ifelse(Ready.to.plot$w.test.1se < 0.05  , "p<0.05", "Not Sig"))
        
        Ready.to.plot$w.test.min <- -log10(Ready.to.plot$w.test.min)
        Ready.to.plot$w.test.1se <- -log10(Ready.to.plot$w.test.1se)
        colnames(Ready.to.plot)[c(9,10)] <- c("minus_log10_w.test.min", "minus_log10_w.test.1se")
        saturation <- -log10(0.001)
        Ready.to.plot[which(Ready.to.plot$`-log10_w.test.min` >= saturation),"minus_log10_w.test.min"] <- saturation
        Ready.to.plot[which(Ready.to.plot$`-log10_w.test.1se` >= saturation),"minus_log10_w.test.1se"] <- saturation
        
        return(Ready.to.plot)
      }))
      return(results_heatmap)
    }))
    return(results_heatmap)
  }))
  return(results_heatmap)
}))
  
# ------------------------------------------------------------------------------------------------------------ #
# GGPLOT: 
# Boxplots with performance across 100 iterations
# ------------------------------------------------------------------------------------------------------------ #

#summary_view.lasso <- subset(summary_view, algorithm == "Lasso")
summary_view.consensus <- subset(summary_view, task == "consensus")
plot_performance <- ggplot(summary_view.consensus, aes(x = algorithm, y = perf_min, fill = dataType)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    ylim(c(0.65,1)) +
    #coord_fixed(ratio = 3) +
    facet_grid(metric ~ CancerType , scales = "fixed", space = "fixed")
    #scale_x_discrete(labels = function(x) format(x, width = 5)) +
    #theme(axis.text.x = element_blank(), strip.background = element_blank(), strip.text.y = element_text(size=7,face="bold"),
          #axis.ticks.x = element_blank()) +
    #strip.text.y = element_blank()) +
    #theme(axis.text.x =  element_blank(), axis.text.y = element_text(size = 7, face = "bold"),
          #axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none") +
    #theme(plot.margin=unit(c(-0.5,0.15,-1,0.2), "cm")) +
    #ggtitle(label= "A") + theme(plot.title = element_text(size=14,face = "bold"))

plot_performance
  
# ------------------------------------------------------------------------------------------------------------ #
# Heatmap: 
# Heatmap showing results from significant analysis
# ------------------------------------------------------------------------------------------------------------ #
results_heatmap.lasso.consensus <- subset(results_heatmap, algorithm == "BEMKL" & task == "consensus")

plot_heatmap <- ggplot(results_heatmap.lasso.consensus, aes(x = task, y = metric, fill = minus_log10_w.test.min)) +
  #add border white colour of line thickness 0.25
  geom_tile(colour="white",size=0.1) +
  scale_fill_gradient2(low="darkgreen", midpoint = -log10(0.05), high="darkred") +
  theme_bw() +
  coord_fixed(ratio = 1) +
  theme(panel.border=element_blank(), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.spacing =  unit(0, "lines"),
        strip.background = element_blank(),strip.text.x = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank()) +
  #geom_text_repel(data=Ready.to.plot, aes(label=ifelse(sig!="Not Sig","*",""))) +
  facet_grid(. ~ CancerType, scales = "fixed", space = "fixed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.7)) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(labels = function(x) format(x, width = 5)) +
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        legend.position = "none") 
  #theme(plot.margin=unit(c(-7,0.56,2,0.8), "cm"))

plot_heatmap

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
