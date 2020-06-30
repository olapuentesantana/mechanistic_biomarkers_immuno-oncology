predict_MT_GLM <- function(DataViews.train, DataViews.test, Label.test, View, K=100,
                         Trained.model, View.info, standardize_any=T){

  # ****************
  # scripts
  source("./R/PanCancer_draft_v1/scaling_function.R")
  source("./R/PanCancer_draft_v1/MT_GLM_cluster/MT_GLM_test.R")
  
  # ****************
  # Initialize variables
  P <- length(View.info)
  models <- names(Trained.model[[1]]$performances$MSE)
  drugs <- names(Trained.model[[1]]$performances$MSE$min.mse)
  Ndrug <- length(drugs)
  labels <- labels.all <- list()
  predictions <- predictions.all.tasks <- predictions.all.models <- list()

  # Per task, per view
  labels[[View]] <- matrix(Label.test$label, nrow=nrow(Label.test), ncol=K, 
                           dimnames = list(Label.test$Sample, seq(1,100,1)))
  
  predictions[[View]] <- matrix(NA, nrow=nrow(Label.test), ncol=K, 
                                dimnames = list(Label.test$Sample, seq(1,100,1)))
  
  labels.all <- labels
  predictions.all.models <- do.call(c, lapply(models, function(X){predictions.all.models[[X]] <- predictions; return(predictions.all.models)}))
  predictions.all.tasks <- do.call(c, lapply(drugs, function(X){predictions.all.tasks[[X]] <- predictions.all.models; return(predictions.all.tasks)}))

  for (i in 1:K){
    
    state <- Trained.model[[i]]$model$cv.glmnet.features
    features.learning <- lapply(1:length(View.info), function(x){names(Trained.model[[i]]$mas.mea.learning.X[[x]])})
    prediction.X <- lapply(names(View.info), function(x){DataViews.test[[x]]})
    
    cat("Iteration", i,"\n")
    
    # standardize
    if (standardize_any==T){
      for (m in 1:P){
        
        # Check same features availability
        keep_pos <- na.omit(match(colnames(prediction.X[[m]]), features.learning[[m]]))
        keep_names <- intersect(colnames(prediction.X[[m]]), features.learning[[m]])
        prediction.X[[m]] <- prediction.X[[m]][,keep_names]
        
        # Normalization should be done taking into account the train set. #
        Trained.model[[i]]$mas.mea.learning.X[[m]] <- Trained.model[[i]]$mas.mea.learning.X[[m]][keep_pos]
        Trained.model[[i]]$mas.std.learning.X[[m]] <- Trained.model[[i]]$mas.std.learning.X[[m]][keep_pos]
        
        prediction.X[[m]] <- standarization(prediction.X[[m]], Trained.model[[i]]$mas.mea.learning.X[[m]], 
                                            Trained.model[[i]]$mas.std.learning.X[[m]])
      }
    }
    
    # perform prediction
    prediction_cv <- lapply(state, function(X){mgaussian_test(prediction.X, X)})
    
    # save predictions
    for (X in drugs){
      predictions.all.tasks[[X]][["1se.mse"]][[View]][,i] <- prediction_cv$`1se.mse`[,X]
      predictions.all.tasks[[X]][["min.mse"]][[View]][,i] <- prediction_cv$min.mse[,X]
    }
  }
  summary_pred <- list(pred = predictions.all.tasks, lab = labels.all)
  return(summary_pred)
}