predict_with_rmtlr <- function(DataViews.train, DataViews.test, Label.test, View, K=100,
                        Trained.model, Algorithm, View.info, standardize_any=T){
  
  # scripts
  source("./tcga_training/scaling_function.R")
  source("./tcga_training/RMTLR/rmtlr_test.R")
  
  # Initialize variables
  P <- length(View.info)
  Ndrug <- length(colnames(Trained.model[[Algorithm]][[1]]$model$cv.glmnet.features$min.mse))
  drugs <- colnames(Trained.model[[Algorithm]][[1]]$model$cv.glmnet.features$min.mse)
  labels <- list()
  predictions <- list()
  predictions.all <- list()
  labels.all <- list()
  
  # Per task, per view
  labels[[View]] <- matrix(Label.test$label, nrow=nrow(Label.test), ncol=K, 
                           dimnames = list(Label.test$Sample, seq(1,100,1)))
  
  predictions[[View]] <- matrix(NA, nrow=nrow(Label.test), ncol=K, 
                                dimnames = list(Label.test$Sample, seq(1,100,1)))
  
  labels.all <- do.call(c, lapply(drugs, function(X){labels.all[[X]] <- labels; return(labels.all)}))
  
  predictions.all <- do.call(c, lapply(drugs, function(X){
    predictions.all[[X]]$min.mse <- predictions
    predictions.all[[X]]$`1se.mse` <- predictions
    return(predictions.all)}))

  for (i in 1:K){
    state <- Trained.model[[Algorithm]][[i]]$model$cv.glmnet.features
    
    features.learning <- lapply(1:length(View.info),function(x){tmp = names(Trained.model[[Algorithm]][[i]]$mas.mea.learning.X[[x]])})
    prediction.X <- lapply(names(View.info), function(x){tmp = DataViews.test[[x]]})
    
    cat("Iteration", i,"\n")
    
    # standardize
    if (standardize_any==T){
      for (m in 1:P){
        
        # Check same features availability
        keep_pos <- na.omit(match(colnames(prediction.X[[m]]), features.learning[[m]]))
        keep_names <- intersect(colnames(prediction.X[[m]]), features.learning[[m]])
        prediction.X[[m]] <- prediction.X[[m]][,keep_names]
        
        # Normalization should be done taking into account the train set. #
        Trained.model[[Algorithm]][[i]]$mas.mea.learning.X[[m]] <- Trained.model[[Algorithm]][[i]]$mas.mea.learning.X[[m]][keep_pos]
        Trained.model[[Algorithm]][[i]]$mas.std.learning.X[[m]] <- Trained.model[[Algorithm]][[i]]$mas.std.learning.X[[m]][keep_pos]
        
        prediction.X[[m]] <- standarization(prediction.X[[m]], Trained.model[[Algorithm]][[i]]$mas.mea.learning.X[[m]], 
                                            Trained.model[[Algorithm]][[i]]$mas.std.learning.X[[m]])
        
      }
    }
    
    # perform prediction
    predictions.model <- lapply(state, function(X) rmtlr_test(prediction.X, X))
    
    # save predictions
    predictions.all <- lapply(drugs, function(X){
      predictions.all[[X]][["min.mse"]][[View]][,i] <- predictions.model[["min.mse"]][,X] 
      predictions.all[[X]][["1se.mse"]][[View]][,i] <- predictions.model[["1se.mse"]][,X]
      return(predictions.all[[X]])})
    names(predictions.all) <- drugs
  }
  
  summary_pred <- list(pred = predictions.all, lab = labels.all)
  
  return(summary_pred)
  
}
