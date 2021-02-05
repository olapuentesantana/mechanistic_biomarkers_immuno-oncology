predict_BEMKL <- function(DataViews.train, DataViews.test, Label.test, View, K=100,
                          Trained.model, Algorithm, View.info, standardize_any=T){
  
  # packages
  library(reshape2)
  library(pdist)
  
  # scripts
  source("./tcga_training/scaling_function.R")
  source("./tcga_training/BEMKL/bemkl_supervised_multioutput_regression_variational_test.R")
  
  P <- length(View.info)

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
  predictions.all <- do.call(c, lapply(drugs, function(X){predictions.all[[X]] <- predictions; return(predictions.all)}))
  
  for (i in 1:K){
    state <- Trained.model[[Algorithm]][[i]]$model
    learning.X <- Trained.model[[Algorithm]][[i]]$training_set
    prediction.X <- lapply(names(View.info), function(x){DataViews.test[[x]]})

    cat("Iteration", i,"\n")
    # standardize
    if (standardize_any==T){
      for (m in 1:P){
        
        if (View.info[m] != "jaccard"){
          # Check same features availability
          keep_pos <- na.omit(match(colnames(prediction.X[[m]]), colnames(learning.X[[m]])))
          keep_names <- intersect(colnames(prediction.X[[m]]), colnames(learning.X[[m]]))
          prediction.X[[m]] <- prediction.X[[m]][,keep_names]
          learning.X[[m]] <- learning.X[[m]][,keep_names]

          # Normalization should be done taking into account the train set. #
          Trained.model[[Algorithm]][[i]]$mas.mea.learning.X[[m]] <- Trained.model[[Algorithm]][[i]]$mas.mea.learning.X[[m]][keep_pos]
          Trained.model[[Algorithm]][[i]]$mas.std.learning.X[[m]] <- Trained.model[[Algorithm]][[i]]$mas.std.learning.X[[m]][keep_pos]
          
          prediction.X[[m]] <- standarization(prediction.X[[m]], Trained.model[[Algorithm]][[i]]$mas.mea.learning.X[[m]], 
                                              Trained.model[[Algorithm]][[i]]$mas.std.learning.X[[m]])
          
        }
      }
    }
    
    # compute prediction kernel
    Nlearning <- nrow(learning.X[[1]])
    Nprediction <- nrow(prediction.X[[1]])
    Kx_prediction = array(rep(0, Nlearning*Nprediction*P), c(Nlearning, Nprediction, P))
    
    for (m in 1:P){
      Kx_prediction[, , m] <- exp(-(as.matrix(pdist(learning.X[[m]], prediction.X[[m]])))^2/ncol(learning.X[[m]])/2)
    }
    
    Ktest <- Kx_prediction #should be an Ntra x Ntest x P matrix containing similarity values between training and test samples
    
    #perform prediction
    prediction <- bemkl_supervised_multioutput_regression_variational_test(Ktest, state)
    predictions <- t(prediction$Y$mu)
    
    predictions_IS.min[[j]][,i] <- predictions[,2]
   
    # save predictions
    predictions.all[[View]][,i] <- predictions
  
  }
  
  summary_pred <- list(pred = predictions.all, lab = labels.all)
  return(summary_pred)
}

