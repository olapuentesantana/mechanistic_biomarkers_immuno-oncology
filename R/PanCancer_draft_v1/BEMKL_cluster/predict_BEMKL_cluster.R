predict_BEMKL <- function(DataViews.train, DataViews.test, Label.test, View, K=100,
                          Trained.model, View.info, standardize_any=T){
  
  # ****************
  # packages
  library(reshape2,lib.loc= Sys.getenv("R_LIBS"))
  library(pdist,lib.loc= Sys.getenv("R_LIBS"))  
  
  # ****************
  # scripts
  source("./R/PanCancer_draft_v1/scaling_function.R")
  source("./R/PanCancer_draft_v1/BEMKL_cluster/bemkl_supervised_multioutput_regression_variational_test.R")
  
  P <- length(View.info)
  drugs <- names(Trained.model[[1]]$performances$MSE)
  Ndrug <- length(drugs)
  labels <- labels.all <- list()
  predictions <- predictions.all <- list()
  
  # Per task, per view
  labels[[View]] <- matrix(Label.test$label, nrow=nrow(Label.test), ncol=K, 
                           dimnames = list(Label.test$Sample, seq(1,100,1)))
  
  predictions[[View]] <- matrix(NA, nrow=nrow(Label.test), ncol=K, 
                                dimnames = list(Label.test$Sample, seq(1,100,1)))
  
  labels.all <- labels
  predictions.all <- do.call(c, lapply(drugs, function(X){predictions.all[[X]] <- predictions; return(predictions.all)}))
  
  
  for (i in 1:K){
    state <- Trained.model[[i]]$model
    learning.X <- Trained.model[[i]]$training_set
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
          Trained.model[[i]]$mas.mea.learning.X[[m]] <- Trained.model[[i]]$mas.mea.learning.X[[m]][keep_pos]
          Trained.model[[i]]$mas.std.learning.X[[m]] <- Trained.model[[i]]$mas.std.learning.X[[m]][keep_pos]
          
          prediction.X[[m]] <- standarization(prediction.X[[m]], Trained.model[[i]]$mas.mea.learning.X[[m]], 
                                              Trained.model[[i]]$mas.std.learning.X[[m]])
          
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
    colnames(predictions) <- drugs
    rownames(predictions) <- rownames(prediction.X[[1]])
    
    # save predictions
    for (X in drugs){
      predictions.all[[X]][[View]][,i] <- predictions[,X]
    }
  }
  
  summary_pred <- list(pred = predictions.all, lab = labels.all)
  return(summary_pred)
}

