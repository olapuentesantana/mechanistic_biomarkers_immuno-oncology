predict_BEMKL <- function(DataViews, DataViews.test.pre, ImmuneResponse.test.pre, j, K,
                          all_cv_res, alg, view_combination, standardize_any=T, X){
  
  library(pdist)
  
  # ****************
  # packages
  library(reshape2)
  library(pdist)
  
  # ****************
  # scripts
  source("./R/scaling_function.R")
  source("./R/BEMKL/bemkl_supervised_multioutput_regression_variational_test.R")
  
  P = length(view_combination)
  
  label_IS[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  predictions_IS.min[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  
  label_CYT[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  predictions_CYT.min[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  
  label_IPS[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  predictions_IPS.min[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  
  label_impres[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  predictions_impres.min[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  
  label_consensus[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  predictions_consensus.min[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  
  label_common[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  predictions_common.min[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  
  for (i in 1:K){
    state <- all_cv_res[[alg]][[i]]$model
    learning.X <- all_cv_res[[alg]][[i]]$training_set
    prediction.X <- lapply(names(view_combination), function(x){DataViews.test.pre[[X]][[x]]})

    cat("Iteration", i,"\n")
    # standardize
    if (standardize_any==T){
      for (m in 1:P){
        
        if (view_combination[m] != "jaccard"){
          # Check same features availability
          keep_pos = na.omit(match(colnames(prediction.X[[m]]), colnames(learning.X[[m]])))
          keep_names = intersect(colnames(prediction.X[[m]]), colnames(learning.X[[m]]))
          prediction.X[[m]] <- prediction.X[[m]][,keep_names]
          learning.X[[m]] <- learning.X[[m]][,keep_names]

          all_cv_res[[alg]][[i]]$mas.mea.learning.X[[m]] <- all_cv_res[[alg]][[i]]$mas.mea.learning.X[[m]][keep_pos]
          all_cv_res[[alg]][[i]]$mas.std.learning.X[[m]] <- all_cv_res[[alg]][[i]]$mas.std.learning.X[[m]][keep_pos]
          
          prediction.X[[m]] <- standarization(prediction.X[[m]], all_cv_res[[alg]][[i]]$mas.mea.learning.X[[m]], 
                                              all_cv_res[[alg]][[i]]$mas.std.learning.X[[m]])
          
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
    
    label_IS[[j]][,i] <- as.character(ImmuneResponse.test.pre[[X]][,1]) # Keep them as character
    predictions_IS.min[[j]][,i] <- predictions[,2]
    
    label_CYT[[j]][,i] <- as.character(ImmuneResponse.test.pre[[X]][,1])
    predictions_CYT.min[[j]][,i] <- predictions[,1]
    
    label_IPS[[j]][,i] <- as.character(ImmuneResponse.test.pre[[X]][,1])
    predictions_IPS.min[[j]][,i] <- predictions[,3]
    
    label_impres[[j]][,i] <- as.character(ImmuneResponse.test.pre[[X]][,1])
    predictions_impres.min[[j]][,i] <- predictions[,4]
    
    label_consensus[[j]][,i] <- as.character(ImmuneResponse.test.pre[[X]][,1])
    predictions_consensus.min[[j]][,i] <- predictions[,5]
    
    label_common[[j]][,i] <- as.character(ImmuneResponse.test.pre[[X]][,1])
    predictions_common.min[[j]][,i] <- apply(predictions[,1:2], 1, mean)
  
  }
  return(list(label_IS = label_IS, predictions_IS.min = predictions_IS.min,
              label_CYT = label_CYT , predictions_CYT.min = predictions_CYT.min,
              label_IPS = label_IPS, predictions_IPS.min = predictions_IPS.min, 
              label_impres = label_impres, predictions_impres.min = predictions_impres.min,
              label_consensus = label_consensus, predictions_consensus.min = predictions_consensus.min, 
              label_common = label_common, predictions_common.min = predictions_common.min))
  
}

