predict_GLMs <- function(DataViews, DataViews.test.pre, ImmuneResponse.test.pre, j, K,
                         all_cv_res, alg, view_combination, standardize_any=T, X){
  
  # ****************
  # packages
  
  # ****************
  # scripts
  source("./R/scaling_function.R")
  source("./R/GLMs/elastic_net_test.R")
  
  P = length(view_combination)
  Ndrug = length(names(all_cv_res[[alg]][[1]]$model$Coef))

  label_IS[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  predictions_IS.min[[j]] <- predictions_IS.1se[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  
  label_CYT[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  predictions_CYT.min[[j]] <- predictions_CYT.1se[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  
  label_IPS[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  predictions_IPS.min[[j]] <- predictions_IPS.1se[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  
  label_impres[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  predictions_impres.min[[j]]<- predictions_impres.1se[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  
  label_consensus[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  predictions_consensus.min[[j]] <- predictions_consensus.1se[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  
  label_common[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  predictions_common.min[[j]] <- predictions_common.1se[[j]] <- matrix(NA, nrow=nrow(ImmuneResponse.test.pre[[X]]), ncol=K)
  
  for (y in 1:Ndrug){
    o = names(all_cv_res[[alg]][[1]]$model$Coef)[y]
    cat(y,".drug source: ", o, "\n")
    
    for (i in 1:K){
      state.min <- all_cv_res[[alg]][[i]]$model$Coef[[y]]$min.mse
      state.1se <- all_cv_res[[alg]][[i]]$model$Coef[[y]]$`1se.mse`
      
      features.learning <- names(all_cv_res[[alg]][[i]]$mas.mea.learning.X)
      prediction.X <- do.call(cbind,lapply(names(view_combination), function(x){
        tmp = DataViews.test.pre[[X]][[x]]}))
      
      cat("Iteration", i,"\n")
      
      # standardize
      if (standardize_any==T){
        
        # Check same features availability
        keep_pos = na.omit(match(colnames(prediction.X), features.learning))
        keep_names = intersect(colnames(prediction.X), features.learning)
        prediction.X <- prediction.X[,keep_names]
        
        all_cv_res[[alg]][[i]]$mas.mea.learning.X <- all_cv_res[[alg]][[i]]$mas.mea.learning.X[keep_pos]
        all_cv_res[[alg]][[i]]$mas.std.learning.X <- all_cv_res[[alg]][[i]]$mas.std.learning.X[keep_pos]
        prediction.X = standarization(prediction.X,  all_cv_res[[alg]][[i]]$mas.mea.learning.X, 
                                      all_cv_res[[alg]][[i]]$mas.std.learning.X)
      }
      
      # perform prediction
      predictions.min <- elastic_net_test(prediction.X, state.min)
      predictions.1se <- elastic_net_test(prediction.X, state.1se)
      
      if (y == 1){
        label_CYT[[j]][,i] <- as.character(ImmuneResponse.test.pre[[X]][,1])
        predictions_CYT.min[[j]][,i] <- predictions.min[,1]
        predictions_CYT.1se[[j]][,i] <- predictions.1se[,1]
        
      }else if (y == 2){
        label_IS[[j]][,i] <- as.character(ImmuneResponse.test.pre[[X]][,1]) # Keep them as character
        predictions_IS.min[[j]][,i] <- predictions.min[,1]
        predictions_IS.1se[[j]][,i] <- predictions.1se[,1]
        
      }else if (y == 3){
        label_IPS[[j]][,i] <- as.character(ImmuneResponse.test.pre[[X]][,1])
        predictions_IPS.min[[j]][,i] <- predictions.min[,1]
        predictions_IPS.1se[[j]][,i] <- predictions.1se[,1]
        
      }else if (y == 4){
        label_impres[[j]][,i] <- as.character(ImmuneResponse.test.pre[[X]][,1])
        predictions_impres.min[[j]][,i] <- predictions.min[,1]
        predictions_impres.1se[[j]][,i] <- predictions.1se[,1]
        
      }else if (y == 5){
        label_consensus[[j]][,i] <- as.character(ImmuneResponse.test.pre[[X]][,1])
        predictions_consensus.min[[j]][,i] <- predictions.min[,1]
        predictions_consensus.1se[[j]][,i] <- predictions.1se[,1]
        
        label_common[[j]][,i] <- as.character(ImmuneResponse.test.pre[[X]][,1])
        predictions_common.min[[j]][,i] <- apply(cbind(predictions_CYT.min[[j]][,i], predictions_IS.min[[j]][,i]),1, mean)
        predictions_common.1se[[j]][,i] <- apply(cbind(predictions_CYT.1se[[j]][,i], predictions_IS.1se[[j]][,i]),1, mean)                                           
      }
    }
  }
  
  summary_kfold <- list(label_IS = label_IS, predictions_IPS.min = predictions_IPS.min, predictions_IS.1se = predictions_IS.1se,
                        label_CYT = label_CYT , predictions_CYT.min = predictions_CYT.min, predictions_CYT.1se = predictions_CYT.1se,
                        label_IPS = label_IPS, predictions_IPS.min = predictions_IPS.min, predictions_IPS.1se = predictions_IPS.1se,
                        label_impres = label_impres, predictions_impres.min = predictions_impres.min, predictions_impres.1se = predictions_impres.1se,
                        label_consensus = label_consensus, predictions_consensus.min = predictions_consensus.min, predictions_consensus.1se = predictions_consensus.1se,
                        label_common = label_common, predictions_common.min = predictions_common.min, predictions_common.1se = predictions_common.1se)
  
  return(summary_kfold)
  
}