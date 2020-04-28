predict_GLMs <- function(DataViews.train, DataViews.test, Label.test, View, K=100,
                         Trained.model, Algorithm, View.info, standardize_any=T){
  
  # ****************
  # packages
  
  # ****************
  # scripts
  source("../R/scaling_function.R")
  source("../R/GLMs/elastic_net_test.R")
  
  # ****************
  # Initialize variables
  P <- length(View.info)
  Ndrug <- length(names(Trained.model[[Algorithm]][[1]]$model$Coef))
  drugs <- names(Trained.model[[Algorithm]][[1]]$model$Coef)
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
  
  for (y in drugs){
    
    cat("Drug source: ", y, "\n")
    
    for (i in 1:K){
      state.min <- Trained.model[[Algorithm]][[i]]$model$Coef[[y]]$min.mse
      #state.1se <- Trained.model[[Algorithm]][[i]]$model$Coef[[y]]$`1se.mse`
      
      features.learning <- names(Trained.model[[Algorithm]][[i]]$mas.mea.learning.X)
      prediction.X <- do.call(cbind,lapply(names(View.info), function(x){tmp = DataViews.test[[x]]}))
      
      cat("Iteration", i,"\n")
      
      # standardize
      if (standardize_any==T){
        
        # Check same features availability
        keep_pos <- na.omit(match(colnames(prediction.X), features.learning))
        keep_names <- intersect(colnames(prediction.X), features.learning)
        prediction.X <- prediction.X[,keep_names]
        
        # Normalization should be done taking into account the train set. #
        Trained.model[[Algorithm]][[i]]$mas.mea.learning.X <- Trained.model[[Algorithm]][[i]]$mas.mea.learning.X[keep_pos]
        Trained.model[[Algorithm]][[i]]$mas.std.learning.X <- Trained.model[[Algorithm]][[i]]$mas.std.learning.X[keep_pos]
        
        prediction.X <- standarization(prediction.X, Trained.model[[Algorithm]][[i]]$mas.mea.learning.X, 
                                      Trained.model[[Algorithm]][[i]]$mas.std.learning.X)
      }
      
      # perform prediction
      predictions.min <- elastic_net_test(prediction.X, state.min)
      # predictions.1se <- elastic_net_test(prediction.X, state.1se)
      # save predictions
      predictions.all[[y]][[View]][,i] <- predictions.min[,1]

      }
    }
  
  summary_pred <- list(pred = predictions.all, lab = labels.all)
  
  return(summary_pred)
  
}