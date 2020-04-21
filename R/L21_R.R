##########################################################################################################################
# Script to perform regularized Multi-task Learning using the package RMTL - L21 regularization 
##########################################################################################################################

L21_R <- function(drug_source, views_source, view_combination, learning_indices, prediction_indices,
                  standardize_any=F, standardize_response=F,
                  parameters, iteration) {

  
  # ****************
  # scripts
  source("./R/scaling_function.R")
  source("./R/unscaling_function.R")
  source("./R/civalue.R")
  source("./R/L21_regularization_EN/L21_test.R")
  source("./R/L21_regularization_EN/L21_train.R")
  
  # ****************
  # General variables:
  names_view <- names(view_combination) # Need data type 
  P <- length(view_combination) # number of views
  Ndrug <- ncol(drug_source) # number of tasks (output variables)
  
  # Output variables: 
  performances <- list()
  model <- list()
  training_set <- list()
  mas.mea.learning.X <- mas.std.learning.X <- list()
  mas.mea.learning.Y <- mas.std.learning.Y <- list()
  
  # ****************
  # Initialize variables
  names_view <- names(view_combination) # Need data type 
  #mas.std.learning.X <- mas.mea.learning.X <- vector("list", length = random)
  #mas.std.learning.Y <- mas.mea.learning.Y <- vector("list", length = random)
  
  # separate input data in learing and prediction
  learning.X <- lapply(names_view, function(x){views_source[[x]][learning_indices,]})
  prediction.X <- lapply(names_view, function(x){views_source[[x]][prediction_indices,]})
  
  # separate output data in learing and prediction
  learning.Y <- drug_source[learning_indices,]
  validation.Y <- drug_source[prediction_indices,]
  
  if (standardize_any==T){
    # view normalization
    for (m in 1:P){
      # nan_indices_learning = is.na(learning.X[[m]])
      # nan_indices_prediction = is.na(prediction.X[[m]])
    
      mas.mea.learning.X[[m]] = colMeans(learning.X[[m]], na.rm = T)
      mas.std.learning.X[[m]] = colSds(as.matrix(learning.X[[m]]), na.rm = T)
      
      ## same but slower
      # mas.mea2 = apply(learning.X[[m]], 2, mean, na.rm = T)
      # mas.std2 = apply(learning.X[[m]], 2, sd, na.rm = T)
      
      mas.std.learning.X[[m]][mas.std.learning.X[[m]]==0] = 1
      
      learning.X[[m]]= standarization(learning.X[[m]], mas.mea.learning.X[[m]], mas.std.learning.X[[m]])
      prediction.X[[m]] = standarization(prediction.X[[m]], mas.mea.learning.X[[m]],mas.std.learning.X[[m]])
      
      # learning.X[[m]][nan_indices_learning==T] <- 0
      # prediction.X[[m]][nan_indices_prediction==T] <- 0
    }
    cat("input data normalization done","\n")
    
    # drug response standardization
    # mas.mea = apply(learning.Y, 2, mean, na.rm = T)
    
    mas.mea.learning.Y = colMeans(as.matrix(learning.Y), na.rm = T)
    mas.std.learning.Y =  rep(1, ncol(learning.Y)) #
    
    learning.Y <- sweep(learning.Y, 2,   mas.mea.learning.Y, FUN = "-")
    validation.Y <- sweep(validation.Y, 2, mas.mea.learning.Y, FUN = "-")
    
    if (standardize_response == T){
      
      mas.std.learning.Y = colSds(as.matrix(learning.Y), na.rm = T)
      mas.std.learning.Y[mas.std.learning.Y == 0] <- 1
      
      learning.Y <- sweep(learning.Y, 2, mas.std.learning.Y, FUN = "/")
      validation.Y <- sweep(validation.Y, 2, mas.std.learning.Y, FUN = "/")
    }
    
    cat("output data normalization done","\n")
    
  }else{
    # if standardize_any=F we assume data have been already normalized,
    # NA are set to 0 (i.e. average value)
    # for (m in 1:P){
    #   nan_indices_learning = is.na(learning.X[[m]])
    #   learning.X[[m]][nan_indices_learning==T] <- 0
    # }
    # 
    cat("careful: data are assumed to be already normalized (if that is not the case set standardize_any=T)","\n")
  }
  
  # Hyperparameters estimation
  state <- L21_train(drug_source = learning.Y,
                             views_source = learning.X,
                             view_combination = names_view,
                             parameters = parameters, iteration = iteration)
  
  coef_values <- state$cv.MTL.features
  hyperparameter_values <- state$cv.MTL.hyperparameters
  cat("training performed","\n")
  
  # perform prediction
  prediction_cv <- L21_test(prediction.X, coef_values)
  cat("prediction computed","\n")
  
  MSE <- apply((validation.Y - prediction_cv)^2, 2, mean)
  SpCorr <- diag(cor(validation.Y, prediction_cv, method = "spearman"))
  PeCorr <- diag(cor(validation.Y, prediction_cv, method = "pearson"))
  CI <- sapply(1:Ndrug, function(x) civalue(validation.Y[,x], prediction_cv[,x]))
  
  names(MSE) <- names(SpCorr) <- names(PeCorr) <- names(CI) <- colnames(drug_source)
  
  performances <- list(MSE=MSE, SpCorr=SpCorr, PeCorr=PeCorr, CI=CI)
  model <- state
  
  return(list(performances = performances,
              model = model,
              training_set = NULL,
              mas.std.learning.Y = mas.std.learning.Y,
              mas.mea.learning.Y = mas.mea.learning.Y,
              mas.std.learning.X = mas.std.learning.X,
              mas.mea.learning.X = mas.mea.learning.X))
  
}



