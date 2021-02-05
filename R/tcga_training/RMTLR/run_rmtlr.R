##########################################################################################################################
# Script to run multiresponse (multi-task learning) from glmnet
##########################################################################################################################

run_rmtlr <- function(drug_source, views_source, view_combination, learning_indices, prediction_indices,
                        standardize_any=F, standardize_response=F,
                        parameters, iteration) {
  
  
  # ****************
  # scripts
  source("./tcga_training/scaling_function.R")
  source("./tcga_training/civalue.R")
  source("./tcga_training/RMTLR/rmtlr_test.R")
  source("./tcga_training/RMTLR/rmtlr_train.R")
  
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
  
  # separate input data in learing and prediction
  learning.X <- lapply(names_view, function(x){views_source[[x]][learning_indices,]})
  prediction.X <- lapply(names_view, function(x){views_source[[x]][prediction_indices,]})
  
  # separate output data in learing and prediction
  learning.Y <- drug_source[learning_indices,]
  validation.Y <- drug_source[prediction_indices,]

  # separate output data in learing and prediction
  learning.Y <- drug_source[learning_indices,]
  validation.Y <- drug_source[prediction_indices,]
  
  if (standardize_any==T){
    # view normalization
    for (m in 1:P){
      
      mas.mea.learning.X[[m]] = colMeans(learning.X[[m]], na.rm = T)
      mas.std.learning.X[[m]] = colSds(as.matrix(learning.X[[m]]), na.rm = T)
      
      mas.std.learning.X[[m]][mas.std.learning.X[[m]]==0] = 1
      
      learning.X[[m]]= standarization(learning.X[[m]], mas.mea.learning.X[[m]], mas.std.learning.X[[m]])
      prediction.X[[m]] = standarization(prediction.X[[m]], mas.mea.learning.X[[m]],mas.std.learning.X[[m]])

    }
    cat("input data normalization done","\n")
    
    # drug response standardization
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
  state <- rmtlr_train(drug_source = learning.Y,
                       views_source = learning.X,
                       family = "mgaussian",
                       view_combination = names_view,
                       measure_type = "mse",
                       parameters = parameters,
                       parallelize = T, iteration = iteration)
  
  coef_values <- state$cv.glmnet.features
  hyperparameter_values <- state$cv.glmnet.hyperparameters
  cat("training performed","\n")
  
  # perform prediction
  prediction_cv <- lapply(coef_values, function(X){mgaussian_test(prediction.X, X)})
  cat("prediction computed","\n")
  
  # Four metrics
  MSE <- lapply(prediction_cv, function(X){apply((validation.Y - X)^2, 2, mean)})
  SpCorr <- lapply(prediction_cv, function(X){diag(cor(validation.Y, X, method = "spearman"))})
  PeCorr <- lapply(prediction_cv, function(X){diag(cor(validation.Y, X, method = "pearson"))})
  CI <- lapply(prediction_cv, function(X){civalue(validation.Y, X)})

  # Return performance and model info
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

