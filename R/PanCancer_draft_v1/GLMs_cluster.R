##########################################################################################################################
# Script to perform generalized linear models :
## Linear regression with Lasso regularization (alpha = 1),
## Linear regression with Ridge regularization (alpha = 0), 
## Linear regression with Elastic Net regularization (alpha = 0.5),
## Linear regression with CV-Elastic Net regularization (alpha-lambda)
##########################################################################################################################
GLMs <- function(drug_source, views_source, view_combination, learning_indices, prediction_indices,
                 standardize_any=F, standardize_response=F,
                 parameters, iteration) {
  
  # ****************
  # packages
  library(matrixStats, lib.loc= Sys.getenv("R_LIBS"))
  
  # ****************
  # scripts
  source("./R/PanCancer_draft_v1/scaling_function.R")
  source("./R/PanCancer_draft_v1/civalue.R")
  source("./R/PanCancer_draft_v1/GLMs_cluster/elastic_net_test.R")
  source("./R/PanCancer_draft_v1/GLMs_cluster/elastic_net_train.R")
  source("./R/PanCancer_draft_v1/GLMs_cluster/log_elastic_net_test.R")
  
  # ****************
  # General variables:
  names_view <- names(view_combination) # Need data type 
  Ndrug <- ncol(drug_source) # number of tasks (output variables)
  P <- length(view_combination) # number of views
  
  # Output variables: 
  performances <- list()
  model <- list()
  training_set <- list()
  mas.std.learning.X <- list()
  mas.mea.learning.X <- list()
  
  # ****************
  # Initialize variables
  names_view <- names(view_combination) # Need data type 
  MSE <- SpCorr <- PeCorr <- CI <- MR <- coef_values <- x.train <- hyperparameter_values <- freq_values <- vector("list", length = Ndrug)
  names(MSE) <- names(SpCorr) <- names(PeCorr) <- names(CI) <- names(MR) <- names(coef_values) <- names(x.train) <- names(hyperparameter_values) <- names(freq_values) <- names(drug_source)
  
  # When combining data sets:
  name = NULL
  views_source.comb = NULL
  if (P > 1){
    for (i in 1:P) {
      # names
      tmp <- names_view[[i]]
      # data
      aux <- as.matrix(views_source[[which(names(views_source) %in% tmp)]])
      views_source.comb <- cbind(views_source.comb, aux)
    }
  }else{
    views_source.comb <-  as.matrix(views_source[[names_view]])
  }
  
  mas.std.learning.Y <- mas.mea.learning.Y <- vector("list", length = Ndrug) 
  names(mas.std.learning.Y) <- names(mas.mea.learning.Y) <- names(drug_source)
  
  # ****************
  # Model: NO imputation
  for (ii in 1:Ndrug){
    
    o = names(drug_source)[ii]
    cat(ii,".drug source: ", o, "\n")
    
    # separate input data in learning and prediction
    learning.X <- views_source.comb[learning_indices,]
    prediction.X <- views_source.comb[prediction_indices,]
    
    # separate output data in learning and prediction
    learning.Y <- drug_source[learning_indices, ii, drop = FALSE]
    validation.Y <- drug_source[prediction_indices, ii, drop = FALSE]
    
    if (standardize_any==T){
      
      # view normalization
      nan_indices_patients_learning = which(is.na(rowSums(learning.X)) == T)
      nan_indices_patients_prediction = which(is.na(rowSums(prediction.X)) == T)
      
      if(length(nan_indices_patients_learning) != 0) learning.X <- learning.X[-nan_indices_patients_learning,] 
      if(length(nan_indices_patients_prediction) != 0) prediction.X <- prediction.X[-nan_indices_patients_prediction,] 
      
      mas.mea.learning.X = colMeans(as.matrix(learning.X), na.rm = T)
      mas.std.learning.X = colSds(as.matrix(learning.X), na.rm = T)
      
      mas.std.learning.X [mas.std.learning.X == 0] = 1;
      
      learning.X = standarization(learning.X, mas.mea.learning.X, mas.std.learning.X)
      prediction.X = standarization(prediction.X, mas.mea.learning.X, mas.std.learning.X)
      
      cat("input data normalization done","\n")
      
      # drug response standardization
      if (o != "label"){
        
      if(length(nan_indices_patients_learning) != 0) learning.Y <- learning.Y[-nan_indices_patients_learning]
      if(length(nan_indices_patients_prediction) != 0) validation.Y <- validation.Y[-nan_indices_patients_prediction]
      
      mas.mea.learning.Y[[o]] = colMeans(as.matrix(learning.Y), na.rm = T)
      mas.std.learning.Y[[o]] = 1
      
      learning.Y = sweep(as.matrix(learning.Y), 2, mas.mea.learning.Y[[o]], FUN = "-")
      validation.Y = sweep(as.matrix(validation.Y),2, mas.mea.learning.Y[[o]], FUN = "-")
      
      }
      if (standardize_response == T){
        
        mas.std.learning.Y[[o]]= colSds(as.matrix(learning.Y), na.rm = T)
        mas.std.learning.Y[[o]][mas.std.learning.Y[[o]]  == 0] <- 1
        
        learning.Y = sweep(as.matrix(learning.Y),2, mas.std.learning.Y[[o]], FUN = "/") 
        validation.Y = sweep(as.matrix(validation.Y),2, mas.std.learning.Y[[o]], FUN = "/")
        
        cat("output data normalization done","\n")
      }
      
    }else{
      
      # if standardize_any=F we assume data have been already normalized,
      # NA are set to 0 (i.e. average value)
  
      # nan_indices_learning = is.na(learning.X)
      # learning.X[nan_indices_learning==F] <- 0
      
      # Remove samples with any NA value
      nan_indices_patients_learning = which(is.na(rowSums(learning.X)) == T)
      if(length(nan_indices_patients_learning) != 0) learning.X <- learning.X[-nan_indices_patients_learning,] 
      
      cat("careful: data are assumed to be already normalized (if that is not the case set standardize_any=T)","\n")
    }
    # Hyperparameters estimation
    if (o != "label"){
    state <- elastic_net_train(drug_source = learning.Y,
                               views_source = learning.X,
                               family = view_combination[[1]],
                               view_combination = names_view,
                               measure_type = "mse",
                               parameters = parameters,
                               parallelize = T, iteration = iteration)
    } else{
      
      state <- elastic_net_train(drug_source = learning.Y,
                                 views_source = learning.X,
                                 family = view_combination[[1]],
                                 view_combination = names_view,
                                 measure_type = "class",
                                 parameters = parameters,
                                 parallelize = T, iteration = iteration)
      
    }

    coef_values[[o]]<- state$cv.glmnet.features
    hyperparameter_values[[o]] <- state$cv.glmnet.hyperparameters
    freq_values[[o]] <- state$cv.glmnet.freq
    cat("training performed","\n")
    
    # perform prediction
    if (o != "label"){
      prediction_cv <- lapply(names(coef_values[[o]]), function(X){elastic_net_test(prediction.X, coef_values[[o]][[X]])})
      MSE[[o]] <- lapply(prediction_cv, function(X){apply((validation.Y - X)^2, 2, mean)})
      SpCorr[[o]]  <- lapply(prediction_cv, function(X){diag(cor(validation.Y, X, method = "spearman"))})
      PeCorr[[o]]  <- lapply(prediction_cv, function(X){diag(cor(validation.Y, X, method = "pearson"))})
      CI[[o]] <-   lapply(prediction_cv, function(X){civalue(validation.Y, X)})
      names(MSE[[o]]) <- names(SpCorr[[o]]) <- names(PeCorr[[o]]) <- names(CI[[o]]) <- names(coef_values[[o]])
      performances <- list(MSE=MSE, SpCorr=SpCorr, PeCorr=PeCorr, CI=CI)
      
    }else{
      prediction_cv <- lapply(names(coef_values[[o]]), function(X){log_elastic_net_test(prediction.X, coef_values[[o]][[X]])})
      CI[[o]] <- lapply(prediction_cv, function(X){civalue(validation.Y, X)})
      MR[[o]] <- lapply(prediction_cv, function(X){table(X, validation.Y$label)[1]/length(validation.Y$label)})
      names(MR[[o]]) <- names(CI[[o]]) <- names(coef_values[[o]])
      performances <- list(MR=MR, CI=CI)
      
    }
    cat("prediction computed","\n")

    model <- list(Coef=coef_values, hyperparameters=hyperparameter_values, freq_features=freq_values)
  }
  
  return(list(performances = performances,
              model = model,
              training_set = NULL,
              mas.std.learning.Y = mas.std.learning.Y,
              mas.mea.learning.Y = mas.mea.learning.Y,
              mas.std.learning.X = mas.std.learning.X,
              mas.mea.learning.X = mas.mea.learning.X))
 
}
  