##########################################################################################################################
# Script to perform bayesian efficiente multiple kernel algorithm :
##########################################################################################################################

BEMKL <- function(drug_source, views_source, view_combination, learning_indices, prediction_indices,
                  standardize_any=F, standardize_response=F,
                  parameters, iteration) {
  
  # ****************
  # packages
  library(pdist, lib.loc= Sys.getenv("R_LIBS"))
  
  # ****************
  # scripts
  source("./R/PanCancer_draft_v1/scaling_function.R")
  source("./R/PanCancer_draft_v1/civalue.R")
  source("./R/PanCancer_draft_v1/BEMKL_cluster/bemkl_supervised_multioutput_regression_variational_train.R")
  source("./R/PanCancer_draft_v1/BEMKL_cluster/bemkl_supervised_multioutput_regression_variational_test.R")

  # ****************
  # General variables:
  names_view <- names(view_combination) # Need data type 
  Ndrug <- ncol(drug_source) # number of tasks (output variables)
  P <- length(view_combination) # number of views
  
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
      
      if (view_combination[m] != "jaccard"){
        
        mas.mea.learning.X[[m]] = colMeans(learning.X[[m]], na.rm = T)
        mas.std.learning.X[[m]] = colSds(as.matrix(learning.X[[m]]), na.rm = T)
        
        ## same but slower
        # mas.mea2 = apply(learning.X[[m]], 2, mean, na.rm = T)
        # mas.std2 = apply(learning.X[[m]], 2, sd, na.rm = T)
        
        mas.std.learning.X[[m]][mas.std.learning.X[[m]]==0] = 1
        
        learning.X[[m]]= standarization(learning.X[[m]], mas.mea.learning.X[[m]], mas.std.learning.X[[m]])
        prediction.X[[m]] = standarization(prediction.X[[m]], mas.mea.learning.X[[m]],mas.std.learning.X[[m]])
        
      }
      # 
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
  
  # kernels for views
  Nlearning <- length(learning_indices)
  Nprediction <- length(prediction_indices)
  Kx_learning = array(rep(0, Nlearning*Nlearning*P), c(Nlearning, Nlearning, P))
  Kx_prediction = array(rep(0, Nlearning*Nprediction*P), c(Nlearning, Nprediction, P))
  
  for (m in 1:P){
    
    if (view_combination[m] == "gaussian"){
      Kx_learning[, , m] <- exp(-(as.matrix(dist(learning.X[[m]], method = "euclidean")))^2/ncol(learning.X[[m]])/2)
      Kx_prediction[, , m] <- exp(-(as.matrix(pdist(learning.X[[m]], prediction.X[[m]])))^2/ncol(learning.X[[m]])/2)
      
    }else if(view_combination[m] == "jaccard"){
      Kx_learning[, , m] <- 1 - as.matrix(dist(learning.X[[m]], method = "binary"))
      x.tmp <- as.matrix(learning.X[[m]])
      y.tmp <- as.matrix(prediction.X[[m]])
      Kx_prediction[, , m] <- do.call(cbind, lapply(1:nrow(y.tmp), function(y){
        do.call(rbind, lapply(1:nrow(x.tmp), function(x){
          jaccard(x.tmp[x,], y.tmp[y,])
        }))
      }))
      Kx_learning[, , m][is.na(Kx_learning[, , m])] <- 0
      Kx_prediction[, , m][is.na(Kx_prediction[, , m])] <- 0
    }
    
  }
  cat("kernels computed","\n")
  
  #set the number of outputs
  L <- ncol(learning.Y)
  #set the number of kernels
  P <- dim(Kx_learning)[3]
  
  #initialize the kernels and outputs for training
  Ktrain <- Kx_learning #should be an Ntra x Ntra x P matrix containing similarity values between training samples
  Ytrain <- t(learning.Y) #should be an L X Ntra matrix containing target outputs where L is the number of outputs
  
  #perform training
  state <- bemkl_supervised_multioutput_regression_variational_train(Ktrain, Ytrain, parameters)
  #display the kernel weights
  # print(state$be$mu[(L + 1):(L + P)])
  cat("training performed","\n")
  
  #initialize the kernels for testing
  Ktest <- Kx_prediction #should be an Ntra x Ntest x P matrix containing similarity values between training and test samples
  
  # perform prediction
  prediction <- bemkl_supervised_multioutput_regression_variational_test(Ktest, state)
  cat("prediction computed","\n")
  # #display the predictions
  # print(prediction$Y$mu)
  
  prediction_cv <- t(prediction$Y$mu)
  
  MSE <- apply((validation.Y - prediction_cv)^2, 2, mean)
  SpCorr <- diag(cor(validation.Y, prediction_cv, method = "spearman"))
  PeCorr <- diag(cor(validation.Y, prediction_cv, method = "pearson"))
  CI <- sapply(1:Ndrug, function(x) civalue(validation.Y[,x], prediction_cv[,x]))
  names(MSE) <- names(SpCorr) <- names(PeCorr) <- names(CI) <- colnames(drug_source)
  
  performances <- list(MSE=MSE, SpCorr=SpCorr, PeCorr=PeCorr, CI=CI)
  model <- state
  training_set <- learning.X
  
  return(list(performances = performances,
              model = model,
              training_set = training_set,
              mas.std.learning.Y = mas.std.learning.Y,
              mas.mea.learning.Y = mas.mea.learning.Y,
              mas.std.learning.X = mas.std.learning.X,
              mas.mea.learning.X = mas.mea.learning.X))
}
