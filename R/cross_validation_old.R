###########################################################################################
# Script to perform cross validation regardless of algorithm choice.
###########################################################################################
cross_validation <- function(drug_source, views_source, view_combination, algorithm, standardize_any = F, standardize_response=F,
                             parameters, k_fold=5, random=NULL) {

  # ****************
  # packages
  library(matrixStats)
  library(pdist)
  
  # ****************
  # scripts
  source("./R/scaling_function.R")
  source("./R/unscaling_function.R")
  source("./R/elastic_net_test.R")
  source("./R/elastic_net_train.R")
  source("./R/civalue.R")
  source("./R/bemkl_supervised_multioutput_regression_variational_train.R")
  source("./R/bemkl_supervised_multioutput_regression_variational_test.R")
  source("./R/bemkl_supervised_multioutput_regression_variational_test.R")
  
  # ****************
  # General variables:
  names_view <- names(view_combination) # Need data type 
  N <- nrow(drug_source) # number of observations
  Ndrug <- ncol(drug_source) # number of tasks (output variables)
  P <- length(view_combination) # number of views
  
  # prepare for k-fold CV
  rand_obs_ix <- sample(N) # prepare index of randomized observatios
  folds <- cut(seq(1,N),breaks=k_fold,labels=FALSE) #divide in k folds
  folds_ix <- lapply(1:k_fold, function(i){ 
    rand_obs_ix[which(folds==i,arr.ind=TRUE)]
  })
  
  # Linear regression with regularization
  if (algorithm %in% c("Lasso", "Elastic_Net", "CV_linear_reg_L1&L2","Ridge")){ 
    
    MSE <- SpCorr <- PeCorr <- CI <- coef_values <- x.train <- hyperparameter_values <- freq_values <- vector("list", length = Ndrug)
    names(MSE) <- names(SpCorr) <- names(PeCorr) <- names(CI) <- names(coef_values) <- names(x.train) <- names(hyperparameter_values) <- names(freq_values) <- names(drug_source)
  
    # When combining data sets:
    name = NULL
    views_source.comb = NULL
    if (P > 1){
      for (i in 1:P) {
        # names
        tmp <- names_view[[i]]
        #name <- paste0(name, tmp)
        # data
        aux <- as.matrix(views_source[[which(names(views_source) %in% tmp)]])
        views_source.comb <- cbind(views_source.comb, aux)
      }
      #names_view <- name
    }else{
      views_source.comb <-  as.matrix(views_source[[names_view]])
    }
    
    mas.std.learning.Y <- mas.mea.learning.Y <- vector("list", length = Ndrug) 
    mas.std.learning.X <- mas.mea.learning.X <- vector("list", length = random)
  
    mas.std.learning.Y <- lapply(1:Ndrug, function(X){mas.std.learning.Y[[X]] = vector("list", length = random)})
    mas.mea.learning.Y <- mas.std.learning.Y
    
  }else{ # BEMKL
    
    mas.std.learning.X <- mas.mea.learning.X <- vector("list", length = random)
    mas.std.learning.Y <- mas.mea.learning.Y <- vector("list", length = random)
    
  }
  
  performances <- list()
  model <- list()
  training_set <- list()
  
  #if random is set to a number, a randomised cross-validation, selecting each time
  #20% of the samples is performed and repeated the specified amount of times
  if (!is.null(random)){
    k_fold = random
  }
  
  # ****************
  # k-fold CV:
  #if random is set to a number, a randomised cross-validation, selecting each time
  #20% of the samples is performed and repeated the specified amount of times
  
  for (k in 1:k_fold){
    cat("\n---\niteration", k,":", "\n")
    
    # indices of observations used for learning (training) and prediction (test)
    if (is.null(random)){
      learning_indices = setdiff(1:N, folds_ix[[k]]);
      prediction_indices = folds_ix[[k]];
    }else{
      prediction_indices = sample(1:N, round(N*20/100))
      learning_indices = setdiff(1:N, prediction_indices)
    }
    
    # *******************
    # Elastic Net or BEMKL
    
    if (algorithm %in% c("Lasso", "Elastic_Net", "CV_linear_reg_L1&L2","Ridge")){

      ## Elastic Net ##
     for (ii in 1:Ndrug){
        o = names(drug_source)[ii]
        cat(ii,".drug source: ", o, "\n")
        
        # IMPUTATION NO:

        # separate input data in learning and prediction
        learning.X <- views_source.comb[learning_indices,]
        prediction.X <- views_source.comb[prediction_indices,]
        
        # separate output data in learning and prediction
        learning.Y <- drug_source[learning_indices, ii]
        validation.Y <- drug_source[prediction_indices, ii]

        if (standardize_any==T){
          # view normalization
          nan_indices_patients_learning = which(is.na(rowSums(learning.X)) == T)
          nan_indices_patients_prediction = which(is.na(rowSums(prediction.X)) == T)
          
          if(length(nan_indices_patients_learning) != 0) learning.X <- learning.X[-nan_indices_patients_learning,] 
          if(length(nan_indices_patients_prediction) != 0) prediction.X <- prediction.X[-nan_indices_patients_prediction,] 
          
          mas.mea.learning.X[[k]] = colMeans(as.matrix(learning.X), na.rm = T)
          mas.std.learning.X[[k]] = colSds(as.matrix(learning.X), na.rm = T)

          mas.std.learning.X[[k]] [mas.std.learning.X[[k]] == 0] = 1;

          learning.X = standarization(learning.X, mas.mea.learning.X[[k]], mas.std.learning.X[[k]])
          prediction.X = standarization(prediction.X, mas.mea.learning.X[[k]], mas.std.learning.X[[k]])

          ## same but slower
          # mas.mea2 = apply(learning.X[[m]], 2, mean, na.rm = T)
          # mas.std2 = apply(learning.X[[m]], 2, sd, na.rm = T)

          cat("input data normalization done","\n")
          
          # drug response standardization
          # mas.mea = apply(learning.Y, 2, mean, na.rm = T)
          
          if(length(nan_indices_patients_learning) != 0) learning.Y <- learning.Y[-nan_indices_patients_learning]
          if(length(nan_indices_patients_prediction) != 0) validation.Y <- validation.Y[-nan_indices_patients_prediction]

          mas.mea.learning.Y[[ii]][[k]] = colMeans(as.matrix(learning.Y), na.rm = T)
          mas.std.learning.Y[[ii]][[k]] = 1
          
          learning.Y = sweep(as.matrix(learning.Y), 2, mas.mea.learning.Y[[ii]][[k]], FUN = "-")
          validation.Y = sweep(as.matrix(validation.Y),2, mas.mea.learning.Y[[ii]][[k]], FUN = "-")
          
          if (standardize_response == T){
            
            mas.std.learning.Y[[ii]][[k]] = colSds(as.matrix(learning.Y), na.rm = T)
            mas.std.learning.Y[[ii]][[k]][mas.std.learning.Y[[ii]][[k]]  == 0] <- 1
            
            learning.Y = sweep(as.matrix(learning.Y),2, mas.std.learning.Y[[ii]][[k]], FUN = "/") 
            validation.Y = sweep(as.matrix(validation.Y),2, mas.std.learning.Y[[ii]][[k]], FUN = "/")
          }

          cat("output data normalization done","\n")
      
        }else{
          # if standardize_any=F we assume data have been already normalized,
          # NA are set to 0 (i.e. average value)
          # nan_indices_learning = is.na(learning.X)
          # learning.X[nan_indices_learning==T] <- 0
          
          nan_indices_patients_learning = which(is.na(rowSums(learning.X)) == T)
          if(length(nan_indices_patients_learning) != 0) learning.X <- learning.X[-nan_indices_patients_learning,] 
          
          cat("careful: data are assumed to be already normalized (if that is not the case set standardize_any=T)","\n")
        }
        # Hyperparameters estimation
        state <- elastic_net_train(drug_source = learning.Y,
                                   views_source = learning.X,
                                   view_combination = names_view,
                                   family = view_combination[[1]],
                                   measure_type = "mse",
                                   parameters = parameters[[algorithm]],
                                   parallelize = F, iteration = k)
        
        coef_values[[o]]<- state$cv.glmnet.features
        hyperparameter_values[[o]] <- state$cv.glmnet.hyperparameters
        freq_values[[o]] <- state$cv.glmnet.freq
        cat("training performed","\n")
        
        # perform prediction
        
        prediction_cv <- lapply(names(coef_values[[o]]), function(X){elastic_net_test(prediction.X, coef_values[[o]][[X]])})
        cat("prediction computed","\n")
      
        MSE[[o]] <- lapply(prediction_cv, function(X){apply((validation.Y - X)^2, 2, mean)})
        SpCorr[[o]]  <- lapply(prediction_cv, function(X){diag(cor(validation.Y, X, method = "spearman"))})
        PeCorr[[o]]  <- lapply(prediction_cv, function(X){diag(cor(validation.Y, X, method = "pearson"))})
        CI[[o]] <-   lapply(prediction_cv, function(X){civalue(validation.Y, X)})
        
        names(MSE[[o]]) <- names(SpCorr[[o]]) <- names(PeCorr[[o]]) <- names(CI[[o]]) <- names(coef_values[[o]])
        
        performances[[k]] <- list(MSE=MSE, SpCorr=SpCorr, PeCorr=PeCorr, CI=CI)
        model[[k]] <- list(Coef=coef_values, hyperparameters=hyperparameter_values, freq_features=freq_values )
        training_set[[k]] <- NULL
    }
      
    }else if (algorithm == "BEMKL") {
      ## BEMKL ##
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

            mas.mea.learning.X[[k]][[m]] = colMeans(learning.X[[m]], na.rm = T)
            mas.std.learning.X[[k]][[m]] = colSds(as.matrix(learning.X[[m]]), na.rm = T)
            
            ## same but slower
            # mas.mea2 = apply(learning.X[[m]], 2, mean, na.rm = T)
            # mas.std2 = apply(learning.X[[m]], 2, sd, na.rm = T)
            
            mas.std.learning.X[[k]][[m]][mas.std.learning.X[[k]][[m]]==0] = 1
            
            learning.X[[m]]= standarization(learning.X[[m]], mas.mea.learning.X[[k]][[m]], mas.std.learning.X[[k]][[m]])
            prediction.X[[m]] = standarization(prediction.X[[m]], mas.mea.learning.X[[k]][[m]],mas.std.learning.X[[k]][[m]])
            
          }
          # 
          # learning.X[[m]][nan_indices_learning==T] <- 0
          # prediction.X[[m]][nan_indices_prediction==T] <- 0
        }
        cat("input data normalization done","\n")
        
        # drug response standardization
        # mas.mea = apply(learning.Y, 2, mean, na.rm = T)
      
        mas.mea.learning.Y[[k]] = colMeans(as.matrix(learning.Y), na.rm = T)
        mas.std.learning.Y[[k]] =  rep(1, ncol(learning.Y)) #
        
        learning.Y <- sweep(learning.Y, 2,   mas.mea.learning.Y[[k]], FUN = "-")
        validation.Y <- sweep(validation.Y, 2, mas.mea.learning.Y[[k]], FUN = "-")
        
        if (standardize_response == T){
          
          mas.std.learning.Y[[k]] = colSds(as.matrix(learning.Y), na.rm = T)
          mas.std.learning.Y[[k]][mas.std.learning.Y[[k]] == 0] <- 1
          
          learning.Y <- sweep(learning.Y, 2, mas.std.learning.Y[[k]], FUN = "/")
          validation.Y <- sweep(validation.Y, 2, mas.std.learning.Y[[k]], FUN = "/")
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
      state <- bemkl_supervised_multioutput_regression_variational_train(Ktrain, Ytrain, parameters$BEMKL)
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
      
      performances[[k]] <- list(MSE=MSE, SpCorr=SpCorr, PeCorr=PeCorr, CI=CI)
      model[[k]] <- state
      training_set[[k]] <- learning.X
      
    }
  }
  
  return(list(performances = performances,
              model = model,
              training_set = training_set,
              mas.std.learning.Y = mas.std.learning.Y,
              mas.mea.learning.Y = mas.mea.learning.Y,
              mas.std.learning.X = mas.std.learning.X,
              mas.mea.learning.X = mas.mea.learning.X))
}

# ****************
