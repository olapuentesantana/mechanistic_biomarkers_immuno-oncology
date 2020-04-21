###########################################################################################
# Script to estimate the hyperparameters, specific for L21 regularization
# * Reminder: 
## Alpha corresponds to parameter lam2
## Lambda corresponds to parameter lam1
###########################################################################################

L21_train <- function(drug_source, views_source, view_combination, standardize_any = F, standardize_response = F,
                              family, parameters, k_fold = 5, iteration) {
  
  # ****************
  # packages
  library(RMTL)
  
  # ****************
  # General variables:
  names_view <- view_combination # Need data type 
  alpha_range <- parameters$alpha
  data <- list()
  view_source.comb <- do.call(cbind, lapply(1:length(views_source), function(X) return(views_source[[X]])))
  features <- colnames(view_source.comb)
  drugs <- colnames(drug_source)
  
  data$X <- list(view_source.comb)
  data$Y <- list(as.matrix(drug_source[,1]))
  
  if (length(drugs) >= 2){
    for (ii in 2:length(drugs)){
      data$X[[ii]] <-  as.matrix(view_source.comb)
      data$Y[[ii]] <- as.matrix(drug_source[,ii])
    }
  }

  # L21 built-in cross validation
  assign("fit",lapply(alpha_range, function(X){
    
    # Same foldid for each alpha, but different foldid per iteration.
    set.seed(3000 + 10 * iteration)
    
    # Applying internal cross-validation:
    cvfit <- cvMTL(data$X , data$Y, type="Regression", Regularization="L21", Lam2 = as.numeric(X) , nfolds = k_fold)
  
    # Train
    model <- MTL(data$X, data$Y, type="Regression", Regularization="L21",
                 Lam1=cvfit$Lam1.min, Lam1_seq=cvfit$Lam1_seq, Lam2 = as.numeric(X))
    
    return(list(cv = list(Lam1_seq = cvfit$Lam1_seq, cvm = cvfit$cvm), model = model))
    
  }))
  
  names(fit) <- paste0(alpha_range,"-")
  
  
  # For each iteration, MSE, and coefficients for each alpha across range of lambdas. 
  # MSE
  MSE_for_each_alpha <- lapply(fit, function(X) {
    tmp <- X$cv$cvm
    names(tmp) <- X$cv$Lam1_seq
    return(tmp)})

  # COEF
  coef_for_each_alpha <- lapply(fit, function(X) {
    tmp <- X$model$W ; rownames(tmp) <- features ; colnames(tmp) <- drugs
    intercept <- X$model$C
    tmp <- rbind(intercept, tmp)
    tmp_list <- list(tmp) ; names(tmp_list) <- X$model$Lam1
    return(tmp_list)})
  
  # Saving values -->
  data_task_cv.MTL <- list(MSE = MSE_for_each_alpha, coef = coef_for_each_alpha)

  # Pair of lambda-alpha values:
  value_MSE_alphas <- NULL
  
  # Find alpha-lambda: 1. min MSE, 2. min MSE + 1SE
  value_MSE_alphas <- unlist(data_task_cv.MTL$MSE)
  all_alphas <- sapply(strsplit(names(value_MSE_alphas), split = "-.", fixed = TRUE), head, 1)
  all_lambdas <- sapply(strsplit(names(value_MSE_alphas), split = "-.", fixed = TRUE), function(X) return(X[2]))
  
  ## 1. Lambda-alpha value leading to min MSE
  min_MSE <- value_MSE_alphas[which(value_MSE_alphas == min(value_MSE_alphas,na.rm = TRUE), arr.ind=TRUE)]
  alpha_min_MSE <- sapply(strsplit(names(min_MSE), split = "-.", fixed = TRUE), head, 1)
  lambda_min_MSE <- sapply(strsplit(names(min_MSE), split = "-.", fixed = TRUE), function(X) return(X[2]))
  
  # * strange case where two lambdas obtain the same mse, therefore we have 2 pair alpha-lambda
  if (length(alpha_min_MSE) > 1 & length(lambda_min_MSE) > 1){
    alpha_min_MSE <- unique(alpha_min_MSE)
    lambda_min_MSE <- max(lambda_min_MSE)
  }
  
  # 2. Selecting coefficients for both pairs of alpha-lambda
  # Alpha-lambda min
  coef_min <- as.matrix(data_task_cv.MTL$coef[[paste0(alpha_min_MSE,"-")]][[lambda_min_MSE]])
  features.MTL.min <- coef_min 
  
  # All tasks
  Coef_cv.MTL <- list(min.mse = features.MTL.min)
  MSE_cv.MTL <- list(min.mse = min_MSE) 
  Hyperparameters_cv.MTL <-  list(min.mse = list(Lam2 = alpha_min_MSE, Lam1 = lambda_min_MSE))
  
  return(list(cv.MTL.data = data_task_cv.MTL,
              cv.MTL.hyperparameters = Hyperparameters_cv.MTL,
              cv.MTL.features = Coef_cv.MTL,
              cv.MTL.mse = MSE_cv.MTL))

}
