# #########################################################################################################
# Script to understand how glmnet multi-gaussian works
# #########################################################################################################

load("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/PanCancer/GBM/new_remove_all_genes/DataViews_no_filter_GBM.RData")
X <- DataViews.no_filter$pathways
load("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/PanCancer/GBM/new/ImmuneResponse_no_filter_GBM_matrix_format.RData")
Y <- ImmuneResponse.no_filter

load("~/Desktop/PhD_TU:e/Research/mechanistic_biomarkers_immuno-oncology/data/parameters_4_all.RData")
parameters$mgaussian$alpha <- parameters[["L21"]][["alpha"]]
# ****************
# packages
library(glmnet)
library(parallel)

mas.mea.learning.X = colMeans(as.matrix(X), na.rm = T)
mas.std.learning.X = colSds(as.matrix(X), na.rm = T)

X <- sweep(X, 2,   mas.mea.learning.X, FUN = "-")
X <- sweep(X, 2,   mas.std.learning.X, FUN = "/")

mas.mea.learning.Y = colMeans(as.matrix(Y), na.rm = T)
mas.std.learning.Y = colSds(as.matrix(Y), na.rm = T)

Y <- sweep(Y, 2,   mas.mea.learning.Y, FUN = "-")
Y <- sweep(Y, 2,   mas.std.learning.Y, FUN = "/")


# ****************
# General variables:
names_view <- view_combination # Need data type 
N <- nrow(Y) # number of observations
Ndrug <- ncol(Y) # number of tasks (output variables)
P <- 1 # number of views
alpha_range <- parameters$mgaussian$alpha

# Normalization training set: Z score
## View:
views_source.Z <- X
## Response: 
drug_source.Z <- Y

# ****************
# Hyperparameters estimation: alpha and lambda

# initialize variables:
# Coefficients and mse values for each iteration
features.glmnet.1se <- features.glmnet.min <- lowest.mse.cv.1se <- lowest.mse.cv.min <- NULL

# data (mse, sd, coef) for each iteration, for every alpha
data_task_cv.glmnet <- NULL

# Parameters for each iteration 
parameters.glmnet.1se <- parameters.glmnet.min <- NULL

coef_for_each_alpha <- NULL
MSE_for_each_alpha <- NULL
sd_for_each_alpha <- NULL

iteration = 1
# Internal cross validation: number of folds by default is 10
assign("fit",lapply(alpha_range, function(X){
  # Same foldid for each alpha, but different foldid per iteration.
  set.seed(3000 + 10 * iteration)
  # Applying internal cross-validation:
  cv.glmnet(x = as.matrix(views_source.Z), 
            y = as.matrix(drug_source.Z), 
            family = "mgaussian",
            alpha = X, 
            type.measure = "mse", 
            standardize = FALSE, 
            standardize.response = FALSE,
            intercept = TRUE, 
            keep = TRUE,
            nfolds = 5,
            parallel = TRUE)}))

names(fit) <- paste0(alpha_range,"-")

# For each iteration, MSE, sd and coefficients for each alpha across range of lambdas. 
# MSE
MSE_for_each_alpha <- lapply(fit, function(X) {
  tmp = X$cvm
  names(tmp) = X$lambda
  return(tmp)})
# SD
sd_for_each_alpha <- lapply(fit, function(X) {
  tmp = X$cvsd
  names(tmp) = X$lambda
  return(tmp)})
# COEF
coef_for_each_alpha <- lapply(fit, function(X) {
  
  coef_for_each_alpha <- lapply(colnames(drug_source.Z), function(Y) {
    tmp = coef(X, s = X$lambda)[[Y]]
    tmp@Dimnames[[2]] = as.character(X$lambda)
  return(tmp)})
  names(coef_for_each_alpha) <- colnames(drug_source.Z)
  
return(coef_for_each_alpha)})

# Saving values -->
data_task_cv.glmnet <- list(MSE = MSE_for_each_alpha, SD = sd_for_each_alpha, coef = coef_for_each_alpha)

# Pair of lambda-alpha values:
value_MSE_alphas <- NULL

# Find alpha-lambda: 1. min MSE, 2. min MSE + 1SE
value_MSE_alphas <- unlist(data_task_cv.glmnet$MSE)
all_alphas = sapply(strsplit(names(value_MSE_alphas), split = "-.", fixed = TRUE), head, 1)
all_lambdas = sapply(strsplit(names(value_MSE_alphas), split = "-.", fixed = TRUE), function(X) return(X[2]))

## 1. Lambda-alpha value leading to min MSE
min_MSE <- value_MSE_alphas[which(value_MSE_alphas == min(value_MSE_alphas,na.rm = TRUE), arr.ind=TRUE)]
alpha_min_MSE = sapply(strsplit(names(min_MSE), split = "-.", fixed = TRUE), head, 1)
lambda_min_MSE = sapply(strsplit(names(min_MSE), split = "-.", fixed = TRUE), function(X) return(X[2]))

# * strange case where two lambdas obtain the same mse, therefore we have 2 pair alpha-lambda
if (length(alpha_min_MSE) > 1 & length(lambda_min_MSE) > 1){
  alpha_min_MSE <- unique(alpha_min_MSE)
  lambda_min_MSE <- max(lambda_min_MSE)
}

## 2.  Lambda-alpha value leading leading to 1 SE + min MSE 
one.se <- data_task_cv.glmnet$SD[[paste0(alpha_min_MSE,"-")]][which(names(data_task_cv.glmnet$SD[[paste0(alpha_min_MSE,"-")]]) == lambda_min_MSE)]
high <- min_MSE + one.se
interval_mse <- value_MSE_alphas[value_MSE_alphas < high]
interval_mse <- interval_mse[which(as.numeric(sapply(strsplit(names(interval_mse), split = "-.", fixed = TRUE), function(X) return(X[2])))
                                   >= as.numeric(lambda_min_MSE))]
one.se_MSE <- interval_mse[which(interval_mse == max(interval_mse))]
lambda_1se_MSE <- sapply(strsplit(names(one.se_MSE), split = "-.", fixed = TRUE), function(X) return(X[2]))
alpha_1se_MSE <- sapply(strsplit(names(one.se_MSE), split = "-.", fixed = TRUE), head, 1)

# Saving parameters for each iteration:
parameters.glmnet.1se <- list(alpha_1se_MSE, lambda_1se_MSE)  
parameters.glmnet.min <- list(alpha_min_MSE, lambda_min_MSE)
names(parameters.glmnet.1se) <- names(parameters.glmnet.min) <- c("alpha", "lambda")

# 3. Selecting coefficients for both pairs of alpha-lambda
# Alpha-lambda min
coef_min_per_tasks <- do.call(cbind, lapply(colnames(drug_source.Z), function (Z){
  coef_min <- as.matrix(data_task_cv.glmnet$coef[[paste0(alpha_min_MSE,"-")]][[Z]][,lambda_min_MSE])
  return(coef_min)
}))
colnames(coef_min_per_tasks) <- colnames(drug_source.Z)
features.glmnet.min <- coef_min_per_tasks 

# Alpha-lambda 1se
coef_1se_per_tasks <- do.call(cbind, lapply(colnames(drug_source.Z), function (Z){
  coef_1se <- as.matrix(data_task_cv.glmnet$coef[[paste0(alpha_1se_MSE,"-")]][[Z]][,lambda_1se_MSE])
  return(coef_1se)
}))
features.glmnet.1se <- coef_1se_per_tasks 

# mse.1se and mse.min from internal cross validation for iteration k
lowest.mse.cv.1se = one.se_MSE
lowest.mse.cv.min = min_MSE

# Per task
Coef_cv.glmnet <- list(features.glmnet.1se, features.glmnet.min) ; names(Coef_cv.glmnet) <- c("1se.mse","min.mse")
MSE_cv.glmnet <- list(lowest.mse.cv.1se, lowest.mse.cv.min) ; names(MSE_cv.glmnet) <- c("1se.mse","min.mse")
Hyperparameters_cv.glmnet <-  list(parameters.glmnet.1se, parameters.glmnet.min) ; names(Hyperparameters_cv.glmnet) <- c("1se.mse","min.mse")
