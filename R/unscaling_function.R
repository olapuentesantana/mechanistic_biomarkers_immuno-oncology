####################################################################################
## Removing Scaling function: Z score normalization
####################################################################################
remove_standarization <- function(X, mean, sd){
  
  X.scale = matrix(0,nrow(X), ncol(X), dimnames = list(rownames(X),colnames(X)))
  X.scale <- sweep(X, 2, sd, FUN = "*")
  X.scale <- sweep(X.scale, 2, mean, FUN = "+")
  
  return(as.matrix(X.scale))
}
