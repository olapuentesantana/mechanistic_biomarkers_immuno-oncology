####################################################################################
## Scaling function: Z score normalization
####################################################################################
standarization <- function(X, mean, sd){
  
    X.scale = matrix(0,nrow(X), ncol(X), dimnames = list(rownames(X),colnames(X)))
    
  if (missing(mean) & missing(sd)) {
     mean.X = colMeans(X,na.rm = T)
     sd.X = colSds(as.matrix(X),na.rm = T)
     X.scale <- sweep(X, 2, mean.X, FUN = "-")
     X.scale <- sweep(X.scale, 2, sd.X, FUN = "/")
  } else {
    X.scale <- sweep(X, 2, mean, FUN = "-")
    X.scale <- sweep(X.scale, 2, sd, FUN = "/")
    }
  #out <- list(X.scale = X.scale, mean.X = mean.X, sd.X = sd.X)
  return(as.matrix(X.scale))
}

