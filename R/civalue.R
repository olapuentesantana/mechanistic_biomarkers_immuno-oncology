civalue<-function(y,yhat) {
  
  y <- as.matrix(y) ; yhat <- as.matrix(yhat)
  ifelse(dim(y)[2] > 1, Ly<-nrow(y), Ly<-length(y))
  ci<-rep(0,Ly)
  
  s<-0 # score sum
  n<-0 # pairs of scores
  
  ci <- sapply(1:ncol(y), function(X){
    
    for (i in 1:(Ly-1)) {
      
      for (j in (i+1):Ly) {
        
          if (y[i,X]>y[j,X]) {
            
            s<-s+(yhat[i,X]>yhat[j,X])+0.5*(yhat[i,X]==yhat[j,X]);
            n<-n+1;
            
          } else if (y[i,X]<y[j,X]) {
            
            s<-s+(yhat[i,X]<yhat[j,X])+0.5*(yhat[i,X]==yhat[j,X]);
            n<-n+1;
            
          }
    
      }
      
    }
    ci<-s/n
    return(ci)
  })
  names(ci) <- colnames(y)
  return(ci)
  
}