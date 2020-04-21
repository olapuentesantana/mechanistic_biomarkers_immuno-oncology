civalue<-function(y,yhat) {
  
  Ly<-length(y)
  ci<-rep(0,Ly)
  
  s<-0 # score sum
  n<-0 # pairs of scores
  
  for (i in 1:(Ly-1)) {
    
    for (j in (i+1):Ly) {
      
      if (y[i]>y[j]) {
        
        s<-s+(yhat[i] > yhat[j])+0.5*(yhat[i]==yhat[j]);
        n<-n + 1;
        
      } else if (y[i]<y[j]) {
        
        s<-s+(yhat[i]<yhat[j])+0.5*(yhat[i]==yhat[j]);
        n<-n+1;
        
      }
      
    }
    
  }
  
  ci<-s/n
  return(ci)
  
}