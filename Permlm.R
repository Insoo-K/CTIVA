Tstat=function(X){
  n = nrow(X)
  beta1hat = (sum(X[,1]*X[,2])-(sum(X[,2])*sum(X[,1])/n))/sum((X[,1]-mean(X[,1]))^2)
  beta0hat = mean(X[,2])-beta1hat*mean(X[,1])
  yhat = beta0hat + beta1hat*X[,1]
  err = X[,2] - yhat
  sebeta1hat = sqrt((sum(err^2)/(n-2))/sum((X[,1]-mean(X[,1]))^2))
  Tstat = beta1hat/sebeta1hat

  return(Tstat)
}

permlm=function(X, B=10000){
  Tvalues = NULL
  Tstat = Tstat(X)
  N=nrow(X)
  for(i in 1:B){
    X_ = X
    X_[,1] = X[sample(1:N,N,replace = F),1]
    Tvalues = c(Tvalues, Tstat(X_))
  }
  if(Tstat>=0){  P = 2*(sum(Tvalues>=Tstat)+1)/(B+1)}
  else{  P = 2*(sum(Tvalues<=Tstat)+1)/(B+1)}
  return(P)
}



