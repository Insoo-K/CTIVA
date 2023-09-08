Fstats=function(X){
  cls = unique(X[,2])
  between_group_variance = 0
  within_group_variance = 0
  N = nrow(X)
  K = length(cls)
  m = mean(X[,1])
  for(i in cls){
    idx = X[,2]==i
    m2 = mean(X[idx,1])
    between_group_variance = between_group_variance + sum(idx)*(m-m2)^2
    within_group_variance = within_group_variance + sum((X[idx,1]-m2)^2)
  }
  between_group_variance = between_group_variance/(K-1)
  within_group_variance = within_group_variance/(N-K)
  Fstats = between_group_variance/within_group_variance
  return(Fstats)
}

permanova=function(X, B=10000){
  FF = NULL
  Fstats = Fstats(X)
  N=nrow(X)
  for(i in 1:B){
    X_ = X
    X_[,1] = X[sample(1:N,N,replace = F),1]
    FF = c(FF, Fstats(X_))
  }
  P = (sum(FF>=Fstats)+1)/(B+1)
  return(P)
}



