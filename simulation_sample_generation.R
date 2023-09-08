# library
library(survival);
library(pROC);
library(mnormt)
source("Permanova.R")
source("Permlm.R")

# functions
rowVars = function(E,na.rm=FALSE) {
  n = ncol(E);
  m = rowMeans(E,na.rm=na.rm);
  mm = rowMeans(E^2,na.rm=na.rm);
  v = (n*mm - n*m^2)/(n-1);
  return( v );
}

random_data_surv = function(N,lt=1,lc=1/2,of=0) {
  T1 = rexp(N,lt); T2 = T1 + rexp(N,lt);
  C1 = rexp(N,lc)+of; C2 = rexp(N,lc)+of;
  X1 = Surv(apply(cbind(T1,C1),1,min),(T1<=C1) );
  X2 = Surv(apply(cbind(T2,C2),1,min),(T2<=C2) );
  X = cbind(as.data.frame(X1),X2); colnames(X) = c('X1','X2');
  A = c(0,max(X1[,1])*1.1,0,max(X2[,1])*1.1);
  # X: bivariate survival times
  # A: region of interest
  # T: true event times
  return( list(X=X, A=A, T=cbind(T1,T2)) );
}

# log-normal
random_data_lognorm = function(N,TR=0.5,CR=0,of=0) {
  if( TR == 1 ) {
    T1 = T2 = exp(rnorm(N,0,1));
  } else {
    T = rmnorm(N,c(0,0.5),matrix(c(1,TR,TR,1),nrow=2));
    T1 = exp(T[,1]); T2 = exp(T[,2]);
  }
  if( CR == 1 ) {
    C1 = C2 = exp(rnorm(N,0,1))+of;
  } else {
    C = rmnorm(N,c(0,0),matrix(c(1,CR,CR,1),nrow=2));
    C1 = exp(C[,1]+of); C2 = exp(C[,2]+of);
  }
  X1 = Surv(apply(cbind(T1,C1),1,min),(T1<=C1) );
  X2 = Surv(apply(cbind(T2,C2),1,min),(T2<=C2) );
  X = cbind(as.data.frame(X1),X2); colnames(X) = c('X1','X2');
  A = c(0,max(X1[,1])*1.1,0,max(X2[,1])*1.1);
  return( list(X=X, A=A, T=cbind(T1,T2)) );
}

# clayton q=1
random_data_clayton = function(N,lambda=1/2,q=1,of=0,A=NULL,common.censor=FALSE) {
  U1 = runif(N); U2 = runif(N);
  a = (1-U2)^(-1/q);
  
  T1 = q*log( (1-a) + a*(1-U1)^(-1/(1+q)) );
  T2 = -log(1-U2);
  
  C1 = rexp(N,lambda)+of;
  C2 = rexp(N,lambda)+of;
  if( common.censor ) { C2 = C1; }
  
  X1 = Surv(apply(cbind(T1,C1),1,min),(T1<=C1) );
  X2 = Surv(apply(cbind(T2,C2),1,min),(T2<=C2) );
  X = cbind(as.data.frame(X1),X2); colnames(X) = c('X1','X2');
  if( is.null(A) ) { A = c(0,max(X1[,1])*1.1,0,max(X2[,1])*1.1); }
  
  return( list(X=X, A=A, T=cbind(T1,T2)) );
}

random_data_cov_combined = function(T,NV=c(100,100,100,700),MAXC=5, rate=0.5, SNOISE) {
  
  T1 = T[,1]; T2 = T[,2]; N = nrow(T);
  NC = sample(2:MAXC, sum(NV), replace=TRUE);
  C1 = t(matrix(rep(T1,NV[1])+rnorm(N*NV[1],0,SNOISE^2),nrow=N));
  C2 = t(matrix(rep(T2,NV[2])+rnorm(N*NV[2],0,SNOISE^2),nrow=N));
  Cd = t(matrix(rep(T2-T1,NV[3])+rnorm(N*NV[3],0,SNOISE^2),nrow=N));
  Cr = matrix(rnorm(N*NV[4],0,1),ncol=N);
  
  idx_cont_1 = 1:round(NV[1]*rate)
  idx_cont_2 = 1:round(NV[2]*rate)
  idx_cont_d = 1:round(NV[3]*rate)
  idx_cont_r = 1:round(NV[4]*rate)
  
  idx_cate_1 = (round(NV[1]*rate)+1):NV[1]
  idx_cate_2 = (round(NV[2]*rate)+1):NV[2]
  idx_cate_d = (round(NV[3]*rate)+1):NV[3]
  idx_cate_r = (round(NV[4]*rate)+1):NV[4]
  
  C1_cate = factorize(C1[idx_cate_1,], NC[idx_cate_1])
  C2_cate = factorize(C2[idx_cate_2,], NC[NV[1]+idx_cate_2])
  Cd_cate = factorize(Cd[idx_cate_d,], NC[NV[1]+NV[2]+idx_cate_d])
  Cr_cate = factorize(Cr[idx_cate_r,], NC[NV[1]+NV[2]+NV[3]+idx_cate_r])
  
  COV = cbind(t(C1[idx_cont_1,]), C1_cate,
              t(C2[idx_cont_2,]), C2_cate,
              t(Cd[idx_cont_d,]), Cd_cate,
              t(Cr[idx_cont_r,]), Cr_cate);
  
  colnames(COV) = sprintf('V%04d',1:sum(NV));
  
  return( COV );
}

factorize = function(M, NC){
  nc = ncol(M); nr = nrow(M);
  df = data.frame(matrix(rep(0,nc*nr), nrow=nc, ncol=nr))
  for(i in 1:nr){
    boundary = runif(NC[i]-1, min=0.1, max=0.9)
    df[,i] = separating(M[i,], boundary)
  }
  return(df)
}

separating = function(v, bound){
  b = quantile(v, sort(bound))
  cls = sample(1:(length(b)+1),length(b)+1,replace=FALSE)
  cv = v
  cv[v<b[1]]=cls[1]
  cv[b[length(b)]<=v]=cls[length(cls)]
  for(i in 1:(length(b)-1)){
    idx = b[i]<=v & v<b[i+1]
    cv[idx] = cls[i+1]
  }
  return(as.factor(cv))
}


NV = c(100,100,100,700)
Cond = c(rep(0,NV[1]+NV[2]), rep(1,NV[3]), rep(0,NV[4]))

dir_d = sprintf("./data")

if(!dir.exists(dir_d)) dir.create(dir_d)

for(scenario in 1:3){
  for(iter in 1:1){
    print(sprintf("(%d/%d)",(scenario-1)*100+iter, 300))
    # generate survival times and gene expression matrix
    set.seed(123456*iter);
    if(scenario==1) {D = random_data_surv(500);}
    if(scenario==2) {D = random_data_lognorm(500);}
    if(scenario==3) {D = random_data_clayton(500);}
    
    File_True = sprintf('%s/trueT_%d_%d.txt',dir_d,scenario, iter)
    File_S = sprintf('%s/surv_%d_%d.txt',dir_d,scenario, iter)
    File_C = sprintf('%s/covariate_%d_%d.txt',dir_d,scenario, iter)
    
    # making two input files: input.surv.txt and input.expr.txt
    x.out  = cbind(D$X[,1][,1],D$X[,1][,2],D$X[,2][,1],D$X[,2][,2]);
    write.table(round(x.out*100000)/100000,file=File_S,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t');
    
    t.out  = cbind(D$T[,1],D$T[,2]);
    write.table(round(t.out*100000)/100000,file=File_True,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t');
    
    # read data from files
    TMP = read.table(File_S,sep='\t',header=FALSE);
    # censored survival time
    X1 = Surv(TMP[,1],TMP[,2]); X2 = Surv(TMP[,3],TMP[,4]);
    X = cbind(as.data.frame(X1),X2);
    # covariates
    COV = random_data_cov_combined(D$T, NV=NV,SNOISE=1)
    write.table(COV, file=File_C, sep="\t",row.names=F, col.names=F)
  }
}






