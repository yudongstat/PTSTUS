rm(list=ls())
##### load the data
transactions<-read.csv(file = "G:/My Drive/Topic7_missing labels/missing_label_revision/datasets/european supermarket data/data/Supermarket.csv",header = TRUE)
myX<-transactions[which((transactions$BeginDate=="16/02/2019")&(transactions$WorkstationGroupID==1)),10]
myY<-transactions[which((transactions$BeginDate=="23/02/2019")&(transactions$WorkstationGroupID==1)),10]

##### this function implements the naive unpaired t test
test.naive.oracle <- function(dat)
{
  X <- dat$X
  Y <- dat$Y
  n <- length(X)
  
  #naive two sample test, pretending independent X and Y
  XY.dif.mean = mean(X)-mean(Y)
  XY.dif.sd = (var(X)+var(Y))^0.5/n^0.5
  T_n.naive.two = XY.dif.mean/XY.dif.sd
  p.value.naive.two = 2*(1-pnorm(abs(T_n.naive.two)))
  
  return(result<-list(p.value.naive.two=p.value.naive.two))
}
dat<-list(X=myX,Y=myY)
n<-length(myX)
n1<-length(myX)
n2<-length(myY)
test.naive.oracle(dat)

##### this function implements the proposed sample-splitting-based test T_n
test.our <- function(dat,m_n1,m_n2)
{
  X <- dat$X
  Y <- dat$Y
  n <- length(X)
  m_n <- m_n1
  groupsize<-rep(n%/%m_n,m_n)
  if(n%%m_n>0){
    for(i in 1:(n%%m_n)){
      groupsize[i]<-groupsize[i]+1
    }
  }
  cum_groupsize<-rep(0,m_n)
  for(i in 1:m_n){
    cum_groupsize[i]<-sum(groupsize[1:i])
  }
  
  #Split X
  Ybar = mean(Y)
  Y.sd = sd(Y)
  X.sd = sd(X)
  Un <- rep(0,m_n)
  for (k in 1:m_n)
  {
    if(k==1){
      X.k<-X[1:cum_groupsize[1]]
    }else{
      X.k<-X[(cum_groupsize[k-1]+1):cum_groupsize[k]]
    }
    Un[k] <- groupsize[k]*(mean(X.k)-Ybar)^2/(X.sd)^2
  }
  
  T_n1 <- (sum(Un)-m_n)/(2*m_n)^0.5
  
  n <- length(Y)
  m_n <- m_n2
  groupsize<-rep(n%/%m_n,m_n)
  if(n%%m_n>0){
    for(i in 1:(n%%m_n)){
      groupsize[i]<-groupsize[i]+1
    }
  }
  cum_groupsize<-rep(0,m_n)
  for(i in 1:m_n){
    cum_groupsize[i]<-sum(groupsize[1:i])
  }
  #Split Y
  Xbar = mean(X)
  X.sd = sd(X)
  Y.sd = sd(Y)
  Sn <- rep(0,m_n)
  for (k in 1:m_n)
  {
    if(k==1){
      Y.k<-Y[1:cum_groupsize[1]]
    }else{
      Y.k<-Y[(cum_groupsize[k-1]+1):cum_groupsize[k]]
    }
    Sn[k] <- groupsize[k]*(mean(Y.k)-Xbar)^2/(Y.sd)^2
  }
  
  T_n2 <- (sum(Sn)-m_n)/(2*m_n)^0.5
  
  #combining two samples together
  T_n3 <- (T_n1+T_n2)/sqrt(2)
  p.value1 <- 1-pnorm(T_n3)
  
  return(result<-list(p.value1=p.value1,T_n3=T_n3))
}

###compute Tn with no anticipation
result <- test.our(dat=dat,m_n1 = 125,m_n2 = 120)
result


##### this function implements the proposed subsampling-based test W_n
test.our <- function(dat,m_n=numeric(0),n_m=numeric(0))
{
  X <- dat$X
  Y <- dat$Y
  X.all <- dat$X.all
  n <- length(X)
  #Split X
  Ybar = mean(Y)
  Xbar = mean(X)
  Y.sd = sd(Y)
  X.sd = sd(X)
  Un <- rep(0,m_n)
  Sn <- rep(0,m_n)
  for (k in 1:m_n)
  {
    candidate<-sort(sample(x=1:n1,size = n_m))
    X.k<-X[candidate]
    Un[k] <- sqrt(n_m)*(mean(X.k)-Ybar)/(X.sd)
    candidate<-sort(sample(x=1:n2,size = n_m))
    Y.k<-Y[candidate]
    Sn[k] <- sqrt(n_m)*(mean(Y.k)-Xbar)/(Y.sd)
  }
  
  
  T_n1 <- (sum(Un^2)+sum(Sn^2)-2*m_n)/(4*m_n)^0.5
  adp_var<-(2+(mean((X-Xbar)^4)-3*X.sd^4)/(2*n_m*X.sd^4)+(mean((Y-Ybar)^4)-3*Y.sd^4)/(2*n_m*Y.sd^4))/2
  T_n1 <- T_n1/sqrt(adp_var)
  
  # p1<-2*(1-pnorm(abs(T_n1)))
  p1<-1-pnorm(T_n1)
  return(list(T_n1=T_n1,p1=p1))
}
myseed<-55040
set.seed(myseed)

###compute Wn with no anticipation
test.our(dat = dat,m_n=190,n_m=20)


