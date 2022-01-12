rm(list=ls())
drop<-read.csv("C:/Users/e0265296/My Drive/Topic7_missing labels/2020May JRSSB codes/real data analysis/Thermal_test.csv")
myX<-as.matrix(drop[1:22,2:7])
myY<-as.matrix(drop[23:44,2:7])
mydat<-list(X=myX,Y=myY,Y.match=myY)
myres_ag<-rep(0,7)
myres_mul<-rep(0,7)



test.our <- function(dat,m_n)
{
  X <- dat$X
  Y <- dat$Y
  n <- length(X)
  m_n <- round(m_n)
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
n<-22
m_n_candidate<-2:11
length.candidate<-length(m_n_candidate)
###generate the data
dat<-mydat
my_X<-apply(dat$X,MARGIN = 1,sum)
my_Y<-apply(dat$Y,MARGIN = 1,sum)
# dat_aggre1<-list(X=my_X,Y=my_Y,Y.match=my_Ymatch)
dat_aggre2<-list(X=my_X,Y=my_Y)
set.seed(150)
dat_aggre3<-list(X=dat_aggre2$X[sample(1:22,22)],Y=dat_aggre2$Y)
test.our(dat_aggre3,m_n = 11)
test.our(dat_aggre3,m_n = 9)


test.our.multiple <- function(dat,m_n)
{
  X <- dat$X
  Y <- dat$Y
  n <- nrow(X)
  D <- ncol(X)
  
  m_n <- round(m_n)
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
  Y.bar = apply(Y,2,mean)
  X.cov = cov(X)
  X.cov.inv = solve(X.cov)
  Un <- rep(0,m_n)
  for (k in 1:m_n)
  {
    if(k==1){
      X.k<-X[1:cum_groupsize[1],]
    }else{
      X.k<-X[(cum_groupsize[k-1]+1):cum_groupsize[k],]
    }
    if(is.matrix(X.k)){
      Un[k] <- groupsize[k]*t(apply(X.k,2,mean)-Y.bar)%*%X.cov.inv%*%(apply(X.k,2,mean)-Y.bar)
    }else{
      Un[k] <- groupsize[k]*t(X.k-Y.bar)%*%X.cov.inv%*%(X.k-Y.bar)
    }
  }
  
  T_n1 <- (sum(Un)-m_n*D)/(2*m_n*D)^0.5
  
  #Split Y
  X.bar = apply(X,2,mean)
  Y.cov = cov(Y)
  Y.cov.inv <- solve(Y.cov)
  Sn <- rep(0,m_n)
  for (k in 1:m_n)
  {
    if(k==1){
      Y.k<-Y[1:cum_groupsize[1],]
    }else{
      Y.k<-Y[(cum_groupsize[k-1]+1):cum_groupsize[k],]
    }
    if(is.matrix(Y.k)){
      Sn[k] <- groupsize[k]*t(apply(Y.k,2,mean)-X.bar)%*%Y.cov.inv%*%(apply(Y.k,2,mean)-X.bar)
    }else{
      Sn[k] <- groupsize[k]*t(Y.k-X.bar)%*%Y.cov.inv%*%(Y.k-X.bar)
    }
    
  }
  T_n2 <- (sum(Sn)-m_n*D)/(2*m_n*D)^0.5
  T_n3 <- (T_n1+T_n2)/sqrt(2)
  p.value1 <- 1-pnorm(T_n3)
  
  return(result<-list(p.value1=p.value1,T_n3=T_n3))
}
dat<-mydat
set.seed(1234)
test.our.multiple(dat=list(X=dat$X[sample(1:22,22),],Y=dat$Y[sample(1:22,22),],Y.match=dat$Y.match),m_n = 11)


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
    candidate<-sort(sample(x=1:n,size = n_m))
    X.k<-X[candidate]
    Un[k] <- sqrt(n_m)*(mean(X.k)-Ybar)/(X.sd)
    candidate<-sort(sample(x=1:n,size = n_m))
    Y.k<-Y[candidate]
    Sn[k] <- sqrt(n_m)*(mean(Y.k)-Xbar)/(Y.sd)
  }
  
  
  # T_n1 <- (sum(Un)-sum(Sn))/(2*m_n)^0.5
  T_n1 <- (sum(Un^2)+sum(Sn^2)-2*m_n)/(4*m_n)^0.5
  adp_var<-(2+(mean((X-Xbar)^4)-3*X.sd^4)/(2*n_m*X.sd^4)+(mean((Y-Ybar)^4)-3*Y.sd^4)/(2*n_m*Y.sd^4))/2
  T_n1 <- T_n1/sqrt(adp_var)
  
  # p1<-2*(1-pnorm(abs(T_n1)))
  p1<-1-pnorm(T_n1)
  return(list(T_n1=T_n1,p.value1=p1))
}
dat<-mydat
my_X<-apply(dat$X,MARGIN = 1,sum)
my_Y<-apply(dat$Y,MARGIN = 1,sum)
dat_aggre2<-list(X=my_X,Y=my_Y)
set.seed(4481)
test.our(dat = dat_aggre2,m_n = 12,n_m = 2)
test.our(dat = dat_aggre2,m_n = 8,n_m = 3)



test.our.multiple <- function(dat,m_n=numeric(0),n_m=numeric(0))
{
  X <- dat$X
  Y <- dat$Y
  n <- nrow(X)
  D <- ncol(X)
  
  X.bar = apply(X,2,mean)
  Y.bar = apply(Y,2,mean)
  X.cov = cov(X)
  Y.cov = cov(Y)
  X.cov.inv = solve(X.cov)
  Y.cov.inv = solve(Y.cov)
  Un <- rep(0,m_n)
  Sn <- rep(0,m_n)
  for (k in 1:m_n)
  {
    candidate<-sort(sample(x=1:n,size = n_m))
    X.k<-X[candidate,]
    candidate<-sort(sample(x=1:n,size = n_m))
    Y.k<-Y[candidate,]
    if(n_m==1){
      Un[k] <- t(X.k-Y.bar)%*%X.cov.inv%*%(X.k-Y.bar)
      Sn[k] <- t(Y.k-X.bar)%*%Y.cov.inv%*%(Y.k-X.bar)
    }else{
      Un[k] <- n_m*t(apply(X.k,2,mean)-Y.bar)%*%X.cov.inv%*%(apply(X.k,2,mean)-Y.bar)
      Sn[k] <- n_m*t(apply(Y.k,2,mean)-X.bar)%*%Y.cov.inv%*%(apply(Y.k,2,mean)-X.bar)
    }
  }
  # T_n1 <- (sum(Un)-sum(Sn))/(2*m_n)^0.5
  T_n1 <- (sum(Un)+sum(Sn)-2*m_n*D)/(4*m_n*D)^0.5
  kappa_X<-0
  kappa_Y<-0
  for(ii in 1:n){
    kappa_X<-kappa_X+(t(X[ii,]-X.bar)%*%X.cov.inv%*%(X[ii,]-X.bar))^2
    kappa_Y<-kappa_Y+(t(Y[ii,]-Y.bar)%*%Y.cov.inv%*%(Y[ii,]-Y.bar))^2
  }
  kappa_X<-kappa_X/n
  kappa_Y<-kappa_Y/n
  adp_var<-(kappa_X+kappa_Y-4*D-2*D^2)/(4*n_m*D)+1
  T_n1 <- T_n1/sqrt(adp_var)
  # p1<-2*(1-pnorm(abs(T_n1)))
  p.value1<-1-pnorm(T_n1)
  return(list(p.value1=p.value1,T_n1=T_n1))
}
dat<-mydat
set.seed(4250)
test.our.multiple(dat = mydat,m_n = 12,n_m = 2)

test.naive.oracle.aggregate <- function(dat)
{
  X <- dat$X
  Y <- dat$Y
  Y.match <- dat$Y.match
  n <- nrow(X)
  
  eta <- apply(X,1,sum)
  zeta <- apply(Y,1,sum)
  zeta.match <- apply(Y.match,1,sum)
  
  #naive paiered test, pretending X and Y matches
  naive.dif = eta-zeta
  T_n.naive = n^0.5*mean(naive.dif)/sd(naive.dif)
  p.value.naive = 1-pnorm(T_n.naive)
  
  #naive two sample test, pretending independent X and Y
  eta.zeta.dif.mean = mean(eta-zeta)
  eta.zeta.dif.sd = (var(eta)+var(zeta))^0.5/n^0.5
  T_n.naive.two = eta.zeta.dif.mean/eta.zeta.dif.sd
  p.value.naive.two = 1-pnorm(T_n.naive.two)
  
  return(result<-list(p.value.naive=p.value.naive,
                      p.value.naive.two=p.value.naive.two))
}


test.naive.oracle.multiple <- function(dat)
{
  X <- dat$X
  Y <- dat$Y
  Y.match <- dat$Y.match
  n <- nrow(X)
  D <- ncol(X)
  
  #naive paiered test, pretending X and Y matches
  naive.dif = X-Y
  T_n.naive = n*t(apply(naive.dif,2,mean))%*%solve(cov(naive.dif))%*%apply(naive.dif,2,mean)
  p.value.naive = 1-pchisq(T_n.naive,D)
  
  #naive two sample test, pretending independent X and Y
  X.Y.dif.mean = apply(X-Y,2,mean)
  X.Y.dif.cov = ((n-1)*cov(X)+(n-1)*cov(Y))/(2*n-2)
  T_n.naive.two =  n/2*t(X.Y.dif.mean)%*%solve(X.Y.dif.cov)%*%X.Y.dif.mean
  T_n.naive.two <- T_n.naive.two*(2*n-D-1)/(D*(2*n-2))
  p.value.naive.two = 1-pf(T_n.naive.two,D,df1 = D, df2 = (2*n-D-1))
  
  return(result<-list(p.value.naive=p.value.naive,
                      p.value.naive.two=p.value.naive.two))
}
test.naive.oracle.aggregate(mydat)

test.naive.oracle.multiple(mydat)

