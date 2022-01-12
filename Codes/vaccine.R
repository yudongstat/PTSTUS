rm(list=ls())
vaxdata<-read.csv(file = "C:/Users/e0265296/Desktop/VAERS_data.csv",header = TRUE)
library("mvtnorm")
vaxdata$VAX_DATE<-as.Date(vaxdata$VAX_DATE,"%m/%d/%y")
vaxdata$RECVDATE<-as.Date(vaxdata$RECVDATE,"%m/%d/%y")
start_age<-c(0,16,55,65,75)
end_age<-c(15,54,64,74,Inf)
start_age<-rep(start_age,2)
end_age<-rep(end_age,2)
SEX_seq<-c(rep("F",5),rep("M",5))
PFIZER_mat<-matrix(0,10,15)
PFIZER_list<-as.list(rep(0,10))
tuning_mat<-matrix(c(476,139,2293,1883,1269,2026,1687,1414,1001,418,650,243,877,909,315,275,270,253,369,472),10,2,byrow = TRUE)
for(ttt in 1:10){
  PFIZER_mat[ttt,1]<-start_age[ttt]
  aaa<-(vaxdata$VAX_MANU%in%c("PFIZER\\BIONTECH"))&(vaxdata$VAX_DOSE_SERIES=="1")&(vaxdata$AGE_YRS>=start_age[ttt])&(vaxdata$AGE_YRS<=end_age[ttt])&(vaxdata$SEX==SEX_seq[ttt])&(vaxdata$RECVDATE<=as.Date("2020-09-16"))
  bbb<-(vaxdata$VAX_MANU%in%c("PFIZER\\BIONTECH"))&(vaxdata$VAX_DOSE_SERIES=="2")&(vaxdata$AGE_YRS>=start_age[ttt])&(vaxdata$AGE_YRS<=end_age[ttt])&(vaxdata$SEX==SEX_seq[ttt])&(vaxdata$RECVDATE<=as.Date("2020-09-16"))
  aaa<-which(aaa)
  bbb<-which(bbb)
  X_sample<-vaxdata$NUM_SYPTMS[aaa]
  Y_sample<-vaxdata$NUM_SYPTMS[bbb]
  PFIZER_mat[ttt,2]<-length(X_sample)
  PFIZER_mat[ttt,3]<-length(Y_sample)
  PFIZER_mat[ttt,4]<-mean(X_sample)-mean(Y_sample)
  if((as.numeric(PFIZER_mat[ttt,2])>20)&(as.numeric(PFIZER_mat[ttt,3])>20)){
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
    dat<-list(X=X_sample,Y=Y_sample)
    dat<-list(X=X_sample[sample(1:length(X_sample),length(X_sample))],Y=Y_sample[sample(1:length(Y_sample),length(Y_sample))])
    PFIZER_mat[ttt,8]<-t.test(X_sample,Y_sample)$p.value
    PFIZER_mat[ttt,5]<-t.test(X_sample,Y_sample)$statistic
    
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
    myX<-X_sample
    myY<-Y_sample
    n<-length(myX)
    n1<-length(myX)
    n2<-length(myY)
    ###compute the proposed test stat
    my_seed<-222
    set.seed(my_seed)
    dat<-list(X=X_sample[sample(1:length(X_sample),length(X_sample))],Y=Y_sample[sample(1:length(Y_sample),length(Y_sample))])
    result <- test.our(dat=dat,m_n1=tuning_mat[ttt,1],m_n2=tuning_mat[ttt,2])
    PFIZER_mat[ttt,6]<-result$T_n3
    PFIZER_mat[ttt,9]<-result$p.value1
    myres1<-0
    for(qqq in 1:200){
      dat<-list(X=X_sample[sample(1:length(X_sample),length(X_sample))],Y=Y_sample[sample(1:length(Y_sample),length(Y_sample))])
      result <- test.our(dat=dat,m_n1=tuning_mat[ttt,1],m_n2=tuning_mat[ttt,2])
      myres1<-myres1+as.numeric(result$p.value1<0.05)
    }
    PFIZER_mat[ttt,14]<-myres1/200
    print(c(ttt,PFIZER_mat[ttt,]))
  }
}

tuning_mat<-matrix(c(22,65,70,546,10,1058,5,1046,10,358,16,135,25,701,10,295,20,169,10,184),10,2,byrow = TRUE)
for(ttt in 1:10){
  PFIZER_mat[ttt,1]<-start_age[ttt]
  aaa<-(vaxdata$VAX_MANU%in%c("PFIZER\\BIONTECH"))&(vaxdata$VAX_DOSE_SERIES=="1")&(vaxdata$AGE_YRS>=start_age[ttt])&(vaxdata$AGE_YRS<=end_age[ttt])&(vaxdata$SEX==SEX_seq[ttt])&(vaxdata$RECVDATE<=as.Date("2020-09-16"))
  bbb<-(vaxdata$VAX_MANU%in%c("PFIZER\\BIONTECH"))&(vaxdata$VAX_DOSE_SERIES=="2")&(vaxdata$AGE_YRS>=start_age[ttt])&(vaxdata$AGE_YRS<=end_age[ttt])&(vaxdata$SEX==SEX_seq[ttt])&(vaxdata$RECVDATE<=as.Date("2020-09-16"))
  aaa<-which(aaa)
  bbb<-which(bbb)
  X_sample<-vaxdata$NUM_SYPTMS[aaa]
  Y_sample<-vaxdata$NUM_SYPTMS[bbb]
  if((as.numeric(PFIZER_mat[ttt,2])>20)&(as.numeric(PFIZER_mat[ttt,3])>20)){
    library("mvtnorm")
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
      p1<-1-pnorm(T_n1)
      return(list(T_n1=T_n1,p1=p1))
    }
    myX<-X_sample
    myY<-Y_sample
    n<-length(myX)
    n1<-length(myX)
    n2<-length(myY)
    
    ###compute the proposed test stat
    my_seed<-777
    set.seed(my_seed)
    dat<-list(X=X_sample[sample(1:length(X_sample),length(X_sample))],Y=Y_sample[sample(1:length(Y_sample),length(Y_sample))])
    result <- test.our(dat=dat,m_n = tuning_mat[ttt,2],n_m=tuning_mat[ttt,1])
    PFIZER_mat[ttt,10]<-result$p1
    PFIZER_mat[ttt,7]<-result$T_n1
    myres1<-0
    for(qqq in 1:200){
      dat<-list(X=X_sample[sample(1:length(X_sample),length(X_sample))],Y=Y_sample[sample(1:length(Y_sample),length(Y_sample))])
      result <- test.our(dat=dat,m_n = tuning_mat[ttt,2],n_m=tuning_mat[ttt,1])
      myres1<-myres1+as.numeric(result$p1<0.05) 
      
    }
    PFIZER_mat[ttt,15]<-myres1/200
    print(c(ttt,PFIZER_mat[ttt,]))
  }
}


output_mat<-PFIZER_mat[,2:10]
output_mat<-round(output_mat,digits = 4)
library("xtable")
xtable(output_mat,digits = 4)
old_output_mat<-output_mat
