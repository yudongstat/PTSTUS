#' This function implements the sample-splitting-based test statistic
#'
#' @param dat a data frame with two variables: X and Y, which contain the unordered samples before and after the treatment, respectively.
#' @param m_n the number of blocks
#'
#' @return A list of length 2：
#          test_stat is the computed test statistic;
#          p_value is the computed p-value.
#' @export
#'
#' @examples
sample_split <- function(dat,m_n)
{
  X <- dat$X
  Y <- dat$Y
  n <- nrow(X)
  D <- ncol(X)
  m_n <- round(m_n)
  groupsize<-rep(n%/%m_n,m_n)
  temp<-n%%m_n
  if(temp>0){
    groupsize[1:temp]<-groupsize[1:temp]+1
  }
  cum_groupsize<-Reduce("sum", groupsize, accumulate = TRUE)
  cum_groupsize<-c(0,cum_groupsize)
  #Split X
  X.bar = apply(X,2,mean)
  Y.bar = apply(Y,2,mean)
  X.cov = cov(X)
  Y.cov = cov(Y)
  X.cov.inv = solve(X.cov)
  Y.cov.inv <- solve(Y.cov)
  Un <- rep(0,m_n)
  Sn <- rep(0,m_n)
  for (k in 1:m_n)
  {
    X.k<-X[(cum_groupsize[k]+1):cum_groupsize[k+1],]
    Y.k<-Y[(cum_groupsize[k]+1):cum_groupsize[k+1],]
    if(is.matrix(X.k)){
      temp1<-apply(X.k,2,mean)-Y.bar
      temp2<-apply(Y.k,2,mean)-X.bar
      Un[k] <- groupsize[k]*t(temp1)%*%X.cov.inv%*%temp1
      Sn[k] <- groupsize[k]*t(temp2)%*%Y.cov.inv%*%temp2
    }else{
      temp1<-X.k-Y.bar
      temp2<-Y.k-X.bar
      Un[k] <- groupsize[k]*t(temp1)%*%X.cov.inv%*%temp1
      Sn[k] <- groupsize[k]*t(temp2)%*%Y.cov.inv%*%temp2
    }
  }

  T_n1 <- (sum(Un)-m_n*D)/(2*m_n*D)^0.5
  T_n2 <- (sum(Sn)-m_n*D)/(2*m_n*D)^0.5

  temp1<-X-rep(1,n)%*%t(X.bar)
  temp2<-Y-rep(1,n)%*%t(Y.bar)
  term1<-mean((diag(temp1%*%X.cov.inv%*%t(temp1)))^2)
  term2<-mean((diag(temp2%*%Y.cov.inv%*%t(temp2)))^2)

  adp_var<-1+m_n*(term1+term2-4*D-2*D^2)/(4*n*D)+2/m_n+tr(Y.cov%*%X.cov.inv%*%Y.cov%*%X.cov.inv)/(2*m_n*D)+tr(X.cov%*%Y.cov.inv%*%X.cov%*%Y.cov.inv)/(2*m_n*D)+2*tr(X.cov%*%Y.cov.inv)/(m_n*D)+2*tr(Y.cov%*%X.cov.inv)/(m_n*D)#-8/(m_n*D)#-m_n*D/(2*n)#(2+m_n*(term2-3*Y.sd^4)/(2*n*Y.sd^4)+(X.sd^4+2*X.sd^2*Y.sd^2)/(m_n*Y.sd^4))/2+(X.sd^4+Y.sd^4+2*X.sd^2*Y.sd^2)/(m_n*X.sd^2*Y.sd^2)-0.5*(X.sd^2+Y.sd^2)^2/(m_n*X.sd^4)-0.5*(X.sd^2+Y.sd^2)^2/(m_n*Y.sd^4)-(3/4)*(X.sd^2+Y.sd^2)^2/(m_n*X.sd^2*Y.sd^2)-m_n/(2*n)
  if(adp_var<0){
    stop()
  }
  T_n3 <- (T_n1+T_n2)/sqrt(2)
  T_n3<-T_n3/sqrt(adp_var)
  p.value1 <- 1-pnorm(T_n3)

  return(result<-list(test_stat=T_n3,p_value=p.value1))
}
