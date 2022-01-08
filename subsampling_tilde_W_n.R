#' This function implements the subsampling-based test statistic tilde W_n
#'
#' @param dat a data frame with two variables: X and Y, which contain the unordered samples before and after the treatment, respectively.
#' @param m_n the number of blocks
#' @param n_m the block size
#'
#' @return A list of length 2：
#          test_stat is the computed test statistic;
#          p_value is the computed p-value.
#' @export
#'
#' @examples
subsampling_tilde_W_n <- function(dat,m_n=numeric(0),n_m=numeric(0))
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
  X.cov.inv.sqrt = sqrtm(X.cov.inv)
  Y.cov.inv.sqrt = sqrtm(Y.cov.inv)
  Un <- matrix(0,D,m_n)
  Sn <- matrix(0,D,m_n)
  if(n_m==1){
    for (k in 1:m_n)
    {
      candidate<-sort(sample(x=1:n,size = n_m))
      X.k<-X[candidate,]
      Un[,k] <- sqrt(n_m)*X.cov.inv.sqrt%*%(X.k-Y.bar)
      candidate<-sort(sample(x=1:n,size = n_m))
      Y.k<-Y[candidate,]
      Sn[,k] <- sqrt(n_m)*Y.cov.inv.sqrt%*%(Y.k-X.bar)
    }
  }else{
    for (k in 1:m_n)
    {
      candidate<-sort(sample(x=1:n,size = n_m))
      X.k<-X[candidate,]
      Un[,k] <- sqrt(n_m)*X.cov.inv.sqrt%*%(apply(X.k,2,mean)-Y.bar)
      candidate<-sort(sample(x=1:n,size = n_m))
      Y.k<-Y[candidate,]
      Sn[,k] <- sqrt(n_m)*Y.cov.inv.sqrt%*%(apply(Y.k,2,mean)-X.bar)
    }
  }



  T_n1 <- (apply(Un,MARGIN = 1,sum)-apply(Sn,MARGIN = 1,sum))/(2*m_n)^0.5
  T_n0 <- T_n1
  adp_var<-diag(rep(1,D))*(n-n_m)/(n-1)+m_n*n_m/(2*n)*(X.cov.inv.sqrt+Y.cov.inv.sqrt)%*%(X.cov+Y.cov)%*%t(X.cov.inv.sqrt+Y.cov.inv.sqrt)
  T_n1 <- sqrtm(solve(adp_var))%*%T_n1

  T_n1 <- sum(T_n1^2)
  p.value1<-1-pchisq(T_n1,df=D)
  return(list(test_stat=T_n0,p_value=p.value1))
}
