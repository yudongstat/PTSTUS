#' This function implements the subsampling-based test statistic W_n
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
subsampling_W_n <- function(dat,m_n=numeric(0),n_m=numeric(0))
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
    Un[k] <- n_m*t(apply(X.k,2,mean)-Y.bar)%*%X.cov.inv%*%(apply(X.k,2,mean)-Y.bar)
    candidate<-sort(sample(x=1:n,size = n_m))
    Y.k<-Y[candidate,]
    Sn[k] <- n_m*t(apply(Y.k,2,mean)-X.bar)%*%Y.cov.inv%*%(apply(Y.k,2,mean)-X.bar)
  }


  # T_n1 <- (sum(Un)-sum(Sn))/(2*m_n)^0.5
  T_n1 <- (sum(Un)+sum(Sn)-2*m_n*D)/(4*m_n*D)^0.5
  T_n1 <- T_n1-0.5*sqrt(m_n/D)*(n/(n-1)*n_m*(D+tr(X.cov.inv%*%Y.cov))/n)-0.5*sqrt(m_n/D)*(n/(n-1)*n_m*(D+tr(Y.cov.inv%*%X.cov))/n)
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
  return(list(test_stat=T_n1,p_value=p.value1))
}
