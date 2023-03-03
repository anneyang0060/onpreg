#' function to estimate regression function with one-pass local smoothing method OPLS method) 
#' @param x,y available sample at present 
#' @param K current timestamp
#' @param L length of candidate bandwidth sequence
#' @param domain the domain of x
#' @return a list containing RMISE, time of maintaining and querying
onpreg <- function(x,y,K,L,domain)
{
  # common parameters
  {
    a <- domain[1]+0.05; b <- domain[2]-0.05
    EV <- 50; 
    eval <- seq(a,b,length.out = EV)
    G <- 0.7
    if(Kmax>0){K1 <- 300}else{K1 <- 1000}
  }
  
  # initialize
  if(K==1)
  {
    N <<- 0
    h0 <<- 0.7
    theta <<- 0; sigma_est <<- 0
    res_theta <<- list(); res_sigma <<- list(); res <<- list()
    res_theta$centroids <<- rep(0, L)
    res_theta$P <<- array(0, dim = c(4,4,EV,L))
    res_theta$q <<- array(0, dim = c(4,EV,L))
    res_sigma$centroids <<- rep(0, L)
    res_sigma$P <<- array(0, dim = c(2,2,EV,L))
    res_sigma$q <<- array(0, dim = c(2,EV,L))
    res$centroids <<- rep(0, L)
    res$P <<- array(0, dim = c(2,2,EV,L))
    res$q <<- array(0, dim = c(2,EV,L))
  }
  # generate data
  { 
    NK <- length(y)
    N <<- N + length(y)
  }
  # estimate bandwidth
  if(K<=K1)
  {
    h_theta <- G * N^(-1/7)
    res_theta <<- online_LCub(x, y, eval, h_theta, L,res_theta, N, NK, 1,K)
    m_sec_deri <- sapply(1:EV, function(i){2*(solve(res_theta$P[,,i,1]+diag(1e-12,4)) %*% matrix(res_theta$q[,i,1],4,1))[3]})
    theta <<-  sum(res_theta$P[1,1,,1] * m_sec_deri^2)/sum(res_theta$P[1,1, ,1])
    h_sigma <- G * N^(-1/5)
    res_sigma <<- online_LL(x, y,eval, h_sigma, L, res_sigma, N, NK, 1,K)
    m_est <- sapply(1:EV, function(i){(solve(res_sigma$P[,,i,1]+diag(1e-12,2)) %*% matrix(res_sigma$q[,i,1],2,1))[1]})
    sp <- smooth.spline(eval,m_est)
    m_est<-predict(sp,x)$y
    sigma_est <<- (N-NK)/N * sigma_est + NK/N * mean((y-m_est)^2)
  }
  
  # estimate
  {
    h <- min((15 * sigma_est / theta)^(1/5) * N^(-1/5),h0)
    h0 <<- h
    print(h)
    res <<- online_LL(x, y,eval, h, L, res, N, NK, 1,K)
    m_est <- sapply(1:EV, function(i){(solve(res$P[,,i,1]+diag(1e-12,2)) %*% matrix(res$q[,i,1],2,1))[1]})
  }
  
  return(m_est)
}
