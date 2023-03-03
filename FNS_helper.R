Epan <- function(z)
{ 
  return( 3/4 * (1-z^2) * (abs(z)<1) ) 
} 

online_LCub <- function(x, y, eval, h, L, res_list, N, n, d,K)
{
  eta <- sapply(1:L, function(l){ ((L-l+1) / L) ^ (1/(6+d)) * h}) 
  if(K>1){
    idx <- sapply(1:L,function(l){which.min(abs(eta[l] - res_list$centroids))})}
  else{idx <- 1:L}
  res_list$centroids <- (res_list$centroids[idx] * (N-n) + eta * n) / N
  { 
    EV <- length(eval)
    for(l in 1:L){
      Pnew <- array(0, dim = c(4,4,EV)); qnew <- matrix(0,4,EV)
      for(i in 1:EV){
        side <- cbind(1, x - eval[i], (x - eval[i])^2, (x - eval[i])^3)
        K_vec <- Epan((x - eval[i])/eta[l])/eta[l]
        for(nr in 1:4){
          for(nc in 1:4){
            Pnew[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/n }
          qnew[nr,i] <-  sum(K_vec*side[,nr]*y)/n
        }
      }
      res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
      res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
    }
  }
  return(res_list)
}


online_LL <- function(x, y, eval, h, L, res_list, N, n, d,K)
{
  eta <- sapply(1:L, function(l){ ((L-l+1) / L) ^ (1/(d+4)) * h})

  if(K>1){
    idx <- sapply(1:L,function(l){which.min(abs(eta[l] - res_list$centroids))})}
  else{
    idx <- 1:L
  }
  res_list$centroids <- (res_list$centroids[idx] * (N-n) + eta * n) / N

  {
    EV <- length(eval) 
    for(l in 1:L){
      Pnew <- array(0, dim = c(2,2,EV)); qnew <- matrix(0,2,EV)
      for(i in 1:EV){
        side <- cbind(1, x - eval[i])
        K_vec <- Epan((x - eval[i])/eta[l])/eta[l]
        Pnew[,,i] <- matrix(c(
          sum(K_vec*side[,1]^2), sum(K_vec*side[,1]*side[,2]),
          sum(K_vec*side[,2]*side[,1]),  sum(K_vec*side[,2]^2)
        ),2,2) / n
        qnew[,i] <- matrix(c(
          sum(K_vec*side[,1]*y), sum(K_vec*side[,2]*y)
        ),2,1) / n
      }
      res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
      res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
    }
  }
  
  return(res_list)
}



