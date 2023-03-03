source('data_gene.R')
source('FNS_helper.R')
source('FNS_main.R')
################ define true regression function and parameters #########################
basisnum=8
acoe=rep(0,basisnum)
for(i in 1:basisnum)
{
  acoe[i]=i^{-1.5}
}
f=function(s,k){
  if(k==1) return(1)
  else if(k%%2==0) return(sqrt(2)*cos(k*pi*s))
  else return(sqrt(2)*sin((k-1)*pi*s))
}
m=function(s){Reduce("+",lapply(1:basisnum,function(i){acoe[i]*f(s,i)}))}

n=100000
domain=c(0,1)
snr=2
n0=1000
C=1/2
q0=5
ext=c(0,0.1,0)

################### main process ###################
set.seed(1)
data <- datagene(m=m,snr=snr,n=n,domain=domain)
Kmax <- 100
databatch <- datatobatch(data,Kmax)
rm(data)

m_est_ <- c()
for(K in 1:Kmax){
  
  x <- databatch[[K]][,1]
  y <- databatch[[K]][,2]
  m_est_ <- rbind(m_est_, onpreg(x,y,K,L=10,domain))
  if(K%%10==0)
  {
    plot(m_est_[K,], type='l', main=paste0('K=',K))
    lines(m(seq(0.05,0.95,length.out=50)), lty=2, col='gray')
  }
}
