check_bimodality<-function(a)
{
  library(mixtools)
  mu1<-NULL
  mu2<-NULL
  sig1<-NULL
  sig2<-NULL
  sig2<-NULL
  prop<-NULL
  for(k in 1:(dim(a)[1]))
  {
    if(sum(a[k,]>0)<3)
    {
      next
    }
    x<-as.numeric(a[k,])
    tryCatch({
    y<-log(x[x>0])
  fit<-normalmixEM(y)
  mu1[k]<-fit$mu[1]
  mu2[k]<-fit$mu[2]
  sig1[k]<-fit$sigma[1]
  sig2[k]<-fit$sigma[2]
  prop[k]<-fit$lambda[1]
  print(k)
    },
  error = function(e){ }
  )
  }
  d<-abs(mu1-mu2)/(2*sqrt(sig1*sig2))
  stat1<-abs(log(1-prop)-log(prop))
  stat2<-2*log(d-sqrt(d^2-1))+2*d*sqrt(d^2-1)
  indd<-rep(NA,length(d))
  indd[d<1]<-0
  indd[(d>1)&(stat1>stat2)]<-0
  indd[(d>1)&(stat1<stat2)]<-1
  indd[is.na(indd)] <- 0
  return(indd)
}







