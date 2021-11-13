likelihood_bimodal<-function(data,data1,data2)
{
  mu1<-rep(NA,dim(data)[2])
  mu2<-rep(NA,dim(data)[2])
  sig1<-rep(NA,dim(data)[2])
  sig2<-rep(NA,dim(data)[2])
  pi<-rep(NA,dim(data)[2])
  
  mu11<-rep(NA,dim(data)[2])
  mu12<-rep(NA,dim(data)[2])
  mu21<-rep(NA,dim(data)[2])
  mu22<-rep(NA,dim(data)[2])
  sig11<-rep(NA,dim(data)[2])
  sig12<-rep(NA,dim(data)[2])
  sig21<-rep(NA,dim(data)[2])
  sig22<-rep(NA,dim(data)[2])
  pi1<-rep(NA,dim(data)[2])
  pi2<-rep(NA,dim(data)[2])
  
  l<-NULL
  library(mixtools)
  for(j in 1:(dim(data)[2]))
  {
    tryCatch({
      a<-normalmixEM(log(data[(data[,j]>0),j]))
      if(a$mu[1]<a$mu[2])
      {
        s1<-a$sigma[1]
        s2<-a$sigma[2]
        p<-a$lambda[2]
        indd<-a$posterior[,1]
      }else
      {
        s1<-a$sigma[2]
        s2<-a$sigma[1]
        p<-a$lambda[1]
        indd<-a$posterior[,2]
      }
      mu1[j]<-min(a$mu)
      sig1[j]<-s1
      mu2[j]<-max(a$mu)
      sig2[j]<-s2
      pi[j]<-p
    },
    error = function(e){ }
    )
    if(is.na(mu1[j]))
    {
      yy<-log(data[(data[,j]>0),j])
      Me<-median(yy)
      mu1[j]<-mean(yy[yy<=Me])
      mu2[j]<-mean(yy[yy>Me])
      sig1[j]<-sd(yy[yy<=Me])
      sig2[j]<-sd(yy[yy>=Me])
      pi[j]<-0.5
    }
    if((sig1[j])==0)
    {
      yy<-log(data[(data[,j]>0),j])
      Me<-median(yy)
      mu1[j]<-mean(yy[yy<=Me])
      mu2[j]<-mean(yy[yy>Me])
      sig1[j]<-sd(yy[yy<=Me])
      sig2[j]<-sd(yy[yy>=Me])
      pi[j]<-0.5
    }
    
  }

  mu21<-mu2
  mu22<-mu2
  sig21<-sig2
  sig22<-sig2
  mu11<-mu1
  mu12<-mu1
  sig11<-sig1
  sig12<-sig1
  pi1<-pi
  pi2<-pi
  
  for(j in 1:(dim(data)[2]))
  {
    l[j]<-0
    
    ind<-(pi[j]*dnorm(log(data[(data[,j]>0),j]),mu2[j],sig2[j]))/((pi[j]*dnorm(log(data[(data[,j]>0),j]),mu2[j],sig2[j]))+((1-pi[j])*dnorm(log(data[(data[,j]>0),j]),mu1[j],sig1[j])))
    ind1<-(pi1[j]*dnorm(log(data1[(data1[,j]>0),j]),mu21[j],sig21[j]))/((pi1[j]*dnorm(log(data1[(data1[,j]>0),j]),mu21[j],sig21[j]))+((1-pi1[j])*dnorm(log(data1[(data1[,j]>0),j]),mu11[j],sig11[j])))
    ind2<-(pi2[j]*dnorm(log(data2[(data2[,j]>0),j]),mu22[j],sig22[j]))/((pi2[j]*dnorm(log(data2[(data2[,j]>0),j]),mu22[j],sig22[j]))+((1-pi2[j])*dnorm(log(data2[(data2[,j]>0),j]),mu12[j],sig12[j])))
    
    p1<-mean(data1[,j]>0)
    p2<-mean(data2[,j]>0)
    p<-mean(data[,j]>0)
    
    pi[j]<-mean(ind)
    pi1[j]<-mean(ind1)
    pi2[j]<-mean(ind2)
    
    l[j]<-0
    l[j]<-l[j]+sum(ind1)*log(pi1[j])+sum(1-ind1)*log(1-pi1[j])+sum(ind2)*log(pi2[j])+sum(1-ind2)*log(1-pi2[j])-sum(ind)*log(pi[j])-sum(1-ind)*log(1-pi[j])
    
  }
  
 
   
  return(l)
  
}






