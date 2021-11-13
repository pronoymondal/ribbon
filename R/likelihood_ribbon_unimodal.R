likelihood_unimodal<-function(data,data1,data2)
{
  l<-NULL
  
  for(j in 1:(dim(data)[2]))
  {
    n<-sum(data[,j]>0)
    mu<-mean(log(data[(data[,j]>0),j]))
    sig<-sd(log(data[(data[,j]>0),j]))*sqrt((n-1)/n)
    p<-mean(data[,j]>0)
    n1<-sum(data1[,j]>0)
    mu1<-mean(log(data1[(data1[,j]>0),j]))
    sig1<-sd(log(data1[(data1[,j]>0),j]))*sqrt((n1-1)/n1)
    p1<-mean(data1[,j]>0)
    n2<-sum(data2[,j]>0)
    mu2<-mean(log(data2[(data2[,j]>0),j]))
    sig2<-sd(log(data2[(data2[,j]>0),j]))*sqrt((n2-1)/n2)
    p2<-mean(data2[,j]>0)
    l[j]<-sum(dnorm(log(data1[(data1[,j]>0),j]),mu1,sig1,log=TRUE))+sum(dnorm(log(data2[(data2[,j]>0),j]),mu2,sig2,log=TRUE))-sum(dnorm(log(data[(data[,j]>0),j]),mu,sig,log=TRUE))
    
  }
  return(l)
  
}






