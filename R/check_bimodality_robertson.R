
check_bimodality_robertson<-function(data)
{
	a<-data
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

d<-dim(data)
D<-array(0,dim<-dim(data))
mu1<-NULL
mu2<-NULL
sig1<-NULL
sig2<-NULL
pi<-NULL

for(j in 1:(d[2]))
{
  if(sum(data[,j]>0)<0.0001)
  {
    if(sum(data[,j]>0)>0)
    {
      mu1[j]<-mean(log(data[(data[,j]>0),j]))
      sig1[j]<-sd(log(data[(data[,j]>0),j]))
      mu2[j]<-mean(log(data[(data[,j]>0),j]))
      sig2[j]<-sd(log(data[(data[,j]>0),j]))
      pi[j]<-0.5
      for(i in 1:(d[1]))
      {
        D[i,j]<-0.5
      }
    }else
    {
      mu1[j]<-0
      sig1[j]<-1
      mu2[j]<-0
      sig2[j]<-1
      pi[j]<-0.5
      for(i in 1:(d[1]))
      {
        D[i,j]<-0.5
      }
    }
  }else
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
    if(is.na(mu1[j]+mu2[j]+sig1[j]+sig2[j]))
    {
      yy<-log(data[(data[,j]>0),j])
      Me<-median(yy)
      mu1[j]<-mean(yy[yy<=Me])
      mu2[j]<-mean(yy[yy>Me])
      sig1[j]<-sd(yy[yy<=Me])
      sig2[j]<-sd(yy[yy>Me])
      pi[j]<-0.5
    }
    if((sig1[j])==0)
    {
      yy<-log(data[(data[,j]>0),j])
      Me<-median(yy)
      mu1[j]<-mean(yy[yy<=Me])
      mu2[j]<-mean(yy[yy>Me])
      sig1[j]<-sd(yy[yy<=Me])
      sig2[j]<-sd(yy[yy>Me])
      pi[j]<-0.5
    }
    if(mu1[j]<(-10))
    {

    }
    if(mu2[j]<(-10))
    {

    }
    indices<-indd
    x<-data[data[,j]>0,j]
    for(i in 1:(d[1]))
    {
      if((data[i,j])>0)
      {
        D[i,j]<-(pi[j]*dnorm(log(data[i,j]),mu2[j],sig2[j]))/(pi[j]*dnorm(log(data[i,j]),mu2[j],sig2[j])+(1-pi[j])*dnorm(log(data[i,j]),mu1[j],sig1[j]))
      }
    }
  }
  print(j)
}

prop<-pi
d<-abs(mu1-mu2)/(2*sqrt(sig1*sig2))
stat1<-abs(log(1-prop)-log(prop))
stat2<-2*log(d-sqrt(d^2-1))+2*d*sqrt(d^2-1)
indd<-rep(NA,length(d))

d<-abs(mu1-mu2)/(2*sqrt(sig1*sig2))

stat1<-abs(log(prop)-log(1-prop))
stat2<-2*log(d-sqrt(d^2-1))+2*d*sqrt(d^2-1)

sigg<-sig2/sig1
x1<-abs(mu1-mu2)
x2<-sig1*((2*(sigg^4-sigg^2+1)^(3/2)-(2*sigg^6-3*sigg^4-3*sigg^2+2))/(sigg^2))^(1/2)

muu1<-NULL
muu2<-NULL
sigg1<-NULL
sigg2<-NULL
muu0<-NULL
sigg<-NULL
muu<-NULL
root1<-NULL
root2<-NULL
ci1<-NULL
ci2<-NULL
index<-NULL
y1<-NULL
y2<-NULL
pp<-NULL

library(RConics)
for(i in 1:(length(mu1)))
{
  if(is.na(mu1[i]))
  {
    next
  }
  if(mu1[i]<mu2[i])
  {
    muu1[i]<-mu1[i]
    muu2[i]<-mu2[i]
    sigg1[i]<-sig1[i]
    sigg2[i]<-sig2[i]
    pp[i]<-prop[i]
  }else
  {
    muu1[i]<-mu2[i]
    muu2[i]<-mu1[i]
    sigg1[i]<-sig2[i]
    sigg2[i]<-sig1[i]
    pp[i]<-(1-prop[i])
  }
  muu[i]<-(muu2[i]-muu1[i])/sigg1[i]
  sigg[i]<-sigg2[i]/sigg1[i]
  muu0[i]<-sqrt((2*(sigg[i]^4-sigg[i]^2+1)^(3/2)-(2*sigg[i]^6-3*sigg[i]^4-3*sigg[i]^2+2))/(sigg[i]^2))
  if(muu[i]<muu0[i])
  {
    index[i]<-TRUE
  }else
  {
    index[i]<-FALSE
    sss<-cubic(c((sigg[i]^2-1),-(sigg[i]^2-2)*muu[i],-muu[i]^2,muu[i]*sigg[i]^2))
    y1[i]<-min(sss[(sss>0)&(sss<muu[i])])
    y2[i]<-max(sss[(sss>0)&(sss<muu[i])])
    ci1[i]<-1/(1+(sigg[i]^3*y1[i])/(muu[i]-y1[i])*exp(-0.5*y1[i]^2+0.5*(y1[i]-muu[i])^2/(sigg[i]^2)))
    ci2[i]<-1/(1+(sigg[i]^3*y2[i])/(muu[i]-y2[i])*exp(-0.5*y2[i]^2+0.5*(y2[i]-muu[i])^2/(sigg[i]^2)))
  }

}

indd[((Re(ci1)<pp)&(Re(ci2)>pp)&(x1>x2))]<-1
indd[is.na(indd)]<-0

lst<-list(indd=indd,mu1=mu1,mu2=mu2,sig1=sig1,sig2=sig2,pi=pi,D=D)

return(lst)

}
