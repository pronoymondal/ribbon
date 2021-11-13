

ribbon_estimate<-function(a4)
{
data<-a4
library(mixtools)

data<-data[,(colSums(data>0)>5)]

if(is.null(colnames(data)))
{
  colnames(data)<-1:(dim(data)[2])
}

if(is.null(rownames(data)))
{
  rownames(data)<-1:(dim(data)[1])
}

genes<-colnames(data)

#source("check_bimodality.R")
#indd<-check_bimodality(t(data))

lst<-check_bimodality_robertson(data)

indd<-lst$indd
mu1<-lst$mu1
mu2<-lst$mu2
sig1<-lst$sig1
sig2<-lst$sig2
pi<-lst$pi
D<-lst$D

mu1<-mu1[indd==1]
mu2<-mu2[indd==1]
sig1<-sig1[indd==1]
sig2<-sig2[indd==1]
pi<-pi[indd==1]
D<-D[,(indd==1)]

data1<-data[,(indd==1)]
data2<-data[,(indd==0)]

#source("estimation_model22.R")
xxx1<-estimate_bimodal(data1,mu1,mu2,sig1,sig2,pi,D)

#source("estimation_model21.R")
xxx2<-estimate_unimodal(data2)

lst<-list(indd=indd,xxx1=xxx1,xxx2=xxx2)

return(lst)
}
















