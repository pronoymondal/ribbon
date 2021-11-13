

ribbon_simulate<-function(a4,nCells=-100)
{
if(nCells==(-100))
{
  nCells<-dim(a4)[1]
}

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

data12<-data1
data22<-data2

#matt4<-array(0,dim<-dim(data12))
#matt5<-array(0,dim<-dim(data22))
matt4<-array(0,dim<-c(nCells,(dim(data12)[2])))
matt5<-array(0,dim<-c(nCells,(dim(data22)[2])))



matt<-array(0,dim<-dim(matt5))

alpha<-rep_len(sample(xxx2$alpha),nCells)
c<-xxx2$c
mu<-xxx2$mu
sig<-xxx2$sig
n<-length(alpha)

for(i in 1:n)
{
  for(j in 1:(length(mu)))
  {
    matt[i,j]<-rnorm(1,(mu[j]+alpha[i]),sig[j])
  }
}
E<-array(0,dim<-c(n,length(mu)))
for(i in 1:n)
{
  for(j in 1:(length(mu)))
  {
    E[i,j]<-rbinom(1,1,(pnorm(c[j])))
  }
}
matt5<-E*exp(matt)



matt<-array(0,dim<-dim(matt4))

alpha1<-rep_len(sample(xxx1$alpha1),nCells)
alpha2<-rep_len(sample(xxx1$alpha2),nCells)
c<-xxx1$c
mu1<-xxx1$mu1
sig1<-xxx1$sig1
mu2<-xxx1$mu2
sig2<-xxx1$sig2
pi<-xxx1$pi
n<-length(alpha1)

for(i in 1:n)
{
  for(j in 1:(length(mu1)))
  {
    u<-runif(1)
    if(u<pnorm(c[j]))
    {
      u1<-runif(1)
      if(u1<(pi[j]))
      {
        matt[i,j]<-rnorm(1,(mu2[j]),sig2[j])
      }else
      {
        u2<-runif(1)
        if(1)
        {
          matt[i,j]<-rnorm(1,(mu1[j]+alpha2[i]),sig1[j])
        }

      }
    }else
    {
      matt[i,j]<-0
    }

  }
}
E<-(matt!=0)
matt4<-exp(matt)*E

simulated3<-matt4
simulated4<-matt5




rpc_matrix3<-data12
rpc_matrix4<-simulated3
#source("ks_statistic.R")

ks_stat<-NULL
ks_stat1<-NULL
ks_stat2<-NULL
innddex<-1
for(j in 1:(dim(rpc_matrix3)[2]))
{
  x1<-rpc_matrix3[,j]
  x2<-rpc_matrix4[,j]
  if(!is.na(sum(as.numeric(x2))))
  {
    ks_stat[innddex]<-ks_statistic(x1,x2)
    ks_stat1[innddex]<-ks_statistic(x1[x1>0],x2[x2>0])*sqrt(length(x1[x1>0])*length(x2[x2>0])/(length(x1[x1>0])+length(x2[x2>0])))
    ks_stat2[innddex]<-abs(length(x1[x1==0])-length(x2[x2==0]))/length(x1)
  }
  innddex<-innddex+1
  print(innddex)
}

ks_stat_model12<-ks_stat
ks_stat1_model12<-ks_stat1
ks_stat2_model12<-ks_stat2



rpc_matrix3<-data22
rpc_matrix4<-simulated4
#source("ks_statistic.R")

ks_stat<-NULL
ks_stat1<-NULL
ks_stat2<-NULL
innddex<-1
for(j in 1:(dim(rpc_matrix3)[2]))
{
  x1<-rpc_matrix3[,j]
  x2<-rpc_matrix4[,j]
  if(!is.na(sum(as.numeric(x2))))
  {
    ks_stat[innddex]<-ks_statistic(x1,x2)
    ks_stat1[innddex]<-ks_statistic(x1[x1>0],x2[x2>0])*sqrt(length(x1[x1>0])*length(x2[x2>0])/(length(x1[x1>0])+length(x2[x2>0])))
    ks_stat2[innddex]<-abs(length(x1[x1==0])-length(x2[x2==0]))/length(x1)
  }
  innddex<-innddex+1
  print(innddex)
}

ks_stat_model22<-ks_stat
ks_stat1_model22<-ks_stat1
ks_stat2_model22<-ks_stat2

dataa<-array(0,dim<-c(nCells,dim(a4)[2]))
dataa[,(indd==1)]<-simulated3
dataa[,(indd==0)]<-simulated4

return(dataa)

}
















