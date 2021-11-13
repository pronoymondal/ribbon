ribbon_de_bimodal<-function(dataa1,dataa2)
{

names11<-NULL
names12<-NULL
names21<-NULL
names22<-NULL
names31<-NULL
names32<-NULL

data<-rbind(dataa1,dataa2)

if(is.null(colnames(data)))
{
  colnames(data)<-1:(dim(data)[2])
  colnames(dataa1)<-1:(dim(dataa1)[2])
  colnames(dataa2)<-1:(dim(dataa2)[2])
}

if(is.null(rownames(data)))
{
  rownames(data)<-1:(dim(data)[1])
  rownames(dataa1)<-1:(dim(dataa1)[1])
  rownames(dataa2)<-1:(dim(dataa2)[1])
}

genes<-colnames(data)
genes1<-genes
genes2<-genes




sig<-NULL
mu<-NULL
for(j in 1:(dim(data)[2]))
{
  sig[j]<-sd(log(data[(data[,j]>0),j]))
  mu[j]<-mean(log(data[(data[,j]>0),j]))
}


prob<-colMeans(data>0)
data1<-data[,(prob<0.1)]
data2<-data[,((prob>=0.1)&(prob<=0.95))]
data3<-data[,(prob>0.95)]
prob1<-prob[(prob<0.1)]
prob2<-prob[((prob>=0.1)&(prob<=0.95))]
prob3<-prob[(prob>0.95)]

prob1<-colMeans(dataa1>0)
dataa11<-dataa1[,(prob<0.1)]
dataa21<-dataa1[,((prob>=0.1)&(prob<=0.95))]
dataa31<-dataa1[,(prob>0.95)]
prob11<-prob1[(prob<0.1)]
prob21<-prob1[((prob>=0.1)&(prob<=0.95))]
prob31<-prob1[(prob>0.95)]

prob2<-colMeans(dataa2>0)
dataa12<-dataa2[,(prob<0.1)]
dataa22<-dataa2[,((prob>=0.1)&(prob<=0.95))]
dataa32<-dataa2[,(prob>0.95)]
prob12<-prob1[(prob<0.1)]
prob22<-prob1[((prob>=0.1)&(prob<=0.95))]
prob32<-prob1[(prob>0.95)]

if((dim(data2)[2])>0)
{
  indd<-check_bimodality(t(data2))

  indd[is.na(indd)]<-0

  data21<-t(t(data2[,(indd==0)]))
  data22<-t(t(data2[,(indd==1)]))

  data211<-t(t(dataa21[,(indd==0)]))
  data212<-t(t(dataa21[,(indd==1)]))

  data221<-t(t(dataa22[,(indd==0)]))
  data222<-t(t(dataa22[,(indd==1)]))

  prob21<-colMeans(data21>0)
  prob22<-colMeans(data22>0)

  prob211<-colMeans(data211>0)
  prob212<-colMeans(data212>0)

  prob221<-colMeans(data221>0)
  prob222<-colMeans(data222>0)

  names21<-colnames(data21)
  names22<-colnames(data22)
}


if((dim(t(t(data1)))[2])>0)
{
  indd<-check_bimodality(t(data1))

  indd[is.na(indd)]<-0

  data11<-t(t(t(t(data1))[,(indd==0)]))
  data12<-t(t(t(t(data1))[,(indd==1)]))

  data111<-t(t(t(t(dataa11))[,(indd==0)]))
  data112<-t(t(t(t(dataa11))[,(indd==1)]))

  data121<-t(t(t(t(dataa12))[,(indd==0)]))
  data122<-t(t(t(t(dataa12))[,(indd==1)]))

  names11<-colnames(data11)
  names12<-colnames(data12)

}

if((dim(t(t(data3)))[2])>0)
{
  indd<-check_bimodality(t(data3))

  indd[is.na(indd)]<-0

  data31<-t(t(t(t(data3))[,(indd==0)]))
  data32<-t(t(t(t(data3))[,(indd==1)]))

  data311<-t(t(t(t(dataa31))[,(indd==0)]))
  data312<-t(t(t(t(dataa31))[,(indd==1)]))

  data321<-t(t(t(t(dataa32))[,(indd==0)]))
  data322<-t(t(t(t(dataa32))[,(indd==1)]))

  names31<-colnames(data31)
  names32<-colnames(data32)

}

stat11<-NULL
stat12<-NULL
stat21<-NULL
stat22<-NULL
stat31<-NULL
stat32<-NULL

statt11<-NULL
statt12<-NULL
statt21<-NULL
statt22<-NULL
statt31<-NULL
statt32<-NULL


if("data31"%in%(ls()))
{
if(dim(t(t(data31)))[2]>0)
{
stat31<-likelihood_bimodal(t(t(data31)),t(t(data311)),t(t(data321)))
}
}
if("data11"%in%(ls()))
{
if(dim(t(t(data11)))[2]>0)
{
  stat11<-likelihood_bimodal(t(t(data11)),t(t(data111)),t(t(data121)))
}
}


if("data32"%in%(ls()))
{
if(dim(t(t(data32)))[2]>0)
{
stat32<-likelihood_bimodal(t(t(data32)),t(t(data312)),t(t(data322)))
}
}
if("data12"%in%(ls()))
{
if(dim(t(t(data12)))[2]>0)
{
  stat12<-likelihood_bimodal(t(t(data12)),t(t(data112)),t(t(data122)))
}
}

if("data21"%in%(ls()))
{
if(dim(t(t(data21)))[2]>0)
{
  stat21<-likelihood_bimodal(t(t(data21)),t(t(data211)),t(t(data221)))
}
}

if("data22"%in%(ls()))
{
if(dim(t(t(data22)))[2]>0)
{
  stat22<-likelihood_bimodal(t(t(data22)),t(t(data212)),t(t(data222)))
}
}


if("data21"%in%(ls()))
{
if(dim(t(t(data21)))[2]>0)
{
  statt21<-zero_statistic(t(t(data21>0)),t(t(data211>0)),t(t(data221>0)))
}
}

if("data22"%in%(ls()))
{
if(dim(t(t(data22)))[2]>0)
{
  statt22<-zero_statistic(t(t(data22>0)),t(t(data212>0)),t(t(data222>0)))
}
}


if("data11"%in%(ls()))
{
  if(dim(t(t(data11)))[2]>0)
  {
    statt11<-zero_statistic(t(t(data11>0)),t(t(data111>0)),t(t(data121>0)))
  }
}

if("data12"%in%(ls()))
{
  if(dim(t(t(data12)))[2]>0)
  {
    statt12<-zero_statistic(t(t(data12>0)),t(t(data112>0)),t(t(data122>0)))
  }
}

if("data31"%in%(ls()))
{
  if(dim(t(t(data31)))[2]>0)
  {
    statt31<-zero_statistic(t(t(data31>0)),t(t(data311>0)),t(t(data321>0)))
  }
}

if("data32"%in%(ls()))
{
  if(dim(t(t(data32)))[2]>0)
  {
    statt32<-zero_statistic(t(t(data32>0)),t(t(data312>0)),t(t(data322>0)))
  }
}

pval<-rep(NA,dim(data)[2])
statistic<-rep(NA,dim(data)[2])

statistic[names11]<-stat11+statt11
statistic[names12]<-stat12+statt12
statistic[names21]<-stat21+statt21
statistic[names22]<-stat22+statt22
statistic[names31]<-stat31+statt31
statistic[names32]<-stat32+statt32

pval<-1-pchisq(statistic,2)

lst<-list(statistic=statistic,pval=pval)

return(lst)
}


