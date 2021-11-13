
zero_statistic<-function(E,E1,E2)
{
theta<-rep(0,(dim(E)[1]))
mu<-rep(0,(dim(E)[2]))
c<-0

eps<-1
while(eps>0.0001)
{
  thetaold<-theta
  muold<-mu
  cold<-c
  theta1<-theta[1:(dim(E1)[1])]
  theta2<-theta[(dim(E1)[1]+1):(dim(E)[1])]
  
  matt<-array(c,dim<-dim(E))+rep(theta,dim(E)[2])+rep(mu,each=dim(E)[1])
  matt[matt>7]<-7
  matt[matt<(-7)]<-(-7)
  matt1<-pnorm(matt)
  matt2<-dnorm(matt)
  
  theta<-theta+rowSums((E-matt1)*matt2/(matt1*(1-matt1)))/rowSums(matt2^2/(matt1*(1-matt1)))
  mu<-mu+colSums((E-matt1)*matt2/(matt1*(1-matt1)))/colSums(matt2^2/(matt1*(1-matt1)))
  c<-c+sum((E-matt1)*matt2/(matt1*(1-matt1)))/sum(matt2^2/(matt1*(1-matt1)))
  
  
  theta1<-theta[1:(dim(E1)[1])]
  theta2<-theta[(dim(E1)[1]+1):(dim(E)[1])]
  
  theta1<-theta1-mean(theta1)
  theta2<-theta2-mean(theta2)
  theta<-c(theta1,theta2)
  
  #theta<-theta-mean(theta)
  mu<-mu-mean(mu)
  mu[mu>7]<-7
  mu[mu<(-7)]<-(-7)
  c[c>7]<-7
  c[c<(-7)]<-(-7)
  theta[theta>7]<-7
  theta[theta<(-7)]<-(-7)
  
  eps<-sum(abs(theta-thetaold))+sum(abs(mu-muold))+sum(abs(c-cold))
  print(eps)
  
}

l<-colSums(dbinom(E,1,matt1,log=TRUE))


theta1<-theta[1:(dim(E1)[1])]
theta2<-theta[(dim(E1)[1]+1):(dim(E)[1])]


#theta<-rep(0,(dim(E)[1]))
#theta1<-rep(0,(dim(E1)[1]))
#theta2<-rep(0,(dim(E2)[1]))
#mu1<-rep(0,(dim(E1)[2]))
#mu2<-rep(0,(dim(E2)[2]))
#c<-0
#c2<-0
mu1<-mu
mu2<-mu

eps<-1
iter<-1
while(eps>0.0001)
{
  iter<-iter+1
  thetaold<-theta
  mu1old<-mu1
  mu2old<-mu2
  cold<-c
  #c2old<-c2
  #theta1<-theta[1:(dim(E1)[1])]
  #theta2<-theta[(dim(E1)[1]+1):(dim(E)[1])]
  
  matt<-rbind((array(c,dim<-dim(E1))+rep(theta1,dim(E1)[2])+rep(c(mu1),each=dim(E1)[1])),(array(c,dim<-dim(E2))+rep(theta2,dim(E2)[2])+rep(c(mu2),each=dim(E2)[1])))
  matt[matt>7]<-7
  matt[matt<(-7)]<-(-7)
  
  matt1<-t(t(pnorm(matt)))
  matt2<-t(t(dnorm(matt)))
  matt11<-t(t(matt1[(1:(dim(E1)[1])),]))
  matt12<-t(t(matt1[(((dim(E1)[1])+1):(dim(E)[1])),]))
  matt21<-t(t(matt2[(1:(dim(E1)[1])),]))
  matt22<-t(t(matt2[(((dim(E1)[1])+1):(dim(E)[1])),]))
  
  #theta<-theta+rowSums((E-matt1)*matt2/(matt1*(1-matt1)))/rowSums(matt2^2/(matt1*(1-matt1)))
  mu1<-mu1+colSums((E1-matt11)*matt21/(matt11*(1-matt11)))/colSums(matt21^2/(matt11*(1-matt11)))
  mu2<-mu2+colSums((E2-matt12)*matt22/(matt12*(1-matt12)))/colSums(matt22^2/(matt12*(1-matt12)))
  mu1[mu1>7]<-7
  mu1[mu1<(-7)]<-(-7)
  mu2[mu2>7]<-7
  mu2[mu2<(-7)]<-(-7)
  #c<-c+sum((E-matt1)*matt2/(matt1*(1-matt1)))/sum(matt2^2/(matt1*(1-matt1)))
  #c2<-c2+sum((E2-matt12)*matt22/(matt12*(1-matt12)))/sum(matt22^2/(matt12*(1-matt12)))
  
  #theta1<-theta[1:(dim(E1)[1])]
  #theta2<-theta[(dim(E1)[1]+1):(dim(E)[1])]
  
  #theta1<-theta1-mean(theta1)
  #theta2<-theta2-mean(theta2)
  #theta<-c(theta1,theta2)
  #mu1<-mu1-mean(mu1)
  #mu2<-mu2-mean(mu2)
  
  
  eps<-sum(abs(theta-thetaold))+sum(abs(mu1-mu1old))+sum(abs(mu2-mu2old))+sum(abs(c-cold))
  print(eps)
  plot(theta)
  
}


lalt<-colSums(dbinom(E,1,matt1,log=TRUE))

return(lalt-l)

}








