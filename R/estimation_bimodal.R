
estimate_bimodal<-function(rpc_matrix,mu1,mu2,sig1,sig2,pi,D)
{
  library(mixtools)
  tau1vec<-rep(1,dim(rpc_matrix)[1])
  tau2vec<-rep(1,dim(rpc_matrix)[1])
  data<-rpc_matrix
  d<-dim(data)

  
  mu1old<-mu1
  mu2old<-mu2
  sig1old<-sig1
  sig2old<-sig2
  piold<-pi
  
  
  mu1new<-mu1
  mu2new<-mu2
  cnew<-rep(0,length(mu1))
  c1new<-0
  tau<-1
  tau1<-1
  tau2<-1
  alpha1<-rep(0,(dim(rpc_matrix)[1]))
  alpha2<-rep(0,(dim(rpc_matrix)[1]))
  d<-dim(rpc_matrix)
  E<-rpc_matrix>0
  ones<-array(1,dim<-dim(E))
  boundary<-abs(colMeans(pmax(log(rpc_matrix),0.00000001)*E)/colMeans(E))
  
  c<-NULL
  for(j in 1:(dim(rpc_matrix)[2]))
  {
    #mu[j]<-mean(log(rpc_matrix[(rpc_matrix[,j]>0),j]))
    #sig[j]<-sd(log(rpc_matrix[(rpc_matrix[,j]>0),j]))
    c[j]<-qnorm(mean(E[,j]))
    print(j)
  }
  
  eps<-2
  iter<-1
  
  i<-0
  xxx<-which(rpc_matrix>0)
  jx<-floor(xxx/d[1])+1
  jx[jx>(d[2])]<-(d[2])
  ix<-xxx-d[1]*(jx-1)+1
  ix[ix>(d[1])]<-(d[1])
  xxx1<-which(data==0)
  jx1<-floor(xxx1/d[1])+1
  jx1[(jx1>(d[2]))]<-(d[2])
  ix1<-xxx1-d[1]*(jx1-1)
  
  sig1[sig1==0]<-1
  
  sig<-sig1
  prob<-colMeans(rpc_matrix>0)
  
  sig1o<-sig1
  
  while((eps>1)&(iter<=50))
  {
    iter<-iter+1
    mu1<-mu1new
    tauold<-tau
    tau1old<-tau1
    tau2old<-tau2
    c1<-c1new
    
    xx2<-pmax(pmin(mu1new+c1new,7),(-7))
    xx4<-pmax(pmin(rep(mu1new,each=(d[1]))+rep(alpha1,(d[2]))+c1new,7),(-7))
    matt2<-colSums((1-D)*E)/sig^2
    aaa<-((1-D)*(E-pnorm(xx4)))*dnorm(xx4)/(pnorm(xx4)*(1-pnorm(xx4)))
    vec2<-colSums((1-D)*(log(pmax(rpc_matrix,10^(-8)))*E-rep(mu1,each=dim(E)[1])*E-rep(alpha2,(d[2]))*E))/sig^2+matt2*mu1new
    mu1new<-(vec2)/matt2
    #mu1new[mu1new>boundary]<-boundary[mu1new>boundary]
    #mu1new[mu1new<(-boundary)]<-(-boundary[mu1new<(-boundary)])
    mu1new[mu1new<(-10)]<-mu1old[mu1new<(-10)]
    
    
    matt4<-rowSums((1-D)*E/rep(sig^2,each=(d[1])))+1/tau2
    vec4<-rowSums((1-D)*(log(pmax(rpc_matrix,10^(-8)))*E-rep(mu1new,each=dim(E)[1])*E)/rep(sig^2,each=(d[1])))
    vec4<-vec4
    alpha2<-vec4/matt4
    alpha2<-alpha2-mean(alpha2)
    mu1new<-mu1new+mean(alpha2)
    #alpha22<-alpha2^2+1/matt4
    alpha22<-1/matt4
    tau2<-mean(alpha22)
    
    
    mat1<-(1-D)*((log(pmax(rpc_matrix,10^(-8)))*E-rep(mu1new,each=dim(E)[1]))*E)^2
    mat2<-2*(1-D)*((log(pmax(rpc_matrix,10^(-8)))*E-rep(mu1new,each=dim(E)[1]))*E)*rep(alpha2,(d[2]))*E
    mat3<-(1-D)*rep(alpha22,(d[2]))*E
    sig1<-sqrt(colSums(pmax(mat1-mat2+mat3,0)*E)/colSums((1-D)*E))
    sig1[sig1==0]<-sig1o[(sig1==0)]
    sig<-sig1
    eps<-sum(abs(mu1new-mu1))+abs(tau2-tau2old)
    print(eps)
    plot(colMeans(rpc_matrix>0),pnorm(c1new+mu1new+mean(alpha1)))
  }

  return(lst=list(mu1=mu1new,sig1=sig,tau1=tau1,tau2=tau2,c=c,alpha1=alpha1,alpha2=alpha2,mu1old=mu1old,mu2=mu2,sig2=sig2,pi=pi))

  }


