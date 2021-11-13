
estimate_unimodal<-function(rpc_matrix)
{

E<-(rpc_matrix>0)
eps<-1
d<-dim(E)
c<-rep(0,(dim(E)[2]))
cnew<-c
cold<-c
tau<-1
c1<-0
c1new<-c1
c1old<-c1
#mu<-rep(0,dim(E)[2])
mu<-NULL
sig<-NULL
for(j in 1:(dim(rpc_matrix)[2]))
{
  mu[j]<-mean(log(rpc_matrix[(rpc_matrix[,j]>0),j]))
  sig[j]<-sd(log(rpc_matrix[(rpc_matrix[,j]>0),j]))
  c[j]<-qnorm(mean(E[,j]))
  print(j)
}


#mu<-rep(0,(dim(E)[2]))
#sig<-rep(1,(dim(E)[2]))
alpha<-rep(0,(dim(E)[1]))
tau<-1
eps<-2
sigo<-sig
iter<-1
while((eps>1)&(iter<=50))
{
  iter<-iter+1
  muold<-mu
  sigmaold<-sigma
  tauold<-tau
  alphaold<-alpha
  sigold<-sig
  mu<-colSums((log(pmax(rpc_matrix,0.00000001))*E-E*rep(alpha,(dim(E)[2])))/rep(sig^2,each=(d[1])))/colSums(E/rep(sig^2,each=(d[1])))
  #mu<-colSums((log(pmax(rpc_matrix,0.00000001))*E)/rep(sig^2,each=(d[1])))/colSums(E/rep(sig^2,each=(d[1])))
  
  alpha<-(rowSums(log(pmax(rpc_matrix,0.00000001))*E-E*rep(mu,each=(dim(E)[1]))/rep(sig^2,each=(d[1]))))/(rowSums(E/rep(sig^2,each=(d[1])))+1/tau)
  alpha<-alpha-mean(alpha)
  mu<-mu+mean(alpha)
  
  #alpha2<-alpha^2+1/(rowSums(E/rep(sig^2,each=(d[1])))+1/tau)
  alpha2<-1/(rowSums(E/rep(sig^2,each=(d[1])))+1/tau)
  tau<-mean(alpha2)
  mat1<-((log(pmax(rpc_matrix,10^(-8)))*E-rep(mu,each=dim(E)[1]))*E)^2
  mat2<-2*((log(pmax(rpc_matrix,10^(-8)))*E-rep(mu,each=dim(E)[1]))*E)*rep(alpha,(d[2]))*E
  mat3<-rep(alpha2,(d[2]))*E
  sig<-sqrt(colSums(pmax(mat1-mat2+mat3,0)*E)/colSums(E))
  sig[sig<0.00000001]<-sigo[sig<0.00000001]
  eps<-sum(abs(mu-muold))+sum(abs(alpha-alphaold))+sum(abs(sig-sigold))
  print(eps)
}



return(lst=list(mu=mu,alpha=alpha,tau=tau,sig=sig,c=c))
}


