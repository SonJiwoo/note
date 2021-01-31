#### NELS data
load("nelsSES.RData") 

ids<-sort(unique(nels$sch_id)) 
m<-length(ids)
Y<-list() ; X<-list() ; N<-NULL
for(j in 1:m) 
{
  Y[[j]]<-nels[nels$sch_id==ids[j], 4] 
  N[j]<- sum(nels$sch_id==ids[j])
  xj<-nels[nels$sch_id==ids[j], 3] 
  xj<-(xj-mean(xj))
  X[[j]]<-cbind( rep(1,N[j]), xj  )
}



#### OLS fits
S2.LS<-BETA.LS<-NULL
for(j in 1:m) {
  fit<-lm(Y[[j]]~-1+X[[j]] )
  BETA.LS<-rbind(BETA.LS,c(fit$coef)) 
  S2.LS<-c(S2.LS, summary(fit)$sigma^2) 
} 



#### Figure 11.1
pdf("fig11_1.pdf",family="Times",height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))

plot( range(nels[,3]),range(nels[,4]),type="n",xlab="SES", 
      ylab="math score")
for(j in 1:m) {    abline(BETA.LS[j,1],BETA.LS[j,2],col="gray")  }

BETA.MLS<-apply(BETA.LS,2,mean)
abline(BETA.MLS[1],BETA.MLS[2],lwd=2)

plot(N,BETA.LS[,1],xlab="sample size",ylab="intercept")
abline(h= BETA.MLS[1],col="black",lwd=2)
plot(N,BETA.LS[,2],xlab="sample size",ylab="slope")
abline(h= BETA.MLS[2],col="black",lwd=2)

dev.off()



#### Hierarchical regression model

## mvnormal simulation
rmvnorm<-function(n,mu,Sigma)
{ 
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol(Sigma)) +c(mu))
}

## Wishart simulation
rwish<-function(n,nu0,S0)
{
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}

## Setup
p<-dim(X[[1]])[2]
theta<-mu0<-apply(BETA.LS,2,mean)
nu0<-1 ; s2<-s20<-mean(S2.LS)
eta0<-p+2 ; Sigma<-S0<-L0<-cov(BETA.LS) ; BETA<-BETA.LS
THETA.b<-S2.b<-NULL
iL0<-solve(L0) ; iSigma<-solve(Sigma)
Sigma.ps<-matrix(0,p,p)
SIGMA.PS<-NULL
BETA.ps<-BETA*0
BETA.pp<-NULL
set.seed(1)
mu0[2]+c(-1.96,1.96)*sqrt(L0[2,2])

## MCMC
for(s in 1:10000) {
  ##update beta_j 
  for(j in 1:m) 
  {  
    Vj<-solve( iSigma + t(X[[j]])%*%X[[j]]/s2 )
    Ej<-Vj%*%( iSigma%*%theta + t(X[[j]])%*%Y[[j]]/s2 )
    BETA[j,]<-rmvnorm(1,Ej,Vj) 
  } 
  ##
  
  ##update theta
  Lm<-  solve( iL0 +  m*iSigma )
  mum<- Lm%*%( iL0%*%mu0 + iSigma%*%apply(BETA,2,sum))
  theta<-t(rmvnorm(1,mum,Lm))
  ##
  
  ##update Sigma
  mtheta<-matrix(theta,m,p,byrow=TRUE)
  iSigma<-rwish(1, eta0+m, solve( S0+t(BETA-mtheta)%*%(BETA-mtheta) ) )
  ##
  
  ##update s2
  RSS<-0
  for(j in 1:m) { RSS<-RSS+sum( (Y[[j]]-X[[j]]%*%BETA[j,] )^2 ) }
  s2<-1/rgamma(1,(nu0+sum(N))/2, (nu0*s20+RSS)/2 )
  ##
  ##store results
  if(s%%10==0) 
  { 
    cat(s,s2,"\n")
    S2.b<-c(S2.b,s2);THETA.b<-rbind(THETA.b,t(theta))
    Sigma.ps<-Sigma.ps+solve(iSigma) ; BETA.ps<-BETA.ps+BETA
    SIGMA.PS<-rbind(SIGMA.PS,c(solve(iSigma)))
    BETA.pp<-rbind(BETA.pp,rmvnorm(1,theta,solve(iSigma)) )
  }
  ##
}



## MCMC diagnostics
library(coda)
effectiveSize(S2.b)
effectiveSize(THETA.b[,1])
effectiveSize(THETA.b[,2])

apply(SIGMA.PS,2,effectiveSize)

tmp<-NULL;for(j in 1:dim(SIGMA.PS)[2]) { tmp<-c(tmp,acf(SIGMA.PS[,j])$acf[2]) }

acf(S2.b)
acf(THETA.b[,1])
acf(THETA.b[,2])



#### Figure 11.3 
pdf("fig11_3.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

plot(density(THETA.b[,2],adj=2),xlim=range(BETA.pp[,2]), 
     main="",xlab="slope parameter",ylab="posterior density",lwd=2)
lines(density(BETA.pp[,2],adj=2),col="gray",lwd=2)
legend( -3 ,1.0 ,legend=c( expression(theta[2]),expression(tilde(beta)[2])), 
        lwd=c(2,2),col=c("black","gray"),bty="n") 

quantile(THETA.b[,2],prob=c(.025,.5,.975))
mean(BETA.pp[,2]<0) 

BETA.PM<-BETA.ps/1000
plot( range(nels[,3]),range(nels[,4]),type="n",xlab="SES",
      ylab="math score")
for(j in 1:m) {    abline(BETA.PM[j,1],BETA.PM[j,2],col="gray")  }
abline( mean(THETA.b[,1]),mean(THETA.b[,2]),lwd=2 )
dev.off()



#### Tumor location example
load("tumorLocation.RData") 
Y<-tumorLocation 


#### Figure 11.4
pdf("fig11_4.pdf",family="Times",height=3.5,width=7)
xs<-seq(.05,1,by=.05) 
m<-nrow(Y) 
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
plot(c(0,1),range(Y),type="n",xlab="location",ylab="number of tumors")
for(j in 1:m) { lines(xs,Y[j,],col="gray") }
lines( xs,apply(Y,2,mean),lwd=3)

lya<-log(apply(Y,2,mean))

X<-cbind( rep(1,ncol(Y)),poly(xs,deg=4,raw=TRUE))
fit2<- lm(lya~-1+X[,1:3] )
fit3<- lm(lya~-1+X[,1:4] )
fit4<- lm(lya~-1+X[,1:5] )

yh2<-X[,1:3]%*%fit2$coef
yh3<-X[,1:4]%*%fit3$coef
yh4<-X[,1:5]%*%fit4$coef

plot(xs,lya,type="l",lwd=3,xlab="location",ylab="log average number of tumors",
     ylim=range(c(lya,yh2,yh3,yh4)) )

points(xs,yh2,pch="2",col="black")
lines(xs,yh2,col="gray")
points(xs,yh3,pch="3",col="black")
lines(xs,yh3,col="gray")
points(xs,yh4,pch="4",col="black")
lines(xs,yh4,col="gray")
dev.off()



#### MCMC  

## mvnorm log density
ldmvnorm<-function(X,mu,Sigma,iSigma=solve(Sigma),dSigma=det(Sigma))
{
  Y<-t( t(X)-mu)
  sum(diag(-.5*t(Y)%*%Y%*%iSigma))  -
    .5*(  prod(dim(X))*log(2*pi) +     dim(X)[1]*log(dSigma) )
}

m<-nrow(Y) ; p<-ncol(X)


## priors
BETA<-NULL
for(j in 1:m)
{
  BETA<-rbind(BETA,lm(log(Y[j,]+1/20)~-1+X)$coef) 
}

mu0<-apply(BETA,2,mean) 
S0<-cov(BETA) ; eta0<-p+2
iL0<-iSigma<-solve(S0)

## MCMC
THETA.post<<-SIGMA.post<-NULL ; set.seed(1)
for(s in 1:50000) 
{
  
  ##update theta
  Lm<-solve( iL0 +  m*iSigma )
  mum<-Lm%*%( iL0%*%mu0 + iSigma%*%apply(BETA,2,sum) )
  theta<-t(rmvnorm(1,mum,Lm))
  
  ##update Sigma
  mtheta<-matrix(theta,m,p,byrow=TRUE)
  iSigma<-rwish(1,eta0+m, 
                solve( S0+t(BETA-mtheta)%*%(BETA-mtheta)) )
  
  ##update beta
  Sigma<-solve(iSigma) ; dSigma<-det(Sigma)
  for(j in 1:m)
  {
    beta.p<-t(rmvnorm(1,BETA[j,],.5*Sigma))
    
    lr<-sum( dpois(Y[j,],exp(X%*%beta.p),log=TRUE ) -
               dpois(Y[j,],exp(X%*%BETA[j,]),log=TRUE ) ) +
      ldmvnorm( t(beta.p),theta,Sigma,
                iSigma=iSigma,dSigma=dSigma ) -
      ldmvnorm( t(BETA[j,]),theta,Sigma,
                iSigma=iSigma,dSigma=dSigma )
    
    if( log(runif(1))<lr ) { BETA[j,]<-beta.p }
  }
  
  ##store some output
  if(s%%10==0)
  {  
    cat(s,"\n") 
    THETA.post<-rbind(THETA.post,t(theta)) 
    SIGMA.post<-rbind(SIGMA.post,c(Sigma)) 
  }
  
}


## MCMC diagnostics
library(coda)
round(apply(THETA.post,2,effectiveSize),2)



#### Different posterior regions 
eXB.post<-NULL
for(s in 1:dim(THETA.post)[1])
{
  beta<-rmvnorm(1,THETA.post[s,],matrix(SIGMA.post[s,],p,p))
  eXB.post<-rbind(eXB.post,t(exp(X%*%t(beta) )) )
}

qEB<-apply( eXB.post,2,quantile,probs=c(.025,.5,.975))

eXT.post<- exp(t(X%*%t(THETA.post )) )
qET<-apply( eXT.post,2,quantile,probs=c(.025,.5,.975))
yXT.pp<-matrix( rpois(prod(dim(eXB.post)),eXB.post),
                dim(eXB.post)[1],dim(eXB.post)[2] )

qYP<-apply( yXT.pp,2,quantile,probs=c(.025,.5,.975))



#### Figure 11.5 
pdf("fig11_5.pdf",family="Times",height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))

plot( c(0,1),range(c(0,qET,qEB,qYP)),type="n",xlab="location",
      ylab="number of tumors")
lines(xs, qET[1,],col="black",lwd=1)
lines(xs, qET[2,],col="black",lwd=2)
lines(xs, qET[3,],col="black",lwd=1)

plot( c(0,1),range(c(0,qET,qEB,qYP)),type="n",xlab="location",
      ylab="")
lines(xs, qEB[1,],col="black",lwd=1)
lines(xs, qEB[2,],col="black",lwd=2)
lines(xs, qEB[3,],col="black",lwd=1)

plot( c(0,1),range(c(0,qET,qEB,qYP)),type="n",xlab="location",
      ylab="")
lines(xs, qYP[1,],col="black",lwd=1)
lines(xs, qYP[2,],col="black",lwd=2)
lines(xs, qYP[3,],col="black",lwd=1)

dev.off()


