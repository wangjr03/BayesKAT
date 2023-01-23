
## MCMC based bayesian optimum kernel selection and testing.


rm(list=ls())
###################
# load packages ###
###################
library(Matrix)
library(stringr)
library(MASS)
library(SKAT)
library(kernlab)
library(KRLS)
library(emdbook)
library(matrixcalc)
library(LaplacesDemon)
library(mvtnorm)
library(ttutils)
library(methods)
library(psych)
library(EQL)
library(mvtnorm)
library(Matrix)
library(Rlab)
library(VGAM)
library(coda)
library(mcmcse)
library(BayesianTools)
library(truncnorm)
library(invgamma)


fn_iteration<-function(iter){

set.seed(iter)
#define number of features and number of samples
np<-500
nsamp<-500

#create correlation matrix
r=0.6 
sigma<-toeplitz(sapply(1:np,function(i) r^(i-1)))
Z<-mvrnorm(nsamp,mu=rep(0,np),Sigma = sigma, tol = 0)

X11=rnorm(nsamp,2,1)
X22=rbern(nsamp,0.6)
X<-cbind(X11,X22)
beta1<-c(0.03,0.5)
fixed_part=X%*%beta1 

tao=10 # as p=1000

#signal=tao*(0.6*(Z[,1]*Z[,3]))
#signal=tao*(0.1*(Z[,1]-Z[,3])+0.8*cos(Z[,3])*exp(-(Z[,3]^2)/5)) 
#signal=tao*(0.4*(Z[,1]*Z[,3]))
#signal=tao*(0.4*(Z[,1]*Z[,3]))
signal=tao*(0.4*(Z[,1]*Z[,3]))
y=signal+fixed_part+rnorm(nsamp,0,1)

y=scale(y)
X=scale(X)

gen.ker <- function(covx,kernel.index){ 
#covx:X; 
#kernel.index=c("Gau","Lin","Quad","IBS"): index of kernel function;

  n <- nrow(covx)
  p <- ncol(covx)
  ker <- matrix(0,n,n)
  for (i in 1:n)
    for (j in i:n)
    {
      x <- covx[i,]
      y <- covx[j,]
      
      if (kernel.index=="Gau")
      {
        ker[i,j] <- exp(-sum((x-y)^2)/p) #gaussian kernel 
      }
      
      if (kernel.index=="Lin")
      {
        ker[i,j] <- sum(x*y)  #linear kernel        
      }
      
      if (kernel.index=="Quad")
      {
        ker[i,j]=(sum(x*y)+1)^2 # Quadratic kernel
      }
      if (kernel.index=="IBS")
      {
        ker[i,j] <- 1-sum(abs(x-y))/(2*p) #IBS kernel
      }
    }
    ker <- as.matrix(forceSymmetric(ker))
    ker0 <- ker
    diag(ker0) <- rep(0,n)
    J <- matrix(1,n,n)
    ker.cen <- ker-J%*%ker0/(n-1)-ker0%*%J/(n-1)+J%*%ker0%*%J/n/(n-1)
    v1.cen <- tr(ker.cen)/n
    return(list(ker.cen=ker.cen,v1.cen=v1.cen))
}

###########################################################################################################################################

k1=gen.ker(covx=Z,kernel.index="Lin")
k2=gen.ker(covx=Z,kernel.index="Quad")
k3=gen.ker(covx=Z,kernel.index="Gau")

gker3<- k3$ker.cen/k3$v1.cen
gker1<- k1$ker.cen/k1$v1.cen
gker2<- k2$ker.cen/k2$v1.cen

####MCMC###################################################################################

library(BayesianTools)

#H0.

# Create a general prior distribution by specifying an arbitrary density function and a
# corresponding sampling function
density = function(par){
d1 = dinvgamma(par[1], shape=2,scale=2, log =TRUE)
d2 = dmvnorm(par[2:3], mean= c(0,0), sigma = diag(c(10,10)), log =TRUE)
return(d1 + d2)
}

sampler = function(n=1){
d1 = rinvgamma(n, shape=2,scale=2)
d2 = rmvnorm(n, mean= c(0,0), sigma=diag(c(10,10)))
return(cbind(d1,d2))
}
prior <- createPrior(density = density, sampler = sampler,
lower = c(0,-Inf,-Inf), upper = c(100,Inf,Inf), best = c(1,0,0))

generator = createProposalGenerator(covariance=c(1,1,1))

likelihood0 <- function(param){
      V=param[1]*diag(rep(1,nsamp))
	  beta=param[2:3]
	  mu1=X%*%beta
      singlelikelihoods = dmvnorm(c(y),c(mu1),V,log=TRUE)
      return(singlelikelihoods)   
     }

setUp0 <- createBayesianSetup(likelihood0, 
                                   prior=prior)
settings_2 <- list(proposalGenerator = generator,optimize=T,adapt=T,nrChains = 3,adaptationInterval = 500, adaptationNotBefore
= 2000, iterations = 50000)
out0 <- runMCMC(bayesianSetup = setUp0,sampler="Metropolis",settings=settings_2)

#H1
prob1=rep(1,3)
#install.packages("TruncatedNormal")
library("TruncatedNormal")


density = function(par){
#d1 = dinvgamma(par[1], shape=3,scale=0.5, log =TRUE)
d2=  dinvgamma(par[2], shape=2,scale=2, log =TRUE)
d1 = dunif(par[1], 0,2, log =TRUE)
#d2 = dunif(par[2], 0,2, log =TRUE)
d3 = dmvnorm(par[3:4], mean= c(0,0), sigma = diag(c(10,10)), log =TRUE)
#s=sum(par[5:7])
#par_new=par[5:7]/s
d4=dgamma(par[5],shape=1,log=TRUE)
d5=dgamma(par[6],shape=1,log=TRUE)
d6=dgamma(par[7],shape=1,log=TRUE)
#d4 = ddirichlet(par[5:7],prob1,log =TRUE)
#d4=dtmvnorm(par[5:7],mu=c(1/3,1/3,1/3),sigma=diag(c(1,1,1)),lb=c(0,0,0),ub=c(1,1,1),log=TRUE)
#d3 = dnorm(par[3], mean= 0, sd = 10, log =TRUE)

return(d1 + d2 + d3+ d4+d5+d6)
}
# The sampling is optional but recommended because the MCMCs can generate automatic starting
# conditions if this is provided

sampler = function(n=1){
d1=runif(n,0,2)
d2 = rinvgamma(n, shape=2,scale=2)
d3 = rmvnorm(n, mean= c(0,0), sigma=diag(c(10,10)))
d4=rgamma(n,shape=1)
d5=rgamma(n,shape=1)
d6=rgamma(n,shape=1)
return(cbind(d1,d2,d3,d4,d5,d6))
}
prior <- createPrior(density = density, sampler = sampler,
lower = c(0,0,-Inf,-Inf,0,0,0), upper = c(2,100,Inf,Inf,100,100,100), best = c(0.5,1,0,0,1/3,1/3,1/3))

likelihood1 <- function(param){
	  s=sum(param[5:7])
	  par_new=param[5:7]/s
	  Ker=par_new[1]*gker1 + par_new[2]*gker2 + par_new[3]*gker3
      V=param[2]*(param[1]*Ker+diag(rep(1,nsamp))) ####checked parametrization
	  beta=param[3:4]
	  mu1=X%*%beta
      singlelikelihoods = dmvnorm(c(y),c(mu1),V,log=TRUE)
      return(singlelikelihoods)  
     }

setUp1 <- createBayesianSetup(likelihood1, 
                                   prior=prior)

generator = createProposalGenerator(covariance=c(1,1,1,1,1,1,1))
settings_2 <- list(proposalGenerator = generator,optimize=T,adapt=T,adaptationInterval = 500, adaptationNotBefore
= 2000,nrChains = 3, iterations = 50000)
system.time({
out1 <- runMCMC(bayesianSetup = setUp1,sampler="Metropolis",settings=settings_2)
})

print("The summary of out0 is:")
print(summary(out0))
print("The summary of out1 is:")
print(summary(out1))

M1 = marginalLikelihood(out0, start = 1000)
M2 = marginalLikelihood(out1, start = 1000)
### Calculating Bayes factor
bf=exp(M2$ln.ML - M1$ln.ML)
MAP1=MAP(out0)$parametersMAP
MAP2=MAP(out1)$parametersMAP

return(c(bf,MAP1,MAP2,M1$ln.ML,M2$ln.ML))

}
                       
