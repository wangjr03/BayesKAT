
## BayesKAT_MCMC
#General Suggestion: If you have a large input dataset, try to use BayesKAT_MAP instead.
#Here we provided a simulated dataset which can be used as input. If you are using your own input data. please make sure the configuration matches with the description in readme document.

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
library(BayesianTools)
library(truncnorm)
library(invgamma)
library(TruncatedNormal)

BayesKAT_MCMC<-function(inputAddress_y,inputAddress_X,inputAddress_Z,prior_H1=0.5,output_address){
	
  	y=read.table(inputAddress_y)
  	X=read.table(inputAddress_X)
  	Z=read.table(inputAddress_Z)
	y=as.matrix(y); X=as.matrix(X); Z=as.matrix(Z)
  	y=scale(y); X=scale(X)
  	prior1=as.numeric(prior_H1)
  	prior0=1-prior1
	nsamp=dim(Z)[1]
	np=dim(Z)[2]
	
	
	if(prior_H1<=0 | prior_H1>=1 | !is.numeric(y) | !is.numeric(X) | !is.numeric(Z) | dim(y)[2]>1 | dim(X)[2]!=3 | !(dim(y)[1]==dim(X)[1]) | !(dim(y)[1]==dim(Z)[1]) | sum(Z<0)>0){
		if(prior_H1<0 | prior_H1>1){
			print("prior_H1 is the prior probability of H1, should be between 0 and 1!")
		}else if(prior_H1==0 | prior_H1==1){
			print("If you already know what the true hypothesis is, there is no point of testing it!")
		}else if(!is.numeric(y)| !is.numeric(X) | !is.numeric(y)){
			print("Please make sure the input addresses are correct and the input data matrices are numeric.")
		}else if(dim(y)[2]>1){
			print("Please make sure the input data y follows the basic configuration given in the readme file")
		}else if(dim(X)[2]!=3){
			print("Please make sure the input data matrix X follow the basic configuration given in the readme file")
		}else if(!(dim(y)[1]==dim(X)[1])| !(dim(y)[1]==dim(Z)[1])){
			print("Please make sure the input data matrices contain information on same number of individuals")
		}else if(!(dim(y)[1]==dim(X)[1])| !(dim(y)[1]==dim(Z)[1])){
			print("Please make sure the input data matrices contain information on same number of individuals")
		}else if(sum(Z<0)>0){
			print("Please make sure the input data matrix Z follow the basic configuration given in the readme file")
		}else{
			print("Please check your input")
		}
	
	}else{
	
		if(nsamp<100 & np>100){
			print("Please note that the sample size is too small to have an accurate result!")
		}else if(nsamp>=500){
			print("Please note that, because of high sample size it is going to take a long time to complete the run and you might consider using BayesKAT_MAP instead.")
		}else if(np>=500){
			print("Please note that, due to presence of huge number of genetic features it is going to take a long time to complete the run and you might consider using BayesKAT_MAP instead.")
		}
	
	
	#Function for calculating generalized kernel matrix
	gen.ker <- function(covx,kernel.index){ 
		#covx: Genotype matrix/ gene expression matrix
		#kernel.index=c("Gau","Lin","Quad","IBS"): index of kernel function;
		#Any other kernel function can be added here

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

        ######################################################################################
	##centralized and normalized kernels
	###Note: for demonstration reasons, consider Gau/Lin/Quad kernels as candidate kernel set
	
	k1=gen.ker(covx=Z,kernel.index="IBS")
	k2=gen.ker(covx=Z,kernel.index="Quad")
	k3=gen.ker(covx=Z,kernel.index="Gau")

	gker1<- k1$ker.cen/k1$v1.cen
	gker2<- k2$ker.cen/k2$v1.cen
	gker3<- k3$ker.cen/k3$v1.cen
	
        #######################################################################################
	#sampling from posterior distribution
	#Assuming there are only two demographic covariates.
	##H0: tau=0
	# Create a prior distribution by specifying the density functions and corresponding sampling function
	density = function(par){
		d1 = dinvgamma(par[1], shape=2,scale=2, log =TRUE)
		d2 = dmvnorm(par[2:4], mean= c(0,0,0), sigma = diag(c(10,10,10)), log =TRUE)
		return(d1 + d2)
		}
	
	# This sampling is useful because the MCMCs can generate automatic starting conditions if this is provided
	sampler = function(n=1){
		d1 = rinvgamma(n, shape=2,scale=2)
		d2 = rmvnorm(n, mean= c(0,0,0), sigma=diag(c(10,10,10)))
		return(cbind(d1,d2))
		}
	prior <- createPrior(density = density, sampler = sampler,lower = c(0,-Inf,-Inf,-Inf), upper = c(100,Inf,Inf,Inf), best = c(1,0,0,0))

	generator = createProposalGenerator(covariance=c(1,1,1,1))

	likelihood0 <- function(param){
		V=param[1]*diag(rep(1,nsamp))
		beta=param[2:4]
		mu1=X%*%beta
      		singlelikelihoods = dmvnorm(c(y),c(mu1),V,log=TRUE)
     		return(singlelikelihoods)   
       		}

	setUp0 <- createBayesianSetup(likelihood0, prior=prior)
	
	settings_2 <- list(proposalGenerator = generator,optimize=T,adapt=T,nrChains = 3,adaptationInterval = 500, adaptationNotBefore= 2000, iterations = 50000)
	
	out0 <- runMCMC(bayesianSetup = setUp0,sampler="Metropolis",settings=settings_2)

	##H1: tau>0
	
	prob1=rep(1,3)
	density = function(par){
		d1 = dunif(par[1], 0,2, log =TRUE)
		d2=  dinvgamma(par[2], shape=2,scale=2, log =TRUE)
		d3 = dmvnorm(par[3:5], mean= c(0,0,0), sigma = diag(c(10,10,10)), log =TRUE)
		d4=dgamma(par[6],shape=1,log=TRUE)
		d5=dgamma(par[7],shape=1,log=TRUE)
		d6=dgamma(par[8],shape=1,log=TRUE)
		return(d1 + d2 + d3+ d4+d5+d6)
		}

	sampler = function(n=1){
		d1=runif(n,0,2)
		d2 = rinvgamma(n, shape=2,scale=2)
		d3 = rmvnorm(n, mean= c(0,0,0), sigma=diag(c(10,10,10)))
		d4=rgamma(n,shape=1)
		d5=rgamma(n,shape=1)
		d6=rgamma(n,shape=1)
		return(cbind(d1,d2,d3,d4,d5,d6))
		}
	prior <- createPrior(density = density, sampler = sampler, lower = c(0,0,-Inf,-Inf,-Inf,0,0,0), upper = c(2,100,Inf,Inf,Inf,100,100,100), best = c(0.5,1,0,0,0,1/3,1/3,1/3))

	likelihood1 <- function(param){
	  	s=sum(param[6:8])
	  	par_new=param[6:8]/s
	  	Ker=par_new[1]*gker1 + par_new[2]*gker2 + par_new[3]*gker3
      		V=param[2]*(param[1]*Ker+diag(rep(1,nsamp))) 
	  	beta=param[3:5]
	  	mu1=X%*%beta
      		singlelikelihoods = dmvnorm(c(y),c(mu1),V,log=TRUE)
      		return(singlelikelihoods)  
     		}

	setUp1 <- createBayesianSetup(likelihood1, prior=prior)

	generator = createProposalGenerator(covariance=c(1,1,1,1,1,1,1,1))
	
	settings_2 <- list(proposalGenerator = generator,optimize=T,adapt=T,adaptationInterval = 500, adaptationNotBefore= 2000,nrChains = 3, iterations = 50000)

	out1 <- runMCMC(bayesianSetup = setUp1,sampler="Metropolis",settings=settings_2)


	print("The summary of out0 is:")
	print(summary(out0))
	print("The summary of out1 is:")
	print(summary(out1))
	#check if Gelman Rubin multivariate psrf value is around 1 in each case. Otherwise increase #iterations in settings_2

	M1 = marginalLikelihood(out0, start = 1000)
	M2 = marginalLikelihood(out1, start = 1000)
	### Calculating Bayes factor
	bf=exp(M2$ln.ML - M1$ln.ML)
	MAP_para=MAP(out1)$parametersMAP
	kernel_weights=MAP_para[6:8]/sum(MAP_para[6:8])


post_H1=1/(1+(prior0/prior1)*(1/bf))
result_final=c(post_H1,kernel_weights)
result_names=c("Post_H1","IBS_weight","Quadratic_weight","Gaussian_weight")
write.table(cbind(result_names,result_final),output_address)
return(list(Post_H1=post_H1,Kernel_Weights=kernel_weights))

}
}
