## BayesKAT_MAP
#In Demo Data folder we provided a simulated dataset which can be used as input. If you are using your own input data, please make sure the configuration matches with the description in readme document.


rm(list=ls())
###################
# load packages ###
###################
library(Matrix)
library(MASS)
library(LaplacesDemon)
library(mvtnorm)
library(numDeriv)
library(minqa)

BayesKAT_MAP<-function(inputAddress_y,inputAddress_X,inputAddress_Z,prior_H1=0.5,output_address){
	
  y=read.table(inputAddress_y)
  X=read.table(inputAddress_X)
  Z=read.table(inputAddress_Z)
  y=as.matrix(y); X=as.matrix(X); Z=as.matrix(Z)
  y=scale(y); X=scale(X)
  prior1=as.numeric(prior_H1)
  prior0=1-prior1
  nsamp=dim(Z)[1]
  np=dim(Z)[2]
	
	
  if(prior_H1<=0 | prior_H1>=1 | !is.numeric(y) | !is.numeric(X) | !is.numeric(Z) | dim(y)[2]>1 | dim(X)[2]!=3 | !(dim(y)[1]==dim(X)[1]) | !(dim(y)[1]==dim(Z)[1]) | sum(Z< 0)>0 ){
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
	}
	

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

	k1=gen.ker(covx=Z,kernel.index="IBS")
	gker1<- k1$ker.cen/k1$v1.cen

	k2=gen.ker(covx=Z,kernel.index="Quad")
	gker2<- k2$ker.cen/k2$v1.cen
		
	k3=gen.ker(covx=Z,kernel.index="Gau")
	gker3<- k3$ker.cen/k3$v1.cen
	
####MAP###################################################################################
	id_mat<-diag(rep(1,nsamp))
    
	post1<-function(val1){
      	    s=sum(val1[1:3])
            rho1=val1[1]/s
            rho2=val1[2]/s
            rho3=val1[3]/s
            rho=c(rho1,rho2,rho3)
            tau=val1[4]
            sigsq=val1[5]
            beta=val1[6:8]
            Ker=rho[1]*gker1 + rho[2]*gker2 + rho[3]*gker3
            V<-sigsq*(tau*Ker+id_mat)
            mu1=X%*%beta
            res=dmvnorm(c(y),c(mu1),V,log=TRUE)+dunif(tau, min = 0, max = 2, log = TRUE)+dinvgamma(sigsq, shape=2,scale=2, log =TRUE)+dmvnorm(beta, c(0,0,0),diag(c(10,10,10)), log =TRUE)+dgamma(val1[1],shape=1,log=TRUE)+dgamma(val1[2],shape=1,log=TRUE)+dgamma(val1[3],shape=1,log=TRUE)
            res=ifelse(is.na(res),-1e+1000,res) #for V too close or equal to 0, it might return NaN and model might fail.
            return(-res) #'-'ve because we'll be minimizing
	}
    
    
	inits<-c(1/3,1/3,1/3,0.5,1,0,0,0)
	epsilon=1e-10
	sigsq_lower=0.1

	res=minqa::bobyqa(inits,post1,lower = c(epsilon,epsilon,epsilon,0,sigsq_lower,-1e+10,-1e+10,-1e+10),upper=c(100,100,100,2,100,1e+10,1e+10,1e+10),control=list(npt=12))
	ress=res$par
	s=sum(ress[1:3])
	result1=c(ress[1:3]/s,ress[4:8])
	tau=result1[4]
	if(tau==0){
     	     bf=0
	}else{
      	     #under H0
	     id_mat<-diag(rep(1,nsamp))
	     fn_post0<-function(val0){
        	sigsq=val0[1]
        	beta=val0[2:4]
        	V<-sigsq*id_mat
        	mu1=X%*%beta
        	res=dmvnorm(c(y),c(mu1),V,log=TRUE)+dinvgamma(sigsq, shape=2,scale=2, log =TRUE)+dmvnorm(beta,c(0,0,0),diag(c(10,10,10)), log =TRUE)
       		return(res)
		}
      
	      inits=c(1,0,0,0)
	      res1=optim(inits,fn_post0,method="L-BFGS-B",control=list(fnscale=-1),lower = c(sigsq_lower,-1e+10,-1e+10,-1e+10),upper=c(100,1e+10,1e+10,1e+10),hessian=TRUE)
      	      neg_hes=-res1$hessian
      	      psi_hat=solve(neg_hes)
              det_psi_hat=det(psi_hat)
              log_post_hat=res1$value
      
              lap01= ((2*pi)^(4/2))*sqrt(det_psi_hat)
              lap02=log_post_hat
      
       	      #under H1
	      Ker=result1[1]*gker1 + result1[2]*gker2 + result1[3]*gker3
	      id_mat<-diag(rep(1,nsamp))
              fn_post1<-function(val1){
        	 tau=val1[1]
        	 sigsq=val1[2]
                 beta=val1[3:5]
                 V<-sigsq*(tau*Ker+id_mat)
                 mu1=X%*%beta
        	 res=dmvnorm(c(y),c(mu1),V,log=TRUE)+dinvgamma(sigsq, shape=2,scale=2, log =TRUE)+dmvnorm(beta, c(0,0,0),diag(c(10,10,10)), log =TRUE)+dunif(tau,0,2,log=TRUE)
        	 return(res)
		}
      
     	       val=result1[4:8]
      	       h1=hessian(fn_post1,val)
      
      	       if(sum(is.na(as.array(h1)))==0){  ##adding extra condition
      	          neg_hes1=-h1
      		  psi_hat1=solve(neg_hes1)
      		  det_psi_hat1=det(psi_hat1)
      		  sig_hat=result1[5]
      		  tau_hat=result1[4]
      		  log_post_hat1=fn_post1(val)
      
      		  lap11= ((2*pi)^(5/2))*sqrt(det_psi_hat1)
      		  lap12=log_post_hat1
      
      		  #print(c(lap01,lap02,lap11,lap12))
      		  bf=(lap11/lap01)*(exp(lap12-lap02)) 
      		}else{

        	  fn_post1_B<-function(val1){
          		tau=val1[1]
          		sigsq=val1[2]
          		beta=val1[3:5]
          		V<-sigsq*(tau*Ker+id_mat)
          		mu1=X%*%beta
          		res=dmvnorm(c(y),c(mu1),V,log=TRUE)
          		return(res)
        	   }
        	   inits<-c(0.5,1,0,0,0)
        	   epsilon=1e-10
        	   sigsq_lower=0.1
        	   res=optim(inits,fn_post1_B,method="L-BFGS-B",control=list(fnscale=-1),lower = c(0,sigsq_lower,-1e+10,-1e+10,-1e+10),upper=c(2,100,1e+10,1e+10,1e+10),hessian=TRUE)
                   ress=res$par
                   h1=res$hessian
                   neg_hes1=-h1
                   psi_hat1=solve(neg_hes1)
                   det_psi_hat1=det(psi_hat1)
                   res1=optim(inits,fn_post1_B,method="L-BFGS-B",control=list(fnscale=-1),lower = c(0.1,sigsq_lower,-1e+10,-1e+10,-1e+10),upper=c(1.9,100,1e+10,1e+10,1e+10),hessian=TRUE)
                   ress1=res1$par
                   extra_prob=pmvnorm(mean=ress, sigma =psi_hat1, lower=c(0.1,sigsq_lower,-1e+10,-1e+10,-1e+10), upper=c(1.9,100,1e+10,1e+10,1e+10), maxpts = 25000, abseps =0.001, releps = 0)
                   lap11= ((2*pi)^(5/2))*sqrt(det_psi_hat1)*extra_prob[1]
                   sigsq=ress1[2]
                   beta=ress1[3:5]
                   tau=ress1[1]
        	   log_post_hat1=res$value+ dinvgamma(sigsq, shape=2,scale=2, log =TRUE)+dmvnorm(beta, c(0,0,0),diag(c(10,10,10)), log =TRUE)+dunif(tau,0,2,log=TRUE)
        	   lap12=log_post_hat1
        	   #bf=(lap11/lap01)*(exp(lap12-lap02))
		   bf=(exp(log(lap11))/exp(log(lap01)))*(exp(lap12-lap02))
      		}
	}
kernel_weights=result1[1:3]
post_H1=1/(1+(prior0/prior1)*(1/bf))
result_final=c(post_H1,kernel_weights)
result_names=c("Post_H1","IBS_weight","Quadratic_weight","Gaussian_weight")
write.table(cbind(result_names,result_final),output_address)
return(list(Post_H1=post_H1,Kernel_Weights=kernel_weights))
}
}
			
