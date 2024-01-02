# code for simulating genotype matrix from a given real data genotype.
# requires n-by-p real genotype matrix Z. The rows are individuals and columns are SNPs.  
# requires SNPknock R package: https://cran.r-project.org/web/packages/SNPknock/index.html
# requires fastPHASE software (must be downloaded separately by the user. Visit http://scheet.org/software.html for more information on how to obtain fastPHASE.)
# output: Zk is a n-by-p simulated genotype matrix containing the knockoff variables.

library(SNPknock)
Z=as.matrix(Z)
X_file=writeXtoInp(Z, phased = FALSE, out_file = NULL)
fp_path = "path_fastPHASE_installed"
out_path=runFastPhase(fp_path,X_file, K = 12, numit = 25,phased = FALSE, seed = 1)
hmm = loadHMM(paste0(out_path,"_rhat.txt"), paste0(out_path,"_alphahat.txt"), theta_file=paste0(out_path,"_thetahat.txt"), char_file=paste0(out_path,"_origchars"), compact=TRUE, phased=FALSE)
p=dim(Z)[2]
n=dim(Z)[1]
hmm$r = hmm$r[1:p]
hmm$alpha = hmm$alpha[1:p,]
hmm$theta = hmm$theta[1:p,]
Zk = knockoffGenotypes(Z, hmm$r, hmm$alpha, hmm$theta,seed)

####################################################################################################
# code for simulating datasets using simulated genotype matrix Zk.
# This code can be used to reproduce simulation settings in the "Simulation with discrete features" subsection.
# requires n-by-k covariate matrix given by user or simulated by user. Here k=3 used.

Z=Zk 
Z=as.matrix(Z) 
p<-dim(Z)[2]
n<-dim(Z)[1]
X=read.csv("User_defined_X_path.CSV") #the input_X data in Demo Data folder can be used
X=as.matrix(X)
beta1<-c(0.7,0.01,0.0008) 
fixed_part=X%*%beta1 
signal=0.4*(Z[,1]-Z[,3])+0.4*cos(Z[,3])*exp(-Z[,3]^2/5) #choose any scenario: D,E,F specific form
y=signal+fixed_part+rnorm(n,0,1) 
