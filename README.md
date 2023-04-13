# BayesKAT: A Bayesian Optimal Kernel-based Test for Genetic Association studies

## Summary
We developed the BayesKAT methodology to detect the presence of any genetic association between a trait and a group of preassigned genetic features. BayesKAT can implement a kernel-based association test that automatically selects the optimum kernel using the Bayesian approach and returns an interpretable result by incorporating prior knowledge based on biological experiments or known biological mechanisms.   

## Introduction
Multiple genetic variants often regulate complex diseases or phenotypes. In that case, an individual variant might be weakly associated with a phenotype, but a group of SNPs with a weak marginal signal might jointly regulate some crucial biological process. In that situation, set-based tests using the kernel-based testing (KBT) framework have proven helpful in practice. KBT framework jointly models the group effect of genetic variants controlling for other demographic covariate effects on a phenotype. It measures the similarity between genetic variants through a kernel function and then compares it with the phenotype similarity. Generally, a specific kernel function would be optimum for specific association types between the genetic variants and the phenotype. However, knowing the true functional relationship between the variants and phenotype is impossible, so choosing the optimum kernel beforehand is also tricky. Our method, "BayesKAT," overcomes this kernel specification issue by selecting the optimum kernel based on the dataset while testing for the association. There are two alternative strategies (BayesKAT-MCMC and BayesKAT-MAP) to implement the BayesKAT methodology. BayesKAT-MCMC is a MCMC (Markov Chain Monte Carlo) sampling-based strategy and BayesKAT-MAP is a more flexible and fast optimization based strategy. Here we provide the code for each of the strategy as functions which takes input datasets and return the posterior probability of association. 

## Dependencies
The implementation of the algorithm is based on R. The code for BayesKAT_MCMC depends on BayesianTools. the strategy BayesKAT_MAP depends on these packages: Matrix, MASS, LaplacesDemon, mvtnorm, numDeriv, minqa.

## Input Data
The main function of BayesKAT-MAP or BayesKAT-MCMC can be utilized as follow:
```
BayesKAT_MAP<-function(
              inputAddress_y,
              inputAddress_X,
              inputAddress_Z,
              prior_H1 =0.5,
              output_address
              )
              
BayesKAT_MCMC<-function(
              inputAddress_y,
              inputAddress_X,
              inputAddress_Z,
              prior_H1 =0.5,
              output_address
              )              
```

Arguments:

inputAddress_y: The address for the input data y containing individual level trait or phenotype value. y should be a nx1 matrix where n is the number of individuals.

inputAddress_X: The address for the input data X containing individual level demographic variables. X should be a nxp matrix where n is the number of individuals and p is the number of covariates. Here p=3 for this example function. 

inputAddress_Z: The address for the input data Z containing individual level genetic data. Z should be a nxp matrix where n is the number of individuals and p is the number of genetic features.

prior_H1: The prior probability of alternative hypothesis, i.e the prior belief of association. If there is no prior information, the default value is 0.5.

output_address: The addres where the output file should be saved.

### Demo data:
We have attached example input datasets in demo data folder. The Z matrix is simulated from real individual level genotype data for a group of biologically related variants. X is individual level covariate dataset and y is generated from X and Z using this fuction: 
$y= 2 \times (Z[,1] \times Z[,3]) + X\beta + \epsilon$, where $\epsilon \sim N(0,1)$ and $\beta= c(0.7,0.01,0.0008)$.


## Output data format
The output file contains these informations:

(1) Posterior probability of H1,
(2) Kernel weights in optimum kernel.

## Remark:
The example functions BayesKAT_MCMC and BayesKAT_MAP are given here for demonstration purposes only; they have used only these three kernels: IBS, Quadratic, and Gaussian as candidate kernels, and only three covariates in X. User can change the candidate kernels or the number of parameters inside the function in case that is required. 
