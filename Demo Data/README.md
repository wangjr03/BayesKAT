# Demo Data
We have attached example input datasets in this folder. The Z matrix is simulated from real individual level genotype data for a group of biologically related variants. X is individual level covariate dataset and y is generated from X and Z using this fuction: 
$y= 0.5 \times (Z[,1] \times Z[,3]) + 0.6 \times Z[,1]^2 + Z[,3]+ X\beta + \epsilon$, where $\epsilon \sim N(0,1)$ and $\beta= c(0.7,0.01,0.0008)$.

Here is more information on the moin functions given in Code folder.

## Input Data
The main function for BayesKAT can be utilized as follow:
```
Main_function<-function(
              inputAddress_y,
              inputAddress_X,
              inputAddress_Z,
              prior_H1 =0.5,
              output_address,
              Strategy="BayesKAT-MAP"
              )
```

The individual function for BayesKAT-MAP or BayesKAT-MCMC can be utilized as follow:
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

inputAddress_y: The address for the input data y containing individual-level trait or phenotype value. y should be a nx1 matrix where n is the number of individuals.

inputAddress_X: The address for the input data X containing individual-level demographic variables. X should be a nxp matrix where n is the number of individuals and p is the number of covariates. Here p=3 for this example function. 

inputAddress_Z: The address for the input data Z containing individual-level genetic data. Z should be a nxp matrix where n is the number of individuals and p is the number of genetic features.

prior_H1: The prior probability of the alternative hypothesis, i.e., the prior belief of association. If there is no prior information, the default value is 0.5.

output_address: The address where the output file should be saved.

Strategy: The strategy("BayesKAT-MAP" or "BayesKAT-MCMC") would be used for Bayes factor calculation. By default, BayesKAT-MAP will be used if no strategy is chosen. 


## Expected output
The output file contains these pieces of information:

(1) Posterior probability of H1,
(2) Kernel weights in the optimum kernel.

## Remark:
1. The example functions BayesKAT_MCMC and BayesKAT_MAP are given here for demonstration purposes only; they have used only these three kernels: IBS, Quadratic, and Gaussian as candidate kernels, and only three covariates in X. User can change the candidate kernels or the number of parameters inside the function in case that is required. 

2. In order to use the Main_function, at first, the Sub_functions.R script should be saved on your computer so that it can be called using Main_function.
