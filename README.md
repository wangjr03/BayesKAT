# BayesKAT: A Bayesian Optimal Kernel-based Test for Genetic Association studies

## Summary
We developed the BayesKAT methodology to detect the presence of any genetic association between a trait and a group of preassigned genetic features. BayesKAT can implement a kernel-based association test that automatically selects the optimum kernel using the Bayesian approach and returns an interpretable result by incorporating prior knowledge based on biological experiments or known biological mechanisms.   

## Introduction
Multiple genetic variants often regulate complex diseases or phenotypes. In that case, an individual variant might be weakly associated with a phenotype, but a group of SNPs with a weak marginal signal might jointly regulate some crucial biological process. In that situation, set-based tests using the kernel-based testing (KBT) framework have proven helpful in practice. KBT framework jointly models the group effect of genetic variants controlling for other demographic covariate effects on a phenotype. It measures the similarity between genetic variants through a kernel function and then compares it with the phenotype similarity. Generally, a specific kernel function would be optimum for specific association types between the genetic variants and the phenotype. However, knowing the true functional relationship between the variants and phenotype is impossible, so choosing the optimum kernel beforehand is also tricky. Our method, "BayesKAT," overcomes this kernel specification issue by selecting the optimum kernel based on the dataset while testing for the association. There are two alternative strategies (BayesKAT-MCMC and BayesKAT-MAP) to implement the BayesKAT methodology. BayesKAT-MCMC is a MCMC (Markov Chain Monte Carlo) sampling-based strategy and BayesKAT-MAP is a more flexible and fast optimization based strategy. Here we provide the code for each of the strategy as functions which takes input datasets and return the posterior probability of association. 

## Dependencies
The implementation of the algorithm is based on R. The code for BayesKAT_MCMC depends on BayesianTools. the strategy BayesKAT_MAP depends on these packages: Matrix, MASS, LaplacesDemon, mvtnorm, numDeriv, minqa.

## Input Data

The main function of BayesKAT-MAP can be utilized as follow:
```
BayesKAT_MAP<-function(inputAddress_y,
                        inputAddress_X,
                        inputAddress_Z,
                        prior_H1,
                        output_address
                        )

```

Arguments:

inputAddress_y:

inputAddress_X:

inputAddress_Z:

prior_H1:

output_address:
