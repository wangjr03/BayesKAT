#this is the main function BayesKAT which can be used to call any of the strategies: BayesKAT-MCMC or BayesKAT-MAP

rm(list=ls())
source("../Sub_functions.R") 

main_function<-function(inputAddress_y,inputAddress_X,inputAddress_Z,prior_H1=0.5,output_address,Strategy="BayesKAT-MAP"){
  if(Strategy=="BayesKAT-MAP"){
       result=BayesKAT_MAP(inputAddress_y,inputAddress_X,inputAddress_Z,prior_H1=0.5,output_address)
  }else if(Strategy=="BayesKAT-MCMC"){
       result=BayesKAT_MCMC(inputAddress_y,inputAddress_X,inputAddress_Z,prior_H1=0.5,output_address)
  }else{
        result=c("Please choose one of the two strategies: BayesKAT-MAP or BayesKAT-MAP")
  }
  return(result)    
}

