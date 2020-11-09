#This puts together the simulation, running MCMC and estimating the shannon index.
#Returns the parameters from the all_data() function created

library(nimble)
library(coda)
library(mcmcplots)
library(MCMCglmm)
#library(compareMCMCs)

source('sim_rep.R')
source('nimble.R')
source('estimation.R')
source('shanestimate.R')


#############################################################
 #                 Data Simulation                           #
#############################################################
set.seed(2020)
includecov= FALSE
n.sites=20
n.species = 5
n.visit=3


all_data <- function(n.species, n.visit, n.sites){
covte = covariates(n.sites)
data = sim(n.sites, n.species,n.visit,n.sites ,n.sites, cov=covte, includecov)

#############################################################
#                 Diagnostics Plots                         #
#############################################################
mcmc.out <- est(data, code)
 coda_samples <- mcmc(mcmc.out$mcmc.out$samples)
mcmcplot(coda_samples)

true_shan <- data$shan.index
if(includecov==TRUE){
  est_beta <- matrix(mcmc.out$mode.out[1:(n.species*3)],
                     nrow=n.species, ncol=3)
  est_shan = shanest(est_beta,data$covariates, TRUE)
  mse_shan = sum(((true_shan-est_shan)^2)/n.sites)
  mbias_shan = mean(true_shan-est_shan)
   p_est <- mcmc.out$mode.out[(n.species*3)+1]
   sig_alpha <- NA
   sig_beta <- NA
   sig_gamma <- NA
  }else{
    est_shan <- mcmc.out$mode.out[2:(n.sites+1)]
    mse_shan = sum(((true_shan-est_shan)^2)/n.sites)
    mbias_shan = mean(true_shan-est_shan)
  p_est <- mcmc.out$mode.out[1]
  sig_alpha <- tail(mcmc.out$mode.out,3)[1]
  sig_beta <- tail(mcmc.out$mode.out,3)[2]
  sig_gamma <- tail(mcmc.out$mode.out,3)[3]
  }
data = data.frame(p_est, mse_shan, mbias_shan, sig_alpha, sig_beta, sig_gamma)
return(data)
}

all_data(5,3,20)


