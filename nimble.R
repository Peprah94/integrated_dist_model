###########################################
#   NIMBLE CODE
###########################################

library(nimble)
library(coda)
library(mcmcplots)

##############################################################################################
#                         NOTES                                                              #
##############################################################################################
### We assume n.genus  g = 1,2,...,n.genus
# For each g genera, we assume n.species in that geunus
# We assume the data was collected at n.sites
#We also assume the total number of visits to each site is n.visit.
# We also assume that some of the data are collected with species ID
# Most of them are collected with genus ID
#n.id = number of sites that have ID at the species level
#n.gen = number of sites that have ID at the genus level
#shan.index is the shannon's index
#gamma is the detection probability for the abundance model (genus data)
#p.tag is the detection probability for the occurrence model (species data)
#beta is a n.species* nsites matrix
#includecov (T/F) indicates whether it is a model with a covariate (TRUE) or no covariate (FALSE)
##############################################################################################

code <- nimbleCode({
  #############################################################
  #                 PRIOR DISTRIBUTIONS                       #
  #############################################################
  
  if(includecov){
    #set-up priors for beta
    for(spe.tag in 1:n.species){
      for(pred.tag in 1:n.cov){
        #Covariate effects
        beta[spe.tag, pred.tag] ~ dnorm(0, sd=100)
      }
    }
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        #Linear predictor
        mu[site.tag, spe.tag] <- beta[spe.tag, 1] + beta[spe.tag, 2]* cov[site.tag, 2] + beta[spe.tag, 3]* cov[site.tag, 3]
      }
    }   
  }else{
    sigma.site ~ dgamma(1, 1) # site specific variance
    sigma.species ~ dgamma(1, 1) #species specific variance
      sigma.alpha ~ dgamma(1, 1) #intercept variance
      
      #Site-specific effect
    for(site.tag in 1:n.sites){
      beta.site[site.tag] ~ dnorm(0, var=sigma.site)
    }
      #Species-specific effect
    for(spe.tag in 1:n.species){
      beta.species[spe.tag] ~ dnorm(0, var=sigma.species)
    }
      #Intercept
    alpha ~ dnorm(0, var=sigma.alpha)
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        mu[site.tag,spe.tag] <- alpha +  beta.site[site.tag] + beta.species[spe.tag]
      }
    } 
  }
  
  #Link between the abundance and occupancy
  for(site.tag in 1:n.sites){
    for(spe.tag in 1:n.species){
      log(lambda[site.tag, spe.tag]) <-  mu[site.tag,spe.tag]
      cloglog(psi[site.tag, spe.tag]) <- log(lambda[site.tag, spe.tag])
    }
  } 
  
  #Detection probability
  p.tag ~ dunif(0.001,1) 

  # Likelihood: key definitions in the likelihood
  #############################################################
  #                 Abundance Model                           #
  #############################################################
  
  for(site.tag in 1:n.sites){
    lambda.g[site.tag] <- sum(lambda[site.tag,1:n.species])
    for(k in 1:n.visit){
      y[site.tag,k] ~ dpois(lambda.g[site.tag])
    }
  }
  
  #############################################################
  #                 Occurence Model                           #
  #############################################################
  for(site.tag in 1:n.sites){
    for(spe.tag in 1:n.species){
      z[site.tag,spe.tag] ~ dbern(psi[site.tag, spe.tag])
      for(k in 1:n.visit){
        x[site.tag,spe.tag,k] ~ dbin(z[site.tag,spe.tag]*p.tag, n.replicates)
      }
    }
  }
  
  #############################################################
  #                 Shannon Index                             #
  #############################################################
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        shan[site.tag, spe.tag] <- (lambda[site.tag,spe.tag]/lambda.g[site.tag]) * log(lambda[site.tag,spe.tag]/lambda.g[site.tag])
      }
      shan.index[site.tag] <-  -sum(shan[site.tag,1:n.species])
    }    
})