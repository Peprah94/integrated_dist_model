#Simulation of the data with replication at sites

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
#p.tag is the detectiion probability for the occurence model (species data)
#beta is a n.species* nsites matrix
#includecov (T/F) indicates whether it is a model with a covariate (TRUE) or no covariate (FALSE)
##############################################################################################
set.seed(2020)
myclog <- function(psi){
  return(-log(1-psi))
}

myinvclog <- function(mu){
  return(1-exp(-exp(mu)))
}

sim <- function(n.sites, n.species,n.visit, n.id, n.gen, cov, includecov){
  N <-log_lambda <- vector("numeric", n.sites)
  z <-mu <- lambda.s <- epsilon <- matrix(NA, nrow = n.sites, ncol = n.species)
 x <- array(NA, dim = c(n.sites, n.species, n.visit))
  y <- matrix(NA, nrow=n.sites, ncol=n.visit)
  
  if(includecov== TRUE){
    #covariate effect
    beta.cov = matrix(rnorm(n.species, 1,0.5), nrow=n.species, ncol= 3, byrow=TRUE) 
    mu <- cov %*% t(beta.cov) #linear predictor 
  }else{
    alpha <- rnorm(1,0, sd= 0.4) #intercept
    sigma.species <- 2; sigma.site=1 #variance for species and site effect
    beta.site <- rnorm(n.sites, 0, sd=sqrt(sigma.site)) # site effect
    beta.species <- rnorm(n.species,0,sd=sqrt(sigma.species)) #species effect
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        mu[site.tag, spe.tag] <- alpha + beta.site[site.tag] + beta.species[spe.tag] #linear predictor
      }
    }
  }
  lambda.s <- exp(mu) #mean abundance
  psi.s <- 1-exp(-lambda.s) # estimation of occupancy probability
  lambda.g <- rowSums(lambda.s) #mean abundance for genus
  p.tag<- 0.7 #detection probability.


for(site.tag in 1:n.sites){
  for(spe.tag in 1:n.species){
    #True presence or absence
      z[site.tag, spe.tag] <- rbinom(1,1, psi.s[site.tag, spe.tag]) 
      for(k in 1:n.visit){
        # Genus counts
        y[site.tag,k] <- rpois(1, lambda.g[site.tag])
        # Species Presence and absence observations
        x[site.tag, spe.tag,k] <- rbinom(1,5, z[site.tag, spe.tag]*p.tag)
  }
}
}

  #Shannon Index
    shan.index <- -rowSums(log(lambda.s/lambda.g)*(lambda.s/lambda.g), na.rm = TRUE)   

 # Returning the results
  data <- list(mat.species=x, mat.genus = y, covariates = cov, shan.index=shan.index) 
return(data)
}

# Simulating of covariates.
covariates <- function(n.sites){
  elev.s <- runif(n.sites, -1,1) # The covariate effect (elavation)
  rain.s <- runif(n.sites, -1,1) # The covariate effect (rainfall)
  covte <- cbind(1,elev.s,rain.s)
  return(covte)
}

