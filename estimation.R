library(nimble)
library(coda)

est <- function(data, code){
  dim_data <- dim(data[[1]]) #dimensions of the data
  data_dim <- dim_data[2] 
  
  # Constants for the model 
  const <- list(n.sites = dim_data[1], n.species= (dim_data[2]), n.replicates = 5,
                n.visit=dim_data[3], n.cov=(dim(data[[3]]))[2])
  
  #Data for the model
  idm_data <- list(y = data[[2]], x = data[[1]], cov= data[[3]])
   
  #Retrieving the true Presence/Absence
   zst <- apply(idm_data$x, c(1,2), max)
   zst[zst!= 0] <- 1
  if(includecov){
    idm_inits <- function(){list(beta= matrix(rnorm((const$n.cov *const$n.species), 0, 1), 
                                      nrow=const$n.species, ncol=const$n.cov, byrow = TRUE), #covariate effect
                                 p.tag = runif(1,0,1), #detection probability
                                 z=zst) #True Presence/Absence
      }  
  }else{
    idm_inits <- function(){list(z=zst, #True Presence/Absence
                                 sigma.site = 2, 
                                 sigma.species = 1,
                                sigma.alpha=0.4,
                                 p.tag = runif(1,0,1))}
  }
  initsList <- idm_inits()
  
  modelInfo <- list(
    code = code,
    constants = const,
    data = idm_data,
    inits = idm_inits
  )
  
  #############################################################
  #                 Compile Model                             #
  #############################################################
  
  mwtc <- nimbleModel(code,data = idm_data, constants = const, inits = initsList)
  Cmwtc <- compileNimble(mwtc)
  
  #Configuring the MCMC
  if(includecov){
  mcmcconf <- configureMCMC(Cmwtc, monitors = c("p.tag", "shan.index", "beta"), autoBlock = TRUE)
  }else{
    mcmcconf <- configureMCMC(Cmwtc,monitors = c("p.tag","shan.index", "sigma.alpha", "sigma.site", "sigma.species"))
    }
  Rmcmc <- buildMCMC(mcmcconf, enableWAIC =TRUE)
  cmcmc <- compileNimble(Rmcmc, project = Cmwtc,resetFunctions = TRUE)

  #############################################################
  #                 Run the MCMC                              #
  #############################################################
  mcmc.out <- runMCMC(cmcmc, niter = 500000,nchains = 2,nburnin = 250000,inits = initsList, 
              setSeed = TRUE, thin=10,samples=TRUE, samplesAsCodaMCMC = TRUE, summary = TRUE, WAIC = TRUE)
  output <- mcmc.out$summary$all.chains
  coda_samples <- mcmc(mcmc.out$samples)
  mode <- posterior.mode(mcmc(mcmc.out$samples$chain2))
  ret <- list(mcmc.out= mcmc.out, mode.out = mode, mean.out = output[,1], median.out = output[,2], sd = output[,3], lower = output[,4],upper = output[,5], modelInfo= modelInfo, compile.model = mwtc)
return(ret)
}

