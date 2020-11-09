#estimates the shannon index from the mcmc

source('sim_rep.R')
source('nimble.R')
source('estimation.R')

shanest <- function(beta.cov, cov, includecov){
if(includecov== TRUE){
lambda.s <- exp(cov%*%t(beta.cov))
}else{
lambda.s <- beta.cov
}
lambda.g <- rowSums(lambda.s)

shan.index <- -rowSums(log(lambda.s/lambda.g)*(lambda.s/lambda.g), na.rm = TRUE)
return(shan.index)
}

