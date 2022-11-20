## Scalable Approximate MCMC algorithms for the Horseshoe Prior, Johndrow et al. (2020)
library(mvtnorm)
library(dplyr)

# setting: z = W\beta + e

########### ########### ###########
N <- 10^3
p <- 10^4

########### data ###########
# true beta
beta_gen_ftn <- function(p_input){
  non_null_index <- c(1:23)
  beta <- c(2^(-non_null_index/4 + 9/4), rep(0,(p_input-23)))
  return(beta)
}

# generate z & W
data_gen_ftn <- function(N_input, p_input, Sigma_input){
  mat <- matrix(0, N_input, p_input)
  z <- numeric(N_input)
  for(i in 1:N_input){
    w_i <- rmvnorm(1, rep(0,p_input), Sigma_input)
    mat[i,] <- w_i
    z[i] <- rnorm(1, as.numeric(w_i %*% beta_gen_ftn(p_input)),4) # true sigma^2 = 4
  }
  
  return(list(W=mat,z=z))
}
########### ########### ###########





# exact algorithm
# parameters to be estimated: (eta, xi, sigma^2, beta)





