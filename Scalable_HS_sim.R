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

# eta update function- input: eta_j^t, output: eta_j^(t+1)
eta_update_ftn <- function(eta_j_input,epsilon_input){
  
  a <- 1/5
  b <- 10
  
  target_density <- function(eta_j_input){(exp(-epsilon_input*eta_j_input)/(1+eta_j_input))}
  f_x <- function(x_input){epsilon_input*x_input + log(1+x_input)}
  
  
  A <- f_x(a/epsilon_input)
  I <- f_x(1/epsilon_input)
  B <- f_x(b/epsilon_input)
  
  lambda_2 <- (I-A)/(1/epsilon_input-a/epsilon_input)
  lambda_3 <- (B-I)/(b/epsilon_input-1/epsilon_input)
  


  
  if(eta_j_input < a/epsilon_input){
    case <- 1
  }else if(eta_j_input >= a/epsilon_input & eta_j_input < 1/epsilon_input){
    case <- 2
  }else if(eta_j_input >= 1/epsilon_input & eta_j_input < b/epsilon_input){
    case <- 3
  }else{
    case <- 4
  }
  
 f_L_x <- function(x_input){
   
   value <-    switch (case,
                       log(1+x_input),
                       A + lambda_2*(x_input-a/epsilon_input),
                       I + lambda_3*(x_input-b/epsilon_input),
                       B + epsilon_input*(x_input-b/epsilon_input)
                      )

   return(value)
 }
  
 # normalizing constant
 nu_1 <- log(1+a/epsilon_input)
 nu_2 <- lambda_2^(-1) * exp(-A) * (1- exp(-I+A))
 nu_3 <- lambda_3^(-1) * exp(-I) *(1- exp(-B+I))
 nu_4 <- epsilon_input^(-1)*exp(-B)

 nu <- nu_1 + nu_2 + nu_3 + nu_4
 
 w1 <- nu_1/nu
 w2 <- nu_2/nu
 w3 <- nu_3/nu
 w4 <- nu_4/nu
 
 H2 <- 1-exp(-I+A)
 H3 <- 1-exp(-B+I)
 
 u <- runif(1)
 eta_j_q <- switch(case,
                   ((1+a/epsilon_input)^(u)-1) * w1,
                   (a/epsilon_input + (-log(1-u*H2)/lambda_2)) * w2,
                   (1/epsilon_input + (-log(1-u*H3)/lambda_3)) * w3,
                   (b/epsilon_input + (-log(1-u)/epsilon_input)) * w4
                   )
 eta_j_output <- ifelse(runif(1) < exp(-f_x(eta_j_q)+f_L_x(eta_j_q)), eta_j_q, eta_j_input)
 return(eta_j_output)

}

post_xi <- function(z_input, M_xi_input, w_input, N_input, xi_input){
  det(M_xi_input)^(-1/2) * (w_input/2 + t(z) %*% solve(M_xi_input) %*% z)^(-(N_input + w_input)/2)/(sqrt(xi_input)*(1+xi_input))
}


########### ########### ###########
iter <- 10^3

# exact algorithm
# parameters to be estimated: (eta(1~p), xi, sigma^2, beta(1~p))

########### save ###########
eta_store <- beta_store <- matrix(0, iter, p)
sigma2_store <- xi_store <- numeric(iter)

########### initial values ###########

eta <- runif(p)
beta <- runif(p)
sigma2 <- runif(1)
xi <- runif(1)

eta_store[1,] <- eta
beta_store[1,] <- beta
sigma2_store[1] <- sigma2
xi_store[1] <- xi


for(i in 2:iter){
  
  # eta update
  epsilon <- sapply(beta, function(x){(x^2 * xi) / (2 * sigma2)})
  eta <- sapply(1:p, function(s){eta_update_ftn(eta[s], epsilon[s])})
  
  
  # M_xi
  
  
  # xi update
  log_xi_q <- rnorm(log(xi), 1)
  xi_q <- exp(log_xi_q)
  
  
  xi <- ifelse( runif(1) < post_xi()    )
  
  
  
  
}







