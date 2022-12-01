## Scalable Approximate MCMC algorithms for the Horseshoe Prior, Johndrow et al. (2020)
rm(list=ls())
#install.packages("mvtnorm",repos = "http://cran.us.r-project.org")
#install.packages("dplyr",repos = "http://cran.us.r-project.org")

library(mvtnorm)
library(dplyr)
library(invgamma)
library(truncnorm)
# setting: z = W\beta + e

########### ########### ###########
N <- 100
p <- 200
true_sigma2 <- 4
omega <- 0.001 # noninformative?

########### data ###########
# true beta
beta_gen_ftn <- function(p_input){
  non_null_index <- c(1:23)
  beta <- c(2^(-non_null_index/4 + 9/4)*2, rep(0,(p_input-length(non_null_index))))
  return(beta)
}

true_beta <- beta_gen_ftn(p)

# true_beta <- numeric(p)
# non_null <- 5
# index <- sample(1:p, non_null)
# true_beta[index] <- 4

# generate z & W
data_gen_ftn <- function(N_input, p_input, Sigma_input){
  mat <- matrix(0, N_input, p_input)
  z <- numeric(N_input)
  for(i in 1:N_input){
    w_i <- rmvnorm(1, rep(0,p_input), Sigma_input)
    mat[i,] <- w_i
    z[i] <- rnorm(1, as.numeric(w_i %*% true_beta),2) # true sigma^2 = 4
  }
  
  return(list(W=mat,z=z))
  
}

## define W and z, using identity matrix
Sigma <- matrix(0, p, p)
for(i in 1:p){
  for(j in 1:p){
    Sigma[i,j] <- 0.9^(abs(i-j))
  }
}

W <- data_gen_ftn(N, p, Sigma)$W
z <- data_gen_ftn(N, p, Sigma)$z

# eta update function- input: eta_j^t, output: eta_j^(t+1)
eta_update_ftn <- function(eta_j_input,k_input){
  f_x <- function(x_input){k_input*x_input + log(1+x_input)}
  
  a <- 1/5
  b <- 10^4
  c <- 10^5
  
  A <- f_x(a)
  B <- f_x(b)
  C <- f_x(c)
  
  lambda2 <- (B-A)/(b-a)
  lambda3 <- (C-B)/(c-b)
  
  
  f_L_x <- function(x_input, case_input){
    
    val <- switch(case_input,
                  log(1+x_input),
                  A + lambda2 * (x_input-a),
                  B + lambda3 * (x_input-b),
                  C + k_input * (x_input-c)
    )
    
    return(val)
  }
  
  # normalizing constants
  
  v1 <- log(1+a)
  v2 <- lambda2^(-1) * exp(-A) * (1-exp(-B+A))
  v3 <- lambda3^(-1) * exp(-B) *(1-exp(-C+B))
  v4 <- k_input^(-1) * exp(-C)
  
  v <- v1 + v2 + v3 + v4
  
  p1 <- v1/v
  p2 <- v2/v
  p3 <- v3/v
  p4 <- v4/v
  
  # print(c(p1,p2,p3,p4))
  
  H2 <- 1-exp(-B+A)
  H3 <- 1-exp(-C+B)
  
  u <- runif(1)
  
  case <- sample(c(1:4), 1, prob=c(p1,p2,p3,p4))
  
  eta_j_q <- switch(case,
                    (1+a)^u -1,
                    a - log(1-u*H2)/lambda2,
                    b - log(1-u*H3)/lambda3,
                    c - log(1-u)/k_input
  )
  u <- runif(1)
  eta_j_output <- ifelse(u < exp(-f_x(eta_j_q)+f_L_x(eta_j_q, case)), eta_j_q, eta_j_input)
  return(eta_j_output)
  
}

log_post_xi <- function(z_input, M_xi_input, omega_input, N_input, xi_input){
  as.numeric(determinant(solve(M_xi_input), logarithm = T)$mod)/2 - ((N_input + omega_input)/2) * log(omega_input/2 + (t(z_input) %*% solve(M_xi_input) %*% z_input)/2)-log(sqrt(xi_input)*(1+xi_input))
}

beta_update_ftn <- function(xi_input, eta_input, N_input, z_input, p_input, W_input, M_xi_inv_input, sigma2_input){
  
  u_vector <- t(rmvnorm(1, rep(0, p_input), diag(1/eta_input)/xi_input))
  f_vector <- t(rmvnorm(1, rep(0, N_input), diag(N_input)))
  v_vector <- W_input %*% u_vector + f_vector
  v_star <- M_xi_inv_input %*% (z_input/sqrt(sigma2_input)-v_vector)
  
  beta_vector <- sqrt(sigma2_input) * (u_vector + (diag(1/eta_input) %*% t(W_input) %*% v_star) / xi_input)
  
  return(beta_vector)
  
}




########### ########### ###########
iter <- 10^4

# exact algorithm
# parameters to be estimated: (eta(1~p), xi, sigma^2, beta(1~p))

########### save ###########
eta_store <- beta_store <- matrix(0, iter, p)
sigma2_store <- xi_store <- xi_accept <- numeric(iter)

########### initial values ###########

eta <- runif(p)
beta <- runif(p,0,4)
sigma2 <- runif(1)
xi <- runif(1)

eta_store[1,] <- eta
beta_store[1,] <- beta
sigma2_store[1] <- sigma2
xi_store[1] <- xi

current_time <- Sys.time()
for(i in 2:iter){
  
  # eta update
  k <- sapply(beta, function(x_input){(x_input^2 * xi) / (2 * sigma2)})
  eta <- sapply(1:p, function(s_input){eta_update_ftn(eta[s_input], k[s_input])})
  
  
  # M_xi (needs to be updated after eta and xi updates)
  M_xi <- diag(N) + ( W %*% diag(1/eta) %*% t(W) ) / xi
  
  # xi update

  # truncated version
  # xi_q <- rtruncnorm(1, a=0, b=Inf, mean=xi, sd=prop_var)
  # M_xi_q <- diag(N) + ( W %*% diag(1/eta) %*% t(W) ) / xi_q
  
  # log version
  prop_var <- 5 # sd
  log.xi_q <- rnorm(1, log(xi), prop_var)
  xi_q <- exp(log.xi_q)
  M_xi_q <- diag(N) + ( W %*% diag(1/eta) %*% t(W) ) / xi_q
  
  #log_acc_prob <- log_post_xi(z, M_xi_q, omega, N, xi_q) + dtruncnorm(xi,0,Inf,xi_q,prop_var) - log_post_xi(z, M_xi, omega, N, xi) - dtruncnorm(xi_q,0,Inf,xi,prop_var)
  log_acc_prob <- log_post_xi(z, M_xi_q, omega, N, xi_q) + log(xi_q) - log_post_xi(z, M_xi, omega, N, xi) - log(xi)
  
  if(log(runif(1)) < log_acc_prob){
    xi <- xi_q
    M_xi <- M_xi_q
    xi_accept[i] <- 1
  }
  M_xi_inv <- solve(M_xi)
  
  # sigma2 update
  
  sigma2 <- 1/rgamma(1, shape= (omega + N)/2, rate = 2/(omega + t(z) %*% M_xi_inv %*%z))

  
  # true sigma2
  # sigma2 <- 4
  
  # beta update
  # true beta
  # beta <- true_beta + 0.0001
  
  beta <- beta_update_ftn(xi, eta, N, z, p, W, M_xi_inv, sigma2)
  
  xi_store[i] <- xi
  sigma2_store[i] <- sigma2
  beta_store[i,] <- beta
  eta_store[i,] <- eta
  
  print(i)
  print(beta[1:5])
  
}

end_time <- Sys.time()
end_time - current_time

#save(sigma2_store, beta_store, file="parameters.rda")

apply(beta_store,2,mean)[1:p]
# plot(eta_store[,1],type="l")
# plot(xi_store[1:2000], ylim=c(0,10^2))

## example
# neg_log_p_eta <- function(eta){
#   k <- 16*4/(2*4)
#   return(k*eta+log(1+eta))
# }
# plot(seq(0,10,0.1),neg_log_p_eta(seq(0,10,0.1)),type="l",xlim=c(0,1),ylim=c(0,20))
plot(xi_store,type="l")
plot(sigma2_store,type="l")
# plot(eta_store[,1],type="l")
# plot(eta_store[,2],type="l")
# plot(eta_store[,3],type="l")
# plot(eta_store[,10],type="l") # corresponding to null beta
# range(eta_store[,1])
# range(eta_store[,10])

# plot(sigma2_store,type="l",ylim=c(0,10))
mean(sigma2_store[10001:iter])

apply(eta_store,2,mean) # true beta = 0 when p >5, eta[6:p] must be greater than eta[1:5] -> check
mean(xi_store[10001:iter])
plot(xi_store,type="l")
# plot(sigma2_store,type="l",ylim=c(0,0.00001))
mean(xi_accept)


post <- apply(beta_store[5001:iter,],2,mean)
plot(post,ylim=c(0,10))
points(c(1:p),true_beta,col="red")



plot(apply(eta_store[5001:iter,],2,mean))
points(which(true_beta!=0),true_beta[true_beta!=0],col="red")

# save(beta_store, sigma2_store,xi_store, eta_store, file="beta1.rda")
