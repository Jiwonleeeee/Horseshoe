## generate beta
# W and z are the same
true_eta <- runif(p,0,10)
true_xi <- runif(1,0,10)
mean_vector <- solve(t(W)%*%W + true_xi * solve(diag(1/true_eta))) %*% t(W) %*% z
Sigma_matrix <- solve(t(W)%*%W + true_xi * solve(diag(1/true_eta))) * true_sigma2
true_M_xi <- diag(N) + ( W %*% diag(1/true_eta) %*% t(W) ) / true_xi


temp <- rmvnorm(10^4, mean_vector, Sigma_matrix)
true_beta <- apply(temp,2,mean)



beta_update_ftn <- function(xi_input, eta_input, N_input, z_input, p_input, W_input, M_xi_inv_input, sigma2_input){
  
  u_vector <- t(rmvnorm(1, rep(0, p_input), diag(1/eta_input)/xi_input))
  f_vector <- t(rmvnorm(1, rep(0, N_input), diag(N_input)))
  v_vector <- W_input %*% u_vector + f_vector
  v_star <- M_xi_inv_input %*% (z_input/sqrt(sigma2_input)-v_vector)
  
  beta_vector <- sqrt(sigma2_input) * (u_vector + (diag(1/eta_input) %*% t(W_input) %*% v_star) / xi_input)
  
  return(beta_vector)
  
}

beta_matrix <- matrix(0, 10^4, p)
for(i in 1:10^4){
  beta_matrix[i,] <- beta_update_ftn(true_xi, true_eta, N, z, p, W, solve(true_M_xi),true_sigma2)
}


est_beta <- apply(beta_matrix,2,mean)
plot(true_beta)
points(c(1:p),est_beta,col="red",pch=2,cex=0.3)

# beta update function has no problem


