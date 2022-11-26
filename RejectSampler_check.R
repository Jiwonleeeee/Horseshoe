eta_update_ftn <- function(eta_j_input, k_input){
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
  check <- ifelse(exp(-f_L_x(eta_j_q, case))>=exp(-f_x(eta_j_q)),1,0)
  # print(check)
  acc_prob <- exp(-f_x(eta_j_q)+f_L_x(eta_j_q, case))
  return(eta_j_output)
}


k <- sapply(true_beta, function(x_input){(x_input^2 * 1) / (2 * true_sigma2)})+0.0001


### ex. beta= 2,3,4
y <- function(b_input,x_input){
  
  exp(-b_input*x_input)/(1+x_input)
  
}

x <- seq(0,5,0.01)

plot(x,y(2,x),type="l", ylab="eta",main="eta target function depeding on beta=2,3,4")
lines(x,y(3,x),col="red")
lines(x,y(4,x),col="blue")
legend(4,0.8, legend=c("2","3","4"), col=c("black", "red","blue"), lwd=1, cex=0.8)


k <- sapply(true_beta, function(x_input){(x_input^2 * 1) / (2 * true_sigma2)})+0.0001

## check
rep <- 10^4
eta_matrix <- matrix(0, nrow=rep, ncol=p)
for(j in 1:p){
  eta_matrix[1,j] <- 0.5
  for(i in 2:rep){
    eta_matrix[i,j] <- eta_update_ftn(eta_matrix[i-1,j],k[j])
  }
  print(j)
}
plot(apply(eta_matrix,2,mean), ylab="eta",main="eta check")
points(1:p, true_beta, col="red",cex=0.3)
legend(0, 2000, legend=c("estimated eta", "true beta"), col=c("black","red"),pch=1,cex=0.7)
