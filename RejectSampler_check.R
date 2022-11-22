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
  acc_prob <- exp(-f_x(eta_j_q)+f_L_x(eta_j_q, case))
  return(acc_prob)
}

eta_update_ftn(0.5,1)

xx <- seq(0,100,0.01)
l <- length(xx)
check_v <- numeric(i)
for(i in 1:l){
  check_v[i] <- eta_update_ftn(xx[i],1)
}
hist(check_v)
