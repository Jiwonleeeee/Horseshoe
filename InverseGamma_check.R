IG_pdf <- function(alpha_input, beta_input, sigma2_input){
  beta_input^alpha_input/factorial(alpha_input-1) * sigma2_input^(-alpha_input-1) * exp(-beta_input/sigma2_input)
}

alpha = 3
beta = 4

# mean = 2

sigma2 = seq(0,8,0.001)
n = length(sigma2)
density = IG_pdf(alpha, beta, sigma2)

plot(sigma2, density, type="l")


Sample <- 1/rgamma(n, shape=alpha, rate=beta)
mean(Sample) # 2

