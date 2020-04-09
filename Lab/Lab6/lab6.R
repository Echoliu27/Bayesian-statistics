## Exercise 1: you should have three univariate normal distributions
# p(x|y,z)
mux_1 = 0
mux_2 = c(0,0)
sigmax_11 = 1
sigmax_22 = matrix(c(1,0.1,0.1,1),nrow=2,ncol=2)
sigmax_12 = sigmax_21 = c(0.9,0.1)

mux_1 + sigmax_12%*%solve(sigmax_22)%*%( - mux_2)
sigmax_11 - sigmax_12%*%solve(sigmax_22)%*%sigmax_21

# p(y|x,z)
muy_1 = 0
muy_2 = c(0,0)
sigmay_11 = 1
sigmay_22 = matrix(c(1,0.1,0.1,1),nrow=2,ncol=2)
sigmay_12 = sigmay_21 = c(0.9,0.1)


# p(z|x,y)
muz_1 = 0
muz_2 = c(0,0)
sigmaz_11 = 1
sigmaz_22 = matrix(c(1,0.9,0.9,1),nrow=2,ncol=2)
sigmaz_12 = sigmaz_21 = c(0.1,0.1)


## Exercise 2: write a gibbs sampler that alternates updating each of the variables
n_iter <- 1000; burn_in <- 0.3*n_iter
set.seed(1234)

# innitialize the values
X_n = Y_n = Z_n = 0
X = Y = Z = NULL

for (s in 1:(n_iter+burn_in)){
  #update X_n
  mux = mux_1 + sigmax_12%*%solve(sigmax_22)%*%(c(Y_n,Z_n) - mux_2)
  varx = sigmax_11 - sigmax_12%*%solve(sigmax_22)%*%sigmax_21
  X_n <- rnorm(1, mean = mux, sd= varx)
  
  #update Y_n
  muy = muy_1 + sigmay_12%*%solve(sigmay_22)%*%(c(X_n,Z_n) - muy_2)
  vary = sigmay_11 - sigmay_12%*%solve(sigmay_22)%*%sigmay_21
  Y_n <- rnorm(1, mean = muy, sd= vary)
  
  #update Z_n
  muz = muz_1 + sigmaz_12%*%solve(sigmaz_22)%*%(c(X_n,Y_n) - muz_2)
  varz = sigmaz_11 - sigmaz_12%*%solve(sigmaz_22)%*%sigmaz_21
  Z_n <- rnorm(1, mean = muz, sd= varz)
  
  #save results only past burn-in
  if(s > burn_in){
    X <- rbind(X,X_n)
    Y <- rbind(Y,Y_n)
    Z <- rbind(Z,Z_n)
  }
}

library(mvtnorm)
library(MCMCpack)
X.mcmc <- mcmc(X,start=1); 
summary(X.mcmc)

Y.mcmc <- mcmc(Y,start=1); 
summary(Y.mcmc)

Z.mcmc <- mcmc(X,start=1); 
summary(Z.mcmc)

plot(X.mcmc)
plot(Y.mcmc)
plot(Z.mcmc)

autocorr.plot(X.mcmc)
# comment on the plots

## Exercise 3: Block updates
# p((x,y)|z)
muxy_1 = c(0,0)
muxy_2 = 0
sigmaxy_11 = matrix(c(1,0.9,0.9,1),nrow=2,ncol=2)
sigmaxy_22 = 1
sigmaxy_12 = sigmaxy_21 = c(0.1,0.1)

#p(z|(x,y) is given above

# Gibbs Sampler: one bivaraite normal and one univariate normal
n_iter <- 1000; burn_in <- 0.3*n_iter
set.seed(1234)

# innitialize the values
Z_n = 0
XY = Z = NULL

for (s in 1:(n_iter+burn_in)){
  #update X_n, Y_n
  muxy = muxy_1 + sigmaxy_12%*%solve(sigmaxy_22)%*%(c(Z_n) - muxy_2)
  varxy = sigmaxy_11 - sigmaxy_12%*%solve(sigmaxy_22)%*%sigmaxy_21
  XY_n <- rmvnorm(1, mean = muxy, sigma= varxy)
  
  #update Z_n
  muz = muz_1 + sigmaz_12%*%solve(sigmaz_22)%*%(c(XY_n) - muz_2)
  varz = sigmaz_11 - sigmaz_12%*%solve(sigmaz_22)%*%sigmaz_21
  Z_n <- rnorm(1, mean = muz, sd= varz)
  
  #save results only past burn-in
  if(s > burn_in){
    XY <- rbind(XY,XY_n)
    Z <- rbind(Z,Z_n)
  }
}

colnames(XY) <- c("X","Y")

# diagnostics
XY.mcmc <- mcmc(XY,start=1); 
summary(X.mcmc)

Z.mcmc <- mcmc(X,start=1); 
summary(Z.mcmc)

plot(XY.mcmc[,'X'])
plot(XY.mcmc[,'Y'])
plot(Z.mcmc)

autocorr.plot(XY.mcmc[,'X'])

## what's the diference between the performance of two Gibbs samplers?
