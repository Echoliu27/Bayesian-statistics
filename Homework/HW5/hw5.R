library(mvtnorm)
library(MCMCpack)
require(tidyverse)

## Question 1

# (a) 
# Simulate data
n = 100
set.seed(1234)
theta = c(0,0)
sigma = matrix(c(1,0.8,0.8,1),nrow=2,ncol=2)
Y = rmvnorm(n, mean = theta, sigma= sigma)

## true value
x.points <- seq(-2.5,2.5,length.out=100)
y.points <- x.points
z <- matrix(0,nrow=100,ncol=100)
for (i in 1:100) {
  for (j in 1:100) {
    z[i,j] <- dmvnorm(c(x.points[i],y.points[j]),
                      mean=theta,sigma=sigma)
  }
}
contour(x.points,y.points,z, main= expression(paste('Contour plot using true value of ',theta, ' and ',Sigma)))


## MLE
theta_MLE = apply(Y,2,mean)
sigma_MLE = t(Y-theta_MLE)%*%(Y-theta_MLE)/n

x.points <- seq(-2.5,2.5,length.out=100)
y.points <- x.points
z <- matrix(0,nrow=100,ncol=100)
for (i in 1:100) {
  for (j in 1:100) {
    z[i,j] <- dmvnorm(c(x.points[i],y.points[j]),
                      mean=theta_MLE,sigma=sigma_MLE)
  }
}
contour(x.points,y.points,z, main= expression(paste('Contour plot using MLE of ',theta, ' and ',Sigma)))

# comment on the difference (I don't know what to comment)


# (b) assume independent normal & inverse-Wishart priors for theta and sigma

gibbs_normal = function(mu_0, Lambda_0, nu_0, S_0, n_iter, Y){
  n = nrow(Y)
  ybar = apply(Y,2,mean)
  
  ## Gibbs Sampler
  #Initial values for Gibbs sampler
  Sigma <- cov(Y)
  #Set null matrices to save samples
  THETA <- SIGMA <- NULL
  
  #first set number of iterations and burn-in, then set seed
  burn_in <- 0.3*n_iter
  set.seed(1234)
  
  for (s in 1:(n_iter+burn_in)){
    ##update theta using its full conditional
    Lambda_n <- solve(solve(Lambda_0) + n*solve(Sigma))
    mu_n <- Lambda_n %*% (solve(Lambda_0)%*%mu_0 + n*solve(Sigma)%*%ybar)
    theta <- rmvnorm(1,mu_n,Lambda_n)
    #update Sigma
    S_theta <- (t(Y)-c(theta))%*%t(t(Y)-c(theta))
    S_n <- S_0 + S_theta
    nu_n <- nu_0 + n
    Sigma <- riwish(nu_n, S_n)
    #save results only past burn-in
    if(s > burn_in){
      THETA <- rbind(THETA,theta)
      SIGMA <- rbind(SIGMA,c(Sigma))
    }
  } 
  colnames(THETA) <- c("theta_1","theta_2")
  colnames(SIGMA) <- c("sigma_11","sigma_12","sigma_21","sigma_22") #symmetry in sigma
  
  return(list(theta = THETA, sigma = SIGMA))
}

# Initialize the hyperparameters
mu_0 = c(0,0)
Lambda_0 = matrix(c(1,0.5,0.5,1),nrow=2, ncol=2)
nu_0 = 4
S_0 =  matrix(c(1,0,0,1),nrow=2, ncol=2)

# Performs Gibbs Sampling
normal_mcmc = gibbs_normal(mu_0, Lambda_0, nu_0, S_0, n_iter=1000, Y=Y)

# diagnostics check (autocorrelation plot seems reasonable)
library(coda)
THETA.mcmc <- mcmc(normal_mcmc$theta,start=1);
plot(THETA.mcmc[,"theta_1"])
autocorr.plot(THETA.mcmc[,"theta_1"])


# (c) Compare Bayes estimates of (theta, sigma) to MLE and truth
# comment on the similarities and differences (how much influence prior has)
apply(normal_mcmc$theta, 2, mean)
#theta_1    theta_2 
#0.03985030 0.06716113 

# prior mean is mu_0 = 0
# sample mean is theta_MLE
#Y1         Y2 
#0.03923707 0.06191307 

# so posterior expectation of theta is hugely influenced by sample mean

apply(normal_mcmc$sigma, 2, mean)
#sigma_11  sigma_12  sigma_21  sigma_22 
#1.1899463 0.8687961 0.8687961 0.9786923 

# prior expectation is S_0 = 1,0,0,1
# prior sample estimate is sigma_MLE
#sigma_11  sigma_12  sigma_21  sigma_22
#1.1853870 0.8715629  0.8715629 0.9676537

# so posterior expectation of sigma is hugely influenced by MLE OF sigma (sample estimate)


# (d) Given y_i2 values for 50 new test subjects, how to use Gibbs sampler to predict new y_i1 values from the "conditional posterior predictive distribution"



############################################################################################
## Question 2: Hoff problem 7.4
library(MASS)
age = read.table('agehw.dat', header = TRUE)

#### b) generate a prior predictive dataset of size n = 100
n = 100
m = 15 # 15 datasets?

# Hyperparameters for the priors
mu_0 = c(55,55)
Lambda_0 = matrix(c(225,180,180,225),nrow=2, ncol=2)
nu_0 = 4
S_0 =  Lambda_0

# Set null matrices to save samples
THETA <- SIGMA <- Y.prior <- NULL

for (i in 1:m){
  # generate
  theta = rmvnorm(1, mu_0, Lambda_0)
  Sigma = riwish(nu_0, S_0)
  Y = rmvnorm(n, theta, Sigma)
  
  # update
  THETA <- rbind(THETA,theta)
  SIGMA <- rbind(SIGMA,c(Sigma))
  Y.prior <- rbind(Y.prior,c(Y))
}

# plotting, visual check
ids = sample(1:m, 3)
for (id in ids){
  plot(Y.prior[id, 1:100], Y.prior[id, 101:200], xlab = 'Husband Age', ylab = 'Wife Age')
  abline(a=0, b=1)
}

# The generated prior predictive datasets do represent my prior beliefs on the real dataset.

#### c) Using your prior distribution and the 100 values in the dataset, obtain an MCMC approximation

# Hyperparameters for the priors
mu_0 = c(55,55)
Lambda_0 = matrix(c(225,180,180,225),nrow=2, ncol=2)
nu_0 = 4
S_0 =  Lambda_0

# Gibbs Sampler
normal_age_mcmc = gibbs_normal(mu_0, Lambda_0, nu_0, S_0, n_iter=1000, Y=age)

## Plots
plot(normal_age_mcmc$theta[,"theta_1"], normal_age_mcmc$theta[,"theta_2"], xlab="Husband Ages", ylab="Wife Ages", main = "Joint posterior distribution of THETA")
CORR = normal_age_mcmc$sigma[,2]/sqrt((normal_age_mcmc$sigma[, 1]* normal_age_mcmc$sigma[, 4]))
plot(density(CORR), main = "Posterior Density of Correlation Between Yh and Yw")

##  95% posterior confidence intervals for THETA and correlation coefficient
quantile(normal_age_mcmc$theta[,"theta_1"], probs = c(0.025,0.975))
quantile(normal_age_mcmc$theta[,"theta_2"], probs = c(0.025,0.975))
quantile(CORR, probs = c(0.025,0.975)) # the correlation is even higher

#### d) Obtain 95% posterioe confidence intervals using the following prior distributions

## i. Jeffrey's prior
n = nrow(age)
ybar = apply(age, 2, mean)
Sigma = cov(age)
THETA <- SIGMA <- NULL

n_iter <- 1000; burn_in <- 0.3*n_iter
set.seed(1234)

for (s in 1:(n_iter+burn_in)){
  ##update theta using its full conditional
  theta <- rmvnorm(1,ybar,Sigma/n)
  #update Sigma
  S_theta <- (t(age)-c(ybar))%*%t(t(age)-c(ybar))
  Sigma <- riwish(n+1, S_theta)
  #save results only past burn-in
  if(s > burn_in){
    THETA <- rbind(THETA,theta)
    SIGMA <- rbind(SIGMA,c(Sigma))
  }
} 
colnames(THETA) <- c("theta_1","theta_2")
colnames(SIGMA) <- c("sigma_11","sigma_12","sigma_21","sigma_22") #symmetry in sigma

# 95% posterior confidence interval
CORR = SIGMA[,2]/sqrt((SIGMA[, 1]* SIGMA[, 4]))

quantile(THETA[,"theta_1"], probs = c(0.025,0.975))
quantile(THETA[,"theta_2"], probs = c(0.025,0.975))
quantile(CORR, probs = c(0.025,0.975))

## iii. Diffuse prior
# Hyperparameters for prior
mu_0 = c(0,0)
Lambda_0 = 10^5 * diag(2)
nu_0 = 3
S_0 =  10^3 * diag(2)

# Gibbs sampler
normal_age_mcmc2 = gibbs_normal(mu_0, Lambda_0, nu_0, S_0, n_iter=1000, Y=age)

# 95% posterior confidence interval
CORR = normal_age_mcmc2$sigma[,2]/sqrt((normal_age_mcmc2$sigma[, 1]* normal_age_mcmc2$sigma[, 4]))

quantile(normal_age_mcmc2$theta[,"theta_1"], probs = c(0.025,0.975))
quantile(normal_age_mcmc2$theta[,"theta_2"], probs = c(0.025,0.975))
quantile(CORR, probs = c(0.025,0.975)) # the correlation is even higher

#### e) Compare the confidence intervals from d) to those obtained in c). Discuss whether or not you think that your prior info is helpful
# in estimating THETA and SIGMA

age25 = age[sample(1:100, 25),]

## my prior
mu_0 = c(55,55)
Lambda_0 = matrix(c(225,180,180,225),nrow=2, ncol=2)
nu_0 = 4
S_0 =  Lambda_0

normal_age25_mcmc = gibbs_normal(mu_0, Lambda_0, nu_0, S_0, n_iter=1000, Y=age25)

# post analysis
CORR = normal_age25_mcmc$sigma[,2]/sqrt((normal_age25_mcmc$sigma[, 1]* normal_age25_mcmc$sigma[, 4]))
quantile(normal_age25_mcmc$theta[,"theta_1"], probs = c(0.025,0.975))
quantile(normal_age25_mcmc$theta[,"theta_2"], probs = c(0.025,0.975))
quantile(CORR, probs = c(0.025,0.975))

## diffuse prior
mu_0 = c(0,0)
Lambda_0 = 10^5 * diag(2)
nu_0 = 3
S_0 =  10^3 * diag(2)

normal_age25_mcmc = gibbs_normal(mu_0, Lambda_0, nu_0, S_0, n_iter=1000, Y=age25)

# post analysis
CORR = normal_age25_mcmc$sigma[,2]/sqrt((normal_age25_mcmc$sigma[, 1]* normal_age25_mcmc$sigma[, 4]))
quantile(normal_age25_mcmc$theta[,"theta_1"], probs = c(0.025,0.975))
quantile(normal_age25_mcmc$theta[,"theta_2"], probs = c(0.025,0.975))
quantile(CORR, probs = c(0.025,0.975))

