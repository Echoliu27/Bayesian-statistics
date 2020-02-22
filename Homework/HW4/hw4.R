## Q1: Hoff 6.1

# data and prior
man_A <- c(1, 0, 0, 1, 2, 2, 1, 5, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 2, 1, 3, 
           2, 0, 0, 3, 0, 0, 0, 2, 1, 0, 2, 1, 0, 0, 1, 3, 0, 1, 1, 0, 2, 0, 0, 2, 2, 1, 
           3, 0, 0, 0, 1, 1)
man_B <- c(2, 2, 1, 1, 2, 2, 1, 2, 1, 0, 2, 1, 1, 2, 0, 2, 2, 0, 2, 1, 0, 0, 3, 6, 1, 6,
           4, 0, 3, 2, 0, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 1, 0, 4, 2, 1, 0, 0, 1, 0, 3, 2, 
           5, 0, 1, 1, 2, 1, 2, 1, 2, 0, 0, 0, 2, 1, 0, 2, 0, 2, 4, 1, 1, 1, 2, 0, 1, 1, 
           1, 1, 0, 2, 3, 2, 0, 2, 1, 3, 1, 3, 2, 2, 3, 2, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 
           3, 3, 0, 1, 2, 2, 2, 0, 6, 0, 0, 0, 2, 0, 1, 1, 1, 3, 3, 2, 1, 1, 0, 1, 0, 0, 
           2, 0, 2, 0, 1, 0, 2, 0, 0, 2, 2, 4, 1, 2, 3, 2, 0, 0, 0, 1, 0, 0, 1, 5, 2, 1, 
           3, 2, 0, 2, 1, 1, 3, 0, 5, 0, 0, 2, 4, 3, 4, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 2, 
           0, 0, 1, 1, 0, 2, 1, 3, 3, 2, 2, 0, 0, 2, 3, 2, 4, 3, 3, 4, 0, 3, 0, 1, 0, 1, 
           2, 3, 4, 1, 2, 6, 2, 1, 2, 2)
sum_A <- sum(man_A)
sum_B <- sum(man_B)
n_A <- length(man_A)
n_B <- length(man_B)

a_theta <- 2; b_theta <- 1
ab_gamma <- c(8,16,32,64,128)
m <- length(ab_gamma)
S <- 5000

# starting point
theta_0 = a_theta / b_theta
gamma_0 = 1
gibbs_theta = matrix(nrow = m, ncol = S + 1)
gibbs_gamma = matrix(nrow = m, ncol = S + 1)
gibbs_theta[,1] = theta_0
gibbs_gamma[,1] = gamma_0

# gibbs sampler
for (s in 2:(S+1)) {
  gibbs_theta[,s] = rgamma(m, shape = sum_A + sum_B + a_theta, 
                           rate = n_A + n_B * gibbs_gamma[,s-1] + b_theta)
  gibbs_gamma[,s] = rgamma(m, shape = sum_B + ab_gamma, 
                           rate = n_B * gibbs_theta[,s] + ab_gamma)
}

# analysis
theta_A = gibbs_theta
theta_B = gibbs_theta * gibbs_gamma
E_BmA = rowMeans(theta_B - theta_A)
plot(ab_gamma, E_BmA, "b",
     xlab = "a_gamma = b_gamma",
     ylab = "E[theta_B - theta_A]",
     main = "Effects of Different Priors on the Results")

########################################################################################
## Q2:
# simulate data and prior
lambda = gamma = 1
set.seed(1)
treated = rpois(n=50, lambda=1)
control = rpois(n=50, lambda=1)

sum_tr <- sum(treated)
sum_con <- sum(control)
n_tr <- length(treated)
n_con <- length(control)

S = 10000 # number of samples to draw

# starting point
lambda_0 = 1
gamma_0 = 1
gibbs_lambda = matrix(nrow = 1, ncol = S + 1)
gibbs_gamma = matrix(nrow = 1, ncol = S + 1)
gibbs_lambda[,1] = lambda_0
gibbs_gamma[,1] = gamma_0

# gibbs sampler
set.seed(1234)
for (s in 2:(S+1)){
  gibbs_lambda[,s] = rgamma(1, shape = sum_tr + sum_con + 1, 
                           rate = n_con + n_tr * gibbs_gamma[,s-1] + 1)
  gibbs_gamma[,s] = rgamma(1, shape = sum_tr+1, 
                           rate = n_tr * gibbs_lambda[,s] + 1)
}

# posterior mean for log(lambda)
mean(log(gibbs_lambda[1,])) #-0.0721

# 95% credible interval for log(lambda)
quantile(log(gibbs_lambda[1,]), c(0.025,0.975)) 
#2.5%      97.5% 
 # -0.3423654  0.1977386 


# predictive distribution
#1) for x_i = 0
y.mc.x0 = rpois(1001, gibbs_lambda[1,])
hist(y.mc.x0, breaks = seq(-0.5, max(y.mc.x0)+0.5,1), main="Histogram of predictive distributions of xi=0")
#2) for x_i = 1
y.mc.x1 = rpois(1001, gibbs_lambda[1,]*gibbs_gamma[1,])
hist(y.mc.x1, breaks = seq(-0.5, max(y.mc.x1)+0.5,1), main="Histogram of predictive distributions of xi=1")

## Plotting, convergence diagnostics (is your chain mixing well? what is the effective sample size? does the mixing differ for lambda and gamma?)
library(truncnorm)
library(coda)

# lambda
plot(gibbs_lambda[1,], ylab=expression(lambda),xlab='Iteration',main=expression(paste("Sampled values of ",lambda)))
abline(a=mean(gibbs_lambda[1,]),b=0,col='red4',lwd=2)
traceplot(as.mcmc(gibbs_lambda[1,]), main=expression(paste('Trace plot of ', lambda)))

# gamma
plot(gibbs_gamma[1,], ylab=expression(gamma),xlab='Iteration',main=expression(paste("Sampled values of ",gamma)))
abline(a=mean(gibbs_gamma[1,]),b=0,col='red4',lwd=2)
traceplot(as.mcmc(gibbs_gamma[1,]), main=expression(paste('Trace plot of ', gamma)))

summary(mcmc(gibbs_lambda[1,], start = 1))
summary()
effectiveSize(mcmc(gibbs_lambda[1,], start = 1))
effectiveSize(mcmc(gibbs_gamma[1,], start = 1))

plot(mcmc(gibbs_lambda[1,], start = 1))
plot(mcmc(gibbs_gamma[1,], start = 1))

# autocorrelation plots
autocorr.plot(mcmc(gibbs_lambda[1,]), main=expression(paste("Autocorrelation for ",lambda)))
autocorr.plot(mcmc(gibbs_gamma[1,]), main=expression(paste("Autocorrelation for ",gamma)))
