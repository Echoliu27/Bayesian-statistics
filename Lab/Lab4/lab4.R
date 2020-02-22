library(rstan)
library(tidyverse)
library(rstanarm)
library(magrittr)
library(reshape)

## Prior selection
set.seed(689934)

alpha <- 1   
beta <- -0.25 
sigma <- 1    

N <- 5
x <- array(runif(N, 0, 2), dim=N)                    
y <- array(rnorm(N, beta * x + alpha, sigma), dim=N)


## 1. Flat Priors
stan_dat <- list(y = y, x=x, N=N)
fit.flat <- stan(file = "lab-04-flat_prior.stan", data = stan_dat, chains = 1, refresh = 0, iter = 2000, warmup = 500, seed=48)
alpha.flat <- as.matrix(fit.flat, pars = "alpha")
beta.flat <- as.matrix(fit.flat, pars = "beta")

ggplot(alpha.flat %>% as.data.frame, aes(x = alpha)) +
  geom_histogram(bins = 30)
# the posterior for alpha is quite diffuse (much uncertainty about what alpha is)

print(fit.flat, pars = c("alpha"))

## Exercise 1: write down the posterior means of alpha and beta. Do the results seem surpirsing?
mean(alpha.flat) #0.6982
mean(beta.flat) #0.3537
quantile(alpha.flat, probs = c(0.025, 0.975)) #-3.34, 4.04
quantile(beta.flat, probs = c(0.025, 0.975)) #-2.27, 2.93

#flat priors may actually bias our estimates for a parameter by allowing the posterior to b
#be pulled towards extreme and unlikely values


## 2. Uniform priors
stan_dat <- list(y = y, x=x, N=N, lb = -10, ub = 10)
fit.unif <- stan(file = "lab-04-unif_prior.stan", data = stan_dat, chains = 1, refresh = 0, iter = 2000, warmup = 500, seed=48)
alpha.unif <- as.matrix(fit.unif, pars = c("alpha"))
beta.unif <- as.matrix(fit.unif, pars = c("beta"))

ggplot(alpha.unif %>% as.data.frame, aes(x = alpha)) +
  geom_histogram(bins = 30)

print(fit.unif, pars = c("alpha"))

#the posterior mean under this uniform prior is closer to the true value,
#the posterior is still very spread out. So a diffuse/flat prior is necessarily non-informative
# and in these cases is actually extremely informative.
# --> informative prior (specify how informative we would like our prior beliefs to be)
###############################################################################################

## Weakly informative priors
#( this sort of prior ought to rule out unreasonable parameter values, but is not so strong as
# to rule out possible values which might make sense.)

## 3. light-tailed (N(0,1))
stan_dat <- list(y = y, x=x, N=N)
fit.norm <- stan(file = "lab-04-normal_prior.stan", data = stan_dat, chains = 1, refresh = 0, iter = 2000, warmup = 500, seed=49)
alpha.norm<- as.matrix(fit.norm, pars = c("alpha"))
beta.norm <- as.matrix(fit.norm, pars = c("beta"))

ggplot(alpha.norm %>% as.data.frame, aes(x = alpha)) +
  geom_histogram(bins = 30)
print(fit.norm, pars = c("alpha"))

## Exercise 2: compute the posterior means of alpha and beta
mean(alpha.norm) #0.569
mean(beta.norm) #0.397
quantile(alpha.norm, probs = c(0.025, 0.975)) #-0.548, 1.637
quantile(beta.norm, probs = c(0.025, 0.975)) #-0.426, 1,251

## 4. Hevy-tailed (Cauchy(0,1): t-distribution with df=1, fatter tails than normal distribution)
stan_dat <- list(y = y, x=x, N=N)
fit.cauchy <- stan(file = "lab-04-cauchy_prior.stan",data = stan_dat, chains = 1, refresh = 0, iter = 2000, warmup = 500, seed=55)
alpha.cauchy<- as.matrix(fit.cauchy, pars = c("alpha"))

ggplot(alpha.cauchy %>% as.data.frame, aes(x = alpha)) +
  geom_histogram(bins = 30)
print(fit.cauchy, pars = c("alpha"))

# following plot displays the posteriors for alpha under these two priors
#(utility function)
create_df <- function(post_draws, prior_draws){
  post_draws <- data.frame(post_draws)
  post_draws$distribution <- "posterior"
  
  prior_draws <- data.frame(prior_draws)
  colnames(prior_draws) <- "alpha"
  prior_draws$distribution <- "prior"
  dat <- rbind(post_draws, prior_draws)
  return(dat)
}

plot_dat <- create_df(alpha.norm, alpha.cauchy) %>% 
  mutate(distribution = if_else(distribution == "posterior", "Normal","Cauchy"))

ggplot(plot_dat, aes(alpha, fill = distribution)) + 
  geom_histogram(binwidth = 0.25, alpha = 0.7, position = "identity")+
  geom_vline(xintercept = alpha) +
  scale_fill_brewer()

# the Cauchy prior allocates higher probability mass to extreme values as compared
# to the normal prior (when still concentrating most of the posterior mass for alpha within
# a desired scale)

## Exercise 3: Would you say that a Cauchy prior is more or less informative than a normal
# prior (assume that their interquartile ranges are comparable?)

# Although the Cauchy prior forces most of the posterior within the specified scale, the heavy tail allows for the occasional extreme value 
# that stresses the evaluation of the special functions and can drastically hinder computational performance

################################################################################

## Sensitivity to prior selection
# we happened to choose a good scale for the parameter, what happens when this is not the case

alpha <- 10
N <- 10
x <- runif(N, 0, 2)                    
y <- rnorm(N, beta * x + alpha, sigma)

## 5. Normal prior: α∼N(0,1).
stan_dat <- list(y = y, x=x, N=N)
fit.norm <- stan(file = "lab-04-normal_prior.stan", data = stan_dat, chains = 1, refresh = 0, iter = 2000, warmup = 500, seed=49)
alpha.norm<- as.matrix(fit.norm, pars = c("alpha"))

prior_draws <- rnorm(1000, 0, 1)
plot_dat <- create_df(alpha.norm, prior_draws)

ggplot(plot_dat, aes(alpha, fill = distribution)) + 
  geom_histogram(binwidth = 0.25, alpha = 0.7, position = "identity")+
  geom_vline(xintercept = alpha) +
  scale_fill_brewer()

# how the prior is dominating the likelihood. The posterior is extremely sensitive to 
# the choice of our prior, so much so that we don't observe posterior values close to the true alpha at all.
# Instead, posterior is concentrated around the upper extremes of the prior

## 6. Cauchy prior
stan_dat <- list(y = y, x=x, N=N)
fit.cauchy <- stan(file = "lab-04-cauchy_prior.stan",data = stan_dat, chains = 1, refresh = 0, iter = 2000, warmup = 500, seed=55)
alpha.cauchy<- as.matrix(fit.cauchy, pars = c("alpha"))

prior_draws <- rcauchy(1000, 0, 1)
prior_draws <- prior_draws[abs(prior_draws) < 25]
plot_dat <- create_df(alpha.cauchy, prior_draws)

ggplot(plot_dat, aes(alpha, fill = distribution)) + 
  geom_histogram(binwidth = .5, alpha = 0.7, position = "identity")+
  geom_vline(xintercept = alpha) +
  scale_fill_brewer()

# The posterior is able to concentrate around the true alpha=10. The heavy tails of 
# the Cauchy allow the posterior to move beyond the scale occupied by the prior.
# From this histogram, it's much clearer that the prior we chose was probably inappropriate and
# conflicts with the data

## Exercise 4: what happens as we increase the number of observations?

# data dominates and prior has less influence

## Exercise 5: when might we prefer to use lighter- versus heavier-tailed priors?

# When you're certain about a scale, than using ligher-tailed priors are better.
# Lighter tails induce strong regularization towards the given scale.

## Exercise 6: How might we determine a reasonable scale for the prior?
# If you don't know the scale

#####################################################################################

## Model Reparameterization
theta <- runif(10000,0,1)
hist(theta)

## 7. og-odds ratio: map directly between values theta and li
logit <- function(x){
  ret <- log(x/(1-x))# finish function for log odds
  return(ret)
}
eta <- logit(theta)
hist(eta)

# pretend the true birth rate of females in Paris was 0.3, simulate N=10 obs
set.seed(123);
theta <- 0.3;
N <- 10;
y <- rbinom(N, 1, theta)
  
theta.mle <- sum(y)/N #MLE for theta is the proportion of successes

## 8. Success Probability Parameterization
stan_dat <- list(y = y,N=N)
fit.bayes.prob <- stan(file = "lab-04-prob.stan", data = stan_dat, refresh = 0, iter = 2000)
print(fit.bayes.prob, pars = c("theta", "eta"))

# we get a posterior mean of 0.42

## Exercise 7: what would we expect to be the posterior mode of our samples?
# calculate the posterior mode theoretically and compare it to the estimated mode from the
# posterior samples

# we can take the derivative and set to 0 to find the maximum of the density function
model.prob <- stan_model("lab-04-prob.stan")
fit.pmode.prob <- optimizing(model.prob, data=c("N", "y")) # obtain a point estimate by maximizing the joint posterior
fit.pmode.prob$par
#theta        eta 
#0.4000013 -0.4054598 


## 9. Log-odds parameterization
# what about a uniform prior belief on the log odds li?
# induced prior on theta would be a Beta(0,0)

## Exercise 8: Is this prior proper? Does it result in a proper posterior? If so, under which conditions?
# Not necessarily

fit.logodds <- stan(file = "lab-04-log_odds.stan", data = stan_dat, refresh = 0, iter = 2000)
print(fit.logodds, pars = c("theta", "eta"))

## 10. Odds Parameterization
## Exercise 9: If we set a uniform prior for , what is the induced prior on theta? Is this prior proper?


## 11. Jeffrey's prior
fit.jeffreys <- stan(file = "jeffreys.stan", data = stan_dat, refresh = 0, iter = 1000)
print(fit.jeffreys, pars = c("theta"))

