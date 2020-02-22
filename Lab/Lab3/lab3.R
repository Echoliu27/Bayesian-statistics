library(tidyverse)
library(rstanarm)
library(magrittr)
library(rstan)
library(bayesplot)
library(loo)
library(readxl)

GoT <- read_xlsx("GoT_Deaths.xlsx", col_names = T)
head(GoT)

# Exercise 1: Plot a histogram of the death counts.
hist(GoT$Count)

# Natural to model the data with a Poisson distribution
y <- GoT$Count
n <- length(y)

# Finish: obtain empirical mean
ybar <- mean(y)
sim_dat <- rpois(n, ybar)
qplot(sim_dat, bins = 20, xlab = "Simulated number of deaths", fill = I("#9ecae1"))

# create data frame for side-by-side histograms
df <- rbind(data.frame(y, "observed") %>% rename(count = 1, data = 2),
            data.frame(sim_dat, "simulated") %>% rename(count = 1, data = 2))
ggplot(df, aes(x = count, fill = data))+
  geom_histogram(position = "dodge", bins = 20) +
  scale_fill_brewer()  +
  labs(x = "Number of deaths", y = "Count")

# there are many more 0-valued observations

## Poisson model: lambda ~ Gamma(10,2)
stan_dat <- list(y = y, N = n)
fit <- stan("lab-03-poisson-simple.stan", data = stan_dat, refresh = 0, chains = 2)
lambda_draws <- as.matrix(fit, pars = "lambda") # posterior sampled values of lambda


# Exercise 2: plot smoothed posterior distribuion for lambda with a 90% highest posterior
# density region
mcmc_areas(lambda_draws, prob = 0.9)
print(fit, pars = 'lambda')
# how does the posterior mean for lambda compare to the sample mean?

## Graphical posterior predictive checks (PPC)
# Exercise 3: Generate posterior predictive samples using the posterior values of lambda
# and store the result in a length(lambda_draws) by n matrix

y_rep <- apply(lambda_draws, 1, function(x){rpois(n = n, lambda = x)}) %>% t()
ppc_hist(y, y_rep[1:8, ], binwidth = 1) # compare the empirical distribution of data y to distribution of simulated date
ppc_dens_overlay(y, y_rep[1:60, ])


# calculate the proportions of 0s in y and compare that to the proportion of zeros in 
# posterior predictive samples
prop_zero <- function(x){
  mean(x == 0)
} 
prop_zero(y)
## [1] 0.16
ppc_stat(y, y_rep, stat = "prop_zero")
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
# plot the means and variances for all of simulated datasets
ppc_stat_2d(y, y_rep, stat = c("mean", "var"))
ppc_error_hist(y, y_rep[1:4, ], binwidth = 1) + xlim(-15, 15)

# Exercise 4: Based on these PPCs, does this model appear to be a good fit for the data?
# NO, because it doesn't model number of 0s closely.


## Poisson Hurdle Model
# zero-inflated Poisson model: mixture of a Poisson with a point mass at zero
# hurdle model:
# Exercise 5:
fit2 <- stan("lab-03-poisson-hurdle.stan", data = stan_dat, refresh = 0, chains = 2)

# Extract the sampled vlaues for lambda, and store them in a variable called lambda_draws2:
lambda_draws2 <- as.matrix(fit2, pars = "lambda")

# Compare
lambdas <- cbind(lambda_fit1 = lambda_draws[, 1],
                 lambda_fit2 = lambda_draws2[, 1])

# Shade 90% interval
mcmc_areas(lambdas, prob = 0.9) 

y_rep2 <- as.matrix(fit2, pars = "y_rep")

## Exercise 6: produce the same PPC vizs as before for the new results. Comment on 
# how this second model compares to both the observed data and to the simple Poisson model
ppc_hist(y, y_rep2[1:8, ], binwidth = 1) # compare the empirical distribution of data y to distribution of simulated date
ppc_dens_overlay(y, y_rep2[1:60, ])

prop_zero <- function(x){
  mean(x == 0)
} 
prop_zero(y)
## [1] 0.16
ppc_stat(y, y_rep2, stat = "prop_zero")
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
ppc_stat_2d(y, y_rep2, stat = c("mean", "var"))
ppc_error_hist(y, y_rep2[1:4, ], binwidth = 1) + xlim(-15, 15)

# write comments on this:

## Leave-one-out cross-validation
log_lik1 <- extract_log_lik(fit, merge_chains = FALSE) # convenience function for extracting the pointwise log-likelihood mx from a fitted Stan
r_eff1 <- relative_eff(exp(log_lik1)) #MCMC effective sample size divided by the total sample size
(loo1 <- loo(log_lik1, r_eff = r_eff1)) #

# elpd_loo: point estimate and standard errors of the expected log pointwise predictive density
# p_loo: effective number of aprameters
# looic (-2*elpd_loo): LOO information criterion


log_lik2 <- extract_log_lik(fit2, merge_chains = FALSE)
r_eff2 <- relative_eff(exp(log_lik2)) 
(loo2 <- loo(log_lik2, r_eff = r_eff2))
compare(loo1, loo2)

# Exercise 7: which model performs better in terms of prediction?
#

# Exercise 8: Why are PPCs important?
# We use posterior predictice checks to look for systematic discrepancies between real and simulated data


# Exercise 9: Was the second model a good fit for the data? Why or why not?

# Exercise 10: If someone reported a single LOOCV error to you, would that be useful? Why or why not?

