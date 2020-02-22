install.packages('tidyverse')
install.packages('rstanarm')
install.packages('magrittr')
install.packages('rstan')

library(tidyverse)
library(rstanarm)
library(magrittr)
library(rstan)

tumors <- read.csv(file = url("http://www.stat.columbia.edu/~gelman/book/data/rats.asc"),
                   skip = 2, header = T, sep = " ")[,c(1,2)]
y <- tumors$y
N <- tumors$N
n <- length(y)

plot(seq(0, 1, length.out = 1000), 
     dbeta(seq(0, 1, length.out = 1000), 1, 1),
     type = 'l', 
     xlab = expression(theta), ylab = "Density",
     main = "The Beta(1, 1) density")

stan_dat <- list(n = n, N = N, y =y, a = 1, b = 1)
fit_pool <- stan('lab-02-pool.stan', data = stan_dat, chains = 2, refresh = 0)
pool_output <- rstan::extract(fit_pool)
mean(pool_output$theta)

