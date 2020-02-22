require(tidyverse)
require(magrittr)
require(bayesplot)
require(loo)
require(readxl)
require(plyr)
require(ggrepel)
library(cowplot)
library(truncnorm)

## Visualize the data generation process
n <- 1000
true_Rsquared <- 5
true_sigma <- 1.25
#
u <- runif(n, 0, true_Rsquared)
r <- u + rnorm(n, sd = true_sigma)
theta <- runif(n, 0, 2*pi)
#
ggplot2::ggplot() +
  geom_point(data = data.frame(x = sign(r)*sqrt(abs(r))*cos(theta), y = sign(r)*sqrt(abs(r))*sin(theta)),
             aes(x = x, y = y), shape = 1) +
  geom_path(data = data.frame(R = true_Rsquared) %>%
              plyr::ddply(.(R), function(d){
                data.frame(x = sqrt(d$R)*cos(seq(0, 2*pi, length.out = 100)),
                           y = sqrt(d$R)*sin(seq(0, 2*pi, length.out = 100)))
              }),
            aes(x = x, y = y), alpha = 1, colour = "red") +
  coord_fixed()


#############################################################################################
# hyper-parameters
m <- 3
k <- 1
alpha <- 5/2
beta <- 5/2

#
rpareto <- function(m, k, trunc = NULL){
  p <- m*(1 - runif(1))^(-1/k)
  if(!is.null(trunc)){
    while(p > trunc){
      p <- m*(1 - runif(1))^(-1/k)
    }
  }
  return(p)
}

#
uni_pareto_gibbs <- function(S, r, m, k, alpha, beta, burn_in = min(1000, S / 2), thin = 5){
  # Reparametrize X matrix to squared radius values
  Rsq <- r
  n <- length(Rsq)
  R <- rep(1, S)
  U <- matrix(0, nrow = S, ncol = n)
  U[1, ] <- runif(n, 0, R)
  sigma <- rep(1, S)
  #
  U_curr <- U[1, ]
  R_curr <- R[1]
  sigma_curr <- sigma[1]
  for(s in 1:S){
    # Sample from full conditional of the inner radius
    R_curr <- rpareto(max(c(U_curr, m)), k + n)
    R[s] <- R_curr
    # Sample from full conditional of U values
    U_curr <- truncnorm::rtruncnorm(n, a = 0, b = R_curr, mean = Rsq, sd = sigma_curr)
    U[s, ] <- U_curr
    # Sample from full conditional of sigma
    sigma_curr <-  sqrt(1/rgamma(1,shape=n/2+alpha, rate=0.5*sum((Rsq-U_curr)^2)+beta))#complete this line
    sigma[s] <- sigma_curr
  }
  return(list(R = R[seq(burn_in, S, by = thin)], 
              U = U[seq(burn_in, S, by = thin), ], 
              sigma = sigma[seq(burn_in, S, by = thin)]))
}
#
gibbs_samps <- uni_pareto_gibbs(S = 100000, r, m, k, alpha, beta, burn_in=2000)


## Exercise 1: Complete the Gibbs sampling step for sigma.

# visualize our samples along with the data
ggplot2::ggplot() +
  geom_point(data = data.frame(x = sign(r)*sqrt(abs(r))*cos(theta), y = sign(r)*sqrt(abs(r))*sin(theta)),
             aes(x = x, y = y), shape = 1) +
  geom_path(data = data.frame(R = gibbs_samps$R) %>%
              plyr::ddply(.(R), function(d){
                data.frame(x = sqrt(d$R)*cos(seq(0, 2*pi, length.out = 100)),
                           y = sqrt(d$R)*sin(seq(0, 2*pi, length.out = 100)))
              }),
            aes(x = x, y = y), alpha = 0.005, colour = "blue") +
  geom_path(data = data.frame(R = true_Rsquared) %>%
              plyr::ddply(.(R), function(d){
                data.frame(x = sqrt(d$R)*cos(seq(0, 2*pi, length.out = 100)),
                           y = sqrt(d$R)*sin(seq(0, 2*pi, length.out = 100)))
              }),
            aes(x = x, y = y), alpha = 1, colour = "red") +
  coord_fixed()


## Exercise 2: Make a plot to visualize the marginal posterior density of R^2. How does this
# density compare to the true value of R^2
hist(gibbs_samps$R) # marginal posterior density of R
abline(v = 5, col=c("red"), lty=2, lwd=c(1))
legend("topright", legend=c("Density", "True R square"), col=c("black", "red"), lty=1:2, cex=0.8, bg="transparent", bty = "n")
# true_Rsquared is 5

## Exercise 3: Make a plot to visualize the marginal posterior density of sigma^2. How does this 
# density compare to the true value of sigma^2
hist((gibbs_samps$sigma)^2)
abline(v = 1.25^2, col=c("red"), lty=2, lwd=c(1))
legend("topright", legend=c("Density", "True Sigma square"), col=c("black", "red"), lty=1:2, cex=0.8, bg="transparent", bty = "n")
# true_sigma is 1.25^2 = 1.5625

## Exercise 4: Make a contour plot to visualize the bivariate posterior density of R^2 and sigma^2
# Indicate this truth on the plot
df <- data.frame(rsq=gibbs_samps$R, sigmasq=(gibbs_samps$sigma)^2)
ggplot(df, aes(x=rsq,y=sigmasq)) + 
  #geom_point() +
  #theme_cowplot(12)
  geom_density_2d()+
  geom_point(aes(x=5, y=1.5625), colour="red", size = 3)+
  ggtitle(expression(paste('Contour Plot of Bivaraite Poseterior density ',R^2,' and ', sigma^2)))

# comment on the posterio
