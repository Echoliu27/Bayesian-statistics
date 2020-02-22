library(tidyverse)
library(plotly)
library(reshape2)

### Question 4
# b) sampling distribution
Y = 57
N = 100
theta <- seq(from=0, to=1, by=0.1)
dens <-  choose(N,Y) * theta^Y*(1-theta)^(N-Y) # sum of Y as sufficient statistics is a binomial distribution
df <- data.frame(theta, dens)

#lines(theta,dens,type='l',col='blue')
#plot(x = theta, y = dens, type='l',xlab="theta", ylab="sampling distribution")
ggplot(df, aes(x = theta, y = dens)) +
  geom_bar(stat = 'identity') +
  scale_x_continuous(breaks = theta)

# c) posterior distribution with different constant
#x <- seq(0, 1, length = 20)
prior <- dbeta(theta,1,1)
posterior <- prior * dens/sum(dens)
ggplot(df, aes(x = theta, y = posterior)) +
  geom_bar(stat = 'identity') +
  scale_x_continuous(breaks = theta)

# d) igonre constants
theta.any <- seq(from=0, to=1, by=0.01)
posterior_wlcons <-  choose(N,Y) * theta.any^Y*(1-theta.any)^(N-Y)
plot(theta.any, type='l',posterior_wlcons, xlim=c(0,1))

# e) using formula for beta
posterior_beta <- dbeta(theta.any, 1+Y,1+N-Y)
plot(theta.any, type='l',posterior_beta, xlim=c(0,1))

# Discussion: d is the posterior without normalization by the constant sum(dens) and e is the posterior with normalization
# Their shapes are the same while scales are different due to the lack of influence of prior p(theta).


#############################################################
## Question 5
N = 100
exp.posterior = function(w, theta0, y) {
  (N / (w + N)) * (y / N) + (w / (w + N)) * theta0
}
Theta0 = rev(seq(0.0, 1, by = 0.1))
W = seq(0, 32, by = 0.5)

y = 57
d = outer(Theta0, W, FUN = function(theta0, w) exp.posterior(w, theta0, 57))
rownames(d) = Theta0
colnames(d) = W

df = melt(d)
colnames(df) = c('theta0', 'w', 'theta')

p = ggplot(df, aes(x = w, y = theta0, z = theta)) +
  geom_contour(aes(colour = ..level..)) +
  scale_x_continuous(breaks = c(1, 2, 8, 16, 32), labels = c(1, 2, 8, 16, 32)) +
  scale_y_continuous(breaks = Theta0)
library(directlabels)
direct.label(p, method = 'bottom.pieces')

# One can use this plot to determine whether or not they should believe that θ>0.5 by quantifying 
# their prior beliefs about the proportion with two factors: θ0, an initial estimate of the true proportion, 
# and w, the “sample size” of observed individuals that contributed to the initial estimate.

####################################################
## Question 6
find_n <- function(shape_1, shape_2){
  seq1 <- c()
  for (i in 100:3000){
    cdf <- pbeta(0.0015, shape1 = shape_1,shape2 = shape_2+i)
    seq1 <- c(seq1, cdf)
  }
  n <- Position(function(x) x > 0.95, seq1)
  return(n+100)
}

sample1 <- find_n(1,666) #1231
sample2 <- find_n(0.05,33.33) #146
sample3 <- find_n(1.6,407.4) #2311
sample4 <- find_n(1.05,497) #1564

seq2 <- c()
for (i in 100:3000){
  cdf <- pbeta(0.0015, shape1 = 1,shape2 = 1+i)/pbeta(0.1, shape1 = 1, shape2 = 1+i)
  seq2 <- c(seq2, cdf)
}
n <- Position(function(x) x > 0.95, seq2)
sample5 <- n+100 #199


######################################################
Theta0 = rev(seq(0, 1, by = 0.1))
N0 = seq(0, 32, by = 0.5)

a_formula <- function(n0, theta0){
  n0 * theta0
}

da <- outer(Theta0, N0, function(n0, theta0) a_formula(n0, theta0))
rownames(da) = Theta0
colnames(da) = N0
df_a = melt(da)
colnames(df_a) = c('theta0', 'n0', 'a')
df_a['b'] <- df_a['n0']-df_a['a']

