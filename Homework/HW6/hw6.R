library(mvtnorm)
library(MCMCpack)
require(tidyverse)

## Question 1

# (a)
# Simulate data: train and test
n_train = 100
n_test = 50
theta = c(0,0)
sigma = matrix(c(1,0.8,0.8,1),nrow=2,ncol=2)
set.seed(1234)
Y_train = rmvnorm(n_train, mean = theta, sigma= sigma)
colnames(Y_train) = c('y1','y2')
set.seed(123)
Y_test = rmvnorm(n_test, mean = theta, sigma= sigma)
Y_test_y2 = Y_test[,2]
colnames(Y_test) = c('y1','y2')

############################
# Rerun the gibbs sampler from hw5
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
normal_mcmc = gibbs_normal(mu_0, Lambda_0, nu_0, S_0, n_iter=1000, Y=Y_train)

#############
gibbs_composition = function(normal_mcmc, Y_test_y2){
  # innitialize
  Y1 = NULL
  n_iter <- 750; burn_in <- 0.3*n_iter
  set.seed(1234)
  
  for (s in 1:(n_iter+burn_in)){
    # sample theta and sigma from training set
    theta_1 = normal_mcmc$theta[s,][1]
    theta_2 = normal_mcmc$theta[s,][2]
    sigma_11 = normal_mcmc$sigma[s,][1]
    sigma_22 = normal_mcmc$sigma[s,][4]
    sigma_12 = sigma21 = normal_mcmc$sigma[s,][2]
    
    #update y1 which have 50 new predictions
    y1_vec = NULL
    for (y2 in Y_test_y2){
      theta_y1 = theta_1 + sigma_12%*%solve(sigma_22)%*%(y2 - theta_2)
      var_y1 = sigma_11 - sigma_12%*%solve(sigma_22)%*%sigma_12
      y1 = rnorm(1, mean = theta_y1, sd= var_y1)
      y1_vec = rbind(y1_vec, y1)
    }
    
    #save results only past burn-in
    if(s > burn_in){
      Y1 <- cbind(Y1,y1_vec)
    }
  }
  return(Y1)
}

Y1_predict = gibbs_composition(normal_mcmc, Y_test_y2)


#(b) Make plots
# par(mfrow =c(2,4))
# for (i in 1:50 ){
#   y1_dist = Y1_predict[i,]
#   mean.mc <- mean(y1_dist)
#   probs.mc <- quantile(y1_dist, c(0.025,0.975))
#   plot(density(y1_dist), col="black",lty=1,main=paste("Subject ", i))
#   abline(v = probs.mc, col=c("red"), lty=2, lwd=c(1))
#   abline(v = mean.mc, col="blue", lty=3,lwd=1)
#   legend("topright", legend=c(paste("mean: ", round(mean.mc,2))),
#          bg="transparent", bty = "n")
# }
stats = c()
for (i in 1:50 ){
  y1_dist = Y1_predict[i,]
  mean.mc <- mean(y1_dist)
  probs.mc <- quantile(y1_dist, c(0.025,0.975))
  stats = rbind(stats, c(mean.mc, probs.mc))
}
colnames(stats) <- c("mean", "lwr","upr")
stats = data.frame(stats)
gp = 1:50
plot(gp, stats$mean, ylim=c(-3, 3),xlab="y1 index", ylab="y", pch=16, cex=2)
# Add error bars, I really don't know what pictures they want
arrows(x0=gp, y0=stats$lwr, x1=gp, y1=stats$upr, code=3, angle=90, length=0.1)


# (c) What is the coverage of the 95% predictive intervals out of sample?
y1_true = Y_test[,1]
df = cbind(y1_true,stats)

df = as.data.frame(df) # turn it into a dataframe
count_within = df %>% 
  mutate(within = y1_true >= lwr & y1_true <= upr)%>%
  summarise(sum = sum(within))

count_within[1,1]/50 # 42 of the 95% predictive intervals contain the true y_i1 value

# (d) Fit a (frequentist) linear model to the original 100 train subjects
Y_train_df = as.data.frame(Y_train)
first_model = lm(y1 ~ y2, data = Y_train_df)
summary(first_model)

# predict y1 values given 50 new test y2
Y_test_y2_df = data.frame(Y_test_y2)
colnames(Y_test_y2_df) = 'y2'
y1_predict = predict(first_model, Y_test_y2_df)

# (e) Generate and make a plot of the predictive intervals for the frequentist predictions
y1_predict_interval = predict(first_model, Y_test_y2_df, interval = "predict")
df_fre = data.frame(y1_predict_interval)
names(df_fre) = c("mean", "lwr", "upr")

gp = 1:50
plot(gp, df_fre$mean, ylim=c(-3, 3),xlab="y1 index", ylab="y", pch=16, cex=2)
arrows(x0=gp, y0=df_fre$lwr, x1=gp, y1=df_fre$upr, code=3, angle=90, length=0.1)

## the plots look exactly the same?? why?


## Question 2

## part(c) Implement the Gibbs sampler, present point and interval estimates of the group-specific mean reaction times
## 1. Simulation of the real data
n = c(50, 70, 200, 100, 90)
set.seed(123456)
theta_j_true = rgamma(5, shape=2, rate =3)
#theta_j_true = c(2.2,2.1,2,1.9,1.8)

time = c()
group = c()
iter = 1
for (n_j in n){
  y_j = rgamma(n_j, shape=1, rate=theta_j_true[iter])
  time = c(time, y_j)
  group = c(group, rep(iter, n_j))
  iter = iter + 1
}

Y = cbind(group,time)

## 2. Gibbs sampler

#Data summaries
J <- length(unique(Y[,"group"]))
ysum <- c(by(Y[,"time"],Y[,"group"],sum)) # group-level means
n <- c(table(Y[,"group"])) # n_j how many response time in each group

#Hyperparameters for the priors
b1 = 1
b2 = 1/100
w = 1
#Grid values for sampling a1_grid (assume a1=20 is a reasonably large, distribution approaches 0)
alpha_grid<- 1:5000

#Initial values for Gibbs sampler
theta = theta_mean
alpha = 1
beta = 1/100
  
#First set number of iterations and burn-in, then set seed
n_iter <- 10000; burn_in <- 0.3*n_iter
set.seed(1234)

#Set null matrices to save samples
THETA <- matrix(nrow=n_iter, ncol=J)
OTHER_PAR <- matrix(nrow=n_iter, ncol=2)

#Now, to the Gibbs sampler
for(s in 1:(n_iter+burn_in)){
  #update the theta vector (all the theta_j's)
  alpha_j <- alpha + n
  beta_j <- beta + ysum
  theta <- rgamma(J, alpha_j, beta_j)
  
  #update beta
  bn1 <- b1 + J * alpha
  bn2 <- b2 + sum(theta)
  beta <- rgamma(1,bn1, bn2)
  
  #update alpha
  log_prob_alpha <- alpha_grid*J*log(beta) - J*lgamma(alpha_grid) + (alpha_grid-1)*sum(log(theta)) - w*alpha_grid
  alpha <- sample(alpha_grid,1, prob = exp(log_prob_alpha - max(log_prob_alpha)) )
  
  #save results only past burn-in
  if(s > burn_in){
    THETA[(s-burn_in),] <- theta
    OTHER_PAR[(s-burn_in),] <- c(beta, alpha)
  }
}

colnames(OTHER_PAR) <- c("beta","alpha")

## Diagnostics
library(coda)
OTHER_PAR.mcmc <- mcmc(OTHER_PAR,start=1)
summary(OTHER_PAR.mcmc)
plot(OTHER_PAR.mcmc[,"beta"])
plot(OTHER_PAR.mcmc[,"alpha"])

THETA.mcmc <- mcmc(THETA,start=1)
plot(THETA.mcmc)

autocorr.plot(THETA.mcmc)
autocorr.plot(OTHER_PAR.mcmc) # why this autocorrelation plot is so bad

# (d) Compare results from hierachical specification to the true parameter values that you set. How well does your gibbs sampler perform?
# posterior median and a 95% for theta
 
stats = c()
for (i in 1:5 ){
  theta_i = THETA[,i]
  mean.mc <- mean(theta_i)
  probs.mc <- quantile(theta_i, c(0.025,0.975))
  stats = rbind(stats, c(mean.mc, probs.mc))
}
colnames(stats) <- c("mean", "lwr","upr")
stats = data.frame(stats)
gp = 1:5
plot(gp, stats$mean, ylim=c(0, 3.0),xlab="y1 index", ylab="y", pch=16, cex=2)
# Add error bars and true values
arrows(x0=gp, y0=stats$lwr, x1=gp, y1=stats$upr, code=3, angle=90, length=0.1)
points(x=gp, y=theta_j_true,pch=16, col="red")


# (e) Compare results from hierachical specification to applying the same approach
#Hyperparameters for the priors
alpha = 1 # need to be vectors
beta = 1/100 # need to be vectors

#first set number of iterations and burn-in, then set seed
n_iter <- 10000; burn_in <- 0.3*n_iter
set.seed(1234)

#Set null matrices to save samples
THETA <- matrix(nrow=n_iter, ncol=J)

#Now, to the Gibbs sampler
for(s in 1:(n_iter+burn_in)){
  #update the theta vector (all the theta_j's)
  alpha_j <- alpha + n
  beta_j <- beta + ysum
  theta <- rgamma(J, alpha_j, beta_j)
  
  #save results only past burn-in
  if(s > burn_in){
    THETA[(s-burn_in),] <- theta
  }
}

stats = c()
for (i in 1:5 ){
  theta_i = THETA[,i]
  mean.mc <- mean(theta_i)
  probs.mc <- quantile(theta_i, c(0.025,0.975))
  stats = rbind(stats, c(mean.mc, probs.mc))
}
colnames(stats) <- c("mean", "lwr","upr")
stats = data.frame(stats)
gp = 1:5
plot(gp, stats$mean, ylim=c(0, 3.0),xlab="y1 index", ylab="y", pch=16, cex=2)
# Add error bars and true values
arrows(x0=gp, y0=stats$lwr, x1=gp, y1=stats$upr, code=3, angle=90, length=0.1)
points(x=gp, y=theta_j_true,pch=16, col="red")

