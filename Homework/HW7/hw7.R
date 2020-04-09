library(mvtnorm)
library(BAS)
library(BMA)
## Question 1: 
# (a) using the g-prior with g=n=6, generate samples from the prior predictive distribution for a single swimmer over the 12 weeks
# create a density plot of the predictive draws (one for each week)

Y <- read.table("http://www2.stat.duke.edu/~pdh10/FCBS/Exercises/swim.dat")
Y <- t(Y)

# data summaries
n = nrow(Y)
n_swimmers = ncol(Y)
g = n
W = seq(2,12,length.out=n)
X = cbind(rep(1,n),(W-mean(W)))
p = ncol(X)

# Hyperparameters for the prior
nu_0 = 1
sigma_0_sq = 1/10
beta_0 = c(0,0)


n_iter = 1000
prior_predict = matrix(0, nrow=n_iter, ncol=nrow(X))
for (i in 1:n_iter){
  # Prior predictive distribution
  SIGMA_SQ = 1/rgamma(1,nu_0/2,nu_0*sigma_0_sq/2)
  Sigma_0 = g*SIGMA_SQ*solve(t(X)%*%X)
  BETA = rmvnorm(1, beta_0, Sigma_0)
  prior_predict[i,] = rmvnorm(1, X%*%t(BETA), SIGMA_SQ*diag(1,nrow=n)) # save different weeks
}

colnames(prior_predict) = c('Week 2','Week 4','Week 6','Week 8','Week 10','Week 12')

plot(density(prior_predict[,'Week 2']),col="red3",xlim=c(-3,3),ylim=c(0,1.7),lwd=1.5,
             main="Prior Predictive Distributions",xlab="swimming times")
legend("topright",2,c('Week 2','Week 4','Week 6','Week 8','Week 10','Week 12'),col=c("red3","blue3","orange2","black","green","purple"),lwd=2,bty="n")
lines(density(prior_predict[,"Week 4"]),col="blue3",lwd=1.5)
lines(density(prior_predict[,"Week 6"]),col="orange2",lwd=1.5)
lines(density(prior_predict[,"Week 8"]),col="black",lwd=1.5)
lines(density(prior_predict[,"Week 10"]),col="green",lwd=1.5)
lines(density(prior_predict[,"Week 12"]),col="purple",lwd=1.5)

 ## Are the values plausible? I'm not sure.


# (b) using the data, and the g-prior with g=n=6 for each swimmer, give the posterior distributions of $\beta_0$ and $\beta_1$ and $\sigma^2$ for each swimmer

# Hyperparameters for the prior
# beta_ols <- solve(t(X)%*%X)%*%t(X)%*%Y


# set number of iterations
S <- 10000
BETA = array(0,c(n_swimmers, S, p))
SIGMA_SQ = matrix(0,n_swimmers,S)

# Gibbs Sampler
for (j in 1:n_swimmers){
  beta_ols <- solve(t(X)%*%X)%*%t(X)%*%Y[,j]
  
  # sample sigma_sq
  nu_n <- nu_0 + n
  Hg <- (g/(g+1))* X%*%solve(t(X)%*%X)%*%t(X)
  SSRg <- t(Y[,j])%*%(diag(1,nrow=n) - Hg)%*%Y[,j]
  nu_n_sigma_n_sq <- nu_0*sigma_0_sq + SSRg
  sigma_sq <- 1/rgamma(S,(nu_n/2),(nu_n_sigma_n_sq/2))

  # sample beta
  mu_n <- g*beta_ols/(g+1)
  beta <- matrix(nrow=S,ncol=p)
  for(s in 1:S){
    Sigma_n <- g*sigma_sq[s]*solve(t(X)%*%X)/(g+1)
    beta[s,] <- rmvnorm(1,mu_n,Sigma_n)
  }

  BETA[j,,] = beta
  SIGMA_SQ[j,] = sigma_sq

}

# posterior summaries
beta_postmean <- t(apply(BETA,c(1,3),mean))
colnames(beta_postmean) <- c("Swimmer 1","Swimmer 2","Swimmer 3","Swimmer 4")
rownames(beta_postmean) <- c("beta_0","beta_1")
beta_postmean


## (c) For each swimmer j, plot their posterior predictive distributions for a future time T* two weeks after the last recorded observations
x_new <- matrix(c(1,(14-mean(W))),ncol=1)
post_pred <- matrix(0,nrow=n_iter,ncol=n_swimmers)
for(j in 1:n_swimmers){
  post_pred[,j] <- rnorm(n_iter,BETA[j,,]%*%x_new,sqrt(SIGMA_SQ[j,]))
}
colnames(post_pred) <- c("Swimmer 1","Swimmer 2","Swimmer 3","Swimmer 4")

plot(density(post_pred[,"Swimmer 1"]),col="red3",xlim=c(-100,100),ylim=c(0,0.05),lwd=1.5,
     main="Predictive Distributions",xlab="swimming times")
legend("topleft",2,c("Swimmer1","Swimmer2","Swimmer3","Swimmer4"),col=c("red3","blue3","orange2","black"),lwd=2,bty="n")
lines(density(post_pred[,"Swimmer 2"]),col="blue3",lwd=1.5)
lines(density(post_pred[,"Swimmer 3"]),col="orange2",lwd=1.5)
lines(density(post_pred[,"Swimmer 4"]),lwd=1.5)

## (d) compute P(Yj = max(Y1, Y2, Y3, Y4)) for each swimmer j, and based on this make a recommendation to the coach
post_pred_max <- as.data.frame(apply(post_pred,1,function(x) which(x==max(x))))
colnames(post_pred_max) <- "Swimmers"
post_pred_min$Swimmers <- as.factor(post_pred_min$Swimmers)
levels(post_pred_max$Swimmers) <- c("Swimmer 1","Swimmer 2","Swimmer 3","Swimmer 4")
table(post_pred_max$Swimmers)/n_iter

# We recommend swimmer 1 to the coach since the maximum swimming time for swimmer 1 is lowest among all four swimmers.

##########################################################################################################
## Question 2:
# a) Fit a regression model using the g-prior with g=n, nu_0=2, sigma_sq_0=1
# obtain posterior confidence intervals for all of the parameters

az = read.table("azdiabetes.dat.txt", header = TRUE)[,-8]

# glu ~ 1 + npreg + bp + skin + bmi + ped + age
### data and priors
n = nrow(az)
intercept = as.matrix(rep(1, n), ncol = 1)
colnames(intercept) = c("intercept")
rownames(intercept) = 1:nrow(az)
X = cbind(intercept, as.matrix(az[, -2])) #delete glu col
Y = as.matrix(az[, 2], ncol = 1)
p = ncol(X)
nu_0 = 2
sigma_0_sq = 1
g = n


# Hyperparameters for the prior
beta_ols = solve(t(X)%*%X)%*%t(X)%*%Y

# MC sampling
S = 10000

#sample sigma_sq
nu_n <- nu_0 + n
Hg <- (g/(g+1))* X%*%solve(t(X)%*%X)%*%t(X)
SSRg <- t(Y)%*%(diag(1,nrow=n) - Hg)%*%Y
nu_n_sigma_n_sq <- nu_0*sigma_0_sq + SSRg
sigma_sq <- 1/rgamma(S,(nu_n/2),(nu_n_sigma_n_sq/2))

#sample beta
mu_n <- g*beta_ols/(g+1)
beta <- matrix(nrow=S,ncol=p)
for(s in 1:S){
  Sigma_n <- g*sigma_sq[s]*solve(t(X)%*%X)/(g+1)
  beta[s,] <- rmvnorm(1,mu_n,Sigma_n)
}

#posterior summaries
colnames(beta) <- colnames(X)
mean_beta <- apply(beta,2,mean)
round(mean_beta,4)
round(apply(beta,2,function(x) quantile(x,c(0.025,0.975))),4) 

sigma_sq_mx <- as.matrix(sigma_sq,nrow=S)
colnames(sigma_sq_mx) <- 'sigma^2'
apply(sigma_sq_mx, 2, function(x) quantile(x,c(0.025,0.975)))

# b) Perform the model selection and averaging procedure described in section 9.3. Obtain Pr(beta_j != 0 |y),
# as well as posterior confidence intervals for all of the parameters. Compare to the results in part a).

######## Bayesian Model Selection and Averaging
Data_bas <- bas.lm(glu~npreg+bp+skin+bmi+ped+age, data=az, prior="g-prior",alpha=n,
                   n.models=2^p, update=50, initprobs="Uniform")
plot(Data_bas,which=4)
image(Data_bas)
summary(Data_bas)
coef(Data_bas)
confint(coef(Data_bas))  # confidence intervals
par(mfrow=c(2,2))
plot(coef(Data_bas), subset=2:7,ask=T)

# In part (b) since we have marginal inclusion probability (MIP) of each parameter through model averaging, we can see that the CIs for coefficients with low MIP
# are narrower and closer to 0 compared to CIs in part(a) with only a single model. The coefficients of low MIP parameters are almost 0 (which indicates we should not include those parameters, i.e.npreg, bp, skin).
#





