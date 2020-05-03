qbeta(c(0.025,0.975),8,2)
qbeta(c(0.025,0.975),64,16)
qbeta(c(0.025,0.975),24,6)


#
#Q9/10: Boston dataset
qgamma(c(0.025,0.975),6,12)

library(MASS)
library(mvtnorm)
data(Boston)
head(Boston)

###### g-Prior: with g=n using full model
#Data summaries
n <- nrow(Boston)
p <- 3
g <- n
X <- cbind(1,as.matrix(Boston[,c("lstat","age")]))
Y <- matrix(Boston$medv,ncol=1)

#Hyperparameters for the priors
beta_ols <- solve(t(X)%*%X)%*%t(X)%*%Y
SSR_beta_ols <- (t(Y - (X%*%(solve(t(X)%*%X))%*%t(X)%*%Y)))%*%(Y - (X%*%(solve(t(X)%*%X))%*%t(X)%*%Y))
sigma_ols <- SSR_beta_ols/(n-p)
#sigma_0_sq <- sigma_ols
sigma_0_sq <- 1
nu_0 <- 2

#set number of iterations
S <- 10000

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

#####################################################################
library(BAS)
n = nrow(Boston)
Data_bas <- bas.lm(medv~lstat+age+crim+rm+black, data=Boston, prior="g-prior",alpha=n,
                   n.models=2^5, update=50, initprobs="Uniform")
plot(Data_bas,which=4)
summary(Data_bas)
coef(Data_bas)
#############################################################
sigma = matrix(c(2,1.3,1.3,1),ncol=2)
sig_0 = matrix(c(3,4.3,4.3,8),ncol=2)
t(t(sigma)+t(sig_0))

qbeta(c(0.025,0.975),5,4)
qbeta(c(0.025,0.975),5,6)
