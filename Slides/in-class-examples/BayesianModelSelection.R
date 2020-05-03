
###########################################################################
###########################################################################
############################# Health plans ################################
###########################################################################
###########################################################################


###### Clear environment and load libraries
rm(list = ls())
require(lattice)
library(pls)
library(calibrate)
library(BAS)
library(BMA)
library(mvtnorm)



###### Data
Data <- read.table("costs.txt",header=TRUE)[,-9]



###### Very basic EDA
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
pairs((Data), panel=panel.smooth,diag.panel=panel.hist)
#Check correlation
levelplot(cor(Data))



###### g-Prior: with g=n using full model
#Data summaries
n <- 29
p <- 8
g <- n
X <- cbind(1,as.matrix(Data[,c("RXPM","GS","RI","COPAY","AGE","F","MM")]))
Y <- matrix(Data$COST,ncol=1)

#Hyperparameters for the priors
beta_ols <- solve(t(X)%*%X)%*%t(X)%*%Y
SSR_beta_ols <- (t(Y - (X%*%(solve(t(X)%*%X))%*%t(X)%*%Y)))%*%(Y - (X%*%(solve(t(X)%*%X))%*%t(X)%*%Y))
sigma_ols <- SSR_beta_ols/(n-p)
#sigma_0_sq <- sigma_ols
sigma_0_sq <- 1/100
nu_0 <- 1

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



######## Bayesian Model Selection and Averaging
#library(BAS)
Data_bas <- bas.lm(COST~RXPM+GS+RI+COPAY+AGE+F+MM, data=Data, prior="g-prior",alpha=n,
                    n.models=2^p, update=50, initprobs="Uniform")
plot(Data_bas,which=4)
image(Data_bas)
summary(Data_bas)
coef(Data_bas)
par(mfrow=c(3,3))
plot(coef(Data_bas), subset=2:8,ask=T)


######## Performs Bayesian simultaneous variable selection and outlier identification
#library(BMA)
Data_bma <- MC3.REG(Data[, 1], as.matrix(Data[, -1]),num.its=10000,outliers=TRUE, 
                   M0.out=rep(FALSE, 29), outs.list=1:29, M0.var=rep(TRUE,7))
summary(Data_bma)


# Five outliers identified using BMA under the model with RI included... - 5,10,12,19 and 21
# Observation 19 is always an outlier.
# Remove observation, fit again using g-Prior: with g=n using smaller model





