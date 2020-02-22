
## Question 1
set.seed(123)
u <- runif(1000,0,1)
fa <- pbeta(0.4, shape1=5, shape2 = 10)
fb <- pbeta(0.75, shape1=5, shape2 = 10)
u_trun <- (fb-fa)*u+fa
theta <- qbeta(u_trun, shape1=5, shape2 = 10)

# mean of the random draw
theta_mean <- sum(theta)/1000
theta_mean #0.48
# variance of the random draw
theta_var <- sum((theta-theta_mean)^2)/(1000-1)
theta_var #0.0044

## Question 3
#(a)
aA <- 120; bA <- 10
y_a <- c(12,9,12,14,13,13,15,8,15,6)
y_b <- c(11,11,10,9,9,8,7,10,6,8,8,9,7)
A_quantile <- qgamma(c(0.025,0.975), 237, 20)
B_quantile <- qgamma(c(0.025,0.975), 125, 14)

#(b) posterior expectation of thetaB
n_0 <- seq(1,50,by=1)
exp <- (12*n_0 + 113)/(n_0 + 13)
plot(n_0,exp)
# Since e(theta_a|y) = 11.85, in order for e(theta_b) to be close to 11.85, we need n_0 to be 
# close to or even more than more than 50, which means we need a very strong prior belief of theta_b

#(c)
# It's possible that the parameters of B are quite different from A, so we should view them as independent.


## Question 4
theta1.mc <- rbeta(5000,58,44)
theta2.mc <- rbeta(5000,31,21)
mean(theta1.mc<theta2.mc) #0.6408


## Question 5
#(a)
m = c(10,100,1000)
theta.mc10 <- rgamma(10, 68, 45)
theta.mc100 <- rgamma(100, 68, 45)
theta.mc1000 <- rgamma(1000, 68, 45)

# mark those points
mean.mc10 <- mean(theta.mc10)
mean.mc100 <- mean(theta.mc100)
mean.mc1000 <- mean(theta.mc1000)

probs.mc10 <- quantile(theta.mc10, c(0.025,0.975))
probs.mc100 <- quantile(theta.mc100, c(0.025,0.975))
probs.mc1000 <- quantile(theta.mc1000, c(0.025,0.975))


sample_plot <- function(samples){
  set.seed(123)
  theta.mc <- rgamma(samples, 68, 45)
  mean.mc <- mean(theta.mc)
  probs.mc <- quantile(theta.mc, c(0.025,0.975))
  
  plot(density(theta.mc), xlab="Number of children", ylab="density", col="black",lty=1,main=paste(samples, " MC Samples"))
  abline(v = probs.mc, col=c("red"), lty=2, lwd=c(1))
  abline(v = mean.mc, col="blue", lty=3,lwd=1)
  legend("topright", legend=c("Density", paste("95% Quantile: ", round(probs.mc[[1]],2), round(probs.mc[[2]],2)), paste("mean: ", round(mean.mc,2))),
         col=c("black", "red","blue"), lty=1:3, cex=0.8, bg="transparent", bty = "n")
}

sample_plot(10)
sample_plot(100)
sample_plot(1000)

# true mean and quantile-based CI
true_mean <- 68/45
true_quantile <- qgamma(c(0.025,0.975),68,45)
true_mean #1.51
true_quantile #1.17, 1.89
# As long as number of samples for MC approximation is big enough, the result of mean and 95% CI is very close to the true mean and CI.


#(b) posterior probability theta2<1.5
mean(theta.mc10 < 1.5) #0.4
mean(theta.mc100 < 1.5) #0.48
mean(theta.mc1000 < 1.5) #0.509

#(c)

var.estimate <- 68/(45^2)
s.min = var.estimate/(0.001/2)^2
s.min

## Question 6: Hoff 4.8
#a. obtain 5000 samples from th posterior predictive distribution (negative binomial)
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
a <- 2
b <- 1
sum_A <- sum(man_A)
sum_B <- sum(man_B)
n_A <- length(man_A)
n_B <- length(man_B)

set.seed(123)
Y_A <- rnbinom(5000, size = a+sum_A, prob = 1-1/(b+n_A+1))
Y_B <- rnbinom(5000, size = a+sum_B, prob = 1-1/(b+n_B+1))

# plot the MC approximations to these two posterior predictive dist
hist(Y_A, breaks = seq(-0.5, max(Y_A)+0.5,1), main = 'MC Generated Posterior Predictive Y_A')
hist(Y_B, breaks = seq(-0.5, max(Y_B)+0.5,1), main = 'MC Generated Posterior Predictive Y_B')

#b. find 95% quantile-based posterior CI for theta_B-theta_A and Y_B-Y_A
set.seed(123)
theta_A <- rgamma(5000,a+sum_A, b+n_A)
theta_B <- rgamma(5000,a+sum_B, b+n_B)

quantile(theta_B - theta_A, c(0.025,0.975))
quantile(Y_B-Y_A, c(0.025,0.975))

hist(theta_B-theta_A)
hist(Y_B-Y_A)
# difference between these two populations

#c. 
hist(man_B, breaks = seq(-0.5, max(man_B)+0.5,1), xlab = 'Number of children', main = 'Empirical Distribution of Group B')
set.seed(123)
pos_sam <- rpois(200,1.4)
hist(pos_sam, breaks = seq(-0.5, max(pos_sam)+0.5,1), xlab = 'Number of children', main = 'Poisson Distribution with mean 1.4')

#d.
count_zeros <- c()
count_ones <- c()
for (theta in theta_B){
  pos_rv <- rpois(218, theta)
  count_zeros <- c(count_zeros, sum(pos_rv==0))
  count_ones <- c(count_ones, sum(pos_rv==1))
}
plot(x=count_zeros, y=count_ones)
points(x=sum(man_B == 0), y=sum(man_B==1),pch=16, col="red")
