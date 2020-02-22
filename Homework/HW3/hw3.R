
# Q1: Hoff 5.2
yA <- 75.2; yB <- 77.5
sA <- 7.3; sB <- 8.1
nA <- 16; nB <- 16
mu0 <- 10; sigma20 <- 100

N = 5
Pr = c()
for (i in 0:N){
  k0 = v0 = 2^i
  
  k_nA <- k0 + nA; k_nB <- k0 + nB
  v_nA <- v0 + nA; v_nB <- k0 + nB
  
  mu_nA <- (k0*mu0 + nA*yA)/k_nA
  mu_nB <- (k0*mu0 + nB*yB)/k_nB
  sigma2_nA <- (1/v_nA) * (v0*sigma20 + sA^2*(nA-1) + (nA*k0/k_nA)*(yA-mu0)^2)
  sigma2_nB <- (1/v_nB) * (v0*sigma20 + sB^2*(nB-1) + (nB*k0/k_nB)*(yB-mu0)^2)
  
  # Posterior
  S = 10000
  tauA_postsample <- rgamma(S, k_nA/2, k_nA*sigma2_nA/2)
  thetaA_postsample <- rnorm(S, mu_nA, sqrt(1/(k_nA*tauA_postsample)))
  tauB_postsample <- rgamma(S, k_nB/2, k_nB*sigma2_nB/2)
  thetaB_postsample <- rnorm(S, mu_nB, sqrt(1/(k_nB*tauB_postsample)))
  
  # MC Sampling
  Pr <- c(Pr, mean(thetaA_postsample < thetaB_postsample))
}

plot(2^(0:(N)), Pr, "b", xlab = "k0 = v0", ylab = "Pr(theta_A < theta_B)",
     main="Probability of theta_A < theta_B with different prior using MC")

