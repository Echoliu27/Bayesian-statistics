
## Question 1
theta_list = c(0.25,0.46,0.5, 0.54)
a = b = 1


for (theta in theta_list){
  n = 1000
  y = c()
  bf10_vec = c()
  for(i in 1:n){
    y_new = rbinom(1,1,theta)
    y = c(y, y_new)
    Y = sum(y)
    bf10 = beta(a+Y, b+i-Y)/(0.5^i)
    bf10_vec = c(bf10_vec, bf10)
  }
  plot(x=1:n, y=bf10_vec, main=theta, type="l")
}

## Question 2
