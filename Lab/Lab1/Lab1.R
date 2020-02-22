## 1. Reading and Writing Files
# for the entire help file for read.table(), type ?read.table into your console
data <- read.table(file = url("http://www2.stat.duke.edu/~pdh10/FCBS/Exercises/azdiabetes.dat"), header = TRUE)
head(data)

## 2. Objects
# Unsure about a variable tyle
str(data)
# Vectors can be assigned using the c() function

# create numeric vector of length 5
num <- c(1,2,3,4,5)

#the following are equivalent:
num + 1
num + c(1,1,1,1,1)

num + c(1,1) # why does this throw an error?

# multiplication
num2 <- num*2
num2*c(0,0,1,1,1)

# power
num2^2

# seq(): generate a numerical vector of values from the specified lower and upped bounds
seq1 <- seq(from = 1, to = 100, by = 1) # the 'by' parameters determines the interval
seq1

seq2 <- seq(from = 0, to = 1, length.out = 100)
seq2

# Retrieve the ith element in a vector
seq1[10]
seq1[c(1, 100)]
seq1[-c(25,50)] # except 25th and 50th
seq1[1] <- 0 # reassign element
seq1

# Mean, standard deviation, variance, length
seq3 <- seq(from = -3, to = 3, by = .5)
mean(seq3)
sd(seq3)
var(seq3)
length(seq3)
abs(seq3)
exp(seq3)
sqrt(seq3)
is.na(sqrt(seq3))

# Sort the entries in seq3 from greatest to least
seq3 <- sort(seq3, decreasing = TRUE)
seq3

## 3. Matrices
mat1 <- matrix(data = seq(from = 1,to = 6, by =1), nrow = 3, ncol = 2, byrow = T)
mat2 <- matrix(data = rep(x = 2, times = 6), nrow = 3, ncol = 2)
mat3 <- mat1+mat2

#transpose of matrix: t()
t(mat3)

# matrix multiplication: %*% 
# make sure dimensions agree!
dim(mat1); dim(mat2)
mat1 %*% t(mat3)

# inverse (if non-singular): solve()
mat4 <- matrix(data = c(1,2,3,4), nrow = 2, ncol = 2, byrow = F)
solve(mat4)

# obtain elements of main diagonal: diag()
diag(mat4)

# create a diagonal matrix: diag(x, nrow, ncol). Default is x=1, which creates an identity matrix
diag(4) # creates 4x4 identity matrix
diag(x = 2, nrow = 4) # creates 4x4 diagonal matrix with 2 on diagonal

# indexing matrices
mat4[1,1]
mat4[,2]
mat4[,2] <- c(0,0)
mat4

# generate large matrix
mat5 <- matrix(seq(1,100,1), nrow = 4, ncol = 25, byrow = T)

# apply(X = object, MARGIN = 1 for rows or 2 for columns, FUN = function of choice)
# find mean of every column
apply(X = mat5, MARGIN = 2, FUN = mean)

# Find the variance of each row of mat5
apply(X= mat5, MARGIN = 1, FUN = var)

# You can define you own functions and pass them into FUN
# create a function to calculate log of maximum for arbitrary x, and returns the value stored in ret
log_max <- function(x) {
  ret <- log(max(x))
  return(ret)
}

# find log of maximum for each column 
apply(X = mat5, MARGIN = 2, FUN = log_max)


## 4. Data Frames
# Extract blood pressure from data --> column
bp1 <- data$bp
str(bp1)

# Select blood pressure from data -->data frame
bp2 <- data %>% select(bp)
str(bp2)

## 5. Random Number Generation and Distribution Functions
# Default is such that the first argument specifies the number of random samples you would like
X <- rnorm(10000, mean = 0, sd = sqrt(2))
Y <- rgamma(10000, shape = 1/2, rate = 1/2)
Z <- rpois(10000, lambda = 5)

# Numerical values of the quantiles of these distributions
std_norm_qt <- qnorm(0.95) # For what value x will the CDF function of a N(0,1) R.V. return 0.95?
std_norm_cdf <- pnorm(-2) # What is the value of the CDF function of a N(0,1) R.V. at -2?
std_norm_dens <- dnorm(0.5) # What is the value of the PDF function of a N(0,1) R.V. at 0.5?

# Exercise 5 (beta distribution)
w <- rbeta(500, shape1 = 0.5, shape2 = 0.5)

## 6. Plotting

# 1). Base R Plotting
norm_samples <- rnorm(10000)
#
par(mfrow = c(2, 2)) # Set the number of rows and columns for display panels
#
hist(norm_samples,
     main = "Base R histogram", 
     xlab = "x", ylab = "Count")
#
plot(x = norm_samples, y = pnorm(norm_samples), 
     main = "Base R scatterplot", 
     xlab = "x", ylab = "Phi(x)")
#
boxplot(norm_samples, 
        main = "Base R boxplot", 
        ylab = "x")
#
plot(density(norm_samples),
     main = "Base R density", 
     xlab = "x", ylab = "Density")

# 2) ggplot2
norm_samples %>%
  data.frame(x = .) %>%
  ggplot2::ggplot() +
  geom_histogram(aes(x = x, y = ..density..), 
                 fill = "#756bb1", colour = "white", 
                 alpha = 0.5, bins = 30) +
  geom_density(aes(x = x), colour = "#756bb1") +
  geom_vline(aes(xintercept = std_norm_qt)) +
  geom_point(x = 0.5, y = std_norm_dens) +
  labs(x = "x", y = "Density", title = "ggplot density / histogram")
