require(tidyverse)
require(rstanarm)
require(magrittr)
library(ggplot2)
library(mlmRev)
library(tidybayes)
library(ggstance)
library(modelr)
require(mlmRev)
require(tidybayes)

data(Gcsemv, package = "mlmRev")
summary(Gcsemv)

# Make male the reference category and rename variable
Gcsemv$female <- relevel(Gcsemv$gender, "M")
# Use only total score on coursework paper
GCSE <- subset(x = Gcsemv,
               select = c(school, student, female, course))
set.seed(40)
samps <- sample(1:nrow(GCSE), 20, replace = F)
(GCSE_test <- GCSE[samps,])

GCSE <- GCSE[-samps,]

# Count unique schools and students
m <- length(unique(GCSE$school)) #73 unique schools
N <- nrow(GCSE) #1885 rows

GCSE %>% 
  dplyr::group_by(school) %>% 
  na.omit()%>%
  dplyr::summarise(avg_course = mean(course)) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = avg_course))+
  geom_histogram()

# Pooled and unpooled framework
pooled <- stan_glm(course ~ 1 + female, data = GCSE, refresh = 0)
unpooled <- stan_glm(course ~ -1 + school + female,data=GCSE, refresh = 0)

#################################################################################
## Model 1: varying intercept model with no predictor
mod1 <- stan_lmer(formula = course ~ 1 + (1 | school),
                  data = GCSE,
                  seed = 349,
                  refresh = 0)

prior_summary(object = mod1)
sd(GCSE$course, na.rm = T)

# posterior estimates (median and mean)
print(mod1, digits = 3)
summary(mod1,
        pars = c("(Intercept)", "sigma", "Sigma[school:(Intercept),(Intercept)]"),
        probs = c(0.025, 0.975),
        digits = 3)

# extract posterior draws
mod1_sims <- as.matrix(mod1)
dim(mod1_sims)

par_names <- colnames(mod1_sims)
head(par_names)
tail(par_names)


# obtain draws for mu_theta
mu_theta_sims <- as.matrix(mod1, pars = "(Intercept)")

# obtain draws for each school's contribution to intercept
theta_sims <- as.matrix(mod1,
                        regex_pars ="b\\[\\(Intercept\\) school\\:")


# to finish: obtain draws for sigma and tau^2
sig_sims <- as.matrix(mod1,
                      pars = "sigma")
tau2_sims <- as.matrix(mod1,
                       pars = "Sigma[school:(Intercept),(Intercept)]")

# 73 total varying intercepts and store the posterior mean 95% credible intervals for each intercept
int_sims <- as.numeric(mu_theta_sims) + theta_sims

# posterior mean
int_mean <- apply(int_sims, MARGIN = 2, FUN = mean)

# credible interval
int_ci <- apply(int_sims, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))
int_ci <- data.frame(t(int_ci))

# combine into a single df
int_df <- data.frame(int_mean, int_ci)
names(int_df) <- c("post_mean","Q2.5", "Q97.5")

# sort DF according to posterior mean
int_df <- int_df[order(int_df$post_mean),]

# create variable "index" to represent order
int_df <- int_df %>% mutate(index = row_number())

# plot posterior means of school-varying intercepts, along with 95 CIs
ggplot(data = int_df, aes(x = index, y = post_mean))+
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5))+
  scale_x_continuous("Index", breaks = seq(0,m, 5)) +
  scale_y_continuous(expression(paste("varying intercept ", theta[j])))

##################################################################################
## Model 2: Varying intercept model with a single individual-level predictor



