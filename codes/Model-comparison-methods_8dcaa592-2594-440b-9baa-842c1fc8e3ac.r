library(ggplot2)
library(dplyr)

# Out of 100 trials, correct responses are 87

x <- c(87,80,77,66,90)
n <- 100

# Model 1
# x_i ~ Binomial(n=100,theta)
# theta ~ Beta(80,20)

# Model 2
# x_i ~ Binomial(n=100,theta)
# theta ~ Beta(50,50)

# Cross-validation
x <- c(87,80,77,66,90,88)
n <- 100

elpd_m1 <- c()
elpd_m2 <- c()

for(i in 1:3){
  xtrain <- x[- c(i*2-1,i*2)]
  xtest <- x[c(i*2-1,i*2)]
  # Train model 1,2
  post_samples_m1 <- rbeta(2000,80+sum(xtrain),20+sum(n-xtrain-1))
  post_samples_m2 <- rbeta(2000,50+sum(xtrain),50+sum(n-xtrain-1))
  lpd_m1 <- 0
  lpd_m2 <- 0
  for(j in 1:length(xtest)){
    lpd_m1_j <- log(mean(dbinom(xtest[j],size=n,post_samples_m1)))
    lpd_m2_j <- log(mean(dbinom(xtest[j],size=n,post_samples_m2)))
    lpd_m1 <- lpd_m1+lpd_m1_j
    lpd_m2 <- lpd_m2+lpd_m2_j
  }
  elpd_m1 <- c(elpd_m1,lpd_m1)
  elpd_m2 <- c(elpd_m2,lpd_m2)
}
  
# ELPD for m1
sum(elpd_m1)

# ELPD for m2
sum(elpd_m2)

# Diff of elpds
sum(elpd_m1-elpd_m2)

####################################

# Importance density is Beta (1,1)
hist(rbeta(10000,1,1))

theta_samples <- rbeta(100000,1,1)

# ML for model m1
likelihood <- dbinom(sum(x),size=100*length(x),prob=theta_samples)
prior <- dbeta(theta_samples,80,20)
imp.density <- dbeta(theta_samples,1,1)
ML_m1 <- mean(likelihood*prior/imp.density)

# ML for model m2
likelihood <- dbinom(sum(x),size=100*length(x),prob=theta_samples)
prior <- dbeta(theta_samples,50,50)
imp.density <- dbeta(theta_samples,1,1)
ML_m2 <- mean(likelihood*prior/imp.density)

ML_m1
ML_m2

bf <- ML_m1/ML_m2
bf
