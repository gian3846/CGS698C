library(ggplot2)
library(dplyr)
library(brms)

load("Data/df_pupil.rda")
head(df_pupil)

hist(df_pupil$p_size)

plot(df_pupil$load,df_pupil$p_size)
ggplot(df_pupil,aes(x=load,y=p_size))+geom_point()+
  geom_smooth(method="lm")
lm(p_size~load, df_pupil)
hist(df_pupil$p_size)

# load affects pupil size
# Average pupil size changes as the load changes
# Hypothesis 1: 
# The effect of attentional load on pupil size is [10,100]
# Hypthesis 2: 
# The effect of attentional load on pupil size can be [-200,200]


# Model 1:

# p_size_i ~ Normal(mu_i,sigma)
# mu_i = alpha + beta*load_i
# alpha ~ Normal(500,200)
# beta ~ Normal(10,40) [lb=10,ub=100]
# sigma ~ Normal_+(0,100)

# Model 2:
# beta ~ Normal(0,100) [lb=-200,ub=200]

# Model 1:
fit_pupil <- brm(p_size ~ 1 + load,
                 data = df_pupil,
                 family = gaussian(),
                 prior = c(
                   prior(normal(500, 200), class = Intercept),
                   prior(normal(10,40), class = b, lb=10,ub=100),
                   prior(normal(0,100), class=sigma)
                 ),
                 chains = 4,
                 iter = 2000,
                 warmup = 1000,
                 cores=4)

save(fit_pupil,file="Fit_pupil_M1.Rda")
plot(fit_pupil)
summary(fit_pupil)

# Model 2:

fit_pupil_2 <- brm(p_size ~ 1 + load,
                 data = df_pupil,
                 family = gaussian(),
                 prior = c(
                   prior(normal(500, 200), class = Intercept),
                   prior(normal(0,100), class = b, lb=-200,ub=200),
                   prior(normal(0,100), class=sigma)
                 ),
                 chains = 4,
                 iter = 2000,
                 warmup = 1000,
                 cores=4)

plot(fit_pupil_2)
summary(fit_pupil_2)

#
load("FittedModels/fit_pupil1.rda")
load("FittedModels/fit_pupil2.rda")

# ML_m1
#Importance density: 
#   alpha ~ N(500,100)
#   beta ~ N(0,40)
# Priors:
#   alpha ~ N(500,200)
#   beta ~ N(10,40)
# size ~ Normal(mu_i,sigma) where mu_i = alpha + beta load_i
library(truncnorm)

nsamp <- 10000
alpha_samples <- rnorm(nsamp,500,100)
beta_samples <- rtruncnorm(nsamp,a=-Inf,b=Inf,mean=0,sd=40)
sigma_samples <- rtruncnorm(nsamp,a=0,b=Inf,mean=0,sd=50)

lklpr <- rep(NA,nsamp)
for(i in 1:nsamp){
  mu <- alpha_samples[i]+beta_samples[i]*df_pupil$load
  likelihood <- prod(dnorm(df_pupil$p_size,mu,sigma_samples[i]))
  prior <- dnorm(alpha_samples[i],500,200)*dnorm(beta_samples[i],10,40)*dnorm(sigma_samples,0,100)
  imp <- dnorm(alpha_samples,500,100)*dnorm(beta_samples,0,40)*dnorm(sigma_samples,0,50)
  lklpr[i] <- likelihood*prior/imp
}

ml_m1 <- mean(lklpr)

nsamp <- 10000
alpha_samples <- rnorm(nsamp,500,100)
beta_samples <- rtruncnorm(nsamp,a=-Inf,b=Inf,mean=0,sd=40)
sigma_samples <- rtruncnorm(nsamp,a=0,b=Inf,mean=0,sd=50)

lklpr <- rep(NA,nsamp)
for(i in 1:nsamp){
  mu <- alpha_samples[i]+beta_samples[i]*df_pupil$load
  likelihood <- prod(dnorm(df_pupil$p_size,mu,sigma_samples[i]))
  prior <- dnorm(alpha_samples[i],500,200)*dnorm(beta_samples[i],0,100)*dnorm(sigma_samples,0,100)
  imp <- dnorm(alpha_samples,500,100)*dnorm(beta_samples,0,40)*dnorm(sigma_samples,0,50)
  lklpr[i] <- likelihood*prior/imp
}

ml_m2 <- mean(lklpr)
ml_m2/ml_m1

ml_m1/ml_m2
######################

dat <- read.table("Data/temperature.txt",sep=",",header=F)
head(dat)
colnames(dat) <- c("date_id","year","month","day","day_id","temp")
  
hist(dat$temp)

dat$year.m <- dat$year - mean(dat$year)

plot(dat$year.m,dat$temp)

###### 
# Model
# Likelihood assumption
# temp_i ~ Normal(mu_i, sigma)
# mu_i = alpha + beta*year_i
# Priors
# alpha ~ Normal(0,1.5)
# beta ~ Normal(0, 0.05)
# sigma ~ Normal(0,1)


fit_temp <- brm(temp ~ 1 + year.m^2,
                 data = dat,
                 family = gaussian(),
                 prior = c(
                   prior(normal(0, 1.5), class = Intercept),
                   prior(normal(0,0.05), class = b, coef = year.m),
                   prior(normal(0,1), class=sigma)
                 ),
                 chains = 4,
                 iter = 2000,
                 warmup = 1000,
                 cores=4)

plot(fit_temp)
summary(fit_temp)

save(fit_temp,file="FittedModels/temperature.rda")

load("FittedModels/temperature.rda")
plot(fit_temp)
summary(fit_temp)

################################

# Hypothesis: Temperature increases at a faster rate with respect to time after 1950 compared to before 1950


dat <- read.table("Data/temperature.txt",sep=",",header=F)
head(dat)
colnames(dat) <- c("date_id","year","month","day","day_id","temp")

hist(dat$temp)

dat$year.m <- dat$year - mean(dat$year)
dat$period <- ifelse(dat$year<1950,-1,+1)

library(ggplot2)
ggplot(dat,aes(x=year,y=temp,group=period,color=period))+
  geom_point()

###### 
# Model
# Likelihood assumption
# temp_i ~ Normal(mu_i, sigma)
# mu_i = alpha + beta_1*year_i + beta_2*period_i + beta_3*year_i*period_i
# Priors
# alpha ~ Normal(0,1.5)
# beta_1 ~ Normal(0, 0.05)
# beta_2 ~ Normal(0, 1.5)
# beta_3 ~ Normal(0,0.05)
# sigma ~ Normal(0,1)


fit_temp <- brm(temp ~ 1 + year.m + period + year.m*period,
                data = dat,
                family = gaussian(),
                prior = c(
                  prior(normal(0, 1.5), class = Intercept),
                  prior(normal(0,0.05), class = b, coef = year.m),
                  prior(normal(0,1.5), class = b, coef = period),
                  prior(normal(0,0.05), class = b, coef = year.m:period),
                  prior(normal(0,1), class=sigma)
                ),
                chains = 4,
                iter = 2000,
                warmup = 1000,
                cores=4)

plot(fit_temp)
summary(fit_temp)

save(fit_temp,file="FittedModels/temperature_interaction.rda")

load("FittedModels/temperature_interaction.rda")
plot(fit_temp)
summary(fit_temp)


################################################


x <- rlnorm(10000,meanlog = 6, sdlog=0.5)
hist(x)

hist(log(x))

#######################

















