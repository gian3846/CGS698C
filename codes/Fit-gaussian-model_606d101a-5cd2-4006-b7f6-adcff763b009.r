library(rstan)

y <- rnorm(500,320,23)
dat <- data.frame(type=rep(c(+1,-1),250),y=y)
head(dat)

ls_dat <- list(N=length(dat$y),y=dat$y,type=dat$type)
ls_dat

fit_normal <- stan("Simple-gaussian-mdoel.stan",
                   data=ls_dat,
                   chains=4,cores=4)

samples <- extract(fit_normal)
hist(samples$alpha)
hist(samples$sigma)

