library(ggplot2)


gradient <- function(mu,sigma,y,n,m,s,a,b){
  grad_mu <- (((n*mu)-sum(y))/(sigma^2))+((mu-m)/(s^2))
  grad_sigma <- (n/sigma)-(sum((y-mu)^2)/(sigma^3))+((sigma-a)/(b^2))
  return(c(grad_mu,grad_sigma))
}

V <- function(mu,sigma,y,n,m,s,a,b){
  nlpd <- -(sum(dnorm(y,mu,sigma,log=T))+dnorm(mu,m,s,log=T)+dnorm(sigma,a,b,log=T))
  nlpd
}


y <- rnorm(500,901,25) 
n <- length(y)

#Model parameters
m <- 1000
s <- 100
a <- 20
b <- 5
# HMC sampler
# Internal parameters
# Step size
step <- 0.02
# Number of leapfrog steps
L <- 12

hmc <- function(step=0.02,L=12,nsamp){
  # Markov chain
  nsamp <- 8000
  mu_chain <- rep(NA,nsamp)
  sigma_chain <- rep(NA,nsamp)
  reject <- 0
  #Initialization of Markov chain
  mu_chain[1] <- rnorm(1,1000,10)
  sigma_chain[1] <- rnorm(1,10,1)
  #Evolution of Markov chain
  i <- 1
  
  chain_ev <- data.frame(sample=i,pos=mu_chain[i])
  
  while(i < nsamp){
    q <- c(mu_chain[i],sigma_chain[i]) # Current position of the particle
    p <- rnorm(length(q),0,1) # Generate random momentum at the current position
    current_q <- q
    current_p <- p
    current_V = V(current_q[1],current_q[2],y,n,m,s,a,b) # Current potential energy
    current_T = sum(current_p^2)/2 # Current kinetic energy
    # Take L leapfrog steps
    for(l in 1:L){
      # Change in momentum in 'step/2' time
      p <- p-((step/2)*gradient(q[1],q[2],y,n,m,s,a,b))
      # Change in position in 'step' time
      q <- q + step*p
      # Change in momentum in 'step/2' time
      p <- p-((step/2)*gradient(q[1],q[2],y,n,m,s,a,b))
    }
    proposed_q <- q
    proposed_p <- p
    proposed_V = V(proposed_q[1],proposed_q[2],y,n,m,s,a,b) # Proposed potential energy
    proposed_T = sum(proposed_p^2)/2 # Proposed kinetic energy
    accept.prob <- min(1,exp(current_V+current_T-proposed_V-proposed_T))
    # Accept/reject the proposed position q
    if(accept.prob>runif(1,0,1)){
      mu_chain[i+1] <- proposed_q[1]
      sigma_chain[i+1] <- proposed_q[2]
      i <- i+1
      chain_ev <- rbind(chain_ev,data.frame(sample=i,pos=mu_chain[i]))
      #print(ggplot(chain_ev,aes(x=sample,y=pos))+geom_point(size=1.5,color="red"))
    }else{
      reject <- reject+1
    }
  }
  return(list(mu_chain=mu_chain,sigma_chain=sigma_chain,reject=reject))
}


chain1 <- hmc(step=0.05,L=22)
chain2 <- hmc(step=0.05,L=22)
chain3 <- hmc(step=0.05,L=22)
chain4 <- hmc(step=0.05,L=22)

nsamp <- 8000
rejection_count <- chain1$reject
rate <- (rejection_count/(nsamp+rejection_count))*100
rate

chains <- data.frame(matrix(nrow = nsamp,ncol=4))
colnames(chains) <- c("chain 1","chain 2","chain 3","chain 4")
chains$`chain 1`<- chain1$mu_chain
chains$`chain 2`<- chain2$mu_chain
chains$`chain 3` <- chain3$mu_chain
chains$`chain 4` <- chain4$mu_chain

chains$sample <- 1:nsamp
library(ggplot2)
library(reshape2)
ggplot(subset(melt(chains,id="sample"),sample>1000),aes(x=sample,y=value,group=variable,color=variable))+
  geom_line()

chain.m <- melt(chains,id="sample")
chain.m <- rbind(chain.m,data.frame(sample=1:nsamp*4,
                                    variable=rep("0",nsamp*4),
                                    value=rgamma(nsamp*4,23,20)))
chain.m$type <- factor(c(rep("MCMC",nsamp*4),rep("Analytical",nsamp*4)))

ggplot(chain.m,aes(x=value,group=type,color=type))+
  geom_density()


##########################################

library(brms)

m1 <- brm(y~1,data=data.frame(y=y),
          family = gaussian(),
          prior = c(prior(normal(1000,100),class=Intercept),
                    prior(normal(20,5),class=sigma)),
          chains=4,
          cores=4)

samples <- posterior_samples(m1)
hist(samples$b_Intercept)
hist(samples$sigma)

summary(m1)

save(m1,file="brms_fit_object.Rda")
library(bayesplot)
plot(m1,variable = "b_Intercept")

load("brms_fit_object.Rda")
plot(m1)
