// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  array[N] real y;
  array[N] int type;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

transformed parameters {
  array[N] real mu;
  for(i in 1:N){
    mu[i] = alpha + beta * type[i];
  }
}

model {
  //y ~ normal(mu, sigma);
  target += normal_lpdf(alpha | 300, 50);
  target += normal_lpdf(beta | 0, 25);
  target += normal_lpdf(sigma | 0, 50) - 
  normal_lccdf(0 | 0, 50);
  
  for(i in 1:N){
    target += normal_lpdf(y[i] | mu[i], sigma);
  }
}

