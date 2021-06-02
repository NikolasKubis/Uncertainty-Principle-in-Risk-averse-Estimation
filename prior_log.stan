

data{ // INPUTS
int<lower=0> N; 
real sigma1;
real sigma2;
}

parameters {
real<lower=0> x; 
}


model{
  x ~ lognormal(0,sigma1); 
}

generated quantities {
  real  y_tilde=lognormal_rng(log(x),sigma2);
  }
