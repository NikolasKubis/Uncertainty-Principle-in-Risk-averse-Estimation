

data{ // INPUTS
int<lower=0> N; 
}

parameters {
real x; 
}


model{
  x ~ exponential(0.5); 
}

generated quantities {
  real  y_tilde=normal_rng(x,3*x);
  }


