

data{ // INPUTS
int<lower=0> N; //number of evaluations of the model (declared as always positive integer !)
}

parameters {
real x; // parameter w.r.t. we want to sample.
}


model{
  x ~ exponential(0.5); 
}

generated quantities {
  real  y_tilde=normal_rng(x,3*x);
  }


