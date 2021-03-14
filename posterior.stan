data {
    int<lower=0> N;
    int<lower=0> BATCH;
    vector[BATCH] y;
    real sigma;
}

parameters {
  vector<lower=0>[BATCH] x; 
}

model {
  x ~ exponential(sigma);
  for (n in 1:BATCH) {
    y[n] ~ normal(x[n], 3*x[n]);
  }
}

