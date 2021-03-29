data {
    int<lower=0> N;
    int<lower=0> BATCH;
    vector[BATCH] y;
    real sigma1;
    real sigma2;
}

parameters {
  vector<lower=0>[BATCH] x; // a vector of length BATCH
}

model {
  x ~ lognormal(0,sigma1);
  for (n in 1:BATCH) {
    y[n] ~ lognormal(log(x[n]), sigma2);
  }
}
