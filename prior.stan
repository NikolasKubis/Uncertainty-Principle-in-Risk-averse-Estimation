// Creating a model using stan.
// The first thing is that the stan code is structured with programm blocks.
//The three blocks are the data{}, the parameters{} and the model{}.
//The data block refers to variables that are read in as data (inputs). 
// For example , in my case I need to pass in 1) the value of alpha which is parametrized the model..
// 2) I need some observations from the model i.e to evaluate the model to N points.

//The variables declared as parameters are those that are sampled by stan. might also be vectors !

//Finally the model contains the distribution.

//In this file we create a simple model to examine how things work. 

//starting with the data 

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


