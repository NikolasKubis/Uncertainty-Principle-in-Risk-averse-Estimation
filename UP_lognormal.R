
# The sampling of the prior is made using the prior stan.model.


#Loading the libs..
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(plyr)
library(base)
library(rstan)
library(grid)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


#The maximally_risk_averse_estimator.
maximally_risk <- function(jm) {
  m1=colMeans(jm)
  m2=colMeans(jm^2)
  m3=colMeans(jm^3)
  return ((m3-m2*m1)/(2*(m2-m1^2)))
}


#The risk_averse_estimator (inputs = samples from one posterior and 
# parameter mu) 
#Takes as inputs the Joint Matrix and the parameter mu
# Outputs a vector with length equal to the length of the observables
# -equivalently- to the length of the prior samples.

risk_averse <- function(jm, mu) {
  m1 <- colMeans(jm)
  m2 <- colMeans(jm^2)
  m3 <- colMeans(jm^3)
  
  return( (m1 - 2*mu*(m1*m2 - m3))/(1 + 4*mu*(m2 - m1^2)) )
}

max_eig<-function(jm) {
  
  m1 <- colMeans(jm)
  m2 <- colMeans(jm^2)
  me<-m2-m1^2
  return(me)
}

#Computes the conditinal mean. 
# Takes as input the joint matrix and returns a vector with length 
# equal to the number of observables=number of samples from prior.
conditional_mean<-function(jm){
  return(colMeans(jm)) 
  
}

# This function creates the mse as a function of mu. 
# call with mse=t(t(as.matrix(sapply(mu,FUN=multiplier,joint_matrix,x))))
multiplier <-function(x,A,v)
{
  m1=colMeans(A)
  m2=colMeans(A^2)
  m3=colMeans(A^3)
  X=(m1 - 2*x*(m1*m2 - m3))/(1 + 4*x*(m2 - m1^2)) 
  mv=mean((X-v)^2)
  return (mv)
}



#B=matrix(0,nrow=10000,ncol = 8000)
risky <-function(x,jm)
{
  m1=colMeans(jm)
  
  m2=colMeans(jm^2)
  
  m3=colMeans(jm^3)
  
  XX=(m1 - 2*x*(m1*m2 - m3))/(1 + 4*x*(m2 - m1^2)) 
  
  oros=colMeans((jm-rep(c(as.matrix(XX)),each=nrow(jm)))^2)
  
  MM=(matrix(c(oros), nrow=6000, ncol=length(c(oros)), byrow=TRUE))
  
  h=colMeans((((jm-rep(c(as.matrix(XX)),each=nrow(jm)))^2 - MM ) )^2)
  
  r=mean(h) 
  
  return(r)
}






#DEFINING THE PRIOR 
# This stan model samples from the prior and also generates output samples
# normally distributed with state dependency.
model_prior <- stan_model('prior.stan') #

skewness_prior <- 0.85 #change -> 1.6, 1.75, 1.9
                       # change ->  0.060, 0.085, 0.1

skewness_posterior <- 0.25 #0.25



N=1
stan_data <- list( 
  'N' = N,
  'sigma1'=skewness_prior, 
  'sigma2'=skewness_posterior)


fit_prior <- sampling( model_prior, data = stan_data, chains=4, iter=2500, warmup =500,
                      control = list('adapt_delta' = 0.99, 'max_treedepth' = 10))


#----------------------------------------leave it or it dies...softly ----------------------------------------

samples_prior <- as.data.frame(extract(fit_prior, permuted=TRUE)) 

y_out=samples_prior[,2]

x=samples_prior[,1]

nSamples_prior <- length(y_out) #or length(x)

model_posterior <- stan_model('posterior.stan')



stan_data = list('y' = y_out, # feed in the observations generated randomly from the stan_prior.
                 'N' = N,
                 'BATCH' = nSamples_prior, # length of the observable vector.
                 'sigma1' = skewness_prior,
                 'sigma2' = skewness_posterior)


fit_posterior <- sampling(model_posterior, data = stan_data,
                iter = 3500, warmup = 2000, chains = 4,
                control = list('adapt_delta' = 0.99, 'max_treedepth' = 10))



#-------------------------------------leave it or die....softly------------------------------------------



samples_posterior <- as.data.frame(extract(fit_posterior, permuted=TRUE))


joint_matrix<-data.matrix(samples_posterior[,1:(ncol(samples_posterior)-1)])

# Each column of this matrix corresponds to an observable. 
# and defines f(x | the observable). 
# So each column contains samples from f(x | the observable).



x_max_risk=maximally_risk(joint_matrix)
cond_mean=conditional_mean(joint_matrix)
# Compute the mse width
mean(abs(x_max_risk-cond_mean))
#

x_risk=matrix(0,nSamples_prior)

conditional_risk=matrix(0,nSamples_prior); # The predictive variance

#mse=matrix(0,length(mu))


# Calling the multiplier to compute the mse as a function of mu.
#mse=t(t(as.matrix(sapply(mu,FUN=multiplier,joint_matrix,x))))


n_mu=50
mu=(matrix(0:49,n_mu))*1


#MSE
mse=t(t(as.matrix(sapply(mu,FUN=multiplier,joint_matrix,x)))) #5'
plot(mu,mse,type='l')

#RISK

risk=t(t(as.matrix(sapply(mu,FUN=risky,joint_matrix)))) #12'

plot(mu,risk,type='l')


#PRODUCT
product=mse*risk

plot(mu,product,type='l')
  


#####################################       PLOTS      ###############################################

#MSE as a function of mu (if imported from R.data use mse$V1)
A=cbind(mse,mu) # merging the vectors into a matrix
newdata_mse<-as.data.frame(A) # creating data structure data.frame

# in order to make the artistic plot.
ggplot(newdata_mse, aes(x=mu, y=mse))+
  geom_line(size=0.7)



#RISK as a function of mu
B=cbind(risk,mu) # merging the vectors into a matrix
newdata_risk<-as.data.frame(B)

ggplot(newdata_risk, aes(x=mu, y=risk))+
  geom_line(size=0.7)



#Product as a function of mu

C=cbind(product,mu)
newdata_product<-as.data.frame(C)

ggplot(newdata_product, aes(x=mu, y=product))+
geom_line(size=1) 



