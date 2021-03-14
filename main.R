# setwd(" ")


# The sampling of the posterior is made using the prior stan.model.
# Subsequently the sampler for the posterior is made by utilizing a second stan.model.


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
# Outputs a vector with length equal to the length of the observables\\
# or equivalently to the length of the prior samples.

risk_averse <- function(jm, mu) {
  m1 <- colMeans(jm)
  m2 <- colMeans(jm^2)
  m3 <- colMeans(jm^3)
  
  return( (m1 - 2*mu*(m1*m2 - m3))/(1 + 4*mu*(m2 - m1^2)) )
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



risky <-function(x,jm)
{
  m1=colMeans(jm)
  
  m2=colMeans(jm^2)
  
  m3=colMeans(jm^3)
  
  XX=(m1 - 2*x*(m1*m2 - m3))/(1 + 4*x*(m2 - m1^2)) 
  
  oros=colMeans((jm-rep(c(as.matrix(XX)),each=nrow(jm)))^2)
  
  MM=(matrix(c(oros), nrow=10000, ncol=length(c(oros)), byrow=TRUE))
  
  h=colMeans((((jm-rep(c(as.matrix(XX)),each=nrow(jm)))^2 - MM ) )^2)
  
  r=mean(h) 
  
  return(r)
}







#DEFINING THE PRIOR 
# This stan model samples from the prior and also generates ranrom output samples
# normally distributed with state dependency.
model_prior <- stan_model('UP.stan') #


N=1
stan_data <- list( 
  'N' = N
)


fit_prior <- sampling( model_prior, data = stan_data, chains=4, iter=2500, warmup =500,
                      control = list('adapt_delta' = 0.99, 'max_treedepth' = 10))

## samples=chains*(iter-warmup)

#----------------------------------------leave it or die...softly ----------------------------------------

samples_prior <- as.data.frame(extract(fit_prior, permuted=TRUE)) 

y_out=samples_prior[,2]

x=samples_prior[,1]

nSamples_prior <- length(y_out) #or length(x)

model_posterior <- stan_model('exp2.stan')

expo_rate <- 0.5

stan_data = list('y' = y_out, # feed in the observations generated randomly from the stan_prior.
                 'N' = N,
                 'BATCH' = nSamples_prior, # length of the observable vector.
                 'sigma' = expo_rate)


fit_posterior <- sampling(model_posterior, data = stan_data,
                iter = 4500, warmup = 2000, chains = 4,
                control = list('adapt_delta' = 0.99, 'max_treedepth' = 10))



#-------------------------------------leave it or die..softly------------------------------------------



samples_posterior <- as.data.frame(extract(fit_posterior, permuted=TRUE))


joint_matrix<-data.matrix(samples_posterior[,1:(ncol(samples_posterior)-1)])



x_max_risk=maximally_risk(joint_matrix)
cond_mean=conditional_mean(joint_matrix)
# Compute the mse width
mean(abs(x_max_risk-cond_mean))
#

x_risk=matrix(0,nSamples_prior)

conditional_risk=matrix(0,nSamples_prior); # The predictive variance


n_mu=31
mu=matrix(0:30,n_mu)*0.1
#v=c(1,1,1,1)
mse=t(t(as.matrix(sapply(mu,FUN=multiplier,joint_matrix,x))))
plot(mu,mse)

risk=t(t(as.matrix(sapply(mu,FUN=risky,joint_matrix))))

plot(mu,risk)


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

ggplot(newdata_risk, aes(x=mu, y=risk, alpha=risk))+
  geom_line(size=0.7)



#Product as a function of mu

C=cbind(product,mu)
newdata_product<-as.data.frame(C)

ggplot(newdata_product, aes(x=mu, y=product,alpha=10/product))+
 geom_line(size=1) 



#write.table(mse, file="MSE_mu_UP.csv", sep=",")
#write.table(risk, file="RISK_mu_UP.csv", sep=",")
#write.table(product, file="PRODUCT_mu_UP.csv", sep=",")
#write.table(sample_matrix, file="posterior_mu_UP.csv", sep=",")

#write.table(sample_matrix, file="posterior_sampling.Rdata")
#posterior_sampling<-read.table("posterior_sampling.Rdata") 

#write.table(mse, file="mse.Rdata")
# mse_data<-read.table("mse.Rdata") 

#write.table(risk, file="risk.Rdata")
# risk_data<-read.table("risk.Rdata") 

#write.table(product, file="product.Rdata")
#product_data<-read.table("product.Rdata") 
