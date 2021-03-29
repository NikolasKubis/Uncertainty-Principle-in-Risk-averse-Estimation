setwd(" ")

# This script measures the average skewness of P(x|y).


joint_matrix<-data.matrix(samples_posterior[,1:(ncol(samples_posterior)-1)])

# The conditional mean 
conditional_mean<-function(jm){
  return(colMeans(jm)) 
  
}




#The maximally risk averse

maximally_risk <- function(jm) {
  m1=colMeans(jm)
  m2=colMeans(jm^2)
  m3=colMeans(jm^3)
  return ((m3-m2*m1)/(2*(m2-m1^2)))
}


# max eigenvalue

max_eig<-function(jm) {
  
  m1 <- colMeans(jm)
  m2 <- colMeans(jm^2)
  me<-m2-m1^2
  return(me)
}

dP=sqrt(mean(abs(maximally_risk(joint_matrix)-conditional_mean(joint_matrix))/(max_eig(joint_matrix))))


# Compute the mse of Xo 
mse_conditional_mean<- function(jm)
{
  m1=colMeans(jm)
  conditioned_y=(jm-rep(c(as.matrix(m1)),each=nrow(jm)))^2
  term=colMeans(conditioned_y)
  return(mean(term))
}

msexo=mse_conditional_mean(joint_matrix)

# Compute the sev of Xinf

sev_max_risk_averse <-function(jm,n)
{
  m1=colMeans(jm)
  
  m2=colMeans(jm^2)
  
  m3=colMeans(jm^3)
  
  XX=(m3-m2*m1)/(2*(m2-m1^2))
  
  oros=colMeans((jm-rep(c(as.matrix(XX)),each=nrow(jm)))^2)
  
  MM=(matrix(c(oros), nrow=n, ncol=length(c(oros)), byrow=TRUE))
  
  h=colMeans((((jm-rep(c(as.matrix(XX)),each=nrow(jm)))^2 - MM ) )^2)
  
  r=mean(h) 
  
  return(r)
}

sevxoo=sev_max_risk_averse(joint_matrix)

