#This file demonstrates how to estimate a probit model using R's 'optim' command
library(MASS)
library(foreign)


####################################################################################
# Below I show three different ways to estimate a probit model.  But first, some
# notes about optim.  

# 1. Optim does not find the maximum of the likelihood function.  Instead, we
# feed optime the negative of the function and it finds the minimum.  This is why
# we return the negative of the sum.

# 2. Because the maximum of a function is also the maximum of its log, we use
# the log of the function sent to the normal density.  This is done by either
# specifying the log.p=TRUE option in the pnorm function (Model 1), or by 
# taking the log by hand (Models 2 and 3).

#Define The First Likelihood
llik.probit1 <- function(par, X, y){
  Y <- as.matrix(y)
  X <- cbind(1,X)
  phi <- pnorm(ifelse(Y == 0, -1, 1) * X%*%par, log.p = TRUE)
  out <- return(-sum(phi))
}
# Note that we have optim minizizing the negative of the likelihood
# function instead of maximizing the likelihood



# This shows what the 'log.p = True' option is doing
llik.probit2 <- function(beta,X,y){
  Y <- as.matrix(y)
  X<-cbind(1,X)
  phi <- log(pnorm(ifelse(Y == 0, -1, 1)*(X %*% beta)))
  out <- return(-sum(phi))
}




# Specifying the value of y in a different manner
llik.probit3 <- function(beta,X,y){
  X<-cbind(1,X)
  normal_den <- pnorm(X %*% beta)
  phi <- log(y * (normal_den) + (1-y) * (1-normal_den))
  out <- return(-sum(phi))
}





####################################################################################
# Test the models with simulated data

# Generate data to test the functions with
n<-10000;
e<-rnorm(n);
x<-rnorm(n);
a<- .5
b<- -.4
y<-a+b*x+e;
y<-ifelse(y<=0,0,1)



# Maximize the likelihoods
# Note the par = c(0,1) option spefifys starting values for the two parameters
out1 <- optim(par = c(0,1), llik.probit1 , y = y, X = x, method="BFGS", hessian = T)
out2 <- optim(par = c(0,1), llik.probit2 , y = y, X = x, method="BFGS", hessian = T)
out3 <- optim(par = c(0,1), llik.probit3 , y = y, X = x, method="BFGS", hessian = T)


# Create a function to calculate standard errors
display <- function(estimate){
se <- sqrt(diag(solve(abs(estimate$hessian))))
Zv<-estimate$par/se
pv<-2*(1-pnorm(abs(Zv)))

result<-cbind(estimate$par,se,Zv,pv)
colnames(result)<-c("Coef","Std.err","Zvalue","P-value")
rownames(result)<-c("Intercept","x")
round(result,4)
}

display(out1)
display(out2)
display(out3)

# check with canned function; almost identical
canned_probit <- glm(y ~ x, family=binomial(link="probit"))
summary(canned_probit)


