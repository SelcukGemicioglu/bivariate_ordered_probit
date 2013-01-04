# This file estimates an ordered probit model by hand using R's optim command

library(MASS)
library(foreign)
library(mvtnorm)


# Generate a dataset
n<-10000;
e<-rnorm(n);
x <- matrix(data = rnorm(n),nrow = n, ncol  =1 )
a<- 0.5
b<-0.2
ystar<-a+b*x+e;
y <- ifelse(ystar<=0,0,ystar)
y <- ifelse(y > 0 & y <=1,1,y)
y <- ifelse(y > 1,2,y)





# The first ordered probit function uses the pnorm function to evaluate the normal cdf
# This is much faster than the second and third models, which use the pmvnorm function
# The pmvnorm function is necessary to extend the ordeded probit to a bivariate ordered probit model


# This is a basic ordeded probit model
llk.oprobit1 <- function(param, x, y) {     
  x <- cbind(1, x)                         # adding a constant
  b <- param[1:ncol(x)]                     # generating a placeholder for the parameters
  t2 <- param[(ncol(x)+1)]                  # this is the second cutpoint (note that the first is the constant in the model)

  # probabilities and penalty function
  xb <- x%*%b                               
  p1 <- log(pnorm(-xb))                     # the first cutpoint is zero (because the constant is in the model)
  if (t2<=0)  p2 <- -(abs(t2)*10000)        # make sure t2>0
  else p2 <- log(pnorm(t2-xb)-pnorm(-xb))   # y=1 
  p3 <- log(1-pnorm(t2-xb))                 # y=2
                                             
phi <- (y==0)*p1 + (y==1)*p2 + (y==2)*p3 
return(-sum(phi))     
}

ls.result <- lm(y~x)                    # generate initial values
stval <- c(ls.result$coefficients,2)  

oprobit.result1 <- optim(stval, llk.oprobit1, method="BFGS", x=x, y=y, hessian=T) # Estimate the model

df <- length(y) - length(oprobit.result1$par)
se <- sqrt(diag(abs(solve(oprobit.result1$hessian))))
t <- oprobit.result1$par/se
p <- (1-pt(abs(t),df))*1.96
display1 <- cbind(oprobit.result1$par,se,t,p)
display1 # Displays the coefficients and standard errors of the parameters



# This also estimates an ordered probit model, but used pmvnorm, which takes longer
llk.oprobit2 <- function(param, x, y) {     
  x <- cbind(1, x)                         
  b <- param[1:ncol(x)]                     
  t2 <- param[(ncol(x)+1)]                  
  xb <- x%*%b                               
sigma <- 1
llik <- 0 
for (i in 1:nrow(xb)){ 	
	if (y[i]==0) 
			llik <- llik + log(pmvnorm(lower = -Inf, upper = 0, mean = xb[i], sigma = sigma)[1]) 
	else if (y[i]==1)
			llik <- llik + log(pmvnorm(lower = 0, upper = Inf, mean = (t2-xb[i]), sigma = sigma)[1]  - pmvnorm(lower = -Inf, upper = 0, mean = xb[i], sigma = sigma)[1])  
	else 
			llik <- llik +  log(1 - pmvnorm(lower = 0, upper = Inf, mean = (t2-xb[i]), sigma = sigma)[1])
}								
return(llik) 
}


ls.result <- lm(y~x)                   
stval <- c(ls.result$coefficients,2) 
oprobit.result2 <- optim(stval, llk.oprobit2, method="L-BFGS-B", x=x, y=y, hessian=T, control = list(fnscale = -1), lower=c(-Inf,-Inf,0), upper=c(Inf,Inf,Inf))

df <- length(y) - length(oprobit.result2$par)
se <- sqrt(diag(abs(solve(oprobit.result2$hessian))))
t <- oprobit.result2$par/se
p <- (1-pt(abs(t),df))*1.96
display2 <- cbind(oprobit.result2$par,se,t,p)
display2


# Check to make sure the results are the same as the canned ordered probit routine
y <- as.factor(y)
summary(polr(y~x, method="probit"))
display1
display2