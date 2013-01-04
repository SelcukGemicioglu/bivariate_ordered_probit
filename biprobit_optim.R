# Below are two different ways to estimate a bivariate probit model

library(MASS)
library(foreign)
library(mvtnorm)


rho.convert <- function(r){
rho_out <- (exp(r)-1)/(1+exp(r))	
return(rho_out)
}	


rho.convert2 <- function(r){			#This function is necessary to convert back the estimate of rho
rho_out <- .5*log((1+r)/(1-r))	
return(rho_out)
}	

n <- 1000
x1 <- rnorm(n)
x2 <- rnorm(n)
sig <- diag(2)
sig[1,2] <- sig[2,1] <- .5
mu <- c(0,0)



e <- rmvnorm(n,mu,sig)
a <- .5
b <- .7
c <- -.2
d <- -.9

y1 <- a+b*x1+e[,1]
y2 <- c+d*x2+e[,2]

y1<-ifelse(y1<=0,0,1)
y2<-ifelse(y2<=0,0,1)


y <- cbind(y1,y2)
x <- cbind(x1,x2)

mydata <- cbind(y1,y2,x1,x2)
mydata <- as.data.frame(mydata)
#write.dta(mydata, file="/Users/michaeltiernay/Desktop/biprobit_r.dta")

start.val <- c(0,0,0,0,1)






log.lik <- function(par , X1 , X2 , Y1 , Y2) { 

X1 <- cbind(1,X1)
X2 <- cbind(1,X2)

Beta1 <- par[1:ncol(X1)]		#parameters for the first equation
Beta2 <- par[(ncol(X1)+1):(ncol(X1)+ncol(X2))]	#parameters for the second equation
gamma <- par[(ncol(X1)+ncol(X2)+1)]			#parameter for the correlation coefficient

#multiply gamma by a column of 1's, and then transform: (exp(rho)-1)/(1+exp(rho))
rho <- (exp((matrix(1,nrow(X1),1)) %*% gamma) - 1) / (1 + exp((matrix(1,nrow(X1),1)) %*% gamma)) 
mu1 <- X1 %*% Beta1 
mu2 <- X2 %*% Beta2  



llik <- 0 
for (i in 1:nrow(mu1)){ 
	Sigma <- matrix(c(1, rho[i,], rho[i,], 1), 2, 2) 
	if (Y1[i]==1) 
		if (Y2[i]==1) 
			llik <- llik + log(pmvnorm(lower = c(0, 0), upper = c(Inf, Inf), mean = c(mu1[i,],mu2[i,]), corr = Sigma)) 
		else 
			llik <- llik + log(pmvnorm(lower = c(0, -Inf), upper = c(Inf, 0), mean = c(mu1[i,],mu2[i,]), corr = Sigma)) 
	else 
		if (Y2[i]==1) 
			llik <- llik + log(pmvnorm(lower = c(-Inf, 0), upper = c(0, Inf), mean = c(mu1[i,],mu2[i,]), corr = Sigma)) 
		else 
			llik <- llik + log(pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 0), mean = c(mu1[i,],mu2[i,]), corr = Sigma)) }

return(llik) 

}
 

res1 <- optim(start.val, log.lik, method = "L-BFGS-B", hessian = TRUE, control = list(fnscale = -1), X1 = x1, X2 = x2, Y1 = y1, Y2 = y2) 

res1$par
-solve(res1$hessian)
rho.convert(res1$par[5])	 	




# This function also estimates a bivariate probit model
log.lik2 <- function(par , X1 , X2 , Y1 , Y2) { 

X1 <- cbind(1,X1)
X2 <- cbind(1,X2)

Beta1 <- par[1:ncol(X1)]		#parameters for the first equation
Beta2 <- par[(ncol(X1)+1):(ncol(X1)+ncol(X2))]	#parameters for the second equation
gamma <- par[(ncol(X1)+ncol(X2)+1)]			#parameter for the correlation coefficient

#multiply gamma by a column of 1's, and then transform: (exp(rho)-1)/(1+exp(rho))
#rho <- (exp(gamma) - 1) / (1 + exp(gamma)) 
#rho <- (exp((matrix(1,nrow(X1),1)) %*% gamma) - 1) / (1 + exp((matrix(1,nrow(X1),1)) %*% gamma)) 
rho <- .5*(log((1+gamma)/(1-gamma)))
mu1 <- X1 %*% Beta1 
mu2 <- X2 %*% Beta2  



llik <- 0 
for (i in 1:nrow(mu1)){ 
	Sigma <- matrix(c(1, rho, rho, 1), 2, 2) 
	if (Y1[i]==1) 
		if (Y2[i]==1) 
			llik <- llik + log(pmvnorm(lower = c(-mu1[i,],-mu2[i,]), upper = c(Inf, Inf), mean = c(0,0), corr = Sigma)) 
		else 
			llik <- llik + log(pmvnorm(lower = c(-mu1[i,], -Inf), upper = c(Inf, -mu2[i,]), mean = c(0,0), corr = Sigma)) 
	else 
		if (Y2[i]==1) 
			llik <- llik + log(pmvnorm(lower = c(-Inf, -mu2[i,]), upper = c(-mu1[i,], Inf), mean = c(0,0), corr = Sigma)) 
		else 
			llik <- llik + log(pmvnorm(lower = c(-Inf, -Inf), upper = c(-mu1[i,],-mu2[i,]), mean = c(0,0), corr = Sigma)) }

return(llik) 

}
 


res2 <- optim(start.val, log.lik2, method = "L-BFGS-B", hessian = TRUE, control = list(fnscale = -1), X1 = x1, X2 = x2, Y1 = y1, Y2 = y2, lower = c(-Inf,-Inf,-Inf,-Inf,-1), upper = c(Inf,Inf,Inf,Inf,1)) 

res2$par
-solve(res2$hessian)
#rho.convert(res$par[5])	 	
rho.convert2(res2$par[5])	 	

