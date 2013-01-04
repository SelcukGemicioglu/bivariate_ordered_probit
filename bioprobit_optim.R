#This program estimates a bivariate ordered probit model with each dependent variable having 
# 3 outcomes using R's optim function

library(MASS)
library(foreign)
library(mvtnorm)


rho.convert <- function(r){			#This function is necessary to convert back the estimate of rho
rho_out <- (exp(r)-1)/(1+exp(r))	
return(rho_out)
}	


#Generating data to test the model
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

y1 <- ifelse(y1<=0,0,y1)
y1 <- ifelse(y1 > 0 & y1 <=1,1,y1)
y1 <- ifelse(y1 > 1,2,y1)

y2 <- ifelse(y2<=0,0,y2)
y2 <- ifelse(y2 > 0 & y2 <=1,1,y2)
y2 <- ifelse(y2 > 1,2,y2)


mydata <- cbind(y1,y2,x1,x2)
mydata <- as.data.frame(mydata)
write.dta(mydata, file="/Users/tiernay/Desktop/probit/bioprobit_r.dta")

y1.start <- as.factor(y1)
y2.start <- as.factor(y2)




# The bivariate ordered probit function
log.lik <- function(par , X1 , X2 , Y1 , Y2) {
X1 <- cbind(1,X1)
X2 <- cbind(1,X2)

Beta1 <- par[1:ncol(X1)]		#parameters for the first equation
Beta2 <- par[(ncol(X1)+1):(ncol(X1)+ncol(X2))]   #parameters for the second equation

cut21 <- par[(ncol(X1)+ncol(X2)+1)]	#Second cut point for the first equation
cut22 <- par[(ncol(X1)+ncol(X2)+2)]	#Second cut point for the second equation

gamma <- par[(ncol(X1)+ncol(X2)+3)]	  #parameter for the correlation coefficient

#multiply gamma by a column of 1's, and then transform: (exp(rho)-1)/(1+exp(rho))
rho <- (exp(gamma) - 1) / (1 + exp( gamma)) 

mu1 <- X1 %*% Beta1 
mu2 <- X2 %*% Beta2



llik <- 0 

for (i in 1:nrow(mu1)){ 
	Sigma <- matrix(c(1, rho, rho, 1), 2, 2) 
  
	if (Y1[i]==0) 	#All the calculations if y1 = 0	
		if (Y2[i]==0)  # If y2 = 0
			llik <- llik + log(pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 0), mean = c(mu1[i,],mu2[i,])  , corr = Sigma))		
		else if (Y2[i]==1)  # If y2 = 1
			llik <- llik + log(pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 0), mean = c(mu1[i,],-(cut22 - mu2[i,]))  , corr = Sigma) - pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 0), mean = c(mu1[i,],mu2[i,]), corr = Sigma) )		
		else 			# If y2 = 2
			llik <- llik + log(pmvnorm(lower=c(-Inf,-Inf), upper=c(0,Inf), mean=c(mu1[i,],0),corr=Sigma) - pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 0), mean = c(mu1[i,],-(cut22 - mu2[i,]))  , corr = Sigma) ) 
		
	else if (Y1[i]==1) #All the calculations if y1 = 1	
		if (Y2[i]==0)  # If y2 = 0
			llik <- llik + log(pmvnorm(lower = c(-Inf,-Inf), upper = c( 0,0), mean = c(-(cut21 - mu1[i,]),mu2[i,])  , corr = Sigma) - pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 0), mean = c(mu1[i,],mu2[i,]), corr = Sigma) )
		else if (Y2[i]==1)  # If y2 = 1
			llik <- llik + log(pmvnorm(lower = c( -Inf,-Inf), upper = c( 0,0), mean = c(-(cut21 - mu1[i,]),-(cut22 - mu2[i,])), corr = Sigma) - pmvnorm(lower = c(-Inf,-Inf), upper = c( 0,0), mean = c(-(cut21 - mu1[i,]),mu2[i,])  , corr = Sigma) - pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 0), mean = c(mu1[i,],-(cut22 - mu2[i,]))  , corr = Sigma) + pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 0), mean = c(mu1[i,],mu2[i,])  , corr = Sigma)    )		
		else 			# If y2 = 2
			llik <- llik + log ( pmvnorm(lower=c(c-Inf,-Inf), upper=c(0, Inf),mean=c(-(cut21-mu1[i,]),0),corr=Sigma) -pmvnorm(lower = c( -Inf,-Inf), upper = c( 0,0), mean = c(-(cut21 - mu1[i,]),-(cut22 - mu2[i,])), corr = Sigma) - pmvnorm(lower=c(-Inf,-Inf),upper=c(0, Inf),mean=c(mu1[i,],0),corr=Sigma) + pmvnorm(lower = c(-Inf,-Inf), upper = c( 0,0), mean = c(mu1[i,],-(cut22-mu2[i,]))  , corr = Sigma)    )		
		

	else 				#All the calculations if y1 = 2
		if (Y2[i]==0)  # If y2 = 0		
			llik <- llik +  log(pmvnorm(lower= c(-Inf,-Inf),upper=c(Inf,0),mean=c(0,mu2[i,]),corr=Sigma) - pmvnorm(lower = c(-Inf,-Inf), upper = c(0,0), mean = c(-(cut21 - mu1[i,]), mu2[i,])  , corr = Sigma) ) 
		else if (Y2[i]==1)  # If y2 = 1			
			llik <- llik + log(pmvnorm(lower= c(-Inf,-Inf),upper=c(Inf,0),mean=c(0,-(cut22-mu2[i,])),corr=Sigma) - pmvnorm(lower = c(-Inf,-Inf), upper = c( 0,0), mean = c(-(cut21 - mu1[i,]),-(cut22 - mu2[i,])), corr = Sigma) - pmvnorm(lower=c(-Inf,-Inf),upper=c(Inf,0),mean=c(0,mu2[i,]),corr=Sigma) +pmvnorm(lower = c(-Inf,-Inf), upper = c( 0,0), mean = c(-(cut21 - mu1[i,]),mu2[i,]),corr=Sigma)   )		
		else 			# If y2 = 2
			llik <- llik + log(1   - pmvnorm(lower=c(-Inf,-Inf),upper=c(Inf,0),mean=c(0,-(cut22-mu2[i,])),corr=Sigma) - pmvnorm(lower=c(-Inf,-Inf),upper=c(0,Inf),mean=c(-(cut21-mu1[i,]),0),corr=Sigma) +  pmvnorm(lower = c(-Inf,-Inf), upper = c( 0,0), mean = c(-(cut21 - mu1[i,]),-(cut22 - mu2[i,])),corr=Sigma)) 
  }
return(llik) 
}
 





# Generate starting values with an ordinary probit
op.result1 <- polr(y1.start ~ x1, method=c("probit"))
op.result2 <- polr(y2.start ~ x2, method=c("probit"))
start.val <- c(-op.result1$zeta[1],op.result1$coef[1],-op.result2$zeta[1],op.result2$coef[1],op.result1$zeta[2],op.result2$zeta[2],0) 

res <- optim(start.val, log.lik, method = "L-BFGS-B", hessian = TRUE, control = list(fnscale = -1), X1 = x1, X2 = x2, Y1 = y1, Y2 = y2, lower = c(-Inf,-Inf,-Inf,-Inf,0,0,-Inf), upper = c(Inf,Inf,Inf,Inf,Inf,Inf,Inf)) 
res$par  # Return the parameters
-solve(res$hessian) # Return the covariance matrix
rho.convert(res$par[7])	 	# Return the correlation coefficint
	



