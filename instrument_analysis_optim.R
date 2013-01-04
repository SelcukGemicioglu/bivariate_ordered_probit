# This program estimates a censored bivariate ordered probit model where the dependent
# variable may take on more than one value

# The model is designed to implement the test called for in "Peace Agreements and Peacekeepers: Accounting for Censoring"
# by Michael Tiernay.  See tiernay_instrument_02_13_12.pdf for the results.



rm(list = ls(all = TRUE))

#1st section is qp.result3
#2nd section i try to clean up the likelihood function

library(MASS)
library(foreign)
library(mvtnorm)


rho.convert <- function(r){			#This function is necessary to convert back the estimate of rho
rho_out <- (exp(r)-1)/(1+exp(r))	
return(rho_out)
}	

#instrument <- read.dta("/Volumes/mrt265/backup/instrument/data/instrument_analysis.dta")  #Read in the data
attach(instrument)





#########################
#       #       #       #
#  PK   # PA/PK # PA/PK #
#       #       #       #
######################### Vertical Axis:
#       #       #       # UN's willingness
#   ~   # PA/PK # PA/PK # to get involved
#       #       #       #
#########################
#       #       #       # 
#   ~   #  ~    #  PA   #
#       #       #       #
#########################

#Horizontal Axis: Combatant's willingness to cooperate 


# Note: This same table will appear in the code below



#The likelihood function

log.lik <- function(par , X1 , X2 , Y) { # The function consists of 4 agruments, the parameters to be estimated, 
												# the X's for the UN (X1), the combatants (X2), and the outcome variable


#This section sets up the X's and the parameters that are going to be estimated
X1 <- cbind(1,X1) #adding constants to the colummns of X's
X2 <- cbind(1,X2)

Beta1 <- par[1:ncol(X1)]		#parameters for the UN are Beta1
Beta2 <- par[(ncol(X1)+1):(ncol(X1)+ncol(X2))]   #parameters for the combatants are Beta2

cut21 <- par[(ncol(X1)+ncol(X2)+1)]	#Second cut point for the first equation
cut22 <- par[(ncol(X1)+ncol(X2)+2)]	#Second cut point for the second equation

gamma <- par[(ncol(X1)+ncol(X2)+3)]	  #parameter for the correlation coefficient

#multiply gamma by a column of 1's, and then transform: (exp(rho)-1)/(1+exp(rho))
rho <- (exp((matrix(1,nrow(X1),1)) %*% gamma) - 1) / (1 + exp((matrix(1,nrow(X1),1)) %*% gamma)) 

mu1 <- X1 %*% Beta1  #Generating X'B 
mu2 <- X2 %*% Beta2



# Now we get to the actual function.  We start by giving it a value of zero
# Then we loop over every observation and add it's contibution to the likelihood function	
	
llik <- 0	

for (i in 1:nrow(mu1)){		# Looping over the entire dataset
	
	Sigma <- matrix(c(1, rho[i,], rho[i,], 1), 2, 2)  #Sigma is the correlation matrix
	 
	 
# When Y=1, the combatants signed an agremeent without the UN sending peacekeepers.
# The contribution to the function is the last 2 lines of code.  In order to ensure that 
# pmvnorm is able to return positive probabilities given the parameters, loss functions
# are specified such that if the second cutpoints are less than or equal to zero, then the 
# function is way off.  This ensures we get positive cutpoints.  Finally, to ensure that 
# taking the log of the function does not produce NaN's, a loss function is specified
# such that if taking the log of the function would produce a NaN, then the function is way off.	 
	if (Y[i]==1) 	#Peace Agreement only
			if (cut21<=0)  llik <- -(abs(cut21)*10000) # Loss function for the first cutpoint
			else
			if (cut22<=0)  llik <- -(abs(cut22)*10000) # Loss function for the first cutpoint
			else
			if (pmvnorm(lower=c(-Inf,-Inf), upper=c(-mu1[i,],Inf),corr=Sigma) - #Loss function to prevent NaN's
				pmvnorm(lower=c(-Inf,-Inf), upper=c(-mu1[i,],cut22-mu2[i,]),corr=Sigma) <= 0) llik <- -10000
			else
			llik <- llik + log(pmvnorm(lower=c(-Inf,-Inf), upper=c(-mu1[i,],Inf),corr=Sigma) - 
			pmvnorm(lower=c(-Inf,-Inf), upper=c(-mu1[i,],cut22-mu2[i,]),corr=Sigma))
# The 2 lines above enter the contribution to the likelihood function if Y=1.  
# Using the table below as a guide, these last 2 lines inegrate over the bottom-right cell 
# of the bivarite normal (the cell labeled 'PA').  The first part:
# pmvnorm(lower=c(-Inf,-Inf), upper=c(-mu1[i,],Inf),corr=Sigma)
# integrates over the entire bottom row of the bivariate normal, whil the second part:
# pmvnorm(lower=c(-Inf,-Inf), upper=c(-mu1[i,],cut22-mu2[i,]),corr=Sigma)
# integrates over the two bottom cells labaled '~'.  Subtracting the first part from the second
# produces the area in the 'PA' cell.

# Note that '-mu1[i,]' refers to the UN's first cutpoint, and 'cut21-mu1[i,]' refers to its second
# cutpoint.  The same is true for '-mu2[i,]' and 'cut22-mu2[i,]' for the combatants.


#########################
#       #       #       #
#  PK   # PA/PK # PA/PK #
#       #       #       #
######################### Vertical Axis:
#       #       #       # UN's willingness
#   ~   # PA/PK # PA/PK # to get involved
#       #       #       #
#########################
#       #       #       # 
#   ~   #  ~    #  PA   #
#       #       #       #
#########################

#Horizontal Axis: Combatant's willingness to cooperate 





# When Y=2 there was a peace enforcement mission.  The loss functions are the same as above.

	else if (Y[i]==2) 	#Peace Enforcement Mission
			if (cut21<=0)  llik <- -(abs(cut21)*10000)
			else
			if (cut22<=0)  llik <- -(abs(cut22)*10000)
			else
			if (pmvnorm(lower=c(-Inf,-Inf),upper=c(Inf,-mu2[i,]),corr=Sigma) -
						pmvnorm(lower=c(-Inf,-Inf),upper=c(cut21-mu1[i,],-mu2[i,]),corr=Sigma) <= 0) llik <- -10000
			else
			llik <- llik + log(pmvnorm(lower=c(-Inf,-Inf),upper=c(Inf,-mu2[i,]),corr=Sigma) -
						pmvnorm(lower=c(-Inf,-Inf),upper=c(cut21-mu1[i,],-mu2[i,]),corr=Sigma))

# The 2 lines above enter the contribution to the likelihood function if Y=2  
# Using the table below as a guide, these last 2 lines inegrate over the top-left cell 
# of the bivarite normal (the cell labeled 'PK').  The first part:
# pmvnorm(lower=c(-Inf,-Inf),upper=c(Inf,-mu2[i,]),corr=Sigma)
# integrates over the entire left column of the bivariate normal, whil the second part:
# pmvnorm(lower=c(-Inf,-Inf),upper=c(cut21-mu1[i,],-mu2[i,]),corr=Sigma)
# integrates over the two left cells labaled '~'.  Subtracting the first part from the second
# produces the area in the 'PK' cell.






	else if (Y[i]==3) 	#Peace Agreement with a peacekeeping mission
			if (cut21<=0)  llik <- -(abs(cut21)*10000)
			else
			if (cut22<=0)  llik <- -(abs(cut22)*10000)
			else
			if (1 - pmvnorm(lower=c(-Inf,-Inf),upper=c(Inf,-mu2[i,]),corr=Sigma)
						- pmvnorm(lower=c(-Inf,-Inf),upper=c(-mu1[i,],Inf),corr=Sigma)
						+ pmvnorm(lower=c(-Inf,-Inf),upper=c(-mu1[i,],-mu2[i,]),corr=Sigma) <= 0) llik <- -10000
			else
			llik <- llik + log(1 - pmvnorm(lower=c(-Inf,-Inf),upper=c(Inf,-mu2[i,]),corr=Sigma)
						- pmvnorm(lower=c(-Inf,-Inf),upper=c(-mu1[i,],Inf),corr=Sigma)
						+ pmvnorm(lower=c(-Inf,-Inf),upper=c(-mu1[i,],-mu2[i,]),corr=Sigma))

			#llik <- llik + log(pmvnorm(lower=c(-mu1[i,],-mu2[i,]),upper=c(Inf, Inf),corr=Sigma))


# The 3 lines above enter the contribution to the likelihood function if Y=3  
# Using the table below as a guide, these 3 lines inegrate over the 4 cells 
# of the bivarite normallabeled 'PA/PK'.  The first part:
# pmvnorm(lower=c(-Inf,-Inf),upper=c(Inf,-mu2[i,]),corr=Sigma)
# integrates over the entire left column of the bivariate normal, whil the second part:
# pmvnorm(lower=c(-Inf,-Inf),upper=c(cut21-mu1[i,],-mu2[i,]),corr=Sigma)
# integrates over the two left cells labaled '~'.  Subtracting the first part from the second
# produces the area in the 'PK' cell.	




#########################
#       #       #       #
#  PK   # PA/PK # PA/PK #
#       #       #       #
######################### Vertical Axis:
#       #       #       # UN's willingness
#   ~   # PA/PK # PA/PK # to get involved
#       #       #       #
#########################
#       #       #       # 
#   ~   #  ~    #  PA   #
#       #       #       #
#########################

#Horizontal Axis: Combatant's willingness to cooperate 


	else if (Y[i]==4) #	No Peace Agreement and no peacekeeping mission
			if (cut21<=0)  llik <- -(abs(cut21)*10000)
			else
			if (cut22<=0)  llik <- -(abs(cut22)*10000)
			else
			if (pmvnorm(lower = c(-Inf, -Inf), upper =c(cut21-mu1[i,],-mu2[i,]),corr = Sigma) 
			+	pmvnorm(lower = c(-Inf,-Inf), upper =c(-mu1[i,],cut22 - mu2[i,]),corr = Sigma)
			-	pmvnorm(lower = c(-Inf,-Inf), upper =c(-mu1[i,], -mu2[i,]),corr = Sigma) <= 0) llik <- -10000
			else
			llik <- llik + log(pmvnorm(lower = c(-Inf, -Inf), upper =c(cut21-mu1[i,],-mu2[i,]),corr = Sigma) 
			+	pmvnorm(lower = c(-Inf,-Inf), upper =c(-mu1[i,],cut22 - mu2[i,]),corr = Sigma)
			-	pmvnorm(lower = c(-Inf,-Inf), upper =c(-mu1[i,], -mu2[i,]),corr = Sigma))
	
}
	
return(llik)	
}	


#############################

# X1 refers to the UN's resolve whereas X2 refers to the combatants level of cooperation

# Or use sqrt_black , square_black

X1 <- cbind(black_hawk_down,duration,d2,d3 ,landarea ,lpop ,gdp_per_cap, Oil ,ef ,lmtnest ,cumulative_deaths)
X2 <- cbind(duration,d2,d3 ,landarea ,lpop ,gdp_per_cap, Oil ,ef ,lmtnest ,cumulative_deaths)
Y <- instrument_outcome


#X1 <- cbind(lpop ,gdp_per_cap, Oil ,ef ,lmtnest)
#X2 <- cbind(lpop ,gdp_per_cap,lmtnest, eeurop, lamerica, asia,ef)
#Y <- instrument_outcome


#names <- rbind("constant","population","polity","gdp","moutains","elf","army size","deaths","duration","europe","asia","americas","constant","population","polity","gdp","moutains","elf","army size","deaths","duration","cut1","cut2","rho")



start.val <- c(matrix(0,1,(ncol(X1)+ncol(X2)+2)),1,1,0)

res <- optim(start.val, log.lik, method = "L-BFGS-B", hessian = TRUE, control = list(fnscale = -1, maxit=2500,trace=1), X1 = X1, X2 = X2, Y = Y)


#res <- optim(start.val, log.lik, method = "L-BFGS-B", hessian = TRUE, control = list(fnscale = -1, maxit=2500,trace=1), X1 = X1, X2 = X2, Y = Y,  lower = c(matrix(-Inf,1,(ncol(X1)+ncol(X2)+2)),0,0,-Inf), upper = c(matrix(Inf,1,(ncol(X1)+ncol(X2)+2)),Inf,Inf,Inf)) 

df <- length(Y) - length(res$par)
se <- sqrt(diag(abs(solve(res$hessian))))
t <- res$par/se
p <- (1-pt(abs(t),df))*1.96
display <- cbind(names,res$par,se,t,p)
display
rho <- rho.convert(res$par[(length(res$par))])
rho






























