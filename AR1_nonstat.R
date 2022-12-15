#_____________ MCMC routine for Gaussian AR(1)  __________________

# Data statistics
xsts <- function( x )
{
	T <- length(x)
	xt <- x[2:T]
	xt1 <- x[1:(T-1)]
	xt_bar <- mean(xt)
	xt1_bar <- mean(xt1)
	st12 <- sum(xt1^2)
	stt1 <- sum(xt*xt1)
	nu <- sum((xt-xt_bar)*(xt1-xt1_bar))
	de <- sum((xt1-xt1_bar)^2)
	bhat <- nu/de; ahat <- xt_bar - bhat*xt1_bar
	shat <- mean((xt -ahat -bhat*xt1)^2)
	cat("MLE of slope=",format(bhat,digits=3),"\n")
  cat("MLE of intercept=",format(ahat,digits=3),"\n")
	return(list(xt=xt,xt1=xt1,xtbar=xt_bar,xt1bar=xt1_bar,st12=st12,stt1=stt1,T=T-1, bhat= bhat, ahat = ahat, lhat = 1/shat))
}

# main routine.  Gibbs sampler. 
# aprior are the parameters of a Gaussian prior for the steady state
# bprior are the parameters of a Gaussian prior for the slope
# lprior are the parameters of a Gamma prior for the error precision
# x0 is the initial point
# output is a list with the chains, DIC and data
MCMC_AR1non <- function( M,dats,x0=list(a=0,b=0,l=1),aprior,bprior,lprior )
{
	cat("Gibbs routine AR(1) without restrictions\n")
	# initialise output
	Al <- vector('numeric',M)
	Be <- vector('numeric',M)
	La <- vector('numeric',M)
	dic <- vector('numeric',M)
	# initialise chains
	ca <- x0$a
	cb <- x0$b
	cl <- x0$l
	Al[1] <- ca
	Be[1] <- cb
	La[1] <- cl
	dats <- dats[!is.na(dats)]
	# calculate statistics
	S <- xsts(dats)
	astar <- S$T/2+lprior[1]
	dic[1] <- -(S$T*log(cl) - cl*sum((S$xt- ca - cb*S$xt1)^2))
	for(i in 2:M){
		pstar <- S$T*cl + aprior[2]
		mstar <- (S$T*cl*(S$xtbar - cb*S$xt1bar) + prod(aprior) )/pstar
		ca <- mstar + rnorm(1)/sqrt(pstar)
		bstar <- lprior[2] + sum((S$xt- ca - cb *S$xt1)^2)/2
		cl <- rgamma(1,astar,rate=bstar)
		pstar <- cl*S$st12 + bprior[2]
		mstar <- (cl*(S$stt1 - ca*S$T*S$xtbar) + prod(bprior))/pstar
		cb <- mstar + rnorm(1)/sqrt(pstar)
		La[i] <- cl
		Be[i] <- cb
		Al[i] <- ca
		dic[i] <- -(S$T*log(cl) - cl*sum((S$xt- ca - cb*S$xt1)^2))
		if(i%%(M/10)==0) cat(i/M*100,"% done \n") 
	}
	DIC <- 2*mean(dic) + S$T*log(mean(La)) - mean(La)*sum((S$xt- mean(Al) - mean(Be)*S$xt1)^2)
	return(list(alpha=Al,beta=Be,lambda=La,data=S,DIC=DIC))
}