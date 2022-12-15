#_____________ MCMC routine for Gaussian AR(1) with stationarity __________________

# pdf of a re-scaled beta 
rebe <- function( x,a=2,b=2 )
{
	p <- -(a+b-1)*log(2) - lbeta(a,b) + (a-1)*log(1 +x) + (b-1)*log(1-x)
	return(p)
}

# full conditional of beta
fc_beta <- function( be,al,la,st12,stt1,st1,a,b )
{
	p <- -la/2*(be^2*st12 - 2*be*(stt1-al*st1)) + (a-1)*log(be+1) + (b-1)*log(1-be)
	return(p)
}

# data statistics
datasts <- function( x )
{
	T <- length(x)
	xt <- x[2:T]
	xt1 <- x[1:(T-1)]
	xt_bar <- mean(xt)
	xt1_bar <- mean(xt1)
    st1 <- sum(xt1)
	st12 <- sum(xt1^2)
	stt1 <- sum(xt*xt1)
	nu <- sum((xt-xt_bar)*(xt1-xt1_bar))
	de <- sum((xt1-xt1_bar)^2)
	bhat <- nu/de; ahat <- xt_bar - bhat*xt1_bar
	shat <- mean((xt -ahat -bhat*xt1)^2)
	cat("MLE of slope=",format(bhat,digits=3),"\n")
  cat("MLE of intercept=",format(ahat,digits=3),"\n")
	return(list(xt=xt,xt1=xt1,xtbar=xt_bar,xt1bar=xt1_bar,st1=st1,st12=st12,stt1=stt1,T=T-1, bhat= bhat, ahat = ahat, lhat = 1/shat))
}


# mh step for beta using an independent proposal with mean at the MLE, use v to achieve acceptance rate
move.rw_beta <- function(cb,al,la,S,a,b,v){
  # if(abs(S$mle)<1){
    aux <- (S$bhat^2+v-1)/(2*v)
    pa <- -(S$bhat+1)*aux
    pb <- (S$bhat-1)*aux
	if(pb*pa<=0){
		cat("variance not allowed\n")
		cat("proposal parameters are negative\n")
		cat(c(pa,pb),"\n")
		stop("reduce proposal variance")
	}
	prop <- 2*rbeta(1,pa,pb)-1
	mhr <- fc_beta(prop,al,la,S$st12,S$stt1,S$st1,a,b) - fc_beta(cb,al,la,S$st12,S$stt1,S$st1,a,b)
	u <- log(runif(1))
	alpha <- min(0,mhr)
	dec <- (alpha>=u)*1
	new <- prop*dec + cb*(1-dec)
	return(list(value=new,acc=dec))
}

move.beta <- function(cb,al,la,S,a,b,p){
  b.c <- p*(1-cb)/(1+cb)
	prop <- 2*rbeta(1,p,b.c)-1
  c.b <- p*(1-prop)/(1+prop)
	mhr <- fc_beta(prop,al,la,S$st12,S$stt1,S$st1,a,b) - fc_beta(cb,al,la,S$st12,S$stt1,S$st1,a,b) + rebe(cb,p,c.b) - rebe(prop,p,b.c)
  # print(mhr)
  # print(prop)
	u <- log(runif(1))
	alpha <- min(0,mhr)
	dec <- (alpha>=u)*1
	new <- prop*dec + cb*(1-dec)
	return(list(value=new,acc=dec))
}

# main routine.
# aprior are the parameters of a Gaussian prior for the steady state
# bprior are the parameters of a re-scaled Beta prior for the slope, default are (3,3)
# lprior are the parameters of a Gamma prior for the error precision
# x0 is the initial point
# v is the variance of the proposal distribution, use this to control the acceptance rate
# output is a list with the chains, DIC, acceptance rates and data
MCMC.RW_AR1stat <- function(M, dats, x0=list(a=0,b=0,l=1), aprior, bprior=c(2,2), lprior, v=0.1 )
{
	cat("MCMC routine with independent proposal\n")
	# initialise output
	Al <- vector('numeric',M)
	Be <- vector('numeric',M)
	La <- vector('numeric',M)
	dic <- vector('numeric',M)
	# initialise chains
	ca <- x0$a
	cb <- x0$b
	if(abs(cb)>1) stop("starting point for beta is not stationary")
	cl <- x0$l
	Al[1] <- ca
	Be[1] <- cb
	La[1] <- cl
	# calculate statistics
	dats <- dats[!is.na(dats)]
	S <- datasts(dats)
	astar <- S$T/2+lprior[1]
	dic[1] <- -(S$T*log(cl) - cl*sum((S$xt- ca - cb*S$xt1)^2))
	acb <- 0
	pm <- prod(aprior)
	for(i in 2:M){
		pstar <- S$T*cl + aprior[2]
		mstar <- (S$T*cl*(S$xtbar - cb*S$xt1bar) + pm )/pstar
		ca <- mstar + rnorm(1)/sqrt(pstar)
		bstar <- lprior[2] + sum((S$xt- ca - cb *S$xt1)^2)/2
    # if(astar)
		cl <- rgamma(1,astar,rate=bstar)
		aux <- move.rw_beta(cb,ca,cl,S,bprior[1],bprior[2],v)
		cb <- aux$value
		acb <- acb + aux$acc
		La[i] <- cl
		Be[i] <- cb
		Al[i] <- ca
		dic[i] <- -(S$T*log(cl) - cl*sum((S$xt- ca - cb*S$xt1)^2))
		if(i%%(M/10)==0) 	cat(i/M*100,"% done \n",sep="")
	}
	DIC <- 2*mean(dic) + S$T*log(mean(La)) - mean(La)*sum((S$xt- mean(Al) - mean(Be)*S$xt1)^2)
	cat("acc rate=",format(acb/M,digits=3),"\n")
	return(list(alpha=Al,beta=Be,lambda=La,acc=acb,data=S,DIC=DIC, bhat=S$bhat, ahat=S$ahat, lhat=S$lhat))
}

# Find a mode
moda <- function(x){
   den <- density(x)
   aux <- which(den$y == max(den$y))
   return(den$x[aux])
}

# one step ahead predictive 
AR1.predict<- function(sal, keep){
	xT <- tail(sal$data$xt, n=1)
	M <- length(sal$beta[keep])
	y <- sal$alpha[keep] + sal$beta[keep]*xT + rnorm(M)/sqrt(sal$lambda[keep])
	return(y)
}