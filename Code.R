load("WB_GDP-growth.RData")
source("AR1_nonstat.R")
source("AR1_stationary.R")
source("graph.R")
set.seed(200357876)

#_________________________________________________________
# Choose two countries: Germany and Thailand.
C1 <- "DEU"
Tag1 <- codes$Country[codes$Code==C1]
assign(paste("GrData",C1,sep="."),Growth[,C1])

C2 <- "THA"
Tag2 <- codes$Country[codes$Code==C2]
assign(paste("GrData",C2,sep="."),Growth[,C2])

for (cc in c(C1,C2)) {
  print(paste("Mean of growth in",cc))
  print(mean(get(paste("GrData",cc,sep="."))))
  print(paste("The range of growth in",cc))
  print(range(get(paste("GrData",cc,sep="."))))
}

#_________________________________________________________
# Plot both series
plot(get(paste("GrData",C1,sep=".")), 
     ylab="Growth rate (%)", 
     xlab="Year", main="",
     col="navyblue",
     lwd=2, axes=F, 
     ylim=c(-13,13))
lines(get(paste("GrData",C2,sep=".")), 
      lwd=2, col="OrangeRed")
axis(1, pos=c(-13,0), at = seq(1987,2024,by=3))
axis(2, pos=c(1990,0), at=seq(-27,27,by=2))
grid(ny=6,lwd=2)
legend(1999, -4, c(Tag1, Tag2), 
       col=c("navyblue", "OrangeRed"), 
       text.col=c("navyblue", "OrangeRed"), 
       bty="n", lwd=2, lty=c(1,1), 
       cex=1, y.intersp=1.4)

#_________________________________________________________
# Fitting both models, using my elicitation
a.prior <- c(1,0.2)
l.prior <- c(1,1)
b.prior <- c(0,0.6)
ini.point <- list(a=6, b=0.9, l=0.1) # starting point
assign(paste("v",C1,sep="."),0.25) # proposal variances for country 1
assign(paste("v",C2,sep="."),0.20) # proposal variances for country 2
MM <- 1000 # run length

#_________________________________________________________
# The stationary model
for (c in c(C1,C2)){
  assign(paste("MCout_st",c,sep="."),
         MCMC.RW_AR1stat(
           M=MM,
           dats=get(paste("GrData",c,sep=".")), 
           aprior=a.prior, 
           lprior=l.prior, 
           x0=ini.point, 
           v=get(paste("v",c,sep="."))))
}
# The non-stationary model
for (d in c(C1,C2)) {
  assign(paste("MCout_ns",d,sep="."),
         MCMC_AR1non(
           M=MM,
           dats=get(paste("GrData",d,sep=".")), 
           aprior=a.prior, 
           bprior=b.prior, 
           lprior=l.prior, 
           x0=ini.point))
}

#_________________________________________________________
# Trace plots
graph(C1, C2, stat=TRUE, g.prior=FALSE)
graph(C1, C2, stat=FALSE, g.prior=FALSE)

#_________________________________________________________
# Fitting both models, using given elicitation
a.prior.g <- c(2, 0.1)
l.prior.g <- c(2, 0.6)
b.prior.g <- c(0, 0.5)
ini.point.g <- list(a=2, b=0, l=5)
assign(paste("v",C1,"g",sep="."),0.25)
assign(paste("v",C2,"g",sep="."),0.20)
MM.g <- 1.1e5

#_________________________________________________________
# The stationary model
for (a in c(C1,C2)) {
  assign(paste("MCout_st",a,"g",sep="."),
         MCMC.RW_AR1stat(
           M=MM.g, 
           dats=get(paste("GrData",a,sep=".")), 
           aprior=a.prior.g, 
           lprior=l.prior.g, 
           x0=ini.point.g, 
           v=get(paste("v",a,"g",sep="."))))
}
# The non-stationary model
for (b in c(C1,C2)) {
  assign(paste("MCout_ns",b,"g",sep="."),
         MCMC_AR1non(
           M=MM.g, 
           dats=get(paste("GrData",b,sep=".")), 
           aprior=a.prior.g, 
           bprior=b.prior.g, 
           lprior=l.prior.g, 
           x0=ini.point.g))
}

#_________________________________________________________
# Trace plots
graph(C1, C2, stat=TRUE, g.prior=TRUE)
graph(C1, C2, stat=FALSE, g.prior=TRUE)

#_________________________________________________________
# Posterior odds (st C1, st C2, ns C1, ns C2)
burn <- 10000 
# thin <- 3
keep <- seq(burn+1, MM.g)
allcases <- c(paste("MCout_st",C1,"g",sep="."),
              paste("MCout_st",C2,"g",sep="."),
              paste("MCout_ns",C1,"g",sep="."),
              paste("MCout_ns",C2,"g",sep="."))
for (data in allcases) {
  o <- sum(abs(get(data)$beta[keep]) < 0.02)/
    sum(abs(get(data)$beta[keep]) >= 0.02)
  print(o)
}

#_________________________________________________________
# Growth rates
for (ddata in allcases[c(1,2)]) {
  rho = get(ddata)$ahat/(1-get(ddata)$bhat)
  print(rho)
}
for (dddata in allcases[c(3,4)]) {
  rho = get(dddata)$data$ahat/(1-get(dddata)$data$bhat)
  print(rho)
}
