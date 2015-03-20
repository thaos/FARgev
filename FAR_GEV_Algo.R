packages <- c("boot", "quantreg", "evmix", "ismev", "parallel","foreach","doMC")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
	  install.packages(setdiff(packages, rownames(installed.packages())),repos="http://cran.us.r-project.org")
}
library(boot)
library(quantreg)
library(evmix)
library(ismev)
library(parallel)
library(foreach)
library(doMC)

getP <- function(p2,y.fit,ydat,to.plot=FALSE){
	xcord=p2[1]
	ycord=p2[2]
	x=ydat$year
	covariate=ydat$mua
	if (xcord < min(x) | xcord > max(x))
		stop(" Point outside of data range ")
	xcloser=which((abs(x-xcord))==min(abs(x-xcord)))
	xcloser=min(unique(x[xcloser]))
	ccloser=unique(covariate[which(x==xcloser)])
	mle <- y.fit$results$par
	mu=mle[1]+ccloser*mle[2]
	sigma=mle[3]+ccloser*mle[4]
	shape=mle[5]
	1-evir::pgev(ycord,mu=mu,sigma=sigma,xi=shape)
	1-pevd(ycord,loc=mu,scale=sigma,shape=shape)
}
# getP(c(1920,110),y.fit,ydat,to.plot=TRUE)

getFAR <- function(p1,x2,y.fit,ydat){
	prob1=getP(p1,y.fit,ydat)
	p2=c(x2,p1[2])
	prob2=getP(p2,y.fit,ydat)
	if(prob1==0 & prob2==0)
		FAR=1
	else
	FAR=1-(prob2/prob1)
	FAR
}
# getFAR(c(1920,110),2102,y.fit,ydat)

getFAR.theo=function(xp,t0,t1){
	p0=1-pgev(xp,mu=mu[t0],sigma=sigma[t0],xi=1)
	p1=1-pgev(xp,mu=mu[t1],sigma=sigma[t1],xi=1)
	1-p0/p1
}
FARBoot <- function(ydat,indice,p1,x2){
	data.b=ydat[indice,]
	y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga)
	getFAR(p1,x2,y.fit,ydat)
}
