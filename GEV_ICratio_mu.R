# graphics.off()
# rm(list=ls())

packages <- c("boot", "quantreg", "evmix", "ismev", "gWidgetstcltk", "gWidgets","evir","extRemes")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
	  install.packages(setdiff(packages, rownames(installed.packages())),repos="http://cran.us.r-project.org")
}
library(boot)
library(quantreg)
library(evmix)
library(ismev)
library(evir)
library(extRemes)
require(gWidgetstcltk) #use the gWidgetstcltk package from CRAN
require(gWidgets) #use the gWidgetstcltk package from CRAN
source("FAR_GEV_Algo.R")


gev.ratio.ic.mu=function(xp,t0,t1,y.fit,ydat,ci.p=0.95,like.num = 1000 ,mulink=identity, siglink = identity, show = TRUE,  method="Nelder-Mead",maxit = 10000,to.plot=FALSE, ...){
	require(evir)
	require(foreach)
	require(doMC)
	if (t0 >= t1)
		stop("t0 must be inferior to t1")
	covariate=ydat$mua
	y.fit=fevd(y,ydat,location.fun=~mua,scale.fun=~siga)
	print(y.fit)
	translation=-covariate[t0]
	t.delta=covariate[t1]+translation
	ydat.t=ydat
	ydat.t$mua=ydat$mua+translation
	#         y.fit=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,optim.arg=list(method=method))
	y.fit=fevd(y,ydat.t,location.fun=~mua,scale.fun=~siga,type="GEV")
	#         print(y.fit)
	xdat=ydat$y
	mle <- y.fit$results$par
	init=mle[-2]
	mumat=as.matrix(ydat.t[,c("mub","mua")])
	sigmat=as.matrix(ydat.t[,c("sigb","siga")])
	shmat=as.matrix(ydat.t[,"xi"])
	mu.v <- mumat %*% mle[1:2]
	sig.v <- sigmat %*% mle[3:4]
	sha.v <- shmat * mle[5]
	print(levd(x=xdat,location=mu.v, scale=sig.v, shape=sha.v,type="GEV",npy=1))
	p0=1-evir::pgev(xp,mu=mu.v[t0],sigma=sig.v[t0],xi=sha.v[t0])
	p1=1-evir::pgev(xp,mu=mu.v[t1],sigma=sig.v[t1],xi=sha.v[t1])
	p0=pevd(xp,loc=mu.v[t0],scale=sig.v[t0],shape=sha.v[t0],lower.tail=FALSE)
	p1=pevd(xp,loc=mu.v[t1],scale=sig.v[t1],shape=sha.v[t1],lower.tail=FALSE)
	print(p0/p1)
	if(p0==0 & p1==0)
		p0p1=1
	else
		p0p1=p0/p1
	init.o=mle
	ratio2mua=function(mub,sigb,siga,xi,ratio,xp){
		sc <- siglink(sigmat %*% (c(sigb,siga)))
		sig0=sc[t0]
		sig1=sc[t1]
		T0=1-exp(-(1-(xi/sig0)*(mub-xp))^(-1/xi))
		T0=(1/ratio)*T0
		T0=(-log(1-T0))^(-xi)
		T0=1-T0
		T0=T0-(xi/sig1)*(mub-xp)
		mua=T0*sig1/(xi*t.delta)
		mua
	}
	gev.lik <- function(par,ratio,xpi) {
		xi=par[4]
		mub=par[1]
		sigb=par[2]
		siga=par[3]
		sc <- siglink(sigmat %*% (c(sigb,siga)))
		tc=try({
		mua=ratio2mua(mub,sigb,siga,xi,ratio,xpi)
		})
		if(class(tc)=="try-error" | is.na(mua)|is.infinite(mua))
			return(10^6)
		mu.v <- mulink(mumat %*% (c(mub,mua)))
		res=levd(x=xdat,location=mu.v, scale=sc, shape=xi,type="GEV",npy=1)
		if(is.na(res) | is.infinite(res))
			return(10^6)
		res
	}
	ratio.l <- sort(c(seq(0,2,length.out=like.num),p0p1))[-1]
	overallmax <- -gev.lik(init,p0p1,xp)
	#         parmax=numeric(length(ratio.l))
	registerDoMC(5)
	parmax=foreach(ratio=ratio.l)%dopar%{
		#         for(i in 1:length(ratio.l)) {
		fit <- optim(par=c(init), gev.lik,xpi=xp,ratio=ratio)
		#                 parmax[i] <- fit$value
		fit$value
	}
	parmax <-  - simplify2array(parmax)
	overallmax <- -gev.lik(init,p0p1,xp)
	if(abs(max(parmax)-overallmax)>0.0001){
		print("non equal mle")
		print(max(parmax))
		print(overallmax)
	}
	crit <- overallmax - qchisq(0.999, 1)/2
	cond <- parmax > crit
	ratio.l <- ratio.l[cond]
	parmax <- parmax[cond]
	aalpha <- qchisq(ci.p, 1)
	cond <- !is.na(ratio.l) & !is.na(parmax)
	smth <- spline(ratio.l[cond], parmax[cond], n = 300)
	if(to.plot){
		plot(ratio.l, parmax, type = "l", xlab = "", ylab = "")
		abline(h = overallmax - aalpha/2, lty = 2, col = 2)
		lines(smth$x,smth$y, lty = 2, col = 2)
	}
	ci <- smth$x[smth$y > overallmax - aalpha/2]
	out <- c(min(ci), p0p1, max(ci),p0,mu.v[t0],sig.v[t0],sha.v[t0],p1,mu.v[t1],sig.v[t1],sha.v[t1])
	names(out) <- c("LowerCI", "Estimate", "UpperCI","p0","mu0","sc0","sh0","p1","mu1","sc1","sh1")
	print(out)
        out
}

