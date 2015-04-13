rm(list=ls())
library(ismev)
library(evmix)
library(evir)
library(extRemes)
years=1850:2200
n=length(unique(years))
changement=1
t=seq(1,n)
mu=(t>changement)*.005*(t-changement)+100
sigma=(t>changement)*.001*(t-changement)+1
y=mapply(revd,1,mu,sigma,rep(1,n))
ydat=data.frame(years,y,rep(1,n),years,rep(1,n),years,rep(1,n))
names(ydat)=c("year","y","mub","mua","sigb","siga","xi")


t0=with(ydat,which(year==1860))
t1=with(ydat,which(year==2117))

getFAR.theo=function(xp,t0,t1){
	p0=1-pgev(xp,mu=mu[t0],sigma=sigma[t0],xi=1)
	p1=1-pgev(xp,mu=mu[t1],sigma=sigma[t1],xi=1)
	p0/p1
}


gev.ratio.ic.mu=function(xp,t0,t1,y.fit,ydat,ci.p=0.95,like.num = 1000 ,mulink=identity, siglink = identity, show = TRUE,  method="Nelder-Mead",maxit = 10000,to.plot=FALSE, ...){
	require(evir)
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
	print(y.fit)
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
	browser()
	p0=1-evir::pgev(xp,mu=mu.v[t0],sigma=sig.v[t0],xi=sha.v[t0])
	p1=1-evir::pgev(xp,mu=mu.v[t1],sigma=sig.v[t1],xi=sha.v[t1])
	print(p0/p1)
	init.o=mle
	gev.lik.o <- function(par) {
		xi=par[5]
		mub=par[1]
		mua=par[2]
		sigb=par[3]
		siga=par[4]
		sc <- siglink(sigmat %*% (c(sigb,siga)))
		mu.v <- mulink(mumat %*% (c(mub,mua)))
		res=levd(x=xdat,location=mu.v, scale=sc, shape=xi,type="GEV",npy=1)
		browser()
		res
	}
	#         debug(gev.lik.o)
	cat("Normal")
	print(gev.lik.o(init.o))
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
		if(class(tc)=="try-error" | is.na(mua))
			return(10^6)
		mu.v <- mulink(mumat %*% (c(mub,mua)))
		y <- (xdat - mu.v)/sc
		y <- 1 + xi * y
		res=levd(x=xdat,location=mu.v, scale=sc, shape=xi,type="GEV",npy=1)
		if(is.infinite(res))
			return(10^6)
		res
	}
	#         debug(gev.lik)
#	debug(ratio2mua)
	cat("Paramétrisation")
	print(gev.lik(par=init,ratio=p0/p1,xpi=xp))
	ratio.l <- sort(c(seq(0,2,length.out=like.num),p0/p1))[-1]
	#         ratio.l <- sort(p0/p1)
	#         print(ratio.l)
	overallmax <- -gev.lik(init,p0/p1,xp)
	parmax=numeric(length(ratio.l))
	for(i in 1:length(ratio.l)) {
		#                 print(ratio.l[i])
		#                 fit <- optim(par=c(init), gev.lik, hessian = TRUE, method = method,control = list(maxit = maxit, ...),xpi=xp,ratio=ratio.l[i])
		fit <- optim(par=c(init), gev.lik,xpi=xp,ratio=ratio.l[i])
		parmax[i] <- fit$value
		#                 if(-parmax[i]>overallmax){
		#                         print(fit)
		#                         mua=ratio2mua(fit$par[1],fit$par[2],fit$par[3],fit$par[4],ratio.l[i],xp)
		#                         browser()
		#                         print(fevd(y,ydat.t,location.fun=~mua,scale.fun=~mub,initial=list("mu0"=fit$par[1],"mu1"=mua,"sigma0"=fit$par[2],"sigma1"=fit$par[3],"shape"=fit$par[4])))
		#                 }
	}
	parmax <-  - parmax
	overallmax <- -gev.lik(init,p0/p1,xp)
	if(abs(max(parmax)-overallmax)>0.0001){
		print("non equal mle")
		print(max(parmax))
		print(overallmax)
	}
	crit <- overallmax - qchisq(0.999, 1)/2
	cond <- parmax > crit
	#         cond <- parmax > -10^6
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
	out <- c(min(ci), p0/p1, max(ci))
	names(out) <- c("LowerCI", "Estimate", "UpperCI")
	print(out)
	print(getFAR.theo(xp,t0,t1))
	#         browser()
        out
}
# debug(gev.ratio.ic.mu)
gev.ratio.ic.mu(xp=106,t0=t0,t1=t1,y.fit=y.fit,ydat=ydat,to.plot=TRUE,like.num=1000)
