library(ismev)
library(evmix)
library(evir)
library(extRemes)
years=1850:2200
n=length(unique(years))
changement=1
t=seq(1,n)
mu=(t>changement)*.001*(t-changement)+100
sigma=(t>changement)*.001*(t-changement)+1
y=mapply(revd,1,mu,sigma,rep(1,n))
ydat=data.frame(years,y,rep(1,n),years,rep(1,n),years,rep(1,n))
names(ydat)=c("year","y","mub","mua","sigb","siga","xi")


t0=with(ydat,which(year==1860))
t1=with(ydat,which(year==1904))
gev.ratio.ic(xp=108,t0=t0,t1=t1,y.fit=y.fit,ydat=ydat,to.plot=TRUE)

getFAR.theo=function(xp,t0,t1){
	p0=1-pnorm(xp,mean=mu[t0],sd=sigma[t0])
	p1=1-pnorm(xp,mean=mu[t1],sd=sigma[t1])
	p0/p1
}


gev.ratio.ic.mu=function(xp,t0,t1,y.fit,ydat,ci.p=0.95,like.num = 1000 ,mulink=identity, siglink = identity, show = TRUE,  method="Nelder-Mead",maxit = 10000,to.plot=FALSE, ...){
	require(evir)
	if (t0 >= t1)
		stop("t0 must be inferior to t1")
	covariate=ydat$mua
	translation=-covariate[t0]
	t.delta=covariate[t1]+translation
	ydat.t=ydat
	ydat.t$mua=ydat$mua+translation
	y.fit=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,optim.arg=list(method=method))
	xdat=ydat$y
	mle <- y.fit$results$par
	init=mle[-2]
	mumat=as.matrix(ydat.t[,c("mub","mua")])
	sigmat=as.matrix(ydat.t[,c("sigb","siga")])
	shmat=as.matrix(ydat.t[,"xi"])
	mu <- mumat %*% mle[1:2]
	sig <- sigmat %*% mle[3:4]
	sha <- shmat * mle[5]
	p0=1-evir::pgev(xp,mu=mu[t0],sigma=sig[t0],xi=sha[t0])
	p1=1-evir::pgev(xp,mu=mu[t1],sigma=sig[t1],xi=sha[t1])
	print(p0/p1)
	init.o=mle
	gev.lik.o <- function(par) {
		xi=par[5]
		mub=par[1]
		mua=par[2]
		sigb=par[3]
		siga=par[4]
		sc <- siglink(sigmat %*% (c(sigb,siga)))
		mu <- mulink(mumat %*% (c(mub,mua)))
		res=levd(x=xdat,location=mu, scale=sc, shape=xi,type="GEV",npy=1)
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
	gev.lik <- function(par,mua) {
		xi=par[4]
		mub=par[1]
		sigb=par[2]
		siga=par[3]
		sc <- siglink(sigmat %*% (c(sigb,siga)))
		mu <- mulink(mumat %*% (c(mub,mua)))
		res=levd(x=xdat,location=mu, scale=sc, shape=xi,type="GEV",npy=1)
		if(is.infinite(res))
			return(10^6)
		res
	}
	#         debug(gev.lik)
#	debug(ratio2mua)
	cat("Paramétrisation")
	print(gev.lik(par=init,mua=mle[2]))
	ratio.l <- sort(c(seq(-0.03,0.03,length.out=like.num),mle[2]))
	#         ratio.l <- sort(p0/p1)
	#         print(ratio.l)
	overallmax <- -gev.lik(init,mle[2])
	parmax=numeric(length(ratio.l))
	for(i in 1:length(ratio.l)) {
		#                 fit <- optim(par=c(init), gev.lik, hessian = TRUE, method = method,control = list(maxit = maxit, ...),xpi=xp,ratio=ratio.l[i])
		fit <- optim(par=c(init), gev.lik,mua=ratio.l[i])
		parmax[i] <- fit$value
		#                 if(-parmax[i]>overallmax){
		#                         print(fit)
		#                         mua=ratio2mua(fit$par[1],fit$par[2],fit$par[3],fit$par[4],ratio.l[i],xp)
		#                         browser()
		#                         print(fevd(y,ydat.t,location.fun=~mua,scale.fun=~mub,initial=list("mu0"=fit$par[1],"mu1"=mua,"sigma0"=fit$par[2],"sigma1"=fit$par[3],"shape"=fit$par[4])))
		#                 }
	}
	parmax <-  - parmax
	if(abs(max(parmax)-overallmax)>0.0001){
		print("non equal mle")
		print(max(parmax))
		print(overallmax)
	}
	crit <- overallmax - qchisq(0.999, 1)/2
	cond <- parmax > crit
	cond <- parmax > -10^6
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
	print(ci)
	out <- c(min(ci), mle[2], max(ci))
	names(out) <- c("LowerCI", "Estimate", "UpperCI")
	print(out)
	browser()
        out
}
# debug(gev.ratio.ic.mu)
gev.ratio.ic.mu(xp=104,t0=t0,t1=t1,y.fit=y.fit,ydat=ydat,to.plot=TRUE,like.num=500)
