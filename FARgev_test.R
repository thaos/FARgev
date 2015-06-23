library(ismev)
library(evmix)
library(evir)
library(extRemes)
years=1850:2200
months=1:12
mxy=expand.grid(months,years)
n=nrow(mxy)
changement=12*20+1
changement=1
t=seq(1,n)
mu=(t>changement)*.001*(t-changement)+100
sigma=(t>changement)*.001*(t-changement)+1
y=mapply(rnorm,1,mu,sigma)
database=cbind(mxy,y)
names(database)=c("Mois","Année","y")
d.max=aggregate(y~Année,FUN=max,data=database)
d.max.t=d.max
d.max.t$Année.T=d.max.t$Année-1950
ydat=cbind(1,d.max$Année,1,d.max$Année,1)
ydat.t=ydat
ydat.t[,2]=d.max.t$Année.T
y.fit=gev.fit(d.max$y,ydat,mul=2,sigl=4)
y.fit.t=gev.fit(d.max$y,ydat.t,mul=2,sigl=4)
y.fit.ex=fevd(y,d.max,location.fun=~Année,scale.fun=~Année)
y.fit.t.ex=fevd(y,d.max.t,location.fun=~Année.T,scale.fun=~Année)

t0=with(d.max,which(Année==1860))
t1=with(d.max,which(Année==1878))
gev.ratio.ic(xp=108,t0=t0,t1=t1,y.fit=y.fit,ydat=ydat,to.plot=TRUE)

# gev.ratio.ic=function(xp,t0,t1,y.fit,ydat,ci.p=0.95,like.num = 1000 ,mulink=identity, siglink = identity, show = TRUE, method = "Nelder-Mead", maxit = 10000,to.plot=FALSE, ...){
#         require(evir)
#         if (t0 >= t1)
#                 stop("t0 must be inferior to t1")
#         mle=y.fit$mle
#         muscsh=data.frame(y.fit$vals)
#         names(muscsh)=c("mu","sig","sha")
#         p0=1-evir::pgev(xp,mu=muscsh$mu[t0],sigma=muscsh$sig[t0],xi=muscsh$sha[t0])
#         p1=1-evir::pgev(xp,mu=muscsh$mu[t1],sigma=muscsh$sig[t1],xi=muscsh$sha[t1])
#         p0/p1
#         xdat=y.fit$data
#         covariate=ydat[,4]
#         translation=-covariate[t0]
#         t.delta=covariate[t1]+translation
# print(translation)
#         ydat.t=ydat
#         ydat.t[,2]=ydat[,2]+translation
#         y.fit.t=gev.fit(d.max$y,ydat.t,mul=2,sigl=4)
#         d.max.t$Année=d.max.t$Année+translation
#         y.fit.t.ex=fevd(y,d.max.t,location.fun=~Année,scale.fun=~Année)
#         init <- y.fit.t$mle[-4]
# mle <- y.fit.t.ex$results$par
# init <- y.fit.t.ex$results$par[-4]
#         mumat=ydat.t[,1:2]
#         sigmat=ydat.t[,3:4]
#         shmat=ydat.t[,5]
#         mu <- (mumat %*% (mle[1:2]))
#         sig <- (sigmat %*% (mle[3:4]))
#         sha <- (shmat * (mle[5]))
#         p0=1-evir::pgev(xp,mu=mu[t0],sigma=sig[t0],xi=sha[t0])
#         p1=1-evir::pgev(xp,mu=mu[t1],sigma=sig[t1],xi=sha[t1])
#         print(p0/p1)
# print(init)
#         gev.lik <- function(par,ratio,xpi) {
#                 xi=par[4]
#                 mua=par[2]
#                 mub=par[1]
#                 sigb=par[3]
#                 mu <- mulink(mumat %*% (c(mub,mua)))
#                 mu0=mu[t0]
#                 mu1=mu[t1]
#                 C0=xi*(mu0-xpi)
#                 C1=xi*(mu1-xpi)
#                 M=1-(1/ratio)*(1-exp(-(1-C0/sigb)^(-1/xi)))
#                 lM=((-log(M))^(-xi))-1
#                 siga=(C1-sigb*(lM))/(t.delta*lM)
#                 T0=1-exp(-(1-(xi/sigb)*(mub-xp))^(-1/xi))
#                 T0=(1/ratio)*T0
#                 T0=(-log(1-T0))^(-xi)
#                 T0=1-T0
#                 siga=(C1-sigb*T0)/(t.delta*T0)
#                 siga
#         }
#         debug(gev.lik)
#         gev.lik(par=init,ratio=p0/p1,xpi=xp)
# }
# gev.ratio.ic(xp=108,t0=t0,t1=t1,y.fit=y.fit,ydat=ydat,to.plot=TRUE)


gev.ratio.ic.mu=function(xp,t0,t1,y.fit,ydat,ci.p=0.95,like.num = 1000 ,mulink=identity, siglink = identity, show = TRUE,  method="Nelder-Mead",maxit = 10000,to.plot=FALSE, ...){
	require(evir)
	if (t0 >= t1)
		stop("t0 must be inferior to t1")
	mle=y.fit$mle
	muscsh=data.frame(y.fit$vals)
	names(muscsh)=c("mu","sig","sha")
	p0=1-evir::pgev(xp,mu=muscsh$mu[t0],sigma=muscsh$sig[t0],xi=muscsh$sha[t0])
	p1=1-evir::pgev(xp,mu=muscsh$mu[t1],sigma=muscsh$sig[t1],xi=muscsh$sha[t1])
	p0/p1
	xdat=y.fit$data
	covariate=ydat[,4]
	translation=-covariate[t0]
	t.delta=covariate[t1]+translation
#	print(translation)
	ydat.t=ydat
	ydat.t[,2]=ydat[,2]+translation
	y.fit.t=gev.fit(d.max$y,ydat.t,mul=2,sigl=4)
	y.fit.t.ex=fevd(y,d.max.t,location.fun=~Année.T,scale.fun=~Année,optim.arg=list(method=method))
	xdat=y.fit.t.ex$x
	init <- y.fit.t$mle[-4]
	mle=y.fit.t$mle
	mle <- y.fit.t.ex$results$par
	init <- y.fit.t.ex$results$par[-2]
	mumat=ydat.t[,1:2]
	sigmat=ydat.t[,3:4]
	shmat=ydat.t[,5]
	mu <- (mumat %*% (mle[1:2]))
	sig <- (sigmat %*% (mle[3:4]))
	sha <- (shmat * (mle[5]))
	p0=1-evir::pgev(xp,mu=mu[t0],sigma=sig[t0],xi=sha[t0])
	p1=1-evir::pgev(xp,mu=mu[t1],sigma=sig[t1],xi=sha[t1])
	print(p0/p1)
	browser()
	#         print(init)
	init.o=mle
	gev.lik.o <- function(par) {
		xi=par[5]
		mub=par[1]
		mua=par[2]
		sigb=par[3]
		siga=par[4]
		sc <- siglink(sigmat %*% (c(sigb,siga)))
		mu <- mulink(mumat %*% (c(mub,mua)))
		y <- (xdat - mu)/sc
		y <- 1 + xi * y
		if(any(y <= 0) || any(sc <= 0)) return(10^6)
		sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y) * (1/xi + 1))
		levd(x=xdat,location=mu, scale=sc, shape=xi,type="GEV",npy=1)
	}
	cat("Normal")
	print(gev.lik.o(init.o))
	gev.lik <- function(par,ratio,xpi) {
		xi=par[4]
		mub=par[1]
		sigb=par[2]
		siga=par[3]
		sigma <- siglink(sigmat %*% (c(sigb,siga)))
		sc=sigma
		sig0=sigma[t0]
		sig1=sigma[t1]
		tc=try({
		T0=1-exp(-(1-(xi/sig0)*(mub-xp))^(-1/xi))
		T0=(1/ratio)*T0
		T0=(-log(1-T0))^(-xi)
		T0=1-T0
		T0=T0-(xi/sig1)*(mub-xp)
		})
		if(class(tc)=="try-error" | is.na(T0))
			return(10^6)
		mua=T0*sig1/(xi*t.delta)
		mua
		mu <- mulink(mumat %*% (c(mub,mua)))
		y <- (xdat - mu)/sc
		y <- 1 + xi * y
		cat("mua")
		print(mua)
		if(any(y <= 0) || any(sc <= 0)) return(10^6)
		sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y) * (1/xi + 1))
		levd(x=xdat,location=mu, scale=sc, shape=xi,type="GEV",npy=1)
	}
	#         debug(gev.lik)
	cat("Paramétrisation")
	print(gev.lik(par=init,ratio=p0/p1,xpi=xp))
	ratio.l <- sort(c(seq(0,2,length.out=like.num),p0/p1))[-1]
	#         ratio.l <- sort(p0/p1)
	#         print(ratio.l)
	overallmax <- -gev.lik(init,p0/p1,xp)
	parmax=numeric(length(ratio.l))
	for(i in 1:length(ratio.l)) {
		#                 fit <- optim(par=c(init), gev.lik, hessian = TRUE, method = "BFGS",gr=grlevd, control = list(maxit = maxit, ...),xpi=xp,ratio=ratio.l[i])
		fit <- optim(par=c(init), gev.lik, hessian = TRUE, method = "BFGS",control = list(maxit = maxit, ...),xpi=xp,ratio=ratio.l[i])
		#                 parmax[i]=gev.lik(init,ratio.l[i],xp)
		parmax[i] <- fit$value
		#                 if(ratio.l[i]==p0/p1){
		#                 if(-parmax[i]>overallmax){
		#                         print(fit)
		#                         browser()
		#         }
	}
	parmax <-  - parmax
	#         plot(ratio.l,parmax)
	#         scan(n=1)
	overallmax <- -gev.lik(init,p0/p1,xp)
	#         overallmax=max(parmax)
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
	smth <- spline(ratio.l[cond], parmax[cond], n = 200)
	if(to.plot){
		plot(ratio.l, parmax, type = "l", xlab = "", ylab = "")
		abline(h = overallmax - aalpha/2, lty = 2, col = 2)
		lines(smth$x,smth$y, lty = 2, col = 2)
	}
	ci <- smth$x[smth$y > overallmax - aalpha/2]
	out <- c(min(ci), p0/p1, max(ci))
	names(out) <- c("LowerCI", "Estimate", "UpperCI")
        out
}
# debug(gev.ratio.ic.mu)
gev.ratio.ic.mu(xp=104,t0=t0,t1=t1,y.fit=y.fit,ydat=ydat,to.plot=TRUE)



