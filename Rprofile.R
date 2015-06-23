library(MASS)
library(boot)
library(msm)
library(evmix)
library(ismev)
source("../quantreg/quantreg_Algo.R")
source("evmix_Algo.R")
source("gpd.q.ic.simple.R")


n=350
t=seq(1,n)
mu=(t>100)*.01*(t-100)+100
sigma=(t>100)*.01*(t-100)+1
# sigma=rep(1,n)
y=mapply(rnorm,1,mu,sigma)
# l=1.5
# l=1
# y=(l*y+1)^(1/l)
#  sigma=rep(1,n)



tau=c(seq(10,90,10)/100,0.99,0.9999)
q.theo=t(mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=tau)))
q.theo95=mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=0.95))
q.theo99=mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=0.99))
q.theo995=mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=0.995))
q.theo999=mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=0.999))
q.theo995=c(q.theo995,q.theo995)
mediane=q.theo[,5] 
# paramétrique linéaire
fit.ln=several.rq(y~t,tau=tau)
fit.ln=several.rq(y~mediane,tau=tau)
# ic.ln=several.ic(fit.ln,level=0.95,type="none")
ic.ln=several.ic(fit.ln,level=0.95)
# Non paramétrique avec lissage par spline
df=3
degree=3
fit.spline=several.rq( y~ bs(t, df=df,degree=degree), tau=tau)
fit.spline=several.rq( y~ bs(mediane, df=df,degree=degree), tau=tau)
# ic.spline=several.ic(fit.spline,type="none")
ic.spline=several.ic(fit.spline)
q.hat=list("linear"=ic.ln,"spline"=ic.spline)
q.hat=list("linear"=ic.ln)
p=plot.rq.quantile(t,y,t2e=c(0.10,0.5,0.9,0.9999),tau=tau,q.theo=q.theo,q.hat=q.hat)
plot(p)

gpd.ratio.ic=function(xp,t0,t1,y.fit,ydat,ci.p=0.95,like.num = 1000 ,siglink = identity, show = TRUE, method = "Nelder-Mead", maxit = 10000,to.plot=FALSE, ...){
	if (t0 >= t1)
		stop("t0 must be inferior to t1")
	mle=y.fit$mle
	muscsh=getparam(y.fit,ydat,sigl=c(2))
	p0=evmix::pgpd(xp,muscsh$mu[t0],muscsh$sig[t0],muscsh$sha[t0],lower.tail=FALSE)
	p1=evmix::pgpd(xp,muscsh$mu[t1],muscsh$sig[t1],muscsh$sha[t1],lower.tail=FALSE)
	#         print(p0/p1)
	#         print(p0)
	#         print(p1)
        u=y.fit$threshold
	u0=u[t0]
	u1=u[t1]
	if ( xp<u0 | xp<u1)
		stop(" xp is under threshold")
	rate=y.fit$rate
	xdat=y.fit$xdata
	covariate=ydat[,2]
	translation=-covariate[t0]
	t.delta=covariate[t1]+translation
#	print(translation)
	ydat.t=ydat
	ydat.t[,2]=ydat[,2]+translation
	y.fit.t=gpd.fit(xdat=xdat,threshold=y.fit$threshold,ydat=ydat.t,sigl=c(2),show=FALSE)
	init <- y.fit.t$mle[-2]
	#         print(init)
	xdatu <- xdat[xdat > u]
	xind <- which(xdat > u)
	u <- u[xind]
	sigmat <- cbind(rep(1, length(xdatu)), ydat[xind,2]+translation)
	shmat <- rep(1, length(xdatu))
	gpd.lik <- function(par,ratio,xpi,rate,u0,u1) {
		b=par[1]
		xi=par[2]
		R=ratio^(-xi)
		C0=xi*(xpi-u0)
		C1=xi*(xpi-u1)
		a=b/t.delta*(b+C0-R*(b+C1))/(b*(R-1)-C0)
		if(is.na(a))
			browser()
		sig=c(b,a)
		# calculates gpd neg log lik
		xi <- shmat * xi
		sc <- siglink(sigmat %*% (sig))
		y <- (xdatu - u)/sc
		y <- 1 + xi * y
		if(min(sc) <= 0){
			l <- 10^6
			#                         print(ratio)
		}
		else {
			if(min(y) <= 0)
				l <- 10^6
			else {
				l <- sum(log(sc)) + sum(log(y) * (1/xi + 1))
			}
		}
		l
        }
	#         debug(gpd.lik)
	ratio.l <- sort(c(seq(0,1,length.out=like.num),p0/p1))[-1]
	#         print(ratio.l)
	parmax=numeric(length(ratio.l))
	for(i in 1:length(ratio.l)) {
		fit <- optim(par=c(init), gpd.lik, hessian = TRUE, method = method, control = list(maxit = maxit, ...),xpi=xp,rate=rate,u0=u0,u1=u1,ratio=ratio.l[i])
	    	parmax[i] <- fit$value
		if(ratio.l[i]==p0/p1){
			#                         print(fit)
			#                         browser()
		}
	}
	parmax <-  - parmax
	#         plot(ratio.l,parmax)
	#         scan(n=1)
	overallmax <- -gpd.lik(init,p0/p1,xp,rate,u0,u1)
	overallmax <- max(parmax)
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
getFAR.theo=function(xp,t0,t1){
	p0=1-pnorm(xp,mean=mu[t0],sd=sigma[t0])
	p1=1-pnorm(xp,mean=mu[t1],sd=sigma[t1])
	p0/p1
}

Code2Execute <- function(to.plot=FALSE){
getFAR.theo=function(xp,t0,t1){
	p0=1-pnorm(xp,mean=mu[t0],sd=sigma[t0])
	p1=1-pnorm(xp,mean=mu[t1],sd=sigma[t1])
	p0/p1
}
# sigma=rep(1,n)
y=mapply(rnorm,1,mu,sigma)
qthreshold=0.85
tau=c(seq(10,90,10)/100,0.99,0.9999)
q.theo=t(mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=tau)))
covariate=q.theo[,5]
threshold=predict(rq(y~covariate,tau=qthreshold))
ydat=cbind(1,covariate,1)
y.fit=gpd.fit(xdat=y,threshold=threshold,npy=1,ydat=ydat,sigl=c(2),show=FALSE)
xp=107
t0=325
t1=340
r.ic=gpd.ratio.ic(xp=xp,t0=t0,t1=t1,y.fit=y.fit,ydat=ydat,to.plot=to.plot)
r.theo=getFAR.theo(xp=xp,t0=t0,t1=t1)
# print(r.theo)
# print(r.ic)
res=c(r.theo,r.ic)
names(res)[1]="Theoric"
res
}

Coverage <- function(Code2Execute,n=100){
	function(){
		#                 print(Code2Execute)
		#                 debug(Code2Execute)
		print(n)
	res=t(sapply(1:n,function(x){
			   if(x%%100==0)print(x)
			   if(x==n)
				   Code2Execute(to.plot=TRUE)
			   else
				   Code2Execute()
			}
	))
	res=as.data.frame(res)
	#         print(with(res,mean(Theoric >= LowerCI & Theoric <= UpperCI)))
	invisible(res)
	}
}
# Coverage(Code2Execute)

TestSeveralConf <- function(listn=350,listrepet=1:3,lthreshold=seq(0.70,0.95,length=3),Coverage){
	changement=100
	res=expand.grid(lthreshold,listrepet,listn)
	names(res)=c("threshold","repet","n")
	l.coverage=c()
	xp=107
	t0=230
	t1=300
	i=1
	for (n in listn){
		for( repet in listrepet){
			t=seq(1,n)
			mu=(t>changement)*.01*(t-changement)+100
			sigma=(t>changement)*.01*(t-changement)+1
			t=rep(t,repet)
			mu=rep(mu,repet)
			sigma=rep(sigma,repet)
			for (threshold in lthreshold){
				print(i)
				trysection=try({
					print(Coverage)
					ICs=Coverage()
					print(ICs)
					coverage=with(ICs,mean(Theoric >= LowerCI & Theoric <= UpperCI))
					print(coverage)
				})
				if (class(trysection) == "try-error") 
					coverage=NA
				l.coverage=c(l.coverage,coverage)
				i=i+1
			}
		}
	}
	res=cbind(res,l.coverage)
	names(res)[4]="coverage"
	res
}
tsc=TestSeveralConf(listrepet=3:6,lthreshold=seq(0.75,0.95,length=5))	
tsc=TestSeveralConf(listrepet=3:6,lthreshold=0.95)	


getP <- function(p2,y.fit,x,y,covariate,to.plot=TRUE){
	xcord=p2[1]
	ycord=p2[2]
	if (xcord < min(x) | xcord > max(x))
		stop(" Point outside of data range ")
	xcloser=which((abs(x-xcord))==min(abs(x-xcord)))
	xcloser=min(unique(x[xcloser]))
	ccloser=unique(covariate[which(x==xcloser)])
	mu=unique(y.fit$threshold[which(xcloser==x)])
	sigma=y.fit$mle[1]+ccloser*y.fit$mle[2]
	shape=y.fit$mle[3]
	phi=y.fit$rate
	if (ycord>mu)
		res=evmix::pgpd(ycord,mu,sigma,shape,phi)
	else {
		findP <- function(par,x,y,covariate){
			dframe=data.frame("x"=x,"y"=y,"covariate"=covariate)
			rq.fitted=rq(y~covariate,data=dframe,tau=par)
			ndata=data.frame("covariate"=ccloser)
			predicted=predict.rq(rq.fitted,newdata=ndata)
			#                         print((ycord-predicted)^2)
			print(par)
			abs(ycord-predicted)
		}
		#                 debug(findP)
		res.optim=optimize(findP,interval=c(0,1-phi),x=x,y=y,covariate=covariate,tol=0.001)
		res=res.optim$minimum
		if(to.plot){
			dframe=data.frame("x"=x,"y"=y,"covariate"=covariate)
			dframe=dframe[with(dframe,order(x)),]
			with(dframe,plot(x,y,col="grey"))
			lines(unique(cbind(dframe$x,predict(rq(y~covariate,tau=res,data=dframe)))),col="red")
			points(xcord,ycord,lwd=2,col="blue")
		}
	}		
	res
}
y.fit=res$y.fit
getP(c(1,103),y.fit,x=c(t,t),y=y,covariate=variable)
getP(c(500,115),y.fit,x=c(t,t),y=y,covariate=variable)

getFAR <- function(p1,x2,y.fit,x,y,covariate){
	prob1=1-getP(p1,y.fit,x,y,covariate)
	p2=c(x2,p1[2])
	prob2=1-getP(p2,y.fit,x,y,covariate)
	FAR=1-(prob2/prob1)
	FAR
}
getFAR(c(750,104),500,y.fit=y.fit,x=c(t,t),y=y,covariate=variable)

FARBoot <- function(dat,indice,qthreshold,p1,x2){
	data.b=dat[indice,]
	rq.fitted.b=rq(y~covariable,data=data.b,tau=qthreshold)
	threshold.b=predict(rq.fitted.b)
	ydat.b=cbind(1,data.b$covariable,1)
	y.fit=gpd.fit(xdat=data.b$y,threshold=threshold.b,npy=1,ydat=ydat.b,sigl=c(2),show=FALSE)
	getFAR(p1,x2,y.fit,x=data.b$t,y=data.b$y,covariate=data.b$covariable)
}

xp=107
t0=230
t1=300
n=350
repet=2
qthreshold=0.90
changement=100
t=seq(1,n)
mu=(t>changement)*.01*(t-changement)+100
sigma=(t>changement)*.01*(t-changement)+1
t=rep(t,repet)
mu=rep(mu,repet)
sigma=rep(sigma,repet)
y=mapply(rnorm,1,mu,sigma)
covariate=(mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=0.5)))
data2boot=data.frame("t"=t,"y"=y,"covariable"=covariate)
boot.res=boot(data=data2boot,statistic=FARBoot,R=250,qthreshold=qthreshold,p1=c(t1,xp),x2=t0)
r.boot=mean(boot.res$t)
alpha=0.05
ic.boot=quantile(boot.res$t,p=c(alpha/2,1-alpha/2))
print(ic.boot)
getFAR.theo=function(xp,t0,t1){
	p0=1-pnorm(xp,mean=mu[t0],sd=sigma[t0])
	p1=1-pnorm(xp,mean=mu[t1],sd=sigma[t1])
	1-p0/p1
}
print(getFAR.theo(xp,t0,t1))

Boot2Execute <- function(to.plot=FALSE){
	y=mapply(rnorm,1,mu,sigma)
	getFAR.theo=function(xp,t0,t1){
		p0=1-pnorm(xp,mean=mu[t0],sd=sigma[t0])
		p1=1-pnorm(xp,mean=mu[t1],sd=sigma[t1])
		1-p0/p1
	}
	y=mapply(rnorm,1,mu,sigma)
	covariate=(mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=0.5)))
	data2boot=data.frame("t"=t,"y"=y,"covariable"=covariate)
	boot.res=boot(data=data2boot,statistic=FARBoot,R=250,qthreshold=qthreshold,p1=c(t1,xp),x2=t0)
	r.boot=mean(boot.res$t)
	alpha=0.05
	ic.boot=quantile(boot.res$t,p=c(alpha/2,1-alpha/2))
	r.ic=c(ic.boot[1],r.boot,ic.boot[2])
	names(r.ic) <- c("LowerCI", "Estimate", "UpperCI")
	r.theo=getFAR.theo(xp=xp,t0=t0,t1=t1)
	res=c(r.theo,r.ic)
	names(res)[1]="Theoric"
	res
}
CoverageBoot=Coverage(Boot2Execute)
tsc=TestSeveralConf(listrepet=10,lthreshold=0.95,Coverage=CoverageBoot)	

Prof2Execute <- function(to.plot=FALSE){
	getFAR.theo=function(xp,t0,t1){
		p0=1-pnorm(xp,mean=mu[t0],sd=sigma[t0])
		p1=1-pnorm(xp,mean=mu[t1],sd=sigma[t1])
		p0/p1
	}
	y=mapply(rnorm,1,mu,sigma)
	covariate=(mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=0.5)))
	threshold=predict(rq(y~covariate,tau=qthreshold))
	ydat=cbind(1,covariate,1)
	y.fit=gpd.fit(xdat=y,threshold=threshold,npy=1,ydat=ydat,sigl=c(2),show=FALSE)
	r.ic=gpd.ratio.ic(xp=xp,t0=t0,t1=t1,y.fit=y.fit,ydat=ydat,to.plot=to.plot)
	r.theo=getFAR.theo(xp=xp,t0=t0,t1=t1)
	# print(r.theo)
	# print(r.ic)
	res=c(r.theo,r.ic)
	names(res)[1]="Theoric"
	res
}
CoverageProf=Coverage(Prof2Execute)
tsc.prof=TestSeveralConf(listrepet=3:6,lthreshold=0.95,Coverage=CoverageProf)	
