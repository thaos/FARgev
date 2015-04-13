graphics.off()
rm(list=ls())

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


# repet=4
# years=1850:2200
# n=length(unique(years))
# changement=100
# t=seq(1,n)
# mu=(t>changement)*.020*(t-changement)+100
# sigma=(t>changement)*.001*(t-changement)+1
# shape=rep(-0.1,n)
# mu=rep(mu,repet)
# sigma=rep(sigma,repet)
# t=rep(t,repet)
# shape=rep(shape,repet)
# year=rep(years,repet)
# covariate=(mapply(qgev,mu=mu,sigma=sigma,xi=shape,MoreArgs=list("p"=0.5)))
# y=mapply(revd,1,mu,sigma,rep(1,n))
# ydat=data.frame(year,y,rep(1,n),covariate,rep(1,n),covariate,rep(1,n))
# names(ydat)=c("year","y","mub","mua","sigb","siga","xi")
# t0=1920
# t1=2102
# xp=250
# boot.res=boot(data=ydat,statistic=FARBoot,R=250,p1=c(t1,xp),x2=t0)
# r.boot=mean(boot.res$t)
# alpha=0.05
# ic.boot=quantile(boot.res$t,p=c(alpha/2,1-alpha/2))
# boot.ic=c(ic.boot[1],r.boot,ic.boot[2])
# i0=min(which((abs(ydat$year-t0))==min(abs(ydat$year-t0))))
# i1=min(which((abs(ydat$year-t1))==min(abs(ydat$year-t1))))
# r.theo=getFAR.theo(xp=250,t0=i0,t1=i1,mu,sigma,shape)
# 

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

iFARplot <- function() {
	thisEnv <- environment()
	FarHat=NULL
	repet=4
	years=1850:2200
	n=length(unique(years))
	changement=100
	t=seq(1,n)
	mu=(t>changement)*.020*(t-changement)+100
	sigma=(t>changement)*.001*(t-changement)+1
	shape=rep(-0.1,n)
	mu=rep(mu,repet)
	sigma=rep(sigma,repet)
	t=rep(t,repet)
	shape=rep(shape,repet)
	year=rep(years,repet)
	covariate=(mapply(qgev,mu=mu,sigma=sigma,xi=shape,MoreArgs=list("p"=0.5)))
	#On ajoute un bruit insignifiant pour éviter que la covariable soit exactementy identique sur la première période, sans quoi, ce la fait bugguer le profile de vraisemblance: division par zéro pour une différence en covariable nulle.
	covariate=covariate+rnorm(n*repet,sd=0.000001)
	y=mapply(revd,1,mu,sigma,shape)
	ydat=data.frame(year,y,mu,sigma,shape,rep(1,n*repet),covariate,rep(1,n*repet),covariate,rep(1,n*repet))
	names(ydat)=c("year","y","mu","sigma","shape","mub","mua","sigb","siga","xi")
	df.plot=ydat[order(t),]
	getDensity.theo=function(tx,mu,sigma,xi){
		ymin=min(y)
		ymax=max(y)
		x=seq(ymin,ymax,length.out=100)
		res=sapply(x,devd,loc=mu[tx],scale=sigma[tx],shape=xi[tx])
		res=cbind(x,res)
		res
	}
	getDensity=function(mu,sigma,xi){
		ymin=min(y)
		ymax=max(y)
		x=seq(ymin,ymax,length.out=100)
		res=sapply(x,devd,loc=mu,scale=sigma,shape=xi)
		res=cbind(x,res)
		res
	}
	updatePlot <- function(h,...){
		layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
		with(df.plot,plot(year,y,main="Synthetic Time Serie of Maxima"))
		startAdjust[]<<-seq(min(year),svalue(endAdjust),by=1)
		endAdjust[]<<-seq(svalue(startAdjust),max(year),by=1)
		minEvent=with(df.plot,min(y[year>svalue(startAdjust)&year<svalue(endAdjust)]))
		maxEvent=with(df.plot,max(y[year>svalue(startAdjust)&year<svalue(endAdjust)]))
		EventAdjust[] <<- seq(minEvent,maxEvent,length=100)
		abline(v=svalue(startAdjust),col="blue")
		text(x=svalue(startAdjust),y=109,pos=2,labels="t0",col="red")
		abline(v=svalue(endAdjust),col="blue")
		text(x=svalue(endAdjust),y=109,pos=4,labels="t1",col="red")
		abline(h=svalue(EventAdjust),col="blue")
		text(y=svalue(EventAdjust),x=2200,pos=1,labels="xp",col="red")
		xp=svalue(EventAdjust)
		t0=svalue(startAdjust)
		t1=svalue(endAdjust)
		i0=min(which((abs(df.plot$year-t0))==min(abs(df.plot$year-t0))))
		i1=min(which((abs(df.plot$year-t1))==min(abs(df.plot$year-t1))))
		d0=getDensity.theo(tx=i0,mu=df.plot$mu,sigma=df.plot$sigma,xi=df.plot$shape)
		d1=getDensity.theo(tx=i1,mu=df.plot$mu,sigma=df.plot$sigma,xi=df.plot$shape)
		cat("JeanClaudeRichard")
		if(!is.null(get("FarHat",envir=thisEnv))){
			cat("PierreMarieSolange")
			d0.prof=with(FarHat,getDensity(mu=Prof[2],sigma=Prof[3],xi=Prof[4]))
			d1.prof=with(FarHat,getDensity(mu=Prof[6],sigma=Prof[7],xi=Prof[8]))
			d0.boot=with(FarHat,getDensity(mu=Boot[2],sigma=Boot[3],xi=Boot[4]))
			d1.boot=with(FarHat,getDensity(mu=Boot[6],sigma=Boot[7],xi=Boot[8]))
			print(all(d0.boot==d0.prof))
		}	
		r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu=df.plot$mu,sigma=df.plot$sigma,xi=df.plot$shape)
		p0.theo=attr(r.theo,"p0")
		p1.theo=attr(r.theo,"p1")
		plot(d0[,1],d0[,2],main="Theoric density at time t0",type="l",xlab="",ylab="density",ylim=c(0,0.5),sub=paste("theoric p0:" ,format(p0.theo,digits=2),sep=""))
		if(!is.null(get("FarHat",envir=thisEnv))){
			lines(d0.prof[,1],d0.prof[,2],col="blue")
			lines(d0.boot[,1],d0.boot[,2],col="green")
			legend("topleft",legend=c("Profile","Boot"),lty=1,lwd=2,col=c("blue","green"),bty="n")
		}
		abline(v=xp)
		text(x=xp,y=0.45,pos=4,labels="xp",col="red")
		plot(d1[,1],d1[,2],main="Theoric density at time t1",type="l",xlab="",ylab="density",ylim=c(0,0.5),sub=paste("theoric p1:" ,format(p1.theo,digits=2),sep=""))
		if(!is.null(get("FarHat",envir=thisEnv))){
			lines(d1.prof[,1],d1.prof[,2],col="blue")
			lines(d1.boot[,1],d1.boot[,2],col="green")
			legend("topleft",legend=c("Profile","Boot"),lty=1,lwd=2,col=c("blue","green"),bty="n")
		}
		abline(v=xp)
		text(x=xp,y=0.45,pos=4,labels="xp",col="red")
	}
	ComputeFAR <- function(h,...){
		svalue(resText) <<- paste("Computing...\n",sep="")
		xp=svalue(EventAdjust)
		t0=svalue(startAdjust)
		t1=svalue(endAdjust)
		i0=min(which((abs(df.plot$year-t0))==min(abs(df.plot$year-t0))))
		i1=min(which((abs(df.plot$year-t1))==min(abs(df.plot$year-t1))))
		print(paste(t0,i0))
		print(paste(t1,i1))
		y.fit=fevd(y,ydat,location.fun=~mua,scale.fun=~siga)
		cat("y fitted")
		r.ic=gev.ratio.ic.mu(xp=xp,t0=i0,t1=i1,y.fit=y.fit,ydat=df.plot)
		r.param=r.ic[4:length(r.ic)]
		r.ic[1:3]=1-r.ic[1:3]
		r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu=df.plot$mu,sigma=df.plot$sigma,xi=df.plot$shape)
		p0.theo=attr(r.theo,"p0")
		p1.theo=attr(r.theo,"p1")
		svalue(resText) <<- paste("FAR pour t0=",format(t0,digits=2)," et t1=",format(t1,digits=2)," et seuil=",format(xp,digits=2),"\n\n",sep="")
		insert(resText,paste("theoric FAR:\n", format(r.theo,digits=2),"\n",sep=""))
		insert(resText,paste("theoric p0\t\ttheoric p1\n", format(p0.theo,digits=2),"\t\t",format(p1.theo,digits=2),"\n",sep=""))
		prof.txt=paste("Profile Likelihood FAR:\nEstimate\t\tIC95_inf\t\tIC95_sup\n",format(r.ic[2],digits=2),"\t\t",format(r.ic[3],digits=2),"\t\t",format(r.ic[1],digits=2),"\n",sep="")
		insert(resText,prof.txt)
		insert(resText,paste("p0.esti\t\tp1.esti\n", format(r.ic[4],digits=2),"\t\t",format(r.ic[8],digits=2),"\n",sep=""))
		cat("Computing Bootstrap IC\n")
		boot.res=boot(data=df.plot,statistic=FARBoot,R=250,p1=c(t1,xp),x2=t0,parallel="snow")
		r.boot=colMeans(boot.res$t)
		if(dev.cur()== dev.next())
			dev.new()
		else
			dev.set(dev.next())
		par(mfrow=c(1,3))
		hist(boot.res$t[,1],main="BootSample : FAR",xlab="FAR",breaks="FD",freq=FALSE)
		abline(v=r.theo,col="red")
		abline(v=r.boot[1],col="green")
		hist(boot.res$t[,2],main="BootSample : p0",xlab="p0",breaks="FD",freq=FALSE)
		abline(v=p0.theo,col="red")
		abline(v=r.boot[2],col="green")
		hist(boot.res$t[,6],main="BootSample : p1",xlab="p1",breaks="FD",freq=FALSE)
		abline(v=p1.theo,col="red")
		abline(v=r.boot[6],col="green")
		dev.set(dev.prev()) 
		alpha=0.05
		ic.boot=apply(boot.res$t,2,quantile,p=c(alpha/2,1-alpha/2))
		i0=min(which((abs(ydat$year-t0))==min(abs(ydat$year-t0))))
		i1=min(which((abs(ydat$year-t1))==min(abs(ydat$year-t1))))
		r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu,sigma,shape)
		r.theo=c(r.theo,attr(r.theo,"p0"),attr(r.theo,"p1"))
		boot.ic=as.data.frame(t(rbind(r.boot,ic.boot)[,c(1,2,6)]))
		boot.ic=cbind(r.theo,boot.ic)
		names(boot.ic)=c("Theoric","Estimate","LowerIC","UpperIC")
		rownames(boot.ic)=c("FAR","p0","p1")
		print(boot.ic)
		boot.txt=paste("Resampling Bootstrap FAR:\nEstimate\t\tIC95_inf\t\tIC95_sup\n",format(boot.ic[1,2],digits=2),"\t\t",format(boot.ic[1,3],digits=2),"\t\t",format(boot.ic[1,4],digits=2),"\n",sep="")
		insert(resText,boot.txt)
		insert(resText,paste("p0.esti\t\tp1.esti\n",format(boot.ic[2,2],digits=2),"\t\t",format(boot.ic[3,2],digits=2),"\n",sep=""))
		assign("FarHat",list("Prof"=r.param,"Boot"=r.boot[2:length(r.boot)]),envir=thisEnv)
		updatePlot()
		assign("FarHat",NULL,envir=thisEnv)
	}
	##The widgets
	win <- gwindow("Fraction of Attributable Risk")
	gp <- ggroup(horizontal=FALSE, cont=win)
	tmp <- gframe("Results", container=win, expand=TRUE)
	resText <- gtext("FAR",cont=tmp,width=320,height=340)
	tmp <- gframe("Start-End", container=gp, expand=TRUE)
	startAdjust <- gslider(from=min(year),to=max(year),by=1, value=min(year),
				   cont=tmp, expand=TRUE,
				   handler=updatePlot)
	endAdjust <- gslider(from=min(year),to=max(year),by=1, value=max(year),
				   cont=tmp, expand=TRUE,
				   handler=updatePlot)
	minEvent=min(df.plot$y)
	maxEvent=max(df.plot$y)
	tmp <- gframe("event.threshold", container=gp, expand=TRUE)
	EventAdjust <- gslider(from=minEvent,to=maxEvent,length=100, value=minEvent,
				   cont=tmp, expand=TRUE,
				   handler=updatePlot)
	gbutton("ComputeFAR",container=tmp,handler=ComputeFAR)
	updatePlot()
}
iFARplot()
