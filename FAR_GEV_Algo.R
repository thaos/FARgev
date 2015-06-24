packages <- c("boot", "quantreg", "evmix", "ismev", "parallel","foreach","doMC","snow")
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
library(snow)

getP <- function(p2,y.fit,ydat,to.plot=FALSE){
	xcord=p2[1]
	ycord=p2[2]
	x=ydat$year
	covariate=ydat$mua
	if (xcord < min(x) | xcord > max(x))
		stop(" Point outside of data range ")
	xcloser=which((abs(x-xcord))==min(abs(x-xcord)))
	xcloser=min(unique(x[xcloser]))
	ccloser=mean(unique(covariate[which(x==xcloser)]))
	mle <- y.fit$results$par
	mu=mle[1]+ccloser*mle[2]
	sigma=mle[3]+ccloser*mle[4]
	shape=mle[5]
	1-evir::pgev(ycord,mu=mu,sigma=sigma,xi=shape)
	res=pevd(ycord,loc=mu,scale=sigma,shape=shape,lower.tail=FALSE)
	res=c(res,mu,sigma,shape)
	names(res)=c("p","mu","sigma","shape")
	res
}
# getP(c(1920,110),y.fit,ydat,to.plot=TRUE)

getFAR <- function(p1,x2,y.fit,ydat){
	prob1=getP(p1,y.fit,ydat)
	p2=c(x2,p1[2])
	prob2=getP(p2,y.fit,ydat)
	if(prob1[1]==0 & prob2[1]==0)
		FAR=1
	else
	FAR=1-(prob2[1]/prob1[1])
	res=c(FAR,prob2,prob1)
	names(res)=c("FAR","p0","mu0","sc0","sh0","p1","mu1","sc1","sh1")
	res
}
# getFAR(c(1920,110),2102,y.fit,ydat)

getFAR.theo=function(xp,t0,t1,mu,sigma,xi){
	p0=1-pevd(xp,loc=mu[t0],scale=sigma[t0],shape=xi[t0])
	p1=1-pevd(xp,loc=mu[t1],scale=sigma[t1],shape=xi[t1])
	print(p0/p1)
	if(p0==0 & p1==0)
		p0p1=1
	else
		p0p1=p0/p1
	res=1-p0p1
	attr(res,"p0")=p0
	attr(res,"p1")=p1
	res
}

FARBoot <- function(ydat,indice,p1,x2){
	y.fit.dat=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,method="MLE")
	init=as.list(y.fit.dat$results$par)
	data.b=ydat[indice,]
	#         y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga)
	y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga,method="MLE",initial=init)
	getFAR(p1,x2,y.fit,ydat)
}

FARBoot_explore <- function(ydat,indice,p1,x2){
	data.b=ydat[indice,]
	y.fit.dat=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,method="MLE")
	init=as.list(y.fit.dat$results$par)
	y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga,method="MLE",initial=init)
	res=getFAR(p1,x2,y.fit,ydat)
	i0=min(which((abs(ydat$year-t0))==min(abs(ydat$year-t0))))
	i1=min(which((abs(ydat$year-t1))==min(abs(ydat$year-t1))))
	r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu=ydat$mu,sigma=ydat$sigma,xi=ydat$shape)
	#         with(ydat,plot(year,y,yaxt="n",col="red",ylim=c(100,108)))
	#         par(new=TRUE)
	#         with(data.b,hist(year,breaks=250))
	#         par(new=TRUE)
	#         with(data.b,plot(year,y,yaxt="n",col="green", ylim=c(100,108)))
	if (abs(res[1] - (r.theo[1])) > 0.2){
		print("chelou")
		print(res)
		#                 browser()
	}
	if (abs(res[1] - (r.theo[1])) <= 0.2) {
		print("normal")
		print(res)
		#                 browser()
	}
	res
}

FARBoot.Spline <- function(ydat,indice,p1,x2){
	data.b=ydat[indice,]
	covariate=with(data.b,predict(rq( y~ bs(year, df=3,degree=3), tau=0.5)))
	data.b$mua=covariate
	data.b$siga=covariate
	y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga)
	getFAR(p1,x2,y.fit,ydat)
}

function2draw <- function(condition){
	function(ydat,nrepet,year){
		cond.e=eval(condition)
		sb=which(cond.e)
		sample(sb,size=nrepet,replace=TRUE)
	}
}

cond=quote(ydat$year==year)
draw_same_year <- function2draw(cond)
rbind_list <- function(liste){
	if(is.data.frame(liste)|length(liste)==1)
	   return(liste)
	if(length(liste)>1){
		liste[[2]]= rbind(liste[[1]],liste[[2]])
		return(rbind_list(liste[-1]))
	}
}
		

FARBoot_sy <- function(ydat,indice,p1,x2){
	y.fit.dat=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,method="MLE")
	init=as.list(y.fit.dat$results$par)
	years=as.data.frame(table(ydat$year))
	years$Var1=as.numeric(paste(years$Var1))
	#         years=aggregate(y~year,data=ydat,FUN=length)
	names(years)=c("year","eff")
	indice=c(with(years,mapply(draw_same_year,nrepet=eff,year=year,MoreArgs=list(ydat=ydat),SIMPLIFY=TRUE)))
	data.b=ydat[indice,]
	y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga,method="MLE",initial=init)
	getFAR(p1,x2,y.fit,ydat)
}

tol=5
cond2=function(tol){
	substitute(ydat$year <= year+tol &ydat$year >= year-tol,list(tol=tol))
}
draw_around_year <- function2draw(cond2(5))
FARBoot_ay <- function(ydat,indice,p1,x2){
	y.fit.dat=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,method="MLE")
	init=as.list(y.fit.dat$results$par)
	years=aggregate(y~year,data=ydat,FUN=length)
	names(years)=c("year","eff")
	indice=c(with(years,mapply(draw_around_year,nrepet=eff,year=year,MoreArgs=list(ydat=ydat),SIMPLIFY=TRUE)))
	data.b=ydat[indice,]
	years.v=c(with(years,mapply(rep,x=year,times=eff)))
	data.b$year=years.v
	y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga,method="MLE",initial=init)
	getFAR(p1,x2,y.fit,ydat)
}

tol=0.1
year2mu <- function(year,ydat){
	i=which(year==ydat$year)
	ydat[i[1],"mua"]
}
cond3 <- function(tol){
	substitute(ydat$mua <= year2mu(year,ydat)+tol &ydat$mua >= year2mu(year,ydat)-tol,list(tol=tol))
}
draw_around_mua <- function2draw(cond3(0.1))
FARBoot_ax <- function(ydat,indice,p1,x2){
	y.fit.dat=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,method="MLE")
	init=as.list(y.fit.dat$results$par)
	years=aggregate(y~year,data=ydat,FUN=length)
	names(years)=c("year","eff")
	indice=c(with(years,mapply(draw_around_mua,nrepet=eff,year=year,MoreArgs=list(ydat=ydat),SIMPLIFY=TRUE)))
	data.b=ydat[indice,]
	years.v=c(with(years,mapply(rep,x=year,times=eff)))
	data.b$year=years.v
	y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga,method="MLE",initial=init)
	getFAR(p1,x2,y.fit,ydat)
}

FARBoot_gen  <- function(to_draw){
	function(ydat,indice,p1,x2){
		y.fit.dat=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,method="MLE")
		init=as.list(y.fit.dat$results$par)
		years=aggregate(y~year,data=ydat,FUN=length)
		names(years)=c("year","eff")
		indice=c(with(years,mapply(to_draw,nrepet=eff,year=year,MoreArgs=list(ydat=ydat),SIMPLIFY=TRUE)))
		data.b=ydat[indice,]
		years.v=c(with(years,mapply(rep,x=year,times=eff)))
		data.b$year=years.v
		y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga,method="MLE",initial=init)
		getFAR(p1,x2,y.fit,ydat)
	}
}
FARBoot_ax2=FARBoot_gen(draw_around_mua) 
FARBoot_sy2=FARBoot_gen(draw_same_year) 
FARBoot_ay2=FARBoot_gen(draw_around_year) 
