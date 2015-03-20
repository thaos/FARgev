install.packages("parallel")
library(parallel)
source("iFARplot_GEV.R")
TestSeveralConf <- function(listn=350,listrepet=1:3,Coverage){
	thisEnv=environment()
	i=1
	res=expand.grid(listrepet,listn)
	names(res)=c("repet","n")
	for (ln in listn){
		changement=100
		years=1850
		years=years:(years+ln)
		t0=1920
		t1=2102
		xp=103
		n=length(unique(years))
		t=seq(1,n)
		mu=(t>changement)*.020*(t-changement)+100
		sigma=(t>changement)*.001*(t-changement)+1
		shape=rep(-0.1,n)
		for( repet in listrepet){
			l.coverage=c()
			mu=rep(mu,repet)
			sigma=rep(sigma,repet)
			t=rep(t,repet)
			shape=rep(shape,repet)
			year=rep(years,repet)
			covariate=(mapply(qgev,mu=mu,sigma=sigma,xi=shape,MoreArgs=list("p"=0.5)))
			trysection=try({
				#                                 ICs=Coverage(thisEnv)
				ICs=Coverage(thisEnv)
				#                                         print(ICs)
				coverage=with(ICs,mean(Theoric >= LowerCI & Theoric <= UpperCI))
				#                                         print(coverage)
			})
			if (class(trysection) == "try-error") 
				coverage=NA
			l.coverage=c(l.coverage,coverage)
			i=i+1
		}
	}
	res=cbind(res,l.coverage)
	names(res)[3]="coverage"
	res
}

Coverage <- function(Code2Execute,n=100){
	function(env){
		Code2Execute=get("Code2Execute",parent.env(environment()))
		assign("Code2Execute",Code2Execute,envir=env)
		res=t(sapply(1:n,function(x){
			   if(x==n)
				   Code2Execute(env,to.plot=TRUE)
			   else
				   Code2Execute(env)
			}
	))
	res=as.data.frame(res)
	invisible(res)
	}
}

Boot2Execute <- function(env,to.plot=FALSE){
	#         print(environment())
	#         print(ls(,envir=parent.env(environment())))
	list2env(as.list(env),envir=environment())
	y=mapply(revd,1,mu,sigma,shape)
	ydat=data.frame(year,y,rep(1,n*repet),covariate,rep(1,n*repet),covariate,rep(1,n*repet))
	names(ydat)=c("year","y","mub","mua","sigb","siga","xi")
	boot.res=boot(data=ydat,statistic=FARBoot,R=250,p1=c(t1,xp),x2=t0)
	r.boot=mean(boot.res$t)
	if(is.na(r.boot))
		browser()
	alpha=0.05
	ic.boot=quantile(boot.res$t,p=c(alpha/2,1-alpha/2))
	boot.ic=c(ic.boot[1],r.boot,ic.boot[2])
	i0=min(which((abs(ydat$year-t0))==min(abs(ydat$year-t0))))
	i1=min(which((abs(ydat$year-t1))==min(abs(ydat$year-t1))))
	r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu,sigma,shape)
	names(boot.ic) <- c("LowerCI", "Estimate", "UpperCI")
	res=c(r.theo,boot.ic)
	names(res)[1]="Theoric"
	res
}
CoverageBoot=Coverage(Boot2Execute,n=20)

env=globalenv()
changement=100
repet=4
years=1850:2200
t0=1920
t1=2102
xp=105
n=length(unique(years))
t=seq(1,n)
mu=(t>changement)*.005*(t-changement)+100
sigma=(t>changement)*.001*(t-changement)+1
shape=rep(-0.1,n)
mu=rep(mu,repet)
sigma=rep(sigma,repet)
t=rep(t,repet)
shape=rep(shape,repet)
year=rep(years,repet)
covariate=(mapply(qgev,mu=mu,sigma=sigma,xi=shape,MoreArgs=list("p"=0.5)))
y=mapply(revd,1,mu,sigma,shape)
ydat=data.frame(year,y,rep(1,n*repet),covariate,rep(1,n*repet),covariate,rep(1,n*repet))
names(ydat)=c("year","y","mub","mua","sigb","siga","xi")
debug(CoverageBoot)
CoverageBoot(env)
# 
Prof2Execute <- function(to.plot=FALSE){
	list2env(as.list(env),envir=environment())
	y=mapply(revd,1,mu,sigma,shape)
	ydat=data.frame(year,y,rep(1,n*repet),covariate,rep(1,n*repet),covariate,rep(1,n*repet))
	names(ydat)=c("year","y","mub","mua","sigb","siga","xi")
	i0=min(which((abs(ydat$year-t0))==min(abs(ydat$year-t0))))
	i1=min(which((abs(ydat$year-t1))==min(abs(ydat$year-t1))))
	r.ic=gev.ratio.ic.mu(xp=xp,t0=i0,t1=i1,y.fit=y.fit,ydat=ydat)
	r.ic=1-r.ic
	r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu,sigma,shape)
	# print(r.theo)
	# print(r.ic)
	res=c(r.theo,r.ic)
	names(res)[1]="Theoric"
	res
}
CoverageProf=Coverage(Prof2Execute,n=20)
