source("FAR_GEV_Algo.R")

Coverage <- function(Code2Execute,n=100){
	function(env){
		Code2Execute=get("Code2Execute",parent.env(environment()))
		assign("Code2Execute",Code2Execute,envir=env)
		res=t(sapply(1:n,function(x){
			   if(x==n)
				   Code2Execute(env,to.plot=TRUE)
			   else{
				   print(x)
				   Code2Execute(env)
			   }
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
	theta.boot=colMeans(boot.res$t)
	r.boot=mean(boot.res$t)
	r.boot=theta.boot[1]
	if(is.na(r.boot))
		browser()
	alpha=0.05
	ic.boot=quantile(boot.res$t,p=c(alpha/2,1-alpha/2))
	ic.boot=apply(boot.res$t,2,quantile,p=c(alpha/2,1-alpha/2))
	boot.ic=c(ic.boot[1,1],r.boot,ic.boot[2,1])
	i0=min(which((abs(ydat$year-t0))==min(abs(ydat$year-t0))))
	i1=min(which((abs(ydat$year-t1))==min(abs(ydat$year-t1))))
	r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu,sigma,shape)
	names(boot.ic) <- c("LowerCI", "Estimate", "UpperCI")
	res=c(r.theo,boot.ic)
	names(res)[1]="Theoric"
	res
}

Prof2Execute <- function(env,to.plot=FALSE){
	list2env(as.list(env),envir=environment())
	y=mapply(revd,1,mu,sigma,shape)
	ydat=data.frame(year,y,rep(1,n*repet),covariate,rep(1,n*repet),covariate,rep(1,n*repet))
	names(ydat)=c("year","y","mub","mua","sigb","siga","xi")
	y.fit=fevd(y,ydat,location.fun=~mua,scale.fun=~siga)
	i0=min(which((abs(ydat$year-t0))==min(abs(ydat$year-t0))))
	i1=min(which((abs(ydat$year-t1))==min(abs(ydat$year-t1))))
	r.ic=gev.ratio.ic.mu(xp=xp,t0=i0,t1=i1,y.fit=y.fit,ydat=ydat)
	#         r.ic=r.ic[1:3]
	#         r.ic=r.ic[c(3,2,1)]
	r.ic=r.ic[3:1]
	r.ic=1-r.ic
	names(r.ic) <- c("LowerCI", "Estimate", "UpperCI")
	r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu,sigma,shape)
	# print(r.theo)
	# print(r.ic)
	res=c(r.theo,r.ic)
	names(res)[1]="Theoric"
	res
}
