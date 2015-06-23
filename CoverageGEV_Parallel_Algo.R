source("iFARplot_GEV_Parallel_Algo.R")
source("CoverageGEV_Algo.R")

TestSeveralConf <- function(listxi=c(-0.1,.Machine$double.eps,0.1),listt0=c(1880,1900,1970),listt1=c(2100,2150,2180),listxp=c(100,103,106),listrepet=1:3,Coverage){
			    	#         thisEnv=environment()
	i=1
	res=expand.grid(listrepet,listxp,listt1,listt0,listxi)
	names(res)=c("repet","xp","t1","t0","xi")
	Coverage2Run=function(repet,xp,t1,t0,xi){
		thisEnv=environment()
		changement=100
		years=1850:2200
		#                 years=years:(years+ln)
		#                 t0=1920
		#                 t1=2102
		#                 xp=101
		n=length(unique(years))
		t=seq(1,n)
		mu=(t>changement)*.020*(t-changement)+100
		sigma=(t>changement)*.001*(t-changement)+1
		shape=rep(xi,n)
		mu=rep(mu,repet)
		sigma=rep(sigma,repet)
		t=rep(t,repet)
		shape=rep(shape,repet)
		year=rep(years,repet)
		covariate=(mapply(qevd,loc=mu,scale=sigma,shape=shape,MoreArgs=list("p"=0.5)))
		trysection=try({
				#                                 ICs=Coverage(thisEnv)
				ICs=Coverage(thisEnv)
				#                                 browser()
				#                                         print(ICs)
				coverage=with(ICs,mean(Theoric >= LowerCI & Theoric <= UpperCI))
				#                                         print(coverage)
			})
		if (class(trysection) == "try-error") {
			browser()
			coverage=NA
		}
		coverage
	}
	#         debug(Coverage2Run)
	l.coverage=with(res,mapply(Coverage2Run,repet=repet,xp=xp,t1=t1,t0=t0,xi=xi))
	res=cbind(res,l.coverage)
	names(res)[3]="coverage"
	res
}

Coverage <- function(Code2Execute,n=100){
	function(env){
		require(parallel)
		Code2Execute=get("Code2Execute",parent.env(environment()))
		#                 print(env)
		#                 print(ls(,envir=env))
		#parent.env(environment(Code2Execute))=env
		Code2Run=function(x){
			   if(x==n)
				   Code2Execute(env,to.plot=TRUE)
			   else
				   Code2Execute(env)
		}
	res=t(mapply(Code2Run,x=1:n))
	res=as.data.frame(res)
	#         print(with(res,mean(Theoric >= LowerCI & Theoric <= UpperCI)))
	invisible(res)
	}
}

