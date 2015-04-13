install.packages("parallel")
library(parallel)
source("iFARplot_GEV.R")
source("CoverageGEV_Algo.R")
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
CoverageProf=Coverage(Prof2Execute,n=20)
A=CoverageProf(env)
