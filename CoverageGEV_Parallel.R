source("iFARplot_GEV_Parallel.R")
source("CoverageGEV_Algo.R")
source("CoverageGEV_Parallel_Algo.R")

listxi=c(-0.1,0,0.1);listt0=c(1880,1900,1970);listt1=c(2100,2150,2180);listxp=c(100,103,106);listrepet=1:3


repet=4
changement=100
years=1850:2200
t0=1920
t1=2102
xp=101
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
covariate=(mapply(qevd,loc=mu,scale=sigma,shape,shape,MoreArgs=list("p"=0.5)))
# t0=1920
# t1=2102
# xp=101
# n=length(unique(years))
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
y=mapply(revd,1,mu,sigma,shape)
ydat=data.frame(year,y,rep(1,n*repet),covariate,rep(1,n*repet),covariate,rep(1,n*repet))
names(ydat)=c("year","y","mub","mua","sigb","siga","xi")

CoverageBoot=Coverage(Boot2Execute,n=50)
# debug(CoverageBoot)
A=CoverageBoot(environment())

CoverageProf=Coverage(Prof2Execute,n=10)
A=CoverageProf(environment())

sc.boot=TestSeveralConf(listn=350,listrepet=c(1,3,5),Coverage=CoverageBoot)	
tsc.prof=TestSeveralConf(listn=350,listrepet=c(1,3,5),Coverage=CoverageProf)	
tsc.boot=TestSeveralConf(listn=350,listrepet=c(10),Coverage=CoverageBoot)	
tsc.prof=TestSeveralConf(listn=350,listrepet=c(10),Coverage=CoverageProf)	
tsc.prof=TestSeveralConf(Coverage=CoverageProf)	
tsc.prof.xi0=TestSeveralConf(listxi=0,Coverage=CoverageProf)	

tsc.prof.r1=tsc.prof
save(tsc.prof.r1,file="tsc.prof.r1.RData")
tsc.boot.r1=tsc.boot
save(tsc.boot.r1,file="tsc.boot.r1.RData")
