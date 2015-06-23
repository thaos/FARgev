dir.create("~/R/Anetto",recursive=TRUE)
.libPaths("~/R/Anetto") 
source("iFARplot_GEV_Parallel_Algo.R")
source("CoverageGEV_Algo.R")
source("CoverageGEV_Parallel_Algo.R")
source("GEV_ICratio_mu_alabama_IR2.R")

# CoverageBoot=Coverage(Boot2Execute,n=50)
# debug(CoverageBoot)
# A=CoverageBoot(environment())

CoverageProf=Coverage(Prof2Execute,n=100)
# A=CoverageProf(environment())

# tsc.prof=TestSeveralConf(listxp=c(106),listt0=1900,listt1=2150,listxi=c(-0.1),Coverage=CoverageProf,listrepet=3)
tsc.prof=TestSeveralConf(listxp=seq(100,106,3),listt0=1900,listt1=2150,listxi=c(-0.1,0,0.1),Coverage=CoverageProf,listrepet=3)	
write.table(tsc.prof,file="tsc.prof.txt",sep="\t")
