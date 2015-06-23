dir.create("~/R/Anetto",recursive=TRUE)
.libPaths("~/R/Anetto") 
source("iFARplot_GEV_Parallel_Algo.R")
source("CoverageGEV_Algo.R")
source("CoverageGEV_Parallel_Algo.R")

# CoverageBoot=Coverage(Boot2Execute,n=50)
# debug(CoverageBoot)
# A=CoverageBoot(environment())

CoverageBootSpline=Coverage(Boot2Execute.Spline,n=100)
# A=CoverageProf(environment())

tsc.bootspline=TestSeveralConf(listxp=seq(100,106,3),listt0=1900,listt1=seq(2000,2150,50),listxi=c(-0.1,0,0.1),Coverage=CoverageBootSpline)

write.table(tsc.bootspline,file="tsc.bootspline.txt",sep="\t")
