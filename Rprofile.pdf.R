source("Rprofile_Algo.R")
tsc.boot=TestSeveralConf(listn=350,listrepet=1,Coverage=CoverageBoot)	
fsc.prof=TestSeveralConf(listn=501,listrepet=floor(seq(1,20,length=6)),lthreshold=seq(0.70,0.95,0.05),Coverage=CoverageProf)	
