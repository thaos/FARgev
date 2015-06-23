library(extRemes)
library(parallel)
library(evir)
load("wrong.boot.RData")
		
# wrong.boot=ydat
y.fit.dat=fevd(y,ydat,location.fun=~mua,scale.fun=~siga)
init=y.fit.dat$results$par

gev.lik <- function(par,mua,siga) {
	xi=par[3]
	mub=par[1]
	sigb=par[2]
	sc <- (sigmat %*% (c(sigb,siga)))
	mu.v <- (mumat %*% (c(mub,mua)))
	res=levd(x=xdat,location=mu.v, scale=sc, shape=xi,type="GEV",npy=1)
	if(is.na(res) | is.infinite(res))
		return(10^6)
	res
}

gev.lik.allp <- function(par) {
	xi=par[5]
	mub=par[1]
	mua=par[2]
	sigb=par[3]
	siga=par[4]
	sc <- (sigmat %*% (c(sigb,siga)))
	mu.v <- (mumat %*% (c(mub,mua)))
	res=levd(x=xdat,location=mu.v, scale=sc, shape=xi,type="GEV",npy=1)
	if(is.na(res) | is.infinite(res))
		return(10^6)
	res
}

profiler <- function(wrong.boot, ci.p=0.95, pas=0.01, mua_lim=c(0,1), siga_lim=c(0,1),initial=NULL){
	f_roots=function(mub,mua,sigb,siga,shape){
		init=c(mub,sigb,shape)
		fit=optim(par=c(init), gev.lik,mua=mua,siga=siga)
		parmax=-fit$value
		print(i <<- i+1)
		parmax + aalpha/2 - overallmax
		parmax
	}
	mumat=as.matrix(wrong.boot[,c("mub","mua")])
	sigmat=as.matrix(wrong.boot[,c("sigb","siga")])
	shmat=as.matrix(wrong.boot[,"xi"])
	if(!is.null(initial)){
		y.fit=fevd(y,wrong.boot,location.fun=~mua,scale.fun=~siga, initial=initial)
	} else {
		y.fit=fevd(y,wrong.boot,location.fun=~mua,scale.fun=~siga)
	}
	print(y.fit)
	mle <- y.fit$results$par
	init=mle[-c(2,4)]
	mua=mle[2]
	siga=mle[4]
	xdat=wrong.boot$y
	aalpha <- qchisq(ci.p, 1)
	overallmax <- -gev.lik(init,mua,siga) 
	i=1
	demipas=pas/2
	mua_l=seq(mua_lim[1],mua_lim[2],pas)
	siga_l=seq(siga_lim[1],siga_lim[2],pas)
	mu_siga=expand.grid(mua_l, siga_l)
	names(mu_siga)=c("mua","siga")
	mub_l=101.8-mu_siga$mua*102.2
	sigb_l=1.10-mu_siga$siga*102.2
	xi_l=rep(mle[5],nrow(mu_siga))
	mu_siga <- cbind(mu_siga,mub_l, sigb_l , xi)
	names(mu_siga)=c("mua","siga","mub", "sigb", "xi")
	prof2D <- mcmapply(f_roots,mua=mu_siga$mua,siga=mu_siga$siga,mub=mu_siga$mub,sigb=mu_siga$sigb,shape=mu_siga$xi)
	coord_max=which(prof2D == max(prof2D))
	coord_max=mu_siga[coord_max,]
	print(max(prof2D)== aalpha/2)
	attr(prof2D,"mle") = mle
	attr(prof2D,"cmax") = coord_max
	prof2D
}
debug(profiler)

prof2D_wrong = profiler(wrong.boot)
prof2D_wrong_zoom1 = profiler(wrong.boot,mua_lim=c(0, 0.2), siga_lim=c(0, 0.2), pas=0.001)
prof2D_wrong_zoom2 = profiler(wrong.boot, mua_lim=c(0.8, 1), siga_lim=c(0, 0.2), pas=0.001)
prof2D_normal = profiler(ydat)
prof2D_wrong_init = profiler(wrong.boot, initial=as.list(init))

plot_profile <- function(prof2D,seuil_na=NULL,pas=0.01, mua_lim=c(0, 1), siga_lim=c(0, 1), main=""){
	require(fields)
	mua_l=seq(mua_lim[1],mua_lim[2],pas)
	siga_l=seq(siga_lim[1],siga_lim[2],pas)
	mle=attr(prof2D,"mle")
	cmax=attr(prof2D,"cmax")
	prof2D.na=prof2D
	if (!is.null(seuil_na)){
		prof2D.na[prof2D.na < seuil_na]=NA
	}
	prof_mat=matrix(prof2D.na, nrow=length(mua_l))
	library(fields)
	image.plot(t(prof_mat),y=c(min(mua_l)-demipas,mua_l+demipas),x=c(min(siga_l)-demipas,siga_l+demipas),ylab="MU",xlab="SIG", main=main)
	points(mle[4],mle[2],col="white",lwd=2)
	points(cmax[2],cmax[1],col="black",lwd=2)
}

pdf("Profile2D_v2.pdf")
plot_profile(prof2D_wrong,seuil_na=-3000, mua_lim=c(0, 1), siga_lim=c(0, 1),main="Weird Profile")
plot_profile(prof2D_wrong_zoom1, seuil_na=-3000, mua_lim=c(0, 0.2), siga_lim=c(0, 0.2), pas=0.001, main="WP zoom : Estimated Max")
plot_profile(prof2D_wrong_zoom2, seuil_na=-3000, mua_lim=c(0.8, 1), siga_lim=c(0, 0.2), pas=0.001, main="WP zoom : Actual Max")
plot_profile(prof2D_normal,seuil_na=-3000, mua_lim=c(0, 1), siga_lim=c(0, 1), main="Normal Profile")
plot_profile(prof2D_wrong_init,seuil_na=-3000, mua_lim=c(0, 1), siga_lim=c(0, 1),main="Weird Profile + Init")
dev.off()



# prof2D <- Map(f_roots,mua=mua_l,siga=siga_l)
i=1
