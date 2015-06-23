packages <- c("boot", "quantreg", "evmix", "ismev", "gWidgetstcltk", "gWidgets","evir","extRemes","nloptr","rootSolve")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
	  install.packages(setdiff(packages, rownames(installed.packages())),repos="http://cran.us.r-project.org")
}
library(boot)
library(quantreg)
library(evmix)
library(ismev)
library(evir)
library(extRemes)
require(gWidgetstcltk) #use the gWidgetstcltk package from CRAN
require(gWidgets) #use the gWidgetstcltk package from CRAN
library(nloptr)
library(rootSolve)
source("FAR_GEV_Algo.R")
source("ImproveRoots_v2.R")


gev.ratio.ic.mu=function(xp,t0,t1,y.fit,ydat,ci.p=0.95,like.num = 1000 ,mulink=identity, siglink = identity, show = TRUE,  method="Nelder-Mead",maxit = 10000,to.plot=FALSE, ...){
	require(evir)
	require(foreach)
	require(doMC)
	if (t0 >= t1)
		stop("t0 must be inferior to t1")
	#         covariate=ydat$mua
	mumat=as.matrix(ydat[,c("mub","mua")])
	sigmat=as.matrix(ydat[,c("sigb","siga")])
	shmat=as.matrix(ydat[,"xi"])
	#         y.fit0=fevd(y,ydat,location.fun=~mua,initial=list("mu0"=0,"mu1"=1,"sigma"=1,"shape"=.Machine$double.eps))
	#         print(y.fit0)
	#         sigma0=y.fit$results$par[3]
	#         shape0=y.fit$results$par[4]
	#         y.fit=fevd(y,ydat,location.fun=~mua,scale.fun=~siga)
	tc=try({
		y.fit=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,optim.arg=list(method=method))
		print(y.fit)
	})
	if(class(tc)=="try-error" ){
		print("F********* n°1")
		y.fit=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,initial=list("mu0"=0,"mu1"=1,"sigma0"=0,"sigma1"=1,"shape"=.Machine$double.eps),optim.arg=list(method=method))
		print(y.fit)
	}
	#         covariate=ydat$mua
	#         translation=-covariate[t0]
	#         t.delta=covariate[t1]+translation
	#         ydat.t=ydat
	#         ydat.t$mua=ydat$mua+translation
	#         y.fit=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,optim.arg=list(method=method))
	#         mumat=as.matrix(ydat.t[,c("mub","mua")])
	#         sigmat=as.matrix(ydat.t[,c("sigb","siga")])
	#         shmat=as.matrix(ydat.t[,"xi"])
	#         y.fit=fevd(y,ydat.t,location.fun=~mua,scale.fun=~siga,type="GEV")
	#         overallmax=-y.fit$results$value
	#         print(y.fit)
	xdat=ydat$y
	mle <- y.fit$results$par
	#         mle <- y.fit0$results$par
	#         mle <- c(mle[1:3],0,mle[4])
	print(mle)
	mu.v <- mumat %*% mle[1:2]
	sig.v <- sigmat %*% mle[3:4]
	sha.v <- shmat * mle[5]
	p0=pevd(xp,loc=mu.v[t0],scale=sig.v[t0],shape=sha.v[t0],lower.tail=FALSE)
	p1=pevd(xp,loc=mu.v[t1],scale=sig.v[t1],shape=sha.v[t1],lower.tail=FALSE)
	print(p0)
	print(p1)
	print(p0/p1)
	gev.lik <- function(par) {
		xi=par[5]
		sig.v <- sigmat %*% par[3:4]
		mu.v <- mumat %*% par[1:2]
		levd(x=xdat,location=mu.v, scale=sig.v,
		     shape=xi,type="GEV",npy=1)
	}
	print("FIT2")
	print(gev.lik(mle))
	#         mle=c(0,1,1,0,0)
	tc=try({
		y.fit2=nlminb(start=mle,gev.lik)
		#                 y.fit2=optim(par=mle,gev.lik)
	})
	if(class(tc)=="try-error" ){
		print("F********* n°2")
		browser()
	}
	overallmax=-y.fit2$objective
	mle <- y.fit2$par
	#         overallmax=-y.fit2$value
	#         mle <- y.fit2$par
	print(mle)
	mu.v <- mumat %*% mle[1:2]
	sig.v <- sigmat %*% mle[3:4]
	sha.v <- shmat * mle[5]
	p0=pevd(xp,loc=mu.v[t0],scale=sig.v[t0],shape=sha.v[t0],lower.tail=FALSE)
	p1=pevd(xp,loc=mu.v[t1],scale=sig.v[t1],shape=sha.v[t1],lower.tail=FALSE)
	print(p0/p1)
	if(p0==0 & p1==0)
		p0p1=1
	else
		p0p1=p0/p1
	#         browser()
	#         print(levd(x=xdat,location=mu.v, scale=sig.v, shape=sha.v,type="GEV",npy=1))
	contraindre <- function(R=0.5,xp){
		res=function(par){
			mu.v <- mumat %*% par[1:2]
			sig.v <- sigmat %*% par[3:4]
			xi=par[5]
			tc=try({
				p0=pevd(xp,loc=mu.v[t0],scale=sig.v[t0],shape=xi,lower.tail=FALSE)
				p1=pevd(xp,loc=mu.v[t1],scale=sig.v[t1],shape=xi,lower.tail=FALSE)
			},silent=TRUE)
			if(class(tc)=="try-error" ){
				#                                 print("F*********")
				return(10^6)
			}
			p0/p1-R
		}
		#                 debug(res)
		res
	}
     optimiser_profil <- function(R,xp){
	     require(alabama)
	     tc=try({
		     ans <- auglag(par=mle, fn=gev.lik, heq=contraindre(R=R,xp=xp),control.outer=list(method="nlminb",trace=FALSE))
	     })
	     if(class(tc)=="try-error" ){
		                           print("F*********")
		     return(10^6)
	     }
	     ans
     }
     print(gev.lik(mle))
     tc=try({
	     y.fit3=optimiser_profil(p0p1,xp)
	     overallmax=-y.fit3$value
	     mle <- y.fit3$par
	     print("FIT3")
	     print(overallmax)
	     print(mle)
     })
     if(class(tc)=="try-error" ){
	     print("F********* n°3")
	     browser()
     }
     mu.v <- mumat %*% mle[1:2]
     sig.v <- sigmat %*% mle[3:4]
     sha.v <- shmat * mle[5]
     p0=pevd(xp,loc=mu.v[t0],scale=sig.v[t0],shape=sha.v[t0],lower.tail=FALSE)
     p1=pevd(xp,loc=mu.v[t1],scale=sig.v[t1],shape=sha.v[t1],lower.tail=FALSE)
     print(p0/p1)
     if(p0==0 & p1==0)
	     p0p1=1
     else
	     p0p1=p0/p1
     extract_l<- function(x){
	     if (length(x) == 1)
		     res = x
	     else
		     res = x$value
	     res
     }
     aalpha <- qchisq(ci.p, 1)
     profil.optim <- function(ratio,xp){
	     fit <- optimiser_profil(ratio,xp)
	     extract_l(fit)
     }
     f_roots=function(ratio){
	     #              print("***************************************")
	     print(ratio)
	     parmax = -profil.optim(ratio=ratio,xp=xp)
	     parmax + aalpha/2 - overallmax
     }
     print("BORNE INF-------------------------------------------")
     ic_inf=binf_time(c(-0.1,p0p1),fun=f_roots,nbdiv=2,xmax=p0p1,fmax=aalpha/2)
     print("BORNE SUP-------------------------------------------")
     if(p0p1 <= 1)
	     ic_sup=bsup_time(c(p0p1,1.1),fun=f_roots,nbdiv=2,xmax=p0p1,fmax=aalpha/2)
     else
	     ic_sup=bsup_time(c(p0p1,p0p1*2),fun=f_roots,nbdiv=2,xmax=p0p1,fmax=aalpha/2)
     ratio.l=c(ic_inf$x,ic_sup$x)
     parmax=c(ic_inf$likel,ic_sup$likel)
     summary(parmax)
     if(abs(max(parmax)-aalpha/2)>0.0001){
	     print("non equal mle")
	     print(max(parmax))
	     print(aalpha/2)
     }
     crit <- aalpha/2 - qchisq(0.999, 1)/2
     cond <- parmax > crit
     ratio.l <- ratio.l[cond]
     parmax <- parmax[cond]
     cond <- !is.na(ratio.l) & !is.na(parmax)
     r1cond=ratio.l[cond]
     p1cond=parmax[cond]
     if(sum(cond)<4)
	     browser()
     #      poly3=lm(p1cond~poly(r1cond,degree=3))
     #      bs3=lm(p1cond~bs(r1cond,df=30,degree=3))
     smth <- spline(ratio.l[cond], parmax[cond], n = 500,method="natural")
     to.plot=TRUE
     if(to.plot){
	     plot(ratio.l, parmax, type = "l", xlab = "", ylab = "")
	     abline(h = 0, lty = 2, col = 2)
	     lines(smth$x,smth$y, lty = 2, col = 2)
	     #              lines(ratio.l,predict(poly3),col="blue")
	     #              lines(ratio.l,predict(bs3),col="green")
     }
     #      browser()
     ci <- smth$x[smth$y > 0]
     ci <- ratio.l[parmax > 0]
     if(length(ci)==0)
	     browser()
     out <- c(min(ci), p0p1, max(ci),p0,mu.v[t0],sig.v[t0],sha.v[t0],p1,mu.v[t1],sig.v[t1],sha.v[t1])
     names(out) <- c("LowerCI", "Estimate", "UpperCI","p0","mu0","sc0","sh0","p1","mu1","sc1","sh1")
     print(out)
     print(out[1]<out[2])
     print(out[1]==out[2])
     out
}

