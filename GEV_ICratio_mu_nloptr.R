packages <- c("boot", "quantreg", "evmix", "ismev", "gWidgetstcltk", "gWidgets","evir","extRemes","nloptr")
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
source("FAR_GEV_Algo.R")
source("ImproveRoots.R")


gev.ratio.ic.mu=function(xp,t0,t1,y.fit,ydat,ci.p=0.95,like.num = 1000 ,mulink=identity, siglink = identity, show = TRUE,  method="Nelder-Mead",maxit = 10000,to.plot=FALSE, ...){
	require(evir)
	require(foreach)
	require(doMC)
	if (t0 >= t1)
		stop("t0 must be inferior to t1")
	covariate=ydat$mua
	y.fit=fevd(y,ydat,location.fun=~mua,scale.fun=~siga)
	print(y.fit)
	translation=-covariate[t0]
	t.delta=covariate[t1]+translation
	ydat.t=ydat
	ydat.t$mua=ydat$mua+translation
	#         y.fit=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,optim.arg=list(method=method))
	y.fit=fevd(y,ydat.t,location.fun=~mua,scale.fun=~siga,type="GEV")
	overallmax=-y.fit$results$value
	#         print(y.fit)
	xdat=ydat$y
	mle <- y.fit$results$par
	init=mle[-2]
	mumat=as.matrix(ydat.t[,c("mub","mua")])
	sigmat=as.matrix(ydat.t[,c("sigb","siga")])
	shmat=as.matrix(ydat.t[,"xi"])
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
	browser()
	print(levd(x=xdat,location=mu.v, scale=sig.v, shape=sha.v,type="GEV",npy=1))
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
	gev.lik <- function(par) {
		xi=par[5]
		sig.v <- sigmat %*% par[3:4]
		mu.v <- mumat %*% par[1:2]
		levd(x=xdat,location=mu.v, scale=sig.v, shape=xi,type="GEV",npy=1)
	}
     # Example showing how to solve the problem from the NLopt tutorial.
     #
     # min sqrt( x2 )
     # s.t. x2 >= 0
     #      x2 >= ( a1*x1 + b1 )^3
     #      x2 >= ( a2*x1 + b2 )^3
     # where
     # a1 = 2, b1 = 0, a2 = -1, b2 = 1
     #
     # re-formulate constraints to be of form g(x) <= 0
     #      ( a1*x1 + b1 )^3 - x2 <= 0
     #      ( a2*x1 + b2 )^3 - x2 <= 0
     
     library('nloptr')
     # Solve using NLOPT_LD_MMA with gradient information supplied in separate function
     #      res0 <- nloptr( x0=mle, 
     #                     lb=c(0,-1,-10,-10,-1),ub=c(110,1,110,1,1),
     #                     eval_f=gev.lik,
     #                     eval_g_eq = contraindre(R=0.6,xp=xp),
     #                     opts = list("algorithm"="NLOPT_GN_ISRES",maxeval=10^5)
     #      opts = list("algorithm"="NLOPT_GN_ISRES")
     #                     )
     #      print( res0 )
     library(alabama)
     ratio.l <- sort(c(seq(0,1.2,length.out=like.num),p0p1))[-1]
     #      ratio.l <- p0/p1
     optimiser_profil <- function(R,xp){
	     print(R)
	     tc=try({
		     ans <- auglag(par=mle, fn=gev.lik, heq=contraindre(R=R,xp=xp))
		     print(ans)
		     print(y.fit)
	     },silent=TRUE)
	     if(class(tc)=="try-error" ){
		     #                      print("F*********")
		     return(10^6)
	     }
	     ans
     }
     profil_o=mclapply(ratio.l,optimiser_profil,xp=xp)
     extract_l <- function(x){
	     if (length(x) == 1)
		     res = x
	     else
		     res = x$value
	     res
     }
     parmax=unlist(lapply(profil_o, extract_l))
     parmax=-parmax
     aalpha <- qchisq(ci.p, 1)
     browser()
     if(abs(max(parmax)-aalpha/2)>0.0001){
	     print("non equal mle")
	     print(max(parmax))
	     print(aalpha/2)
     }
     crit <- overallmax- qchisq(0.999, 1)/2
     cond <- parmax > crit
     ratio.l <- ratio.l[cond]
     parmax <- parmax[cond]
     cond <- !is.na(ratio.l) & !is.na(parmax)
     smth <- spline(ratio.l[cond], parmax[cond], n = 500)
     if(to.plot){
	     plot(ratio.l, parmax, type = "l", xlab = "", ylab = "")
	     abline(h = overallmax - aalpha/2, lty = 2, col = 2)
	     lines(smth$x,smth$y, lty = 2, col = 2)
     }
     ci <- smth$x[smth$y > 0]
     out <- c(min(ci), p0p1, max(ci),p0,mu.v[t0],sig.v[t0],sha.v[t0],p1,mu.v[t1],sig.v[t1],sha.v[t1])
     names(out) <- c("LowerCI", "Estimate", "UpperCI","p0","mu0","sc0","sh0","p1","mu1","sc1","sh1")
     print(out)
     print(out[1]<out[2])
     print(out[1]==out[2])
     out
}

