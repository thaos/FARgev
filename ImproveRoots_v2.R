# rm(list=ls())
# x=seq(0,1,length.out=1000)
# fun <- function(x)-100*(x-0.5)^2+2*x+10
# y=fun(x)
# plot(x,y)
# 
# nb.div=20
intervaler <- function(x,w){
	if(length(x) < 2) stop("size of vector inferior to 2")
	if(length(x) != length(w) ) stop("vectors x and w of different size")
	i=1:(length(x)-1)
	res=lapply(i, function(i)c(x[i], x[i+1], w[i], w[i+1]))
	return(res)
}

borner_inf <- function(vect)
	vect[3] <= 0 & vect[4] >= 0
borner_sup <- function(vect)
	vect[3] >= 0 & vect[4] <= 0


binf_time <- function(uninterval,fun,nbdiv,xmax,fmax,nbiter=0,p.res=NULL){
	print("**************************************************")
	print(nbiter)
	print("**************************************************")
	max_iter=20
	z=seq(uninterval[1],uninterval[2],length.out=nbdiv)
	pas=(max(z)-min(z))/nbdiv
	if(!is.null(p.res)){
		z.p=p.res$x
		w.p=p.res$likel
		z=z[!(z%in%z.p)]
	}
	w=unlist(mclapply(z,fun))
	if(nbiter==0){
		iz=which(z == xmax)
		print("checking")
		print(w[iz])
		w[iz]=fmax
	}
	if(!is.null(p.res)){
		z=c(z,z.p)
		w=c(w,w.p)
		w=w[order(z)]
		z=z[order(z)]
	}
	res=list("x"=z,"likel"=w,"pas"=pas)
	if(nbiter >= max_iter)
		return(list("x"=z,"likel"=w,"pas"=pas))
	if(all(w>0)){
		print("all value are positive")
		int_r=uninterval[2]-uninterval[1]	
		extended_interval=c(uninterval[1]-int_r,uninterval[1])
		res=binf_time(extended_interval, fun, nbdiv=nbdiv, xmax=xmax, nbiter=nbiter+1,p.res=res)
		return(res)
	}
	if(all(w<0)){
		stop("all values are negative")
	}
	intervalles=intervaler(z,w)
	bornes_inf=Filter(borner_inf,intervalles)
	if(length(bornes_inf) == 0)
		return(list("x"=z,"likel"=w))
	xinf=min(unlist(lapply(bornes_inf,"[",1)))
	min_binf=which(lapply(bornes_inf,"[",1)==xinf)
	min_binf=bornes_inf[[min_binf]]
	ic_inf=min_binf[1]
	condition <- (xmax-ic_inf)>20*pas
	if(!condition){
		res=binf_time(min_binf[1:2], fun, nbdiv=2*nbdiv, xmax=xmax, nbiter=nbiter+1,p.res=res)
		return(res)
	}else{
		return(list("x"=z,"likel"=w,"pas"=pas))
	}
}

bsup_time <- function(uninterval,fun,nbdiv,xmax,fmax,nbiter=0,p.res=NULL){
	print("**************************************************")
	print(nbiter)
	print("**************************************************")
	max_iter=20
	z=seq(uninterval[1],uninterval[2],length.out=nbdiv)
	pas=(max(z)-min(z))/nbdiv
	if(!is.null(p.res)){
		z.p=p.res$x
		w.p=p.res$likel
		z=z[!(z%in%z.p)]
	}
	w=unlist(mclapply(z,fun))
	if(nbiter==0){
		iz=which(z == xmax)
		print("checking")
		print(w[iz])
		w[iz]=fmax
	}
	if(!is.null(p.res)){
		z=c(z,z.p)
		w=c(w,w.p)
		w=w[order(z)]
		z=z[order(z)]
	}
	res=list("x"=z,"likel"=w,"pas"=pas)
	if(nbiter >= max_iter)
		return(list("x"=z,"likel"=w,"pas"=pas))
	if(all(w>0)){
		print("all value are positive")
		int_r=uninterval[2]-uninterval[1]	
		extended_interval=c(uninterval[2],int_r+uninterval[2])
		res=bsup_time(extended_interval, fun, nbdiv=nbdiv, xmax=xmax, nbiter=nbiter+1,p.res=res)
		return(res)
	}
	if(all(w<0)){
		stop("all values are negative")
	}
	intervalles=intervaler(z,w)
	bornes_sup=Filter(borner_sup,intervalles)
	if(length(bornes_sup) == 0)
		return(list("x"=z,"likel"=w))
	xsup=max(unlist(lapply(bornes_sup,"[",1)))
	max_bsup=which(lapply(bornes_sup,"[",1)==xsup)
	max_bsup=bornes_sup[[max_bsup]]
	ic_sup=max_bsup[2]
	condition <- (ic_sup-xmax)>20*pas
	if(!condition){
		res=bsup_time(max_bsup[1:2], fun, nbdiv=2*nbdiv, xmax=xmax, nbiter=nbiter+1,p.res=res)
		return(res)
	}else{
		return(list("x"=z,"likel"=w,"pas"=pas))
	}
}



get_ic <- function(int_inf,int_sup){
	c("ic_inf"=int_inf[1], "ic_sup"=int_sup[2])
}
	

seek_roots <- function(interval1, interval2=NULL,ic=NULL, fun, nbdiv, iter=1,xmax){
	aggregate_res <- function(x_vrais,y_vrais,res){
			x_vrais = c(x_vrais, res$x_vrais) 
			y_vrais = c(y_vrais, res$y_vrais) 
			#                         x_order=order(x_vrais)
			#                         x_vrais=x_vrais[x_order] 
			#                         y_vrais=y_vrais[x_order]
			return(list("ic"=res$ic,"x_vrais"=x_vrais,"y_vrais"=y_vrais))
	}
	print(iter)
	max.iter=10
	if(iter >= max.iter) {
		print("too many iterations")
		return(list("ic"=ic,"x_vrais"=c(),"y_vrais"=c()))
	}
	cond.lvl=0.01
	if (is.null(interval2)){
		ft=try({
			#                         debug(first_time)
			first_time(interval1,fun,nbdiv,xmax=xmax)
		})
		#                 browser()
		if(class(ft)=="try-error"){
			print("pas trop petit")
			browser()
			#                         return(list("ic"=NA,"x_vrais"=seq(interval1[1],interval1[2],length.out=nbdiv),"y_vrais"=rep(10^6,nbdiv)))
		}
		print(str(ft))
		x_vrais=ft$x
		y_vrais=ft$likel
		print(y_vrais)
		ic = with(ft,get_ic(int_inf,int_sup)) 
		if(is.na(ic[1]))ic[1]=min(x_vrais)
		if(is.na(ic[2]))ic[2]=max(x_vrais)
		pas=ft$pas
		condition=max(abs(ic-xmax)) * cond.lvl >= pas
		#                 condition=(ic[2]-ic[1])*cond.lvl >= pas
		if (condition) return(list("ic"=ic,"x_vrais"=x_vrais,"y_vrais"=y_vrais))
		else {
			interval1 = ft$int_inf[1:2]
			interval2 = ft$int_sup[1:2]
			res = seek_roots(interval1,interval2,ic,fun,nbdiv,iter=iter+1,xmax=xmax)
			return(aggregate_res(x_vrais,y_vrais,res))
		}
	}
	else {
		cat("On continue !")
		x_vrais=c()
		y_vrais=c()
		pas=(ic[2]-ic[1])/nbdiv
		ic1=ic[1]
		ic2=ic[2]
		if(!any(is.na(interval1))){
			i1=binf_time(interval1,fun,nbdiv)
			x_vrais=c(x_vrais,i1$x)
			y_vrais=c(y_vrais,i1$likel)
			ic1 = i1$int_inf[1]
			pas=i1$pas
			interval1 = i1$int_inf[1:2]
		}
		else{
			ic1=ic[1]
		}
		if(!any(is.na(interval2))){
			i2=bsup_time(interval2,fun,nbdiv)
			x_vrais=c(x_vrais,i2$x)
			y_vrais=c(y_vrais,i2$likel)
			ic2 = i2$int_sup[2]
			pas=i2$pas
			interval2 = i2$int_sup[1:2]
		}
		else{
			ic2=ic[2]
		}
		ic=c(ic1,ic2)
		#                 print(i1$pas)
		#                 print(i2$pas)
		#                 if(i1$pas==i2$pas) 
		#                 pas=i1$pas
		#                 condition=(ic[2]-ic[1])*cond.lvl >= pas
		condition=max(abs(ic-xmax)) * cond.lvl >= pas
		print(condition)
		if (condition) return(list("ic"=ic,"x_vrais"=x_vrais,"y_vrais"=y_vrais))
		else {
			res = seek_roots(interval1,interval2,ic,fun,nbdiv,iter=iter+1,xmax=xmax)
			return(aggregate_res(x_vrais,y_vrais,res))
		}
	}
}
	
