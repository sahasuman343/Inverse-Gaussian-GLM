#IRWLS

#diagonal matrix
diag_m=function(k,m,n){
	t=m*n
	m=matrix(rep(0,t),m,n)
	diag(m)=k
	return(m)
}

#Link function
g=function(x){
	return(1/x)
}
#derivative of link function
g_prime=function(x){
	return(-1/(x^2))
}
#Log-Likelihood
llk=function(y,mu,l){
    t=log(dinvGauss(y,mu,l))#(-lambda*(y-mu)^2)/(2*mu^2*y)+log(lambda/(2*pi*y^3))/2
    return(sum(t))
    }
#let y be the dependent variable and X be the design matrix
phi_hat=function(y,mu,n,k){
	b_prime=(-mu^3)/2
	w=diag_m(1/b_prime,n,n)
	phi=as.numeric(t(y-mu)%*%w%*%(y-mu))
	return(phi/(n-k))
}
#weight matrix
update_w=function(lambda,mu){
	k=lambda*mu
	return(diag_m(k,length(mu),length(mu)))
}
#update z
update_z=function(y,mu){
	temp=diag_m(g_prime(mu),length(mu),length(mu))%*%(y-mu)
	z=g(mu)+temp
	return(z)
}
#update betas
update_b=function(X,w,z){
	t=solve(t(X)%*%w%*%X)%*%t(X)%*%w%*%z
	return(t)
}
#update mu
update_mu=function(b_prev,t){
	r=g(t%*%b_prev)
	return(r)
}
#Check Convergence
check_conv=function(b_prev,b,tol=0.0001){
	temp=sum((b-b_prev)^2)/sum(b_prev^2)
	if(temp<tol){
		return(TRUE)
	}else{
		return(FALSE)
	}
}
#Iterative reweighted Least squares
IWRLS=function(y,X,max_iter,lambda=1){
	i=1
	j=0
	#initial values
	mu=y
	w=update_w(lambda,mu)
	z=update_z(y,mu)
	b=update_b(X,w,z)
	while(i<=max_iter){
		print(paste0(i, " out of ", max_iter," iterations"))
		mu=update_mu(b,X)
		w=update_w(lambda,mu)
		z=update_z(y,mu)
		b_new=update_b(X,w,z)
		if(check_conv(b,b_new)){
			j=j+1
			if(j>5){
				print(paste0("stopped at ",i,"th iteration")) 
				break
			}
		}
		i=i+1	
		b=b_new
	}
	return(b)
}
#calculate lambda
lambda=function(mu,n,k){
	b_prime=(-mu^3)/2
	w=diag_m(1/b_prime,n,n)
	phi=as.numeric(t(y-mu)%*%w%*%(y-mu))
	phi=phi/(n-k)
    return(-2/phi)
    }
#Confidence Interval
confint=function(fit,alpha=0.05){
    X=fit$X
    y=fit$y
    l=fit$lambda
    mu=fit$fitted_mu
    b=fit$beta
    v=sqrt(diag(solve(t(X)%*%update_w(l,mu)%*%X)))
    t=qt(1-alpha/2,df=fit$n-fit$nparams)
    L=b-t*v
    R=b+t*v
    print(v)
    return(cbind(L,R))
}
#Deviance
dev=function(fits){
  llk_0=llk(fits$y,fits$y,fits$lambda)
  return(2*(llk_0-fits$llk))
}
utils::globalVariables("n", add = TRUE)
#FIT
GI_fit=function(y,x,max_iter=100){
  if(length(y)==nrow(x)){
    n<-nrow(x)
    covariates=colnames(x)
    fit=list()
    X=as.matrix(cbind(rep(1,n),x))
    k=ncol(X)
    #print(X)
    fit[["n"]]=n
    fit[["X"]]<-X
    fit[["y"]]<-y
    fit[["covariates"]]<-c("Intercept",covariates)
    b=round(IWRLS(y,X,max_iter),3)
    fit[["beta"]]=b
    mu=update_mu(b,X)
    fit[["fitted_mu"]]=mu
    fit[["lambda"]]=lambda(mu,n,k)
    fit[["family"]]="Inverse Gaussian"
    fit[["link"]]="inverse"
    fit[["llk"]]=llk(y,mu,fit[["lambda"]])
    fit[["nparams"]]=ncol(X)
    fit[["w"]]=update_w(fit[["lambda"]],mu)
    return(fit)
  }
  else
    print("Number of observations are different!!")
}

#Testing of Hypothesis
#test=function(fit,score="wald"){
	#if(fit$family=="Inverse Gaussian"){
	#	df=matrix(ncol=4,nrow = nrow(fit$nparams))
	#	if(score="lrt"){
	#	}
 #   else if (score="wald"){
      
#    }

#		}}


#Hat matrix
hat=function(fit){sqrt(fit$w)%*%fit$X%*%solve(t(fit$X)%*%fit$w%*%fit$X)%*%t(fit$X)%*%sqrt(fit$w)}

#Residual

res=function(fit,type="deviance"){
    	H=diag(hat(fit))
	if(type=="pearson"){
		r=(fit$y-fit$fitted_mu)/as.vector(sqrt(var(fit$fitted_mu)%*%(1-H)))

		return(r)
	}
	else{
		d=2*(log(dinvGauss(fit$y,fit$y,fit$lambda))-log(dinvGauss(fit$y,fit$fitted_mu,fit$lambda)))
		e=(sign(fit$y-fit$fitted_mu)*sqrt(d))/sqrt(1-H)
		return(e)

	}
	}
