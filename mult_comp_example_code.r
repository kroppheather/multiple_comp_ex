model{
  for(i in 1:Nobs){
	#likelihood for the data for all gropus 
    Y[i]~dnorm(mu[groupID[i]],tau)
    loglike[i]<-.5*log(tau/(2*3.141593))-((tau/2)*pow(Y[i]-mu[groupID[i]],2))
  }
  #assign noninformative prior to group means and variance
  for(i in 1:Ngroup){
      mu[i]~dnorm(m,tau.s)

    }

    #assign noninformative prior for mean
    m~dnorm(0,0.0001)
	#non informative variance priors
    tau<-pow(sig,-2)
    sig[i]~dunif(0,100)
	#sum log likelihood
    sumLL[i]<-sum(loglike[1:Nobs])
  
  
  
  #non-hierichical 
  for(i in 1:Nobs){
    Y.nh[i]~dnorm(mu.nh[groupID[i]],tau.nh)
    loglike.nh[i]<-.5*log(tau.nh/(2*3.141593))-((tau.nh/2)*pow(Y.nh[i]-mu.nh[groupID[i]],2))
  }
  #assign noninformative prior
  for(i in 1:Ngroup){
      mu.nh[i]~dnorm(0,0.0001)

    }
	
    tau.nh<-pow(sig.nh,-2)
    sig.nh~dunif(0,100)
    sumLL.nh<-sum(loglike.nh[1:Nobs])
  
  #calculate pairwise difference
  for(i in 1:Ncomp){
    diff[i]<-mu[trt1.all[i]]-mu[trt2.all[i]]
    diff.nh[i]<-mu.nh[trt1.all[i]]-mu.nh[trt2.all[i]]
  }	
  #complete pooling 
  for(i in 1:Nobs){
    Y.cp[i]~dnorm(mucp[groupID[i]],tau.cp)
    loglike.cp[i]<-.5*log(tau.cp/(2*3.141593))-((tau.cp/2)*pow(Y.cp[i]-mu.cp[groupID[i]],2))
  }
  #assign noninformative prior

    mu.cp~dnorm(0,.0001)
    tau.cp<-pow(sig.cp,-2)
    sig.cp~dunif(0,100)
    sumLL.cp<-sum(loglike.cp[1:Nobs])
  
  
}
#list(Nscenario=252)
#Nobs=vector of observations for 252
#Ntreat=vector of number of treatments for each 252 scenario
