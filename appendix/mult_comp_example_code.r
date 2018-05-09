########################################################################
####### Code to run example multiple comparisions simulation model######
########################################################################


model{
  # hierarhical model (focal model)
  for(i in 1:Nobs){
	#likelihood for the data for all group levels
    Y[i]~dnorm(mu[groupID[i]],tau)
	#log likelihood for each observation
    loglike[i]<-.5*log(tau/(2*3.141593))-((tau/2)*pow(Y[i]-mu[groupID[i]],2))
  }
  #assign noninformative hierarchical prior to the group-level means 
  for(i in 1:Ngroup){
      mu[i]~dnorm(m,tau.s)

    }

    #assign noninformative prior for mean to global mean
    m~dnorm(0,0.0001)
	#non-informative variance priors to the standard deviation terms
    tau<-pow(sig,-2)
    sig~dunif(0,100)
	tau.s<-pow(sig.s,-2)
    sig.s~dunif(0,100)

  
  
  
  #non-hierarchical version
  for(i in 1:Nobs){
    Y.nh[i]~dnorm(mu.nh[groupID[i]],tau.nh)
	#log likelihood
    loglike.nh[i]<-.5*log(tau.nh/(2*3.141593))-((tau.nh/2)*pow(Y.nh[i]-mu.nh[groupID[i]],2))
  }
  #assign independent noninformative prior to group-level means
  for(i in 1:Ngroup){
      mu.nh[i]~dnorm(0,0.0001)

    }
	#assign non-informative prior to standard deviation:	
    tau.nh<-pow(sig.nh,-2)
    sig.nh~dunif(0,100)

  
  #calculate pairwise difference
  #between group-level means Ncomp is the number of pairwise
	# comparisons, and g.c1 and g.c2 are indexing variables that denote the
	# group-level id’s associate with the first and second mean in the comparison

  for(i in 1:Ncomp){
	#hierarchical pairwise difference
    diff[i]<-mu[g.c1[i]]-mu[g.c2[i]]
	#non-hierarchical pairwise difference (likely not relevant or reported)
    diff.nh[i]<-mu.nh[g.c1[i]]-mu.nh[g.c2[i]]
  }	
  #complete pooling 
  for(i in 1:Nobs){
    Y.cp[i]~dnorm(mu.cp,tau.cp)
	#log likelihood
    loglike.cp[i]<-.5*log(tau.cp/(2*3.141593))-((tau.cp/2)*pow(Y.cp[i]-mu.cp,2))
  }
  #assign noninformative priors to global mean and standard deviation:

    mu.cp~dnorm(0,.0001)
    tau.cp<-pow(sig.cp,-2)
    sig.cp~dunif(0,100)

  
}

