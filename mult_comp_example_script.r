#######################################################################
####### Script to run example multiple comparisions simulation ########
########################################################################

#load jags libraries
library(coda)
library(rjags)
library(xtable)
library(plyr)
library(mcmcplots)


#set the random generator seed used in example
set.seed(14101066) 

#working directory for output
wd.out <- "c:\\Users\\hkropp\\Google Drive\\mult_comp"
wd.code <- "c:\\Users\\hkropp\\Documents\\GitHub\\multiple_comp_ex"

###############################################################################
# Run two example scenarios from the manuscript
#note the number of scenarios and the values in them can all be altered
# all information below needs to be filled out for each scenario
###############################################################################
#number of scenarios to look at:
numsims<-2 # replicate in manuscript (r=1,2,...,100)
N <- c(1000,1000) # observations Nk (Fig. 1A)
grps <- c(20, 20) # group sizes K (Fig. 1B)
h.sig <- c(1,1) # sigma global (Fig. 1C) sigma = 0.1 and 10 found in code below
g.sig <- c(1,1)
h.mu <- c(10,10) # mean global (Fig. 1D, m=0,10,100) set manually
A <- c(4, 4) # Number of groups K to have mu k equal to other groups
num.unbal <- c(0,4) # Number of unbalanced groups
#percent missing in unbalanced groups
per.missing1 <- c(0,0.30)

#################################################################################
# No need for user specification below this point  ##############################
#################################################################################

#generate data for runs
#generate group means. Putting in list format in case group or sample sizes were changed to be different sizes
#initialize lists
g.mu1 <- list()
g.mu2 <- list()
g.mu <- list()
A.same <- list()
y.temp<-list()
y.norm<-list()
groupBal <- list()

#loop through the number of simulations
for(i in 1:numsims){
	#generate means for the group data without the same means
	g.mu1[[i]] <- data.frame(grID=seq(1,grps[i]-A[i]),g.meanID=seq(1,grps[i]-A[i]),g.mu=rnorm(n=grps[i]-A[i], mean=h.mu[i], sd=h.sig[i]))
	#now set the means to be the same for some groups
	#first randomly select the means that will be the same
	#also selecting with replacement so more than two groups can share the same mean
	A.same[[i]] <- sample(seq(1,grps[i]-A[i]), A[i], replace=TRUE )
	g.mu2[[i]] <- data.frame(grID=seq(grps[i]-A[i]+1,grps[i]),g.meanID=g.mu1[[i]]$g.meanID[A.same[[i]]],g.mu=g.mu1[[i]]$g.mu[A.same[[i]]])
	#now combine both together
	g.mu[[i]] <- rbind(g.mu1[[i]],g.mu2[[i]])
		#generate data for  scenarios in groups
		#first need to account for scenarios with unbalanced groups
		groupBal[[i]] <- c(rep(1,num.unbal[i]),rep(0, grps[i]-num.unbal[i]))
		
		#generate data for all other groups
		for(j in 1:grps[i]){
		y.temp[[j]]<-data.frame(grID=rep(j,N[i]*((groupBal[[i]][j]*(1-per.missing1[i]))+(1-groupBal[[i]][j]))),
						y.dat=rnorm(N[i]*((groupBal[[i]][j]*(1-per.missing1[i]))+(1-groupBal[[i]][j])),g.mu[[i]]$g.mu[j],sd=g.sig[i]))
	
	}
	y.norm[[i]] <- ldply(y.temp, data.frame)

}
#need to set up comparision ID
comball <- list()
combTable <- list()
for(i in 1:numsims){
	comball[[i]]<- combn(grps[i],2)
	combTable[[i]] <- data.frame(g.c1=comball[[i]][1,], g.c2=comball[[i]][2,], gc.id=seq(1,dim(comball[[i]])[2]))
}

#set up model runs for each scenario to test
for(i in 1:numsims){
	datalist <- list(Nobs=dim(y.norm[[i]])[1], Y=y.norm[[i]]$y.dat, groupID=y.norm[[i]]$grID, Ngroup=grps[i], 
						Y.nh=y.norm[[i]]$y.dat, Y.cp=y.norm[[i]]$y.dat, Ncomp=dim(combTable[[i]])[1], g.c1=combTable[[i]]$g.c1,
						g.c2=combTable[[i]]$g.c2)



	n.adapt=2000
	iter.tot=3000
	#initialize model
	s.mod.init=jags.model(paste0(wd.code, "\\mult_comp_example_code.r" ),
					  data=datalist,
					 # inits=inits,
					  n.adapt=n.adapt,
					  n.chains=3)

	 s.mod.coda=coda.samples(s.mod.init,variable.names=c("mu","m","sig.s","sig","mu.nh","sig.nh",
											  "mu.cp","sig.cp","sumLL.cp","sumLL.nh",
											  "sumLL","diff","diff.nh"),
					  n.iter=30000, thin=10)
	 


	#pull out model stats
	Mod.out<-summary(s.mod.coda)
	dir.create(paste0(wd.out,"\\simulation",i))

	write.table(Mod.out$statistics, paste0(wd.out,"\\simulation",i,"\\mod_stats.csv"),
			sep=",",row.names=TRUE)
	write.table(Mod.out$quantiles, paste0(wd.out,"\\simulation",i,"\\mod_quant.csv"),
			sep=",",row.names=TRUE)


	#save coda

	chain1<-as.matrix(codaobj.init[[1]])
	write.table(chain1,paste0(wd.out,"\\simulation",i,"chain1_coda.csv"), sep=",")
	chain2<-as.matrix(codaobj.init[[2]])
	write.table(chain2,paste0(wd.out,"\\simulation",i,"chain2_coda.csv"), sep=",")
	chain3<-as.matrix(codaobj.init[[3]])
	write.table(chain3,paste0(wd.out,"\\simulation",i,"chain3_coda.csv"), sep=",")
			

#run mcmc plots on key params
			
mcmcplot(codaobj.init, parms=c("mu","m","sig.s","sig","mu.nh","sig.nh",
											  "mu.cp","sig.cp","sumLL.cp","sumLL.nh",
											  "sumLL","diff","diff.nh"),
			dir=paste0(wd.out,"\\simulation",i))



} #end of simulation models



