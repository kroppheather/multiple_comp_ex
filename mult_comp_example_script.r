#######################################################################
####### Script to run example multiple comparisions simulation ########
########################################################################

#load jags libraries
library(coda)
library(rjags)
library(xtable)
library(plyr)


#set the random generator seed used in example
set.seed(14101066) 

#

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
#need to get the start end end grow 



#set up model runs for each scenario to test
for(i in 1:numsims){
datalist <- list(Nobs=dim(ynorm[[i]])[1], Y=y.norm[[i]]$y.dat, groupID=y.norm[[i]]$grID, Ngroup=grps[i], 


# runtime parameters. Specify total iterations and how many runs you want to split it 
# into for RAM considerations. More splits = more final files but less max RAM use
n.adapt=1000
iter.tot=3000
#initialize model
ptm<-proc.time()
poodel=jags.model("source/jagsmodel indexing.R",
                  data=datalist,
                 # inits=inits,
                  n.adapt=n.adapt,
                  n.chains=3)
#burnin
update(poodel, n.iter=1000)

#ptm<-proc.time() #store proc time

#iters=iter.tot/splits #calculate number of iterations per split
#sampling

  coda=coda.samples(poodel,variable.names=c("mu","m","sig.s","sig","mu.nh","sig.nh",
                                          "mu.cp","sig.cp","sumLL.cp","sumLL.nh",
                                          "sumLL","diff","diff.nh"),
                  n.iter=30000, thin=10)
  save(coda,file=paste0("MC_jagsout_sim_", sim, ".R"))
  rm(coda) #dump and clean memory
  gc()
#}

print(paste0("Run time: ", proc.time() - ptm )) #run time
#  =5 hrs for 3000 iterations.


print(paste0("Model ", sim, " done."))

} #end of function

args<-commandArgs(TRUE)
#step1 <- unlist(strsplit(x=args, split="_", perl=TRUE))
#step1 <- unlist(strsplit(x=step1[3], split=".r", perl=TRUE))
print(as.numeric(args))
print(args)



