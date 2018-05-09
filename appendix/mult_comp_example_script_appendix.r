########################################################################
####### Script to implemen Bayesian models that analyze          ####### 
####### pseudo data and computing the associated PPI             #######
########################################################################

#load jags libraries
library(coda)
library(rjags)
library(plyr)
library(mcmcplots)



#################################
# Users specified data directory#
#################################
#SPECIFY LOCAL DIRECTORY WHERE YOUR GITHUB REPOSITORY IS 
local.path <- #"PATH HERE"
#SPECIFY LOCAL DIRECTORY WHERE YOU WOULD LIKE TO SAVE MODEL OUTPUT
wd.out <- #"PATH HERE"

#the file paths for the script and data within the github repository download

#path for model code
wd.code <- paste0(local.path,"\\appendix"
#path to read in data
wd.data <- paste0(local.path,"\\appendix_data")


#################################
# read in data files            #
#################################
#read in all data
#pseudo data for a single simulation replicate, for a subset of scenarios
datMT <- read.csv(paste0(wd.data,"\\randomDataAppendix.csv"))
#description of each scenario
scnData <- read.csv(paste0(wd.data,"\\scenarioDataAppendix.csv"))
#true means for pairwise comparisons among groups within each scenario
groupCompD <-  read.csv(paste0(wd.data,"\\groupComparisionDataAppendix.csv"))


#################################
# organize psuedo data for      #
# example scenarios             #
#################################

#this appendix focuses on a subset of 4 pseudo datasets with:

#scenario identifier for example pseudo dataset
simN <- scnData$scenarioid

#number of simulations that we will look at (here its 4)
Ndatasets<-length(simN)

#subset the data for simulation so that each simulation
#can be run separately
ranD <- list()
groupD <- list()
for(i in 1:Ndatasets){
	ranD[[i]] <- datMT[datMT$scenarioid==simN[i],]
	groupD[[i]] <- groupCompD[groupCompD$scenarioid==simN[i],]
}


#################################
# Set up Bayesian model run #####
#################################

#run the Bayesian models for each scenario. Implement 3 model variants:
#hierarchical, non-hierarchical, and complete pooling models
#are all run in the same model code for each simulation

#set up model runs for each dataset to test
for(i in 1:Ndatasets){
	datalist <- list(Nobs=dim(ranD[[i]])[1],Y=ranD[[i]]$Y, groupID=ranD[[i]]$groupid,
						Y.nh=ranD[[i]]$Y.nh,Y.cp=ranD[[i]]$Y.cp,
						Ngroup=scnData$grps[i], Ncomp=dim(groupD[[i]])[1],
						g.c1=groupD[[i]]$groupid1, g.c2=groupD[[i]]$groupid2)



	#initialize Bayesian (JAGS) models
	s.mod.init=jags.model(file=paste0(wd.code, "\\mult_comp_example_code.r" ),
					  data=datalist,
					  n.adapt=2000,
					  n.chains=3)
	#update the initialized Bayesian models and monitor quantities (parameters) of interest via JAGS
	 s.mod.coda=coda.samples(s.mod.init,variable.names=c("mu","m","sig.s","sig","mu.nh","sig.nh",
											  "mu.cp","sig.cp","sumLL.cp","sumLL.nh",
											  "sumLL","diff","diff.nh", "loglike","loglike.nh","loglike.cp"),
					  n.iter=30000, thin=10)
	 
	#save Bayesian model run information	

	# extract posterior summary statistics of interest.
	Mod.out<-summary(s.mod.coda)
	dir.create(paste0(wd.out,"\\simulation",i))

	write.table(Mod.out$statistics, paste0(wd.out,"\\simulation",i,"\\mod_stats.csv"),
			sep=",",row.names=TRUE)
	write.table(Mod.out$quantiles, paste0(wd.out,"\\simulation",i,"\\mod_quant.csv"),
			sep=",",row.names=TRUE)


	#save coda

	chain1<-as.matrix(s.mod.coda[[1]])
	write.table(chain1,paste0(wd.out,"\\simulation",i,"chain1_coda.csv"), sep=",")
	chain2<-as.matrix(s.mod.coda[[2]])
	write.table(chain2,paste0(wd.out,"\\simulation",i,"chain2_coda.csv"), sep=",")
	chain3<-as.matrix(s.mod.coda[[3]])
	write.table(chain3,paste0(wd.out,"\\simulation",i,"chain3_coda.csv"), sep=",")
			

	#create mcmc plots for key params
			
mcmcplot(s.mod.coda, parms=c("mu","m","sig.s","sig","mu.nh","sig.nh",
											  "mu.cp","sig.cp","sumLL.cp","sumLL.nh",
											  "sumLL","diff","diff.nh"),
			dir=paste0(wd.out,"\\simulation",i))



} #end of simulation models

#################################
# Calculate partial pooling######
#################################
#partial pooling indices 
#for degree of pooling
#define function to calculate PPI
#HB is the hierarchical model
#CP is the complete pooling model,
#NH is the non-hierarchical model 
Pind<-function(HB,CP,NH){
	(NH-HB)/(NH-CP)
}


#################################
# partial pooling index based on#
# with WAIC effective number of #
# parameters                    #
#################################

#read in coda from model runs
#for each of the four scenario examples
#so that WAIC can be calculated
CODA1 <- list()
CODA2 <- list()
CODA3 <- list()
CODAALL <- list()
for(i in 1:Ndatasets){
	CODA1[[i]] <- read.csv(paste0(wd.out,"\\simulation",i,"chain1_coda.csv"))
	CODA2[[i]] <- read.csv(paste0(wd.out,"\\simulation",i,"chain2_coda.csv"))
	CODA3[[i]] <- read.csv(paste0(wd.out,"\\simulation",i,"chain3_coda.csv"))
	#join all iterations together
	CODAALL[[i]] <- rbind(CODA1[[i]],CODA2[[i]],CODA3[[i]])
}
#subset to pull log likelihood out of columns 
#(note: the log-likelihood of each observations is computed inside of the JAGS model code, for each model variant and each dataset):
dexps3 <- "\\.*[[:digit:]]*\\."

loglike.hh <- list()
loglike.nh <- list()
loglike.cp <- list()


for(i in 1:Ndatasets){
	loglike.hh[[i]]<-CODAALL[[i]][,gsub(dexps3,"",colnames(CODAALL[[i]]))=="loglike"]
	loglike.nh[[i]]<-CODAALL[[i]][,gsub(dexps3,"",colnames(CODAALL[[i]]))=="loglikenh"]
	loglike.cp[[i]]<-CODAALL[[i]][,gsub(dexps3,"",colnames(CODAALL[[i]]))=="loglikecp"]

	
}


#calculate the posterior variatiances of the log likelihood,
# but summing the above posterior variances across all observations within each dataseta
#for each observation across iterations of the model for each dataset
var.ll.hh <- list()
var.ll.nh <- list()
var.ll.cp <- list()

for(i in 1:Ndatasets){
	var.ll.hh[[i]] <- apply(loglike.hh[[i]],2,"var")
	var.ll.nh[[i]] <- apply(loglike.nh[[i]],2,"var")
	var.ll.cp[[i]] <- apply(loglike.cp[[i]],2,"var")
}

#final WAIC calculation
#calculate the sum of the variation across all scenarios
pWAIC.nh <- numeric(0)
pWAIC.hh <- numeric(0)
pWAIC.cp <- numeric(0)
for(i in 1:Ndatasets){
	pWAIC.nh[i] <- sum(var.ll.nh[[i]])
	pWAIC.hh[i] <- sum(var.ll.hh[[i]])
	pWAIC.cp[i] <- sum(var.ll.cp[[i]])
}

#add PPI calculation to scenario information for comparison 
scnData$PPI.waic <-  Pind(pWAIC.hh,pWAIC.cp,pWAIC.nh)

#look at scenarios data
scnData
#save scenarios outcomes
write.table(scnData, paste0(wd.out,"\\PPIforScenarios.csv"), row.names=FALSE,
			sep=",")