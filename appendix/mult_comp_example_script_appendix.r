#######################################################################
####### Script to run example multiple comparisions simulation ########
########################################################################

#load jags libraries
library(coda)
library(rjags)
library(plyr)
library(mcmcplots)



#################################
# Users specified data directory#
#################################

#path to save output to
wd.out <- "c:\\Users\\hkropp\\Google Drive\\mult_comp\\sub_out"
#path for model code
wd.code <- "c:\\Users\\hkropp\\Documents\\GitHub\\multiple_comp_ex"
#path to read in data
wd.data <- "c:\\Users\\hkropp\\Google Drive\\mult_comp\\data"


#################################
# read in data files            #
#################################
#read in all data
datMT <- read.csv(paste0(wd.data,"\\randomDataAppendix.csv"))
scnData <- read.csv(paste0(wd.data,"\\scenarioDataAppendix.csv"))
groupCompD <-  read.csv(paste0(wd.data,"\\groupComparisionDataAppendix.csv"))


#################################
# pull out data for each run    #
#################################
simN <- scnData$scenarioid
#get length of sims
numsims<-length(simN)
#subset the simulations
#set up in a list format so the number
#of simulations can be readily changed

ranD <- list()
groupD <- list()
for(i in 1:numsims){
	ranD[[i]] <- datMT[datMT$scenarioid==simN[i],]
	groupD[[i]] <- groupCompD[groupCompD$scenarioid==simN[i],]
}




#################################
# Set up model run ##############
#################################
#run the model for each scenario

#set up model runs for each scenario to test
for(i in 1:numsims){
	datalist <- list(Nobs=dim(ranD[[i]])[1],Y=ranD[[i]]$Y, groupID=ranD[[i]]$groupid,
						Y.nh=ranD[[i]]$Y.nh,Y.cp=ranD[[i]]$Y.cp,
						Ngroup=scnData$grps[i], Ncomp=dim(groupD[[i]])[1],
						g.c1=groupD[[i]]$groupid1, g.c2=groupD[[i]]$groupid2)



	#initialize model
	s.mod.init=jags.model(file=paste0(wd.code, "\\mult_comp_example_code.r" ),
					  data=datalist,
					  n.adapt=2000,
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

	chain1<-as.matrix(s.mod.coda[[1]])
	write.table(chain1,paste0(wd.out,"\\simulation",i,"chain1_coda.csv"), sep=",")
	chain2<-as.matrix(s.mod.coda[[2]])
	write.table(chain2,paste0(wd.out,"\\simulation",i,"chain2_coda.csv"), sep=",")
	chain3<-as.matrix(s.mod.coda[[3]])
	write.table(chain3,paste0(wd.out,"\\simulation",i,"chain3_coda.csv"), sep=",")
			

#run mcmc plots on key params
			
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

#PP is hierarchical model
# CP is the complete pooling model,
# and NH is the non-hierarchical model 
Pind<-function(HH,CP,NH){
	(NH-HH)/(NH-CP)
}


#################################
# partial pooling with WAIC######
#################################

#read in coda
CODA1 <- list()
CODA2 <- list()
CODA3 <- list()
CODAALL <- list()
for(i in 1:numsims){
	CODA1[[i]] <- read.csv(paste0(wd.out,"\\simulation",i,"chain1_coda.csv"))
	CODA2[[i]] <- read.csv(paste0(wd.out,"\\simulation",i,"chain2_coda.csv"))
	CODA3[[i]] <- read.csv(paste0(wd.out,"\\simulation",i,"chain3_coda.csv"))
	#join all iterations together
	CODAALL[[i]] <- rbind(CODA1[[i]],CODA2[[i]],CODA3[[i]])
}
#subset to pull mean and standard deviation out of columns
dexps3 <- "\\.*[[:digit:]]*\\."

mu.hh <- list()
mu.nh <- list()
mu.cp <- list()
sig.hh <- list()
sig.nh <- list()
sig.cp <- list()

for(i in 1:numsims){
	mu.hh[[i]]<-CODAALL[[i]][,gsub(dexps3,"",colnames(CODAALL[[i]]))=="mu"]
	mu.nh[[i]]<-CODAALL[[i]][,gsub(dexps3,"",colnames(CODAALL[[i]]))=="munh"]
	mu.cp[[i]]<-CODAALL[[i]][,gsub(dexps3,"",colnames(CODAALL[[i]]))=="mucp"]
	sig.hh[[i]]<-CODAALL[[i]][,gsub(dexps3,"",colnames(CODAALL[[i]]))=="sig"]
	sig.nh[[i]]<-CODAALL[[i]][,gsub(dexps3,"",colnames(CODAALL[[i]]))=="signh"]
	sig.cp[[i]]<-CODAALL[[i]][,gsub(dexps3,"",colnames(CODAALL[[i]]))=="sigcp"]	
	
}

#calculate log posterior density for each data point and iteration
#each simulation has a matrix that is simulationsXdatapoints
lpd.nh <- matrix()
Lpd.nh <- list()
lpd.cp <- matrix()
Lpd.cp <- list()
lpd.hh <- matrix()
Lpd.hh <- list()


for(i in 1:numsims){
	#set up empty matrix
	lpd.nh <- matrix(rep(NA,dim(ranD[[i]])[1]*dim(CODAALL[[i]])[1]), ncol=dim(ranD[[i]])[1])
	lpd.hh <- matrix(rep(NA,dim(ranD[[i]])[1]*dim(CODAALL[[i]])[1]), ncol=dim(ranD[[i]])[1])
	lpd.cp <- matrix(rep(NA,dim(ranD[[i]])[1]*dim(CODAALL[[i]])[1]), ncol=dim(ranD[[i]])[1])
	for(j in 1:dim(ranD[[i]])[1]){
		lpd.nh[,j]<- log(dnorm(ranD[[i]]$Y[j],mu.nh[[i]][,ranD[[i]]$groupid[j]],sig.nh[[i]] ))
		lpd.hh[,j]<- log(dnorm(ranD[[i]]$Y[j],mu.hh[[i]][,ranD[[i]]$groupid[j]],sig.hh[[i]] ))
		lpd.cp[,j]<- log(dnorm(ranD[[i]]$Y[j],mu.cp[[i]],sig.cp[[i]] ))
	}
	
	Lpd.nh[[i]] <- lpd.nh
	Lpd.hh[[i]] <- lpd.hh
	Lpd.cp[[i]] <- lpd.cp
}


#find the variation of lpd across iterations for a given datapoint/parameter
var.lpd.hh <- list()
var.lpd.nh <- list()
var.lpd.cp <- list()

for(i in 1:numsims){
	var.lpd.hh[[i]] <- apply(Lpd.hh[[i]],2,"var")
	var.lpd.nh[[i]] <- apply(Lpd.nh[[i]],2,"var")
	var.lpd.cp[[i]] <- apply(Lpd.cp[[i]],2,"var")
}

#final WAIC calculation
#calculate the sum of the variation across all observations now
pWAIC.nh <- numeric(0)
pWAIC.hh <- numeric(0)
pWAIC.cp <- numeric(0)
for(i in 1:numsims){
	pWAIC.nh[i] <- sum(var.lpd.nh[[i]])
	pWAIC.hh[i] <- sum(var.lpd.hh[[i]])
	pWAIC.cp[i] <- sum(var.lpd.cp[[i]])
}

#add PPI to scenario data
scnData$PPI.waic <-  Pind(pWAIC.hh,pWAIC.cp,pWAIC.nh)

#look at scenarios data
scnData
#save scenarios outcomes
write.table(scnData, paste0(wd.out,"\\PPIforScenarios.csv"), row.names=FALSE,
			sep=",")