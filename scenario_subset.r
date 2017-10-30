#######################################################################
####### Subsets data generated for scenarios in simulations    ########
#######################################################################
library(plyr)
#read in data
datMT <- read.csv("c:\\Users\\hkropp\\Google Drive\\mult_comp\\data\\megatable_1.csv")
datD <- read.csv("c:\\Users\\hkropp\\Google Drive\\mult_comp\\data\\designall.csv")
datT <- read.csv("c:\\Users\\hkropp\\Google Drive\\mult_comp\\data\\allcomp_1.csv")

#pull out two scenarios
#scenario with moderate sample size, group size and unbalences
# number 35 and 
#scenario with small group size, small sig, and moderate sample size unbalenced
#scenario 4
simN <- c(10,36)
#get length of sims
numsims<-length(simN)

#subset the simulations
#set up in a list format so the number
#of simulations can be readily changed

ranData <- list()
scenarioD <- list()
truthD <- list()

for(i in 1:numsims){
	ranData[[i]] <- datMT[datMT$scenarioid==simN[i],]
	scenarioD[[i]] <- datD[datD$scn==simN[i],]
	truthD[[i]] <- datT[datT$scn==simN[i],]
}

#turn back into a data frame

rData <- ldply(ranData, data.frame)
scnData <- ldply(scenarioD, data.frame)
trueData <- ldply(truthD, data.frame)

#write subset
write.table(rData, "c:\\Users\\hkropp\\Google Drive\\mult_comp\\data\\rDatasub.csv",
			sep=",", row.names=FALSE)
write.table(scnData, "c:\\Users\\hkropp\\Google Drive\\mult_comp\\data\\scnDatasub.csv",
			sep=",", row.names=FALSE)
write.table(trueData, "c:\\Users\\hkropp\\Google Drive\\mult_comp\\data\\trueDatasub.csv",
			sep=",", row.names=FALSE)			
			