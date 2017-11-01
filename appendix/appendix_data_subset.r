library(plyr)
#read in data
datMT <- read.csv("c:\\Users\\hkropp\\Google Drive\\mult_comp\\data\\randomData.csv")
datD <- read.csv("c:\\Users\\hkropp\\Google Drive\\mult_comp\\data\\scenarioData.csv")
datT <- read.csv("c:\\Users\\hkropp\\Google Drive\\mult_comp\\data\\groupComparisionData.csv")

#pull out two scenarios
#scenario with moderate sample size, group size and unbalences
# number 35 and 
#scenario with small group size, small sig, and moderate sample size unbalenced
#scenario 4
simN <- c(3,5,34,37)
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
	scenarioD[[i]] <- datD[datD$scenarioid==simN[i],]
	truthD[[i]] <- datT[datT$scenarioid==simN[i],]
}

#turn back into a data frame

rData <- ldply(ranData, data.frame)
scnData <- ldply(scenarioD, data.frame)
trueData <- ldply(truthD, data.frame)


#write subset
write.table(rData, "c:\\Users\\hkropp\\Google Drive\\mult_comp\\data\\randomDataAppendix.csv",
			sep=",", row.names=FALSE)
write.table(scnData, "c:\\Users\\hkropp\\Google Drive\\mult_comp\\data\\scenarioDataAppendix.csv",
			sep=",", row.names=FALSE)
write.table(trueData, "c:\\Users\\hkropp\\Google Drive\\mult_comp\\data\\groupComparisionDataAppendix.csv",
			sep=",", row.names=FALSE)	