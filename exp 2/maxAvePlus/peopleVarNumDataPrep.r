library(chMorals)
library(chutils)
library(stringr)

cleanOverlapFile <- FALSE

#set up the new RT variables
fitCol <- "fit.RT"
resCol <- "res.RT"
useTwoParameterModel <- FALSE
overlapDataIsComplete <- TRUE

# read in parameters
params<-ch.readMoralsDBfile("moralsDBfile.txt")

#set up the group and item directories
mainDir <- getwd()
ch.newDir (mainDir, params$gpSubDir)
gpDir <- getwd()
setwd(mainDir)

ch.newDir (mainDir, params$itemSubDir)
itemDir <- getwd()
setwd(mainDir)


statsOutputFile <- file.path(mainDir,paste(params$dt.set, params$statsOutputFilePrefix))
sink(statsOutputFile, append = F)
  cat("\n***** New Run ****\n\n")
sink(NULL)


### read in data
data.raw <-read.table(params$moralsTaskDataFile, header=T, sep="\t")
data.ovrlp <-read.table(params$valueOverlapDataFile, header=T, sep="\t", quote="\"")

######_____REMOVE PRACTICE TRIALS _____######
data.raw <- data.raw[data.raw$trial_type >=1, ]

### do Prep analysis
processedData <- ch.moralsDataPrep(data.raw, "sn", "keybRT", "overlap", "direction", "trial", "keyDef", respChoiceVal = c("Yes", "No"), item1cols = c("IA.1", "IA.2"), item2cols = c("IB.1", "IB.2"), overlapItem1cols = c("IA1", "IA2"), overlapItem2cols = c("IB1", "IB2"), statsOutputFile = statsOutputFile, params = params, overlapDataIsComplete = overlapDataIsComplete)

### get HVO quantity
processedData$HVOq <- ifelse((processedData$QuantOption1 == 2 & processedData$direct.xVy == 1) | (processedData$QuantOption2 == 2 & processedData$direct.xVy == -1) , 2, 1)
### get LVO quantity
processedData$LVOq <- ifelse((processedData$QuantOption1 == 2 & processedData$direct.xVy == -1) | (processedData$QuantOption2 == 2 & processedData$direct.xVy == 1) , 2, 1)

### Filter data
analysisReadyData <- ch.moralsFilterDataQ(processedData, "sn", "keybRT", "overlapRound", "correct",c(1,0), statsOutputFile = statsOutputFile, params = params)

### Do RT and p(Hit Analysis on Group Data - remove learning effects for the group)
analysisReadyData.gp <- ch.moralsGrpRTpHit(analysisReadyData, "trial", "keybRT", fitCol, resCol, "overlapRound", "keyDef",c("Yes", "No"), "correct",c(1,0), useTwoParameterModel = useTwoParameterModel, params = params)
write.table(analysisReadyData.gp, file="analysisReadyData.gp.txt", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

### Do RT and p(Hit Analysis on individual subject Data - remove learning effects for each subject)
analysisReadyData.sn <- ch.moralsSnRTpHit(analysisReadyData, "sn", "trial", "keybRT", fitCol, resCol, "overlap", "correct", c(1,0),  useTwoParameterModel = useTwoParameterModel, params = params, minUniqueOverlaps = 4)
write.table(analysisReadyData.sn, file="analysisReadyData.sn.txt", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#Do d'analysis as a group, but use the data whereby the learning effects were removed by subject
df.dPrime <- ch.moralsDprimeAnalysis(analysisReadyData.sn, "overlapRound", "correct", c(1,0), "targetPresent", c(TRUE,FALSE), resCol, params = params, filenameID = "gp")
write.table(df.dPrime, file="df.dPrime.txt", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


#### For experiments with catagory variable manipulations (e.g., different groups), do an analysis
#### by group
#### this one is by side
    grpFitModels <- ch.moralsPlotsByGrpsAndGetModels(analysisReadyData.gp, c("QuantOption1", "QuantOption2"), resCol, "overlapRound", "keyDef", yesNoVal = c("Yes", "No"), "correct", c(1,0), "targetPresent", c(TRUE,FALSE), useTwoParameterModel = TRUE, params = params, minNperOverlap = params$minOverlapN)
    #### this one is by side
    grpFitModels.2 <- ch.moralsPlotsByGrpsAndGetModels(analysisReadyData.gp, c("HVOq", "LVOq"), resCol, "overlapRound", "keyDef", yesNoVal = c("Yes", "No"), "correct", c(1,0), "targetPresent", c(TRUE,FALSE), useTwoParameterModel = TRUE, params = params, minNperOverlap = params$minOverlapN)
    ### and plot the data
    setwd(gpDir)
    ch.moralsPlotFitsByGrps(grpFitModels, c("QuantOption1", "QuantOption2"), "overlapRound", analysisReadyData.gp, filenameID = params$dt.set)
    ch.moralsPlotFitsByGrps(grpFitModels.2, c("HVOq", "LVOq"), "overlapRound", analysisReadyData.gp, filenameID = params$dt.set)
    setwd(mainDir)
