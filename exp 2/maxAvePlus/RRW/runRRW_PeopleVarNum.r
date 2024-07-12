
require(dplyr)
library(RRW)
library(smartGridSearch)
library(chutils)

ExpName <- "pairs"
test <- FALSE
#use multicore
useMultiCore <- TRUE
#fit free parameters
runSmartGridSearch <- TRUE
#Only run new models that have not been run before?
skipModelsWithExistingDirs <- TRUE

multicorePackages = c('RRW')

mainDir <- getwd()

############## smartGridSearch and simulation parameters
#"BIC" "AIC" "R_Square"
minimizeStat <- "BIC"
#do you want to equalize the influence of RT and pHit on the minimization statistic (TRUE) or have the influence be a function of the number of observations (FALSE)
equalizeRTandPhit <- TRUE
#the minimum number of trials per overlapRound to include in the analysis.
minN <- 40
#number of intervals to split each bounded parameter into in order to get the random values between the hi and low values
numGridIntervals <- 100
#the number of times to loop with grid search before looking for the top 10 best fits.
numGridLoops <- ifelse (test, 10, 2500)
#the number of loops for each RW step.
loopsPerRWstep <- ifelse (test, 20, 500)
#paramtable confirmation loops
optBoundLoops <- ifelse (test, 4, 10)
#The number of best fit runs from which to resize the parameter bounds.
#Here, I set it to the sqrt(total number of runs per grid loop) or 10 at minimum.
  #optParamListN <- ifelse ( trunc(sqrt(numGridLoops)) < 10, 10, trunc(sqrt(numGridLoops)) )
optParamListN <- 10
#the number of simulations to average using the final parameters.
numSimsToAverage <- 40
appendRunStats <- FALSE
######################################################

mainDir <- getwd()
inputDataFile <- "../analysisReadyData.sn.txt"
analysisReadyData<-read.table(inputDataFile, header=T, sep="\t")

### direct.xvy = -1 means the lower valued object will be killed by default
### direct.xvy = 1 means the higher valued object will be killed by default
analysisReadyData$defaultKilled <- ifelse(analysisReadyData$direct.xVy < 0 , "defaultLVO","defaultHVO")

analysisReadyData$peopleQuantity<-paste("HVO_",analysisReadyData$HVOq,"-LVO_", analysisReadyData$LVOq, sep="")

tmp.df<-analysisReadyData

source("modelspecs_PeopleVarNum.r")
modelNames <- names(allModels)

fileTagShortName <- paste(ExpName, sep="_")
fileTagOutStats <- paste(format(Sys.time(), "%b_%d_%Y_%H-%M"), fileTagShortName, sep="_")

### this loop runs the smartGridSearch to optimize free parameters
df.runStats <- NULL
if(runSmartGridSearch) {
  for (i in modelNames) {
    subDirExists <- ch.newDir (mainDir, i)
		if(subDirExists == TRUE) {
			#if a run already exists, then add the output to the previous run
			appendRunStats <- TRUE
		}
    #if the sub directory does NOT exist or you do NOT want to skip models, then run the analysis
    if(subDirExists == FALSE | skipModelsWithExistingDirs == FALSE) {
      modDir <- getwd()
      setwd(modDir)

      tmpModelList <- allModels[[i]]
			
      #create a file tag that will be used to save files.  it will be time and date stamped
      fileTagShortName <- paste(ExpName, i, sep="_")
      fileTag <- paste(format(Sys.time(), "%b_%d_%Y_%H-%M"), fileTagShortName, sep="_")

      run.list <- rrwRunSmartGridSearch(tmp.df, tmpModelList, minN = minN, dataOverlapCol = "overlapRound", RwSamplesCol = "Q50", dataRtCol = "res.RT", correctCol = "correct01", correctVals = c(1,0), loopsPerRWstep = loopsPerRWstep, minimizeStat = minimizeStat, equalizeRTandPhit = equalizeRTandPhit, numLoops = numGridLoops, numIntervals = numGridIntervals, optParamListN = optParamListN, optBoundLoops = optBoundLoops, multicore = useMultiCore, multicorePackages = multicorePackages, numSimsToAverage = numSimsToAverage, fileTag = fileTag)

      df.tmp.stats <- rrwRunStatsToDataframe(run.list$runStats)
      df.tmp.stats$model <- i
      df.runStats <- bind_rows(df.runStats, df.tmp.stats)

      rrwPlotSGSoutput(run.list$df.fitted, tmpModelList, dataRtCol = "rt", dataPhitCol = "pHit", rtFitCol = "rtFit", pHitFitCol = "pCross", correctCol = "correct", overlapCol = "overlap", fileTag = fileTag,numSimsToPlot = numSimsToAverage)
    }

    setwd(mainDir)

    #write the runstats data after each model - just in case the computer bombs, you dont want to loose everything.
    filename <- paste(ExpName, fileTagOutStats, "freeParameterRunStats.txt", sep="_")
    write.table(df.runStats, filename, append = appendRunStats, col.names=T, row.names=F, quote=F, sep="\t")

  }
}

#write the final free parameter runstats dataset
filename <- paste(ExpName, fileTagOutStats, "freeParameterRunStats.txt", sep="_")
write.table(df.runStats, filename, append = appendRunStats, col.names=T, row.names=F, quote=F, sep="\t")


### this loop runs the fits fixed parameters from previous runs
df.runStats <- NULL
if(!is.null(allFixedModels)) {
  modelNames <- names(allFixedModels)

  for (i in modelNames) {
    subDirExists <- ch.newDir (mainDir, i)
    #if the sub directory does NOT exist or you do NOT want to skip models, then run the analysis
    if(subDirExists == FALSE | skipModelsWithExistingDirs == FALSE) {
      modDir <- getwd()
      setwd(modDir)

      tmpModelList <- allFixedModels[[i]]
      
			#create a file tag that will be used to save files.  it will be time and date stamped
      fileTagShortName <- paste(ExpName, i, sep="_")
      fileTag <- paste(format(Sys.time(), "%b_%d_%Y_%H-%M"), fileTagShortName, sep="_")

      #do this to run a set of fixed parameters.  Remember, the rrwModelList should have the upperBound==lowerBound
      run.list <- rrwRunWithFixedParameters(tmp.df, tmpModelList, minN = minN, dataOverlapCol = "overlapRound", RwSamplesCol = "Q50", dataRtCol = "res.RT", correctCol = "correct01", correctVals = c(1,0), loopsPerRWstep = loopsPerRWstep, minimizeStat = minimizeStat, equalizeRTandPhit = equalizeRTandPhit, numSimsToAverage = numSimsToAverage, fileTag = fileTag)

      df.tmp.stats <- rrwRunStatsToDataframe(run.list$runStats)
      df.tmp.stats$model <- i
      df.runStats <- bind_rows(df.runStats, df.tmp.stats)


      rrwPlotSGSoutput(run.list$df.fitted, tmpModelList, dataRtCol = "rt", dataPhitCol = "pHit", rtFitCol = "rtFit", pHitFitCol = "pCross", correctCol = "correct", overlapCol = "overlap", fileTag = fileTag,numSimsToPlot = numSimsToAverage)
    }
    setwd(mainDir)

    #write the runstats data after each model - just in case the computer bombs, you dont want to loose everything.
    filename <- paste(ExpName, fileTagOutStats, "fixedParameterRunStats.txt", sep="_")
    write.table(df.runStats, filename, append = appendRunStats, col.names=T, row.names=F, quote=F, sep="\t")
  }
}

#write the final fixed parameter runstats dataset
filename <- paste(ExpName, fileTagOutStats, "fixedParameterRunStats.txt", sep="_")
write.table(df.runStats, filename, append = appendRunStats, col.names=T, row.names=F, quote=F, sep="\t")
