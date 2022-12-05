ch.valuesCombineDatasets <- function(datafilename) {

    require(reshape2)
    require(dplyr)
    require(chValues)
    require(chutils)

      ### Load dbfile with parameters
      params<-read.table(datafilename,row.names=1,header = F,sep=",")

      #### Filenames
      bootstrapOutdataFile <- as.character(params['bootstrapOutdataFile',1])
      bootstrapOutdataFile2 <- as.character(params['bootstrapOutdataFile2',1])
      combinedDataOutFile <- as.character(params['combinedDataOutFile',1])
      sameItemOutFile <- as.character(params['sameItemOutFile',1])
      dataFile1 <- as.character(params['dataFile1',1])
      dataFile2 <- as.character(params['dataFile2',1])

      ######_____SET SWITCHES AND THRESHOLDS______######
      transformDataFile2 <- as.logical(as.character(params['transformDataFile2',1]))
      dataFile1Tag <- as.character(params['dataFile1Tag',1])
      dataFile2Tag <- as.character(params['dataFile2Tag',1])
      outFileTag <- as.character(params['outFileTag',1])
      transformation <- as.numeric(as.character(params['transformation',1]))
      trim1 <- as.numeric(as.character(params['trim1',1]))
      trim2 <- as.numeric(as.character(params['trim2',1]))
      fitType <- as.character(params['fitType',1])
      minNumPerSide <- as.numeric(as.character(params['minNumPerSide',1]))
      maxNumPerSide <- as.numeric(as.character(params['maxNumPerSide',1]))
      numRuns <- as.numeric(as.character(params['numRuns',1]))
      sdThresh <- as.numeric(as.character(params['sdThresh',1]))

    ###### read in datasets

      DataIn1 <-read.table(dataFile1, header=T, sep="\t", quote="\"")
      DataIn2 <-read.table(dataFile2, header=T, sep="\t", quote="\"")

    if(trim1 > 0) {
#      DataIn1 <- ddply(DataIn1, .(prompt), transform, respScale2 = scale(resp))
      DataIn1 <- data.frame(DataIn1 %>% group_by(prompt) %>% mutate(respScale2 = scale(resp)))
      DataIn1 <- DataIn1[DataIn1$respScale > (-1 * trim1) & DataIn1$respScale < trim1, ]
    }
    if(trim2 > 0) {
#      DataIn2 <- ddply(DataIn2, .(prompt), transform, respScale2 = scale(resp))
      DataIn2 <- data.frame(DataIn2 %>% group_by(prompt) %>% mutate(respScale2 = scale(resp)))
      DataIn2 <- DataIn2[DataIn2$respScale > (-1 * trim2) & DataIn2$respScale < trim2, ]
    }


    ##### transform raw responses
    if (transformation == 0) {
      DataIn1$tResp <- ch.altLogTransform(DataIn1$resp)
      if(transformDataFile2) {
        DataIn2$tResp <- ch.altLogTransform(DataIn2$resp)
      } else {
        DataIn2$tResp <- DataIn2$respS
      }
    } else {
      DataIn1$tResp <- ch.altRootTransform(DataIn1$resp, root=transformation)
      if(transformDataFile2) {
        DataIn2$tResp <- ch.altRootTransform(DataIn2$resp, root=transformation)
      } else {
        DataIn2$tResp <- DataIn2$respS
      }
    }

    #make all the prompts as similar as possible
      Items1 <- unique(DataIn1$prompt)
      Items2 <- unique(DataIn2$prompt)

      Items1 <- tolower(trimws(Items1))
      Items2 <- tolower(trimws(Items2))

      DataIn1$prompt <- tolower(trimws(DataIn1$prompt))
      DataIn2$prompt <- tolower(trimws(DataIn2$prompt))

      #get the identical prompts
      sameItems <- Items1[(Items1 %in% Items2)]

      #set the two datasets for the bootstrap
      data1 <- DataIn1[DataIn1$prompt %in% sameItems,]
      data2 <- DataIn2[DataIn2$prompt %in% sameItems,]


    ##############this is where I will explore the transformations ###############

      # data1S <- ddply(data1,.(prompt),summarise, q25s = quantile(tResp,probs=c(.25), type=3),q50s = quantile(tResp,probs=c(.5), type=3), q75s = quantile(tResp,probs=c(.75), type=3))
      # data2S <- ddply(data2,.(prompt),summarise, q25s = quantile(tResp,probs=c(.25), type=3),q50s = quantile(tResp,probs=c(.5), type=3), q75s = quantile(tResp,probs=c(.75), type=3))

      data1S <- data.frame(data1 %>% group_by(prompt) %>% summarize (q25s = quantile(tResp,probs=c(.25), type=3),q50s = quantile(tResp,probs=c(.5), type=3), q75s = quantile(tResp,probs=c(.75), type=3)))
      data2S <- data.frame(data2 %>% group_by(prompt) %>% summarize (q25s = quantile(tResp,probs=c(.25), type=3),q50s = quantile(tResp,probs=c(.5), type=3), q75s = quantile(tResp,probs=c(.75), type=3)))

      data1Sl <- melt(data1S, id.vars = c("prompt"))
      data2Sl <- melt(data2S, id.vars = c("prompt"))

      dataStatMerge <- merge(data1Sl, data2Sl, by=c("prompt","variable"), suffixes=c(".1", ".2"))
      dataStatMerge$scale.1 <- scale(dataStatMerge$value.1)
      dataStatMerge$scale.2 <- scale(dataStatMerge$value.2)
      dataStatMerge <- dataStatMerge[abs(dataStatMerge$scale.1) < sdThresh & abs(dataStatMerge$scale.2) < sdThresh, ]

    #find the linear function that best predicts dataset 2 as a function of dataset 1.
    #I am fitting the quartiles.
    #I regress on tResp because these transformed values are linearly related
      ### linear fit
      if (fitType=="linear") {
        dataFitT.model <- lm(value.2~value.1, data= dataStatMerge)
        dataStatMerge$fit <- coef(dataFitT.model)[1]+coef(dataFitT.model)[2]*dataStatMerge$value.1
      } else {
        dataFitT.model <- nls(value.2~a*value.1^b,start=c(a=1, b=1), data= dataStatMerge)
        dataStatMerge$fit <- coef(dataFitT.model)["a"]*dataStatMerge$value.1^coef(dataFitT.model)["b"]
      }
      dataStatMerge$res <- dataStatMerge$value.2 - dataStatMerge$fit
      lR2 <- 1-(var(dataStatMerge$res)/var(dataStatMerge$value.2))

      ylimVals <- ch.getPlotAxisMinMax(c(dataStatMerge$fit, dataStatMerge$value.2))
      xlimVals <- ch.getPlotAxisMinMax(dataStatMerge$value.1)
      par(mfrow=c(1,1), bg="white",  bty="n", font=2, family='serif', mar=c(5,6,4,7), las=1, cex=2)
      with(dataStatMerge, plot(value.2~value.1, main="distribution statistics for data 1 and data 2", pch=1, ylim=ylimVals, xlim=xlimVals))
          points(fit~value.1, data= dataStatMerge, pch=16)
          mtext(side=4, expression(paste(r^{2}, "=", sep="")),line=-3,at=c(0), cex=1.8)
          mtext(side=4, round(lR2, d=2), line=-2,at=c(0), cex=1.75)

          filename <- paste(outFileTag, "distStatsFit.pdf", sep="")
          dev.copy(pdf, filename, width=12, height=9)
          dev.off()

    ##############################################################################
    ####  Do transformation
    ####  Here I transform dataset 1 using the coefficients I estimated in the fit analysis above
    #these will be output as the final dataset
    #future overlaps should be calculated on Sresp

    if (fitType=="linear") {
      data1$respS <- coef(dataFitT.model)[1]+coef(dataFitT.model)[2]*data1$tResp
      DataIn1$respS <- coef(dataFitT.model)[1] + coef(dataFitT.model)[2]*DataIn1$tResp
    } else {
      ### to deal with negative numbers: sign(tResp)*abs(tResp)^pow. This retains the ordinal position and
      ### applies the transformation.
      data1$respS <- coef(dataFitT.model)["a"]*(sign(data1$tResp)*abs(data1$tResp)^coef(dataFitT.model)["b"])
      DataIn1$respS <- coef(dataFitT.model)["a"] * (sign(DataIn1$tResp)*abs(DataIn1$tResp)^coef(dataFitT.model)["b"])
    }
      data2$respS <- data2$tResp
      DataIn2$respS <- DataIn2$tResp

      if( ("exp" %in% colnames(DataIn2)) == F)
      {
        DataIn2$exp <- dataFile2Tag
      }
      DataIn1$exp <- dataFile1Tag
      DataOut1 <- DataIn1[,c("sn", "exp",  "trial", "cond", "stand", "sValue", "prompt", "respS", "respTime")]
      DataOut2 <- DataIn2[,c("sn", "exp", "trial", "cond", "stand", "sValue", "prompt", "respS", "respTime")]

      DataOut <- rbind(DataOut1, DataOut2)

      #get summary statistics for output later
      # item.Sresp.1 <-ddply(data1,.(prompt),summarise,SrespMean=mean(respS),SrespSd=sd(respS), q50s=median(respS), q25s = quantile(respS,probs=c(.25), type=3), q75s = quantile(respS,probs=c(.75), type=3))
      # item.Sresp.2 <-ddply(data2,.(prompt),summarise,SrespMean=mean(respS),SrespSd=sd(respS), q50s=median(respS), q25s = quantile(respS,probs=c(.25), type=3), q75s = quantile(respS,probs=c(.75), type=3))

      item.Sresp.1 <-data.frame(data1 %>% group_by(prompt) %>% summarize (SrespMean=mean(respS),SrespSd=sd(respS), q50s=median(respS), q25s = quantile(respS,probs=c(.25), type=3), q75s = quantile(respS,probs=c(.75), type=3)))
      item.Sresp.2 <-data.frame(data2 %>% group_by(prompt) %>% summarize (SrespMean=mean(respS),SrespSd=sd(respS), q50s=median(respS), q25s = quantile(respS,probs=c(.25), type=3), q75s = quantile(respS,probs=c(.75), type=3)))

    ###############################################################################
    #do bootstrap to get overlap of the same items with themselves
      totalNumOverlaps <- numCombinations <- length(sameItems)

      averageP <- vector(mode="numeric", length = totalNumOverlaps)
      sdP <- vector(mode="numeric", length = totalNumOverlaps)
      overlap <- vector(mode="numeric", length = totalNumOverlaps)
      direction <- vector(mode="numeric", length = totalNumOverlaps)

      #if the two datasets are on the same scale, then the items from dataset 1 should overlap perfectly with
      #the same items from dataset 2
      i <- 1
      for(l in 1:length(sameItems)) {

        print(paste(l, "of", (numCombinations - 1), sep=" "))
          xValue <- NULL
          yValue <- NULL

          xValue <- data1$respS[data1$prompt==sameItems[l]]
          yValue <- data2$respS[data2$prompt==sameItems[l]]
          pOut <- ch.distOverlap(xDist = xValue, yDist = yValue, numRuns = numRuns)
          averageP[i] <-pOut["percent"]
          sdP[i] <-pOut["sd"]
          overlap[i] <-pOut["overlap"]
          direction[i] <-pOut["direction"]
          i=i+1
          #put data in dataframe
        }

        alldat.1 <-data.frame(sameItems, sameItems, averageP, sdP, overlap, direction)
        for(nps in 1:maxNumPerSide) {
          colnames(alldat.1)[nps] <- paste("IA", nps, sep="")
          colnames(alldat.1)[nps+maxNumPerSide] <- paste("IB", nps, sep="")
        }

        #merge the overlap data with the new transformed statistics for output later.
        tmp <- merge(alldat.1, item.Sresp.1, by.x=c("IA1"), by.y=c("prompt"))
        alldat.2 <- merge(tmp, item.Sresp.2, by.x=c("IA1"), by.y=c("prompt"), suffixes = c(".1", ".2"))

    #############################################################################

      #now get overlap of same items with other same items
      itemSet <- sameItems
      #drop extraneous columns and levels
      itemSet <- factor (itemSet)
      numItems <- length(itemSet)

      #create a vector of all possible combinations on one side of the dilemma equation
      df.combns <- ch.combnVector (numItems, minNumPerSide,maxNumPerSide)

      #Now use df.comns to do the bootstrap.  That is, use the values to index the items and grab samples, if the indexs are NA then ignore.

      for (i in 1:length(df.combns$I1)) {
      print(paste(itemSet[df.combns[i,1]],itemSet[df.combns[i,2]] ))
      }

      #do bootstrap for each dataset combination (1v1, 1v2, 2v2, 2v1)
      # overlap XvX is dataset1 overlap with iteself.  The correlation is less than 1 because the bootstrap process is stochastic
      # overlap XvY is dataset1 overlap with detaset2.  If the 2 datasets are perfectly scaled, then the XvY correlation should == XvX correlation.
      #if the two datasets are measuring the same thing, then 1v1 and 2v2 overlaps should be correlated
      #if the two datasets are on the same scale, then 1v2 and 2v1 should be correlated
      for(changeSides in 1:4) {

        if(changeSides < 3) {
            datasetA <- 1
            xDat <- data1
          } else {
            datasetA <- 2
            xDat <- data2
        }
        if(changeSides == 1 | changeSides == 4) {
            datasetB <- 1
            yDat <- data1
          } else {
            datasetB <- 2
            yDat <- data2
        }

        alldat <- ch.batchOverlapXY (xDat$respS, xDat$prompt, yDat$respS, yDat$prompt,itemSet, df.combns, numRuns)

        if (changeSides == 1) {
          alldat.out <- alldat
          suf1 <- paste(datasetA,datasetB, sep=".")
        } else {
          if (changeSides > 2) {
            suf1 <- suf2
          }
          suf2 <- paste(datasetA,datasetB, sep=".")
          suf <- c(suf1,suf2)
          alldat.out <- merge (alldat.out, alldat, by=c("IA1", "IB1"), suffixes=suf)
        }
      }

    ##############################################################################

        #see how well the two datasets overlaps correlated now that they are transformed
        overlaps <- with(alldat.out, data.frame(overlap1.1, overlap2.2, overlap1.2, overlap2.1))
        cors <- cor(overlaps)

        #output table of the overlaps of the same items from the different datasets
        filename <- paste(bootstrapOutdataFile, sep="")
        write.table(alldat.2, file=filename, quote=F, sep="\t", row.names = F, append=F)
        #output overlaps from each dataset combination (1.1, 1.2, 2.2, 2.1)
        filename <- paste(bootstrapOutdataFile2, sep="")
        write.table(alldat.out, file=filename, quote=F, sep="\t", row.names = F, append=F)
        #output new combines values dataset
        filename <- paste(combinedDataOutFile, sep="")
        write.table(DataOut, file=filename, quote=F, sep="\t", row.names = F, append=F)
        #output items the two datasets had in common
        filename <- paste(sameItemOutFile, sep="")
        write.table(sameItems, file=filename, quote=F, sep="\t", row.names = F, append=F)

    ##############################################################################
    ####  Make plots

    	par(mfrow=c(1,1), bg="white",  bty="n", font=2, family='serif', mar=c(5,6,4,7), las=1, cex=2)
    	with(alldat.out, plot(overlap1.1, overlap2.2, main="1.1 vs 2.2", pch=16, ylim=c(0,1),xlim=c(0,1) ))
    			lmFit <- with(alldat.out,lm(overlap2.2 ~ 0 + overlap1.1))
          lR2 <- 1-(var(residuals(lmFit))/var(alldat.out$overlap2.2))
    			abline(lmFit, col="black", lwd=3)
    			mtext(side=4, expression(paste(r^{2}, "=", sep="")),line=-3,at=c(0), cex=1.8)
    			mtext(side=4, round(lR2, d=2), line=-2,at=c(0), cex=1.75)

          filename <- paste(outFileTag,"withinDatasetComparison.pdf", sep="")
    			dev.copy(pdf, filename, width=12, height=9)
    			dev.off()

      with(alldat.out,plot(overlap1.2, overlap2.1, main="1.2 vs 2.1", pch=16, ylim=c(0,1),xlim=c(0,1) ))
          lmFit <- with(alldat.out,lm(overlap2.1 ~ 0 + overlap1.2))
          lR2 <- 1-(var(residuals(lmFit))/var(alldat.out$overlap2.1))
          abline(lmFit, col="black", lwd=3)
          mtext(side=4, expression(paste(r^{2}, "=", sep="")),line=-3,at=c(0), cex=1.8)
          mtext(side=4, round(lR2, d=2), line=-2,at=c(0), cex=1.75)

          filename <- paste(outFileTag,"betweenDatasetComparison.pdf", sep="")
          dev.copy(pdf, filename, width=12, height=9)
          dev.off()

      with(alldat.1,boxplot(overlap,  main="same item overlaps"))

          filename <- paste(outFileTag,"sameItemOverlaps.pdf", sep="")
          dev.copy(pdf, filename, width=12, height=9)
          dev.off()

    ######### Sink some important statistics #######################################

        filename <- paste(outFileTag,"compareBootstraps.txt", sep="")
        sink(filename)
          print("mean overlap of same items")
          print(mean(alldat.1$overlap))
          print("sd overlap of same items")
          print(sd(alldat.1$overlap))
          print("correlation of overlap of other items")
          print(cors)
          print("resp to tResp transformation:")
          if(transformation==0) {
            print("altLogTransform")
          } else {
            print("altRootTransform")
            print(transformation)
          }
          print("tranformation")
          print(fitType)
          print(summary(dataFitT.model))
        sink(NULL)
}
