require(chutils)

#Create mode function
  Mode <- function(x) {
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
  }
# Create rounding function
  round_any <-  function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

  #Set rounding values for variables
    numberOfBins <- 50
    round.val<-0.1
    trimPhigh <- .995
    trimPlow <- .005
    ###Set thresholds for data exclusion based on respS and RT
    minMedRT<-2750
    minRT<- 1000
    maxRT<- 45000
    MinMedTResp <- 10

    #Read in raw, merged data
    data.grouping<- data.grouping.tmp[!is.na(data.grouping.tmp$sn),]
    #Exclude practice trials
    analysis_data.tmp<- data.grouping[data.grouping$pract == 1,]

    #Clean data based on median resp
    #remove Ss who give a median response response less than MinMedTResp

    sn_TResp<- data.frame(analysis_data.tmp %>% group_by(sn) %>% summarise(medTResp=median(resp),medRT=median(respTime)))
    count_total_Raw_Sn <- length(sn_TResp$sn)
    remove_TResp<-sn_TResp[sn_TResp$medTResp <= MinMedTResp,]
    count_remove_TResp<-length(remove_TResp$sn)
    percent_remove_TResp<-length(count_remove_TResp)/length(sn_TResp$sn)
    analysis_data<-analysis_data.tmp[!(analysis_data.tmp$sn %in% remove_TResp$sn), ]
    remainingSn <- sn_TResp[sn_TResp$medTResp > MinMedTResp,]
    count_remaining_sn <- nrow(remainingSn)

    ##Clean data based on RT
    if(filterSsBasedOnRT) {
      #Identify and remove subjects based on RT thresholds
      bad_dat_sn<-remainingSn[remainingSn$medRT <= minMedRT,]
      total.sn <- length(remainingSn$sn)
      count_remove.sn <- length(bad_dat_sn$sn)
      percent_remove.sn <- length(bad_dat_sn$sn)/length(total.sn)
      analysis_data<-analysis_data[!(analysis_data$sn %in% bad_dat_sn$sn), ]
    }  else {
      bad_dat_sn<- NULL
      total.sn <- length(remainingSn$sn)
      count_remove.sn <- 0
      percent_remove.sn <- length(bad_dat_sn$sn)/length(total.sn)
    }

    #identify and remove trials faster than 1s or slower than 45s
        bad_trials<- analysis_data[analysis_data$respTime <= minRT | analysis_data$respTime >= maxRT,]
        count_badtrials<-length(bad_trials$trial)
        percent_remove_trials<- length(bad_trials$trial)/length(analysis_data$trial)
        analysis_data<-analysis_data[analysis_data$respTime > minRT & analysis_data$respTime < maxRT,]


        ##Use a for loop to separate probes into individual items which each get their own new column and a column that specifies group size
          maxGroupSize <- 5
          subs <- unique(analysis_data$sn)
          outDF <- NULL
          for(i in subs) {
            sess <- unique(analysis_data[analysis_data$sn==i,"sess"])
            for(j in sess){
              tmp.dat <- analysis_data[analysis_data$sn==i & analysis_data$sess == j,]
              #count rows
              numTrials <- length(tmp.dat$sn)
              groupSize <- vector("numeric")
              probeArray <- array(NA, dim=c(numTrials,maxGroupSize))
              for(h in 1:numTrials) {
                #get groupSize for each row
                groupSize[h] <- length(unlist(strsplit(as.character(tmp.dat$probe[h]), ';')))
                tmpStr <- trimws(unlist(strsplit(as.character(tmp.dat$probe[h]), ';')))
                for(g in 1:groupSize[h]) {
                  probeArray[h,g] <- tmpStr[g]
                }
              }
              df.tmp <- data.frame(groupSize = groupSize, probeArray)
              for (k in 1:maxGroupSize) {
                probeCol <- paste("probe",k,sep="")
                colnames(df.tmp)[c(1+k)] <- probeCol
              }
              df.tmp$sn <- i
              df.tmp$sess <- j
              df.tmp$trial <- tmp.dat$trial
              df.tmp$probe <- tmp.dat$probe
              df.tmp$resp <- tmp.dat$resp
              df.tmp$respTime <- tmp.dat$respTime
              df.tmp$prompt <- tmp.dat$prompt
              df.tmp$cond <- tmp.dat$cond
              df.tmp$stand <- tmp.dat$stand
              df.tmp$sValue <- tmp.dat$sValue
              if(is.null(outDF)) {
                outDF <- df.tmp
              } else {
                outDF <- rbind(outDF, df.tmp)
              }
            }
          }


#Merge data with file which specifies probe category (animal, person, or object) - read in "runGrouping.r"
  analysis_data <-merge(outDF,dat.cat,by="probe1")

##Trim data
  #####*** By group, remove top and bottom 0.005% of data ***###
    df.filtered <- ch.filterGrpByQuantile(analysis_data, "resp", grpCol = "groupSize", lowQuantileThreshold=trimPlow, highQuantileThreshold=trimPhigh)
    dat.trim <- df.filtered$datKept
    pTrimmed <- df.filtered$pRemoved

    ##summarize removed data and write to a file
       sink("removed_dat.txt",append=F)
         cat("\n Total number of raw subjects: ", count_total_Raw_Sn, "\n\t")
         print(sn_TResp)
         cat("\n participants removed for median TResp < ", MinMedTResp, "\n\t")
         print(remove_TResp)
         cat("\n Number of participants removed because of median TResp criterion: ",count_total_Raw_Sn - count_remaining_sn, "\n")
         cat("\n Number of participants remaining in analysis: ", count_remaining_sn, "\n")
         cat("\n percentage of participants removed for TResp criterion: ",percent_remove_TResp, "\n")
         cat("\n participants removed for median RT < 2750ms\n\t")
         print(bad_dat_sn)
         cat("\n percentage of participants removed for min. median RT threshold: ",percent_remove.sn, "\n")
         cat("\n percentage of trials removed (faster than ", minRT, "ms or slower than ", maxRT, "ms): ",percent_remove_trials, "\n")
         cat("\n Trimmed high (quantile is greater than", trimPhigh, ") and low (quantile is lower than", trimPlow, ") outliers by group size. Proportion trimmed: ",pTrimmed, "\n")
       sink(NULL)


    #output individual files for each group sizes
    #these will be used for transformations
    dat.trim.sum <- dat.trim %>% group_by(sn) %>% summarize( groupNum = ceiling(mean(groupSize)))
    dat.trim <- merge(dat.trim, dat.trim.sum)

    for(gps in 2:5) {
      df.tmp <- dat.trim[dat.trim$groupNum == gps,]
      filename <- paste("df_",gps,".txt", sep="")
      write.table(df.tmp,filename,row.names=F,quote=F,sep="\t")
    }

    ch.valuesCombineDatasets("../combValueDBfile2.txt")
    ch.valuesCombineDatasets("../combValueDBfile3.txt")
    ch.valuesCombineDatasets("../combValueDBfile4.txt")
    ch.valuesCombineDatasets("../combValueDBfile5.txt")

    #Read in response post-transformation.  The transformation was completed to
    #link datasets as described in the paper.
    trans.data.1<-read.table("commonProbes2345ValuesData.txt",sep="\t", header=T, quote="\"")
    trans.data.2 <- trans.data.1[trans.data.1$exp !="Minwoo" & trans.data.1$exp !="exp1",]


    #Merge data and transformed data
    analysis_dat.trim<-merge(trans.data.2,dat.trim,by=c("sn","trial", "respTime", "sValue", "stand", "cond"),all.x=T)

## stats data was read in "runGrouping.r"
##Create new columns for each item within each probe which list mean and median resp by merging with probe stats
 dat.temp<-merge(analysis_dat.trim,stats.temp,by.x=c("probe1"),by.y=c("probe1.1"),dat.values.categories=T)
 dat.temp2<-merge(dat.temp,stats.temp2,by="probe2",all.x=T)
 dat.temp3<-merge(dat.temp2,stats.temp3,by="probe3",all.x=T)
 dat.temp4<-merge(dat.temp3,stats.temp4,by="probe4",all.x=T)
 dat.temp5<-merge(dat.temp4,stats.temp5,by="probe5",all.x=T)
 dat.trim<-dat.temp5

##Create a variable that sums all item means and medians for each probe
 dat.trim["Sum_Median"]<-NA
 dat.trim$Sum_Median<-rowSums(cbind(dat.trim$tValueMedian1,dat.trim$tValueMedian2,dat.trim$tValueMedian3,dat.trim$tValueMedian4,dat.trim$tValueMedian5),na.rm = T)
##Create variable that is the mean of the median values of each object within each probe
 dat.trim$Mean_Median<-NA
 dat.trim$Mean_Median<-rowMeans(cbind(dat.trim$tValueMedian1,dat.trim$tValueMedian2,dat.trim$tValueMedian3,dat.trim$tValueMedian4,dat.trim$tValueMedian5),na.rm = T)

##Create variable that is the  maximum median  value of any item within each probe
 dat.trim <- transform(dat.trim, Max_Median = pmax(tValueMedian1,tValueMedian2, tValueMedian3,tValueMedian4,tValueMedian5, na.rm=T))

##Creates columns to round IVs to a value specified at the top of this script
  if(numberOfBins == 0) {
   dat.trim$Sum_Median_Round<-round_any(dat.trim$Sum_Median,round.val)
   dat.trim$Max_Median_Round<-round_any(dat.trim$Max_Median,round.val)
   dat.trim$Mean_Median_Round<-round_any(dat.trim$Mean_Median,round.val)
 } else {
   dat.trim$Sum_Median_Round<-ch.binNumbers(dat.trim$Sum_Median,numberOfBins)
   dat.trim$Max_Median_Round<-ch.binNumbers(dat.trim$Max_Median,numberOfBins)
   dat.trim$Mean_Median_Round<-ch.binNumbers(dat.trim$Mean_Median,numberOfBins)
 }

 dat.trim<- transform(dat.trim, predSDmean= sqrt((rowSums(cbind((tValueSd1^2),(tValueSd2^2),(tValueSd3^2),(tValueSd4^2),(tValueSd5^2)),na.rm=T)))/4)
 dat.trim$meanSD<-rowMeans(cbind(dat.trim$tValueSd1,dat.trim$tValueSd2,dat.trim$tValueSd3,dat.trim$tValueSd4,dat.trim$tValueSd5),na.rm = T)
 dat.trim <- transform(dat.trim, predSDmax = pmax(tValueSd1,tValueSd2,tValueSd3,tValueSd4,tValueSd5, na.rm=T))

#for each row, get the range (75th percentile - 25th percentile) for each item, then get mean
  dat.trim$range_1 <- rowDiffs(cbind(dat.trim$tValueQ251,dat.trim$tValueQ751),na.rm=T)
  dat.trim$range_2 <- rowDiffs(cbind(dat.trim$tValueQ252,dat.trim$tValueQ752),na.rm=T)
  dat.trim$range_3 <- rowDiffs(cbind(dat.trim$tValueQ253,dat.trim$tValueQ753),na.rm=T)
  dat.trim$range_4 <- rowDiffs(cbind(dat.trim$tValueQ254,dat.trim$tValueQ754),na.rm=T)
  dat.trim$range_5 <- rowDiffs(cbind(dat.trim$tValueQ255,dat.trim$tValueQ755),na.rm=T)
  dat.trim$range_mean <-rowMeans(cbind(dat.trim$range_1,dat.trim$range_2,dat.trim$range_3,dat.trim$range_4,dat.trim$range_5),na.rm = T)

#write trimmed data to a file
 write.table(dat.trim,"dat.trim.txt",row.names=F,quote=F,sep="\t")
