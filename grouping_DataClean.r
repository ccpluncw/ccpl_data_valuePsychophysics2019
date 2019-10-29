library(dplyr)
library(matrixStats)
#Create mode function
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
# Create rounding function
  round_any <-  function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

  #Set rounding values for variables
    round.val<-0.1
    sdMult <- 1
    trimP <- .995
    ###Set thresholds for data exclusion based on respS and RT
    minMedRT<-2750
    minRT<- 1000
    maxRT<- 45000

#Read in raw, merged data
  data.grouping.tmp <- read.table("data.grouping.txt",sep="\t", header=T, quote="\"")
  data.grouping<- data.grouping.tmp[!is.na(data.grouping.tmp$sn),]
#Exclude practice trials
  analysis_data.tmp<- data.grouping[data.grouping$pract == 1,]

  # Transformation of raw data
    MinMedTResp <-1
    Temp_MinResp <- min(analysis_data.tmp$resp[analysis_data.tmp$resp>0])
    MinResp <- Temp_MinResp/2
    analysis_data.tmp$resp.tmp <- with(analysis_data.tmp, ifelse(resp==0,MinResp,resp))
    analysis_data.tmp$resp.tmp <- with(analysis_data.tmp, ifelse(resp.tmp<0, MinResp/abs(resp.tmp), resp.tmp))
    analysis_data.tmp$TResp<-log(analysis_data.tmp$resp.tmp)

    #Clean data based on median Tresp
    sn_TResp<- data.frame(analysis_data.tmp %>% group_by(sn) %>% summarise(medTResp=median(TResp),medRT=median(respTime)))
    remove_TResp<-sn_TResp[sn_TResp$medTResp <= MinMedTResp,]
    count_remove_TResp<-length(remove_TResp$sn)
    percent_remove_TResp<-length(count_remove_TResp)/length(sn_TResp$sn)
    tmp.remove<-merge(sn_TResp,analysis_data.tmp,by="sn")
    analysis_data<-tmp.remove[tmp.remove$medTResp>1,]

##Clean data based on RT
    #Identify and remove subjects based on RT thresholds
    bad_dat_sn<-sn_TResp[sn_TResp$medRT <= minMedRT,]
    total.sn <- unique(analysis_data$sn)
    count_remove.sn <- length(bad_dat_sn$sn)
    percent_remove.sn <- length(bad_dat_sn$sn)/length(total.sn)
    dat.tmp<-analysis_data[analysis_data$medRT>minMedRT,]

#identify and remove trials faster than 1s or slower than 45s
    bad_trials<- dat.tmp[dat.tmp$respTime.x <= minRT | dat.tmp$respTime >= maxRT,]
    count_badtrials<-length(bad_trials$trial)
    percent_remove_trials<- length(bad_trials$trial)/length(dat.tmp$trial)
    analysis_data<-dat.tmp[dat.tmp$respTime >= minRT & dat.tmp$respTime <= maxRT,]

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
            if(is.null(outDF)) {
              outDF <- df.tmp
            } else {
              outDF <- rbind(outDF, df.tmp)
            }
          }
        }
#Merge data with file which specifies probe category (animal, person, or object)
  dat.cat<-read.table("probecategories.txt",sep="\t", header=T, quote="\"")
  analysis_data<-merge(outDF,dat.cat,by="probe1")

#Read in response post-transformation.  The transformation was completed to
#link datasets as described in the paper.
  trans.data<-read.table("groupsTransData.txt",sep="\t", header=T, quote="\"")

#Merge data and transformed data
  analysis_data<-merge(trans.data,analysis_data,by=c("sn","trial"),all.x=T)

##Trim data
  #####*** By group, remove top 0.005% of data ***###
    dat.trim <- NULL
    grps <- unique(analysis_data$groupSize)
    grps <- grps[!is.na(grps)]
    grps <- grps[grps != 1]
    for (i in grps) {
      df.tmp1 <- analysis_data[analysis_data$groupSize == i, ]
      high_outliers<-quantile(df.tmp1$respS,probs= trimP)
      df.tmp1<-df.tmp1[df.tmp1$respS < high_outliers,]
      if (is.null(dat.trim)) {
        dat.trim <- df.tmp1
      } else {
        dat.trim <- rbind(dat.trim,df.tmp1)
      }
    }

       summarize.resp<-summary(analysis_data$respS)
     ##isolate top 0.5% of respS, get stats of that subset of data
       high_outliers<-quantile(analysis_data$respS,probs= trimP)
       outlier_dat<-analysis_data[analysis_data$respS>= high_outliers,]
       outlier_summary<-summary(outlier_dat$respS)
     #find trials with negative and zero responses
       trials_0<-analysis_data[analysis_data$respS==0,]
       percent_0<-length(trials_0$respS)/length(analysis_data$respS)
       trials_neg<-analysis_data[analysis_data$respS < 0,]
       percent_neg<-length(trials_neg$respS)/length(analysis_data$respS)
     #take summary of data without negative values
       positive.dat<-analysis_data[analysis_data$respS>=0,]
       pos.dat.sum<-summary(positive.dat$respS)
     #trim top 0.5% of trials, take summary of respS after removing
       dat.trim<-analysis_data[analysis_data$respS < high_outliers,]
       dat.trim_summary<-summary(dat.trim$respS)
#           dat.trim<-dat.trim.tmp[dat.trim.tmp$respS>=0,]
       dat.summary<-summary(dat.trim)

       originalLength <- length(analysis_data$respS)
       finalLength <- length(dat.trim$respS)
       pTrimmed <- 1 - (finalLength/originalLength)


##Read in individual probe stats
 stats.temp<-read.table("ItemStats1.txt",sep="\t", header=T, quote="\"")
 stats.temp2<-read.table("ItemStats2.txt",sep="\t", header=T, quote="\"")
 stats.temp3<-read.table("ItemStats3.txt",sep="\t", header=T, quote="\"")
 stats.temp4<-read.table("ItemStats4.txt",sep="\t", header=T, quote="\"")
 stats.temp5<-read.table("ItemStats5.txt",sep="\t", header=T, quote="\"")
##Create new columns for each item within each probe which list mean and median resp by merging with probe stats
 dat.temp<-merge(dat.trim,stats.temp,by.x=c("probe1"),by.y=c("probe1.1"),dat.values.categories=T)
 dat.temp2<-merge(dat.temp,stats.temp2,by="probe2",all.x=T)
 dat.temp3<-merge(dat.temp2,stats.temp3,by="probe3",all.x=T)
 dat.temp4<-merge(dat.temp3,stats.temp4,by="probe4",all.x=T)
 dat.temp5<-merge(dat.temp4,stats.temp5,by="probe5",all.x=T)
 dat.trim<-dat.temp5

##Create a variable that sums all item means and medians for each probe
 dat.trim["Mean_Sum"]<-NA
 dat.trim$Mean_Sum<-rowSums(cbind(dat.trim$tValueMean1,dat.trim$tValueMean2,dat.trim$tValueMean3,dat.trim$tValueMean4,dat.trim$tValueMean5),na.rm=T)
 dat.trim["Median_Sum"]<-NA
 dat.trim$Median_Sum<-rowSums(cbind(dat.trim$tValueMedian1,dat.trim$tValueMedian2,dat.trim$tValueMedian3,dat.trim$tValueMedian4,dat.trim$tValueMedian5),na.rm = T)
##Create variable that is the mean of the median values of each object within each probe
 dat.trim$Mean_Median<-NA
 dat.trim$Mean_Median<-rowMeans(cbind(dat.trim$tValueMedian1,dat.trim$tValueMedian2,dat.trim$tValueMedian3,dat.trim$tValueMedian4,dat.trim$tValueMedian5),na.rm = T)
 dat.trim$Mean_Mean<-rowMeans(cbind(dat.trim$tValueMean1,dat.trim$tValueMean2,dat.trim$tValueMean3,dat.trim$tValueMean4,dat.trim$tValueMean5),na.rm = T)
##Create four variables that are the minimum and maximum median and mean value of any item within each probe
 dat.trim <- transform(dat.trim, minMedianTvalue = pmin(tValueMedian1,tValueMedian2, tValueMedian3,tValueMedian4,tValueMedian5, na.rm=T))
 dat.trim <- transform(dat.trim, maxMedianTvalue = pmax(tValueMedian1,tValueMedian2, tValueMedian3,tValueMedian4,tValueMedian5, na.rm=T))
 dat.trim <- transform(dat.trim, minMeanTvalue = pmin(tValueMean1,tValueMean2, tValueMean3,tValueMean4,tValueMean5, na.rm=T))
 dat.trim <- transform(dat.trim, maxMeanTvalue = pmax(tValueMean1,tValueMean2, tValueMean3,tValueMean4,tValueMean5, na.rm=T))

##Creates columns to round IVs to a value specified at the top of this script
 dat.trim$Median_Sum_Round<-NA
 dat.trim$Median_Sum_Round<-round_any(dat.trim$Median_Sum,round.val)
 dat.trim$Mean_Sum_Round<-NA
 dat.trim$Mean_Sum_Round<-round_any(dat.trim$Mean_Sum,round.val)
 dat.trim$Min_Median_Round<-NA
 dat.trim$Min_Median_Round<-round_any(dat.trim$minMedianTvalue,round.val)
 dat.trim$Max_Median_Round<-NA
 dat.trim$Max_Median_Round<-round_any(dat.trim$maxMedianTvalue,round.val)
 dat.trim$Max_Mean_Round<-NA
 dat.trim$Max_Mean_Round<-round_any(dat.trim$maxMeanTvalue,round.val)
 dat.trim$Min_Mean_Round<-NA
 dat.trim$Min_Mean_Round<-round_any(dat.trim$minMeanTvalue,round.val)
 dat.trim$Mean_Median_Round<-NA
 dat.trim$Mean_Median_Round<-round_any(dat.trim$Mean_Median,round.val)
 dat.trim$Mean_Mean_Round<-round_any(dat.trim$Mean_Mean,round.val)
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
 ##summarize removed data and write to a file

    sink("removed_dat.txt",append=F)
    cat("\n all subjects\n\t")
    print(sn_TResp)
    cat("\n participants removed for median TResp < 1\n\t")
    print(remove_TResp)
    cat("\n percentage of participants removed for TResp threshold\n\t",percent_remove_TResp)
    cat("\n participants removed for median RT < 2750ms\n\t")
    print(bad_dat_sn)
    cat("\n percentage of participants removed for min. median RT threshold\n\t",percent_remove.sn)
    cat("\n percentage of remaining trials removed (faster than 1s or slower than 45s)\n\t",percent_remove_trials)
    cat("\n trials remaining\n\t",originalLength)
    cat("\n trials remaining after removing outliers\n\t",finalLength)
    cat("\n trim threshold\n\t",trimP)
    cat("\n pTrimmed\n\t",pTrimmed)
     sink(NULL)
