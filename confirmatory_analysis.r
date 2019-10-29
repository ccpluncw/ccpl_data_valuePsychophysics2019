### Packages
library(dplyr)
library(matrixStats)
library(nlme)

### Preset values
minN <- 20
minIV.val<- -1

#r squared function
calc.r2 <- function(predY, y) {
  #if predicted Y is not all NAs, then ...
    if(length(na.omit(predY)) > 1) {
      meanY <- mean(y, na.rm=T)
      ssE <- sum((y - predY)^2, na.rm=T)
      ssT <- sum((y - meanY)^2, na.rm=T)
      r2 <- 1 - (ssE/ssT)
    } else {
      r2 <- 0
    }
      return(r2)
}

  ####### Read in confirmatory data #########
confirm.dat<-read.table("confirm.dat.txt",sep="\t", header=T, quote="\"")
  ##### Exclude data from trials with animal probes
confirm.dat<- confirm.dat[confirm.dat$category!="animal",]
confirm.dat$catE <- ifelse(confirm.dat$category == "object", -1 , 1)

###### Run final, accepted model

# Final, accepted model
  model.f <- lme(respS ~ -1 + Mean_Median_Round + Max_Median_Round , random = ~1|sn/(Mean_Median_Round/Max_Median_Round), data=confirm.dat, control=lmeControl(opt='optim'))
    confirm.dat$fit.f <- predict(model.f)
    r2.f <- calc.r2(confirm.dat$fit.f, confirm.dat$respS)
    bic.f <- BIC(model.f)
    m4.coef <- fixef(model.f)

#### Mean Median Graph
  par(mfrow=c(1,1))
        df.fit <- confirm.dat %>% group_by(Mean_Median_Round) %>% summarise(mRespS = mean(respS), mFit = mean(fit.f))
        df.group <- confirm.dat %>% group_by(groupSize, Mean_Median_Round) %>% summarise(mRespS = mean(respS), mFit = mean(fit.f))

        with(df.fit, plot(mRespS ~ Mean_Median_Round, ylim=c(0,4), xlim=c(0,4),pch=16,las=1,frame=F,tck = 0,ylab="Transformed Response",xlab="Average Item Value"))
        with(df.fit, lines(mFit ~ Mean_Median_Round))
        abline(0,1, col=gray(.5))
        with(df.group[df.group$groupSize==3,], points(mFit ~ Mean_Median_Round,col=gray(0.75), pch=1))
        with(df.group[df.group$groupSize==2,], points(mFit ~ Mean_Median_Round,col=gray(0.75), pch=1))
        with(df.group[df.group$groupSize==4,], points(mFit ~ Mean_Median_Round,col=gray(0.75), pch=1))
        with(df.group[df.group$groupSize==5,], points(mFit ~ Mean_Median_Round,col=gray(0.75), pch=1))

        dev.copy(pdf,"confirmatory mFit x Mean_Median_Round.pdf")
        dev.off()

  #### QR analysis
  ########IQR analysis##############
    confirm.dat$mean_25th<- rowMeans(cbind(confirm.dat$tValueQ251,confirm.dat$tValueQ252,confirm.dat$tValueQ253,confirm.dat$tValueQ254,confirm.dat$tValueQ255),na.rm = T)
    confirm.dat$mean_75th<- rowMeans(cbind(confirm.dat$tValueQ751,confirm.dat$tValueQ752,confirm.dat$tValueQ753,confirm.dat$tValueQ754,confirm.dat$tValueQ755),na.rm = T)
    confirm.dat$IQR <- rowDiffs(cbind(confirm.dat$mean_25th, confirm.dat$mean_75th),na.rm=T)
    confirm.dat <- transform(confirm.dat, Z = mean_25th + 0.5 * IQR)
    confirm.dat <- transform(confirm.dat, QR = (respS - Z)/IQR)

    #Histogram collapsed across group sizes
    par(mfrow=c(1,1))
    b1 <- seq(-19.5, 19.5, 1)
       with(confirm.dat, hist(QR,tck=0,las=1,cex.axis=0.9,main="", xlim=c(-10,10),ylim=c(0,0.6),ylab="Percentage of Responses",xlab="Standardized IQR",breaks = b1, freq=F))
      abline(v= -.5,col="red", lwd=2, lty=3)
      abline(v= 0.5,col="red", lwd=2, lty=3)
      dev.copy(pdf,"confirmatory QR Histogram for All Probes.pdf")
      dev.off()

      QR_df <- data.frame(confirm.dat %>% group_by(groupSize) %>% summarise(
        Q25 = (quantile(QR, .25, na.rm=T)),
        Q75 = (quantile(QR, .75, na.rm=T)),
        Q50 = (quantile(QR, .50, na.rm=T)),
        QR_range = (Q75 - Q25)))

### Sink all model outputs
        sink("Confirmatory Analysis.txt",append=F)
          cat("\n\n **** Final Model **** \n\n")
          print(summary(model.f))
          print(anova(model.f))
          print(r2.f)
          print(bic.f)

          cat("\n\n **** QR Stats **** \n\n")
          print(summary(confirm.dat$QR))
          print(QR_df)
          sink(NULL)
