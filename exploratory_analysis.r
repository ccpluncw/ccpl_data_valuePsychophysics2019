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

  ####### Read in exploratory data #########
    exp.dat <-read.table("exp.dat.txt",sep="\t", header=T, quote="\"")
    ##### Exclude data from trials with animal probes
    exp.dat<- exp.dat[exp.dat$category!="animal",]
    exp.dat$catE <- ifelse(exp.dat$category == "object", -1 , 1)

### Run 4 models
    #### 1. Default Model

    model1.gs <- lme(respS ~groupSize , random = ~1|sn, data=exp.dat, control=lmeControl(opt='optim'))
      exp.dat$fit.1.gs <- predict(model1.gs)
      r2.1.gs <- calc.r2(exp.dat$fit.1.gs, exp.dat$respS)
      bic.1.gs<- BIC(model1.gs)
      m1cg.coef <- fixef(model1.gs)

    #### 2. Sum Model

    model1.s <- lme(respS ~Median_Sum_Round , random = ~1|sn/Median_Sum_Round, data=exp.dat, control=lmeControl(opt='optim'))
      exp.dat$fit.1.s <- predict(model1.s)
      r2.1.s <- calc.r2(exp.dat$fit.1.s, exp.dat$respS)
      bic.1.s<- BIC(model1.s)
      m1s.coef <- fixef(model1.s)
    model2.s <- lme(respS ~Median_Sum_Round + catE, random = ~1|sn/Median_Sum_Round, data=exp.dat, control=lmeControl(opt='optim'))
      exp.dat$fit.2.s <- predict(model2.s)
      r2.2.s <- calc.r2(exp.dat$fit.2.s, exp.dat$respS)
      bic.2.s<- BIC(model2.s)
      m2s.coef <- fixef(model2.s)

    anv.out.s <- anova(model1.s, model2.s)

    #### 3. Averaging Model

    model1.m <- lme(respS ~ Mean_Median_Round , random = ~1|sn/Mean_Median_Round, data=exp.dat, control=lmeControl(opt='optim'))
      exp.dat$fit.1.m <- predict(model1.m)
      r2.1.m <- calc.r2(exp.dat$fit.1.m, exp.dat$respS)
      bic.1.m<- BIC(model1.m)
      m1m.coef <- fixef(model1.m)
    model2.m <- lme(respS ~ Mean_Median_Round + catE, random = ~1|sn/Mean_Median_Round, data=exp.dat, control=lmeControl(opt='optim'))
      exp.dat$fit.2.m <- predict(model2.m)
      r2.2.m <- calc.r2(exp.dat$fit.2.m, exp.dat$respS)
      bic.2.m<- BIC(model2.m)
      m3m.coef <- fixef(model2.m)

    anv.out.m <- anova(model1.m, model2.m)

    #### 4. Averaging Plus Model

    model1 <- lme(respS ~ Mean_Median_Round , random = ~1|sn/Mean_Median_Round, data=exp.dat, control=lmeControl(opt='optim'))
      exp.dat$fit.1 <- predict(model1)
      r2.1 <- calc.r2(exp.dat$fit.1, exp.dat$respS)
      bic.1 <- BIC(model1)
      m1.coef <- fixef(model1)
    model2 <- lme(respS ~ Mean_Median_Round + Max_Median_Round, random = ~1|sn/(Mean_Median_Round/Max_Median_Round), data=exp.dat, control=lmeControl(opt='optim'))
      exp.dat$fit.2 <- predict(model2)
      r2.2 <- calc.r2(exp.dat$fit.2, exp.dat$respS)
      bic.2 <- BIC(model2)
      m2.coef <- fixef(model2)
    model3 <- lme(respS ~ Mean_Median_Round + Max_Median_Round + catE, random = ~1|sn/(Mean_Median_Round/Max_Median_Round), data=exp.dat, control=lmeControl(opt='optim'))
      exp.dat$fit.3 <- predict(model3)
      r2.3 <- calc.r2(exp.dat$fit.3, exp.dat$respS)
      bic.3 <- BIC(model3)
      m3.coef <- fixef(model3)

      anv.out <- anova(model1, model2, model3)

### Sink all model outputs
      sink("Exploratory Analysis.txt",append=F)
          cat("\n\n **** Default Model **** \n\n")
          print(summary(model1.gs))
          print(anova(model1.gs))
          print(r2.1.gs)
          print(bic.1.gs)

          cat("\n\n **** Sum Model **** \n\n")
          print(summary(model1.s))
          print(anova(model1.s))
          print(r2.1.s)
          print(bic.1.s)
          cat("\n\n **** Sum + Category **** \n\n")
          print(summary(model2.s))
          print(anova(model2.s))
          print(r2.2.s)
          print(bic.2.s)

          cat("\n\n **** Averaging Model **** \n\n")
          print(summary(model1.m))
          print(anova(model1.m))
          print(r2.1.m)
          print(bic.1.m)
            cat("\n\n **** Average + Category **** \n\n")
          print(summary(model2.m))
          print(anova(model2.m))
          print(r2.2.m)
          print(bic.2.m)

          cat("\n\n **** Averaging Plus Model **** \n\n")
          print(summary(model1))
          print(anova(model1))
          print(r2.1)
          print(bic.1)
          cat("\n\n **** Average + Max **** \n\n")
          print(summary(model2))
          print(anova(model2))
          print(r2.2)
          print(bic.2)
          cat("\n\n **** Average + Max + Category **** \n\n")
          print(summary(model3))
          print(anova(model3))
          print(r2.3)
          print(bic.3)

          sink(NULL)

#### Model Graphs
  par(mfrow=c(1,1))

        df.fit <- exp.dat %>% group_by(Mean_Median_Round) %>% summarise(mRespS = mean(respS), mFit = mean(fit.2))
        df.group <- exp.dat %>% group_by(groupSize, Mean_Median_Round) %>% summarise(mRespS = mean(respS), mFit = mean(fit.2))

        with(df.fit, plot(mRespS ~ Mean_Median_Round, ylim=c(0,4), xlim=c(0,4),pch=16,las=1,frame=F,tck = 0,ylab="mean RespS",xlab="Mean Median Rounded"))
        with(df.fit, lines(mFit ~ Mean_Median_Round))
        abline(0,1, col=gray(.5))
        with(df.group[df.group$groupSize==3,], points(mFit ~ Mean_Median_Round,col=gray(0.75), pch=1))
        with(df.group[df.group$groupSize==2,], points(mFit ~ Mean_Median_Round,col=gray(0.75), pch=1))
        with(df.group[df.group$groupSize==4,], points(mFit ~ Mean_Median_Round,col=gray(0.75), pch=1))
        with(df.group[df.group$groupSize==5,], points(mFit ~ Mean_Median_Round,col=gray(0.75), pch=1))

        dev.copy(pdf,"mFit x Mean_Median_Round.pdf")
        dev.off()

####################### Standard Deviation ##################
  SD.df <- data.frame(exp.dat %>% group_by(groupSize,category) %>% summarise(
      resp_SD=sd(respS),
      Q25 = (quantile(respS, .25, na.rm=T)),
      Q75 = (quantile(respS, .75, na.rm=T)),
      range = (quantile(respS, .75, na.rm=T) - quantile(respS, .25, na.rm=T)),
      median_predSD=median(predSDmean),
      median_maxSD=max(predSDmax),
      median_meanSD=median(meanSD),
      median_range_mean=median(range_mean),
      median_range_max=median(max_range),
      n=length(respS)))


  ######## IQR analysis ##############
    exp.dat$mean_25th<- rowMeans(cbind(exp.dat$tValueQ251,exp.dat$tValueQ252,exp.dat$tValueQ253,exp.dat$tValueQ254,exp.dat$tValueQ255),na.rm = T)
    exp.dat$mean_75th<- rowMeans(cbind(exp.dat$tValueQ751,exp.dat$tValueQ752,exp.dat$tValueQ753,exp.dat$tValueQ754,exp.dat$tValueQ755),na.rm = T)
    exp.dat$IQR <- rowDiffs(cbind(exp.dat$mean_25th, exp.dat$mean_75th),na.rm=T)
    exp.dat <- transform(exp.dat, Z = mean_25th + 0.5 * IQR)
    exp.dat <- transform(exp.dat, QR = (respS - Z)/IQR)

    #Histogram collapsed across group sizes
    par(mfrow=c(1,1))
    b1 <- seq(-19.5, 19.5, 1)
       with(exp.dat, hist(QR,tck=0,las=1,cex.axis=0.9,main="", xlim=c(-10,10),ylim=c(0,0.6),ylab="Percentage of Responses",xlab="Standardized IQR",breaks = b1, freq=F))
      abline(v= -.5,col="red", lwd=2, lty=3)
      abline(v= 0.5,col="red", lwd=2, lty=3)
      dev.copy(pdf,"QR Histogram for All Probes.pdf")
      dev.off()

      QR_df <- data.frame(exp.dat %>% group_by(groupSize) %>% summarise(
        Q25 = (quantile(QR, .25, na.rm=T)),
        Q75 = (quantile(QR, .75, na.rm=T)),
        Q50 = (quantile(QR, .50, na.rm=T)),
        QR_range = (Q75 - Q25)))
    write.table(QR_df,"QR_quantiles.txt",row.names=F,quote=F,sep="\t")
