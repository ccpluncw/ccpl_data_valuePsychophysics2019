### Packages
require(dplyr)
require(matrixStats)
require(MuMIn)
require(lmerTest)
require(BSDA)


### Preset values
minN <- 20

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
confirm.dat$sn <- as.factor(confirm.dat$sn)

sns <- unique(confirm.dat$sn)

df.subFits <- NULL
for(i in sns) {
 df.tmp <- confirm.dat[confirm.dat$sn == i,]
 df.tmp$maxAve <- 0.5*df.tmp$Mean_Median + 0.5*df.tmp$Max_Median
 lm.maxAve <- with(df.tmp, lm(respS ~ 1 + offset(maxAve)))
 df.tmp$biasedAve <- fixef(model.biasedMean)["Mean_Median"]*df.tmp$Mean_Median + fixef(model.biasedMean)["Max_Median"]*df.tmp$Max_Median
 lm.biasedAve <- with(df.tmp, lm(respS ~ 1 + offset(biasedAve)))
 df.tmp$maxAveFit <- predict(lm.maxAve)
 df.tmp$biasedAveFit <- predict(lm.biasedAve)
 maxAveFit.r2 <- calc.r2(df.tmp$maxAveFit, df.tmp$respS)
 biasedAveFit.r2 <- calc.r2(df.tmp$biasedAveFit, df.tmp$respS)

 df.out <- data.frame(sn = i, maxAveInt = coef(lm.maxAve)[1], biasedAveInt = coef(lm.biasedAve)[1], maxAveFit.r2, biasedAveFit.r2)
 df.subFits <- ch.rbind(df.subFits, df.out)
}

#get which model is the best fit based on r2 for each participant?
df.subFits$bestFit <- colnames(df.subFits[,4:5])[apply(df.subFits[,4:5],1,which.max)]
#Do sign tests to see if the model predicts better thnan chance
maxAve.signTest <- SIGN.test(df.subFits$maxAveFit.r2, md = 0)
biasedAve.signTest <- SIGN.test(df.subFits$biasedAveFit.r2, md = 0)
maxAveVsAve.signTest <- SIGN.test(df.subFits$maxAveFit.r2, df.subFits$biasedAveFit.r2)

###### Run final, accepted model

#strong fit
#biasedAve
    #create the prediction that contains only the fixed effects - use the fixed parameter from the model for catE
    confirm.dat$biasedAveFit <- predict(model.biasedMean, confirm.dat, allow.new.levels = T)

    #now get the random effects plus the effect of category
    model.biasedMean.confirm <- lmer(respS ~  -1 + offset(biasedAveFit) + (1+ biasedAveFit|sn), data=confirm.dat, control=lmerControl(optimizer='Nelder_Mead'))

    #now creat the best fit by adding the model prediction to the predictor variable. This contains random and fixed effects
    confirm.dat$fit.overall.biasedAve.a <- predict(model.biasedMean.confirm)
    #approximate the overall and marginal r2
    model.biasedAve.r2.overall <- calc.r2(confirm.dat$fit.overall.biasedAve.a, confirm.dat$respS)
    model.biasedAve.r2.marg <- calc.r2(confirm.dat$biasedAveFit, confirm.dat$respS)
    model.biasedAve.BIC <- BIC(model.biasedMean.confirm)
    model.biasedAve.ranSlope.tt <- t.test(ranef(model.biasedMean.confirm)$sn["biasedAveFit"], mu=0)
    model.biasedAve.ranInt.tt <- t.test(ranef(model.biasedMean.confirm)$sn["(Intercept)"], mu=0)

#MaxAve
    #create the prediction that contains only the fixed effects - use the fixed parameter from the model for catE
    confirm.dat$maxAve <- .5*confirm.dat$Max_Median + .5*confirm.dat$Mean_Median
    #now get the random effects plus the effect of category
    model.maxAve.confirm <- lmer(respS ~  -1 + offset(maxAve) + (1+ maxAve|sn), data=confirm.dat, control=lmerControl(optimizer='Nelder_Mead'))

    #now creat the best fit by adding the model prediction to the predictor variable. This contains random and fixed effects
    confirm.dat$fit.overall.maxAve.a <- predict(model.maxAve.confirm)
    #approximate the overall and marginal r2
    model.maxAve.r2.overall <- calc.r2(confirm.dat$fit.overall.maxAve.a, confirm.dat$respS)
    model.maxAve.r2.marg <- calc.r2(confirm.dat$maxAve, confirm.dat$respS)
    model.maxAve.BIC <- BIC(model.maxAve.confirm)
    model.maxAve.ranSlope.tt <- t.test(ranef(model.maxAve.confirm)$sn["maxAve"], mu=0)
    model.maxAve.ranInt.tt <- t.test(ranef(model.maxAve.confirm)$sn["(Intercept)"], mu=0)

    #### Mean Median Graph

    df.fit <- confirm.dat %>% group_by(Mean_Median_Round) %>% summarise(mRespS = mean(respS),n = length(respS), weakerFit.biasedAve = mean(biasedAveFit), strongFit.maxAve = mean(maxAve))
    df.group <- confirm.dat %>% group_by(groupSize, Mean_Median_Round) %>% summarise(mRespS = mean(respS),n = length(respS), weakerFit.biasedAve = mean(biasedAveFit), strongFit.maxAve = mean(maxAve))

    pdf("confirmatory Prediction x Response.pdf")
      df.tmp <- df.group[df.group$n>minN,]
      for(i in 5:length(df.group)) {
        maxX <- ifelse(max(df.tmp[[i]]) < 4, 4, 20)
        plot(df.tmp$mRespS ~ df.tmp[[i]], ylim=c(0,5), xlim=c(0,maxX),pch=16,las=1,frame=F,tck = 0,ylab="Group Value Response",xlab=colnames(df.tmp[i]), main = colnames(df.tmp[i]))
        abline(0,1, col=gray(.5))
      }

      df.tmp <- df.group[df.group$n>minN,]
      for(i in 5:length(df.group)) {
          maxX <- ifelse(max(df.tmp[[i]]) < 4, 4, 20)
        plot(NULL, ylim=c(0,5), xlim=c(0,maxX),pch=16,las=1,frame=F,tck = 0,ylab="Group Value Response",xlab=colnames(df.tmp[i]), main = colnames(df.tmp[i]))
        points(df.tmp[df.tmp$groupSize==5,][["mRespS"]] ~ df.tmp[df.tmp$groupSize==5,][[i]],col=gray(0.25), pch=16)
        points(df.tmp[df.tmp$groupSize==4,][["mRespS"]] ~ df.tmp[df.tmp$groupSize==4,][[i]],col=gray(0.40), pch=16)
        points(df.tmp[df.tmp$groupSize==3,][["mRespS"]] ~ df.tmp[df.tmp$groupSize==3,][[i]],col=gray(0.55), pch=16)
        points(df.tmp[df.tmp$groupSize==2,][["mRespS"]] ~ df.tmp[df.tmp$groupSize==2,][[i]],col=gray(0.70), pch=16)
        abline(0,1, col=gray(.5))
      }

      df.fit <- df.fit[order(df.fit$Mean_Median_Round),]
      with(df.group[df.group$n>minN,], plot(NULL, ylim=c(0,15), xlim=c(0,4),pch=16,las=1,frame=F,tck = 0,ylab="Prediction",xlab="Average Group Value", main="Strong Model Predictions"))
      with(df.group[df.group$groupSize == 2 & df.group$n>minN, ], points(mRespS ~ Mean_Median_Round, pch=16, col=gray(0.70)))
      with(df.group[df.group$groupSize == 3 & df.group$n>minN, ], points(mRespS ~ Mean_Median_Round, pch=16, col=gray(0.55)))
      with(df.group[df.group$groupSize == 4 & df.group$n>minN, ], points(mRespS ~ Mean_Median_Round, pch=16, col=gray(0.40)))
      with(df.group[df.group$groupSize == 5 & df.group$n>minN, ], points(mRespS ~ Mean_Median_Round,, pch=16, col=gray(0.25)))
      with(df.fit[df.fit$n>minN, ], lines(weakerFit.biasedAve ~ Mean_Median_Round,, lwd=2, col=gray(0.25)))
      abline(0,1, lty=3, lwd=2)

    dev.off()

  #### QR analysis
  ########IQR analysis##############
    confirm.dat$mean_25th<- rowMeans(cbind(confirm.dat$tValueQ251,confirm.dat$tValueQ252,confirm.dat$tValueQ253,confirm.dat$tValueQ254,confirm.dat$tValueQ255),na.rm = T)
    confirm.dat$mean_75th<- rowMeans(cbind(confirm.dat$tValueQ751,confirm.dat$tValueQ752,confirm.dat$tValueQ753,confirm.dat$tValueQ754,confirm.dat$tValueQ755),na.rm = T)
    confirm.dat$IQR <- rowDiffs(cbind(confirm.dat$mean_25th, confirm.dat$mean_75th),na.rm=T)
    confirm.dat <- transform(confirm.dat, Z = mean_25th + 0.5 * IQR)
    confirm.dat <- transform(confirm.dat, QR = (respS - Z)/IQR)

    #Histogram collapsed across group sizes
    b1 <- seq(-39.5, 39.5, 1)
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
          cat("*** Total Sn:", length(unique(confirm.dat$sn)), "\n\n")
          cat("\n\n **** Confirm Weaker Model **** \n\n")
          cat("\n Approximate Overall R2: ", model.biasedAve.r2.overall)
          cat("\n Approximate Marginal R2: ", model.biasedAve.r2.marg)
          cat("\n BIC: ", model.biasedAve.BIC)
          cat("\n R2 summarized by groups, minN = ", minN, "; R2 = ", with(df.group[df.group$n>minN,], calc.r2(weakerFit.biasedAve, mRespS)))
          cat("\n\n Test whether the Random Effects Intercepts == 0:  T-Test,  Mu = 0 \n")
          cat("If the model fits well, the Random Effects should have a mean around 0 \n")
          print(model.biasedAve.ranInt.tt)
          cat("\n\n Test whether the Random Effects Slopes == 0:  T-Test,  Mu = 0 \n")
          cat("If the model fits well, the Random Effects should have a mean around 0 \n")
          print(model.biasedAve.ranSlope.tt)
          cat("\n\n When we fit the model to each individual participant, do the r2's differ from 0? \n")
          cat("If the model fits well, the r2's should be greater than 0. \n")
          cat("We do a non-parametric Sign Test,  Mu = 0, because there are some outliers \n")
          print(biasedAve.signTest)

          cat("\n\n **** Confirm Stronger Model MaxAve **** \n\n")
          cat("\n Approximate Overall R2: ", model.maxAve.r2.overall)
          cat("\n Approximate Marginal R2: ", model.maxAve.r2.marg)
          cat("\n BIC: ", model.maxAve.BIC)
          cat("\n R2 summarized by groups, minN = ", minN, "; R2 = ", with(df.group[df.group$n>minN,], calc.r2(strongFit.maxAve, mRespS)))
          cat("\n\n Test whether the Random Effects Intercepts == 0:  T-Test,  Mu = 0 \n")
          cat("If the model fits well, the Random Effects should have a mean around 0 \n")
          print(model.maxAve.ranInt.tt)
          cat("\n\n Test whether the Random Effects Slopes == 0:  T-Test,  Mu = 0 \n")
          cat("If the model fits well, the Random Effects should have a mean around 0 \n")
          print(model.maxAve.ranSlope.tt)
          cat("\n\n When we fit the model to each individual participant, do the r2's differ from 0? \n")
          cat("If the model fits well, the r2's should be greater than 0. \n")
          cat("We do a non-parametric Sign Test,  Mu = 0, because there are some outliers \n")
          print(maxAve.signTest)

          cat("\n\n Does the Weak model explain significantly more variance than the Strong model in individual participants? \n")
          cat("We do a non-parametric Sign Test,  Max Ave Vs Biased Ave Sign Test, because there are some outliers \n")
          print(maxAveVsAve.signTest)

          cat("\n\n **** QR Stats **** \n\n")
          print(summary(confirm.dat$QR))
          print(QR_df)
        sink(NULL)

        op <- par(mfrow=c(2,2), bty="n", mar=c(6,4,4,1), oma=c(1,1,1,1))
          with(confirm.dat, boxplot(respS ~ groupSize, main = "Data", ylim=c(0,15), ylab="Psychological Value", xlab="Group Size"))
          with(confirm.dat, boxplot(fit.overall.maxAve.a ~ groupSize, main = "Predicted\nBiased Average", ylim=c(0,15), ylab="Psychological Value", xlab="Group Size"))
          with(confirm.dat, boxplot(Sum_Median ~ groupSize, main = "Predicted\nSum", ylim=c(0,15), ylab="Psychological Value", xlab="Group Size"))
          with(confirm.dat, boxplot(Mean_Median ~ groupSize, main = "Predicted\nAverage", ylim=c(0,15), ylab="Psychological Value", xlab="Group Size"))
        dev.copy(pdf,"confirmatory data groupSize Boxplots.pdf", width=8, height=8)
        dev.off()
        par(op)
