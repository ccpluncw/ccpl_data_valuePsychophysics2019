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

  ####### Read in exploratory data #########
    exp.dat <-read.table("exp.dat.txt",sep="\t", header=T, quote="\"")
    ##### Exclude data from trials with animal probes
    exp.dat<- exp.dat[exp.dat$category!="animal",]
    exp.dat$catE <- ifelse(exp.dat$category == "object", -1 , 1)
    exp.dat$sn <- as.factor(exp.dat$sn)


  #### 1. Default Model

  model.default <- lmer(respS ~ groupSize  +  (-1+groupSize|sn), data=exp.dat, control=lmerControl(optimizer='Nelder_Mead'))
    model.default.r2 <- r.squaredGLMM(model.default)
    model.default.BIC <- BIC(model.default)
    model.default.coef <- fixef(model.default)
    model.default.anova <- anova(model.default)

### Run Strong Models

   strong.exp <- exp.dat
   sns <- unique(strong.exp$sn)

  #fit individual participants with the strong models
  df.subFits <- NULL
  for(i in sns) {
    df.tmp <- strong.exp[strong.exp$sn == i,]
    df.tmp$maxAve <- 0.5*df.tmp$Mean_Median + 0.5*df.tmp$Max_Median
    lm.sum <- with(df.tmp, lm(respS ~ 1 + offset(Sum_Median)))
    lm.ave <- with(df.tmp, lm(respS ~ 1 + offset(Mean_Median)))
    lm.maxAve <- with(df.tmp, lm(respS ~ 1 + offset(maxAve)))
    df.tmp$sumFit <- predict(lm.sum)
    df.tmp$aveFit <- predict(lm.ave)
    df.tmp$maxAveFit <- predict(lm.maxAve)
    sumFit.r2 <- calc.r2(df.tmp$sumFit, df.tmp$respS)
    aveFit.r2 <- calc.r2(df.tmp$aveFit, df.tmp$respS)
    maxAveFit.r2 <- calc.r2(df.tmp$maxAveFit, df.tmp$respS)

    df.out <- data.frame(sn = i, sumInt = coef(lm.sum)[1],aveInt = coef(lm.ave)[1],maxAveInt = coef(lm.maxAve)[1], sumFit.r2, aveFit.r2, maxAveFit.r2)
    df.subFits <- ch.rbind(df.subFits, df.out)
  }

  #get which model is the best fit based on r2 for each participant?
  df.subFits$bestFit <- colnames(df.subFits[,5:7])[apply(df.subFits[,5:7],1,which.max)]
  #Do sign tests to see if the model predicts better thnan chance. Use sign tests because the data is highly skewed
  sum.signTest <- SIGN.test(df.subFits$sumFit.r2, md = 0)
  ave.signTest <- SIGN.test(df.subFits$aveFit.r2, md = 0)
  maxAve.signTest <- SIGN.test(df.subFits$maxAveFit.r2, md = 0)
  maxAveVsAve.signTest <- SIGN.test(df.subFits$maxAveFit.r2, df.subFits$aveFit.r2)


  #Sum model.  Value = sum of the items + effect of category
    #now get the random effects plus the effect of category
    model.sum.a <- lmer(respS ~  -1 + offset(Sum_Median) + (1+ Sum_Median|sn), data=strong.exp, control=lmerControl(optimizer='Nelder_Mead'))
    #now creat the best fit by adding the model prediction to the predictor variable. This contains random and fixed effects
    strong.exp$fit.cond.sum <- predict(model.sum.a)
    #not create the best fit that contains only the fixed effects - use the fixed parameter from the model for catE
    strong.exp$fit.marg.sum <- strong.exp$Sum_Median
    #approximate the conditional and marginal r2
    model.sum.r2.overall <- calc.r2(strong.exp$fit.cond.sum, strong.exp$respS)
    model.sum.r2.marg <- calc.r2(strong.exp$fit.marg.sum, strong.exp$respS)
    model.sum.strong.BIC <- BIC(model.sum.a)
    model.sum.ranSlope.tt <- t.test(ranef(model.sum.a)$sn["Sum_Median"], mu=0)
    model.sum.ranInt.tt <- t.test(ranef(model.sum.a)$sn["(Intercept)"], mu=0)

  #Average model.  Value = average of the items + effect of category
    #now get the random effects plus the effect of category
    model.mean.a <- lmer(respS ~  -1 + offset(Mean_Median) + (1+ Mean_Median|sn), data=strong.exp, control=lmerControl(optimizer='Nelder_Mead'))
    #now creat the best fit by adding the model prediction to the predictor variable. This contains random and fixed effects
    strong.exp$fit.cond.ave <- predict(model.mean.a)
    #not create the best fit that contains only the fixed effects - use the fixed parameter from the model for catE
    strong.exp$fit.marg.ave <- strong.exp$Mean_Median
    #approximate the conditional and marginal r2
    model.ave.r2.overall <- calc.r2(strong.exp$fit.cond.ave, strong.exp$respS)
    model.ave.r2.marg <- calc.r2(strong.exp$fit.marg.ave, strong.exp$respS)
    model.ave.strong.BIC <- BIC(model.mean.a)
    model.ave.ranSlope.tt <- t.test(ranef(model.mean.a)$sn["Mean_Median"], mu=0)
    model.ave.ranInt.tt <- t.test(ranef(model.mean.a)$sn["(Intercept)"], mu=0)

  #MaxAve model.  Value = .5 max + .5 ave
    strong.exp$maxAve <- .5*strong.exp$Max_Median + .5*strong.exp$Mean_Median
    #now get the random effects plus the effect of category
    model.maxAve.a <- lmer(respS ~  -1 + offset(maxAve) + (1 + maxAve|sn), data=strong.exp, control=lmerControl(optimizer='Nelder_Mead'))
    #now creat the best fit by adding the model prediction to the predictor variable. This contains random and fixed effects
    strong.exp$fit.cond.maxAve <- predict(model.maxAve.a)
    #not create the best fit that contains only the fixed effects - use the fixed parameter from the model for catE
    strong.exp$fit.marg.maxAve <- strong.exp$maxAve
    #approximate the conditional and marginal r2
    model.maxAve.r2.overall <- calc.r2(strong.exp$fit.cond.maxAve, strong.exp$respS)
    model.maxAve.r2.marg <- calc.r2(strong.exp$fit.marg.maxAve, strong.exp$respS)
    model.maxAve.strong.BIC <- BIC(model.maxAve.a)
    model.maxAve.ranSlope.tt <- t.test(ranef(model.maxAve.a)$sn["maxAve"], mu=0)
    model.maxAve.ranInt.tt <- t.test(ranef(model.maxAve.a)$sn["(Intercept)"], mu=0)


### Run less strong models

    #### 4. Averaging Plus Model
    model.biasedMean <- lmer(respS ~ -1 + Mean_Median + Max_Median + (1 + Mean_Median + Max_Median|sn), data=exp.dat, control=lmerControl(optimizer='Nelder_Mead'))
      model.biasedMean.r2 <- r.squaredGLMM(model.biasedMean)
      model.biasedMean.BIC <- BIC(model.biasedMean)
      model.biasedMean.coef <- fixef(model.biasedMean)
      model.biasedMean.anova <- anova(model.biasedMean)
      model.biasedMean.r2.overall <- calc.r2(predict(model.biasedMean), exp.dat$respS)

### Sink all model outputs
      sink("Exploratory Analysis.txt",append=F)
        cat("\n******************* Exploratory Analysis *******************\n")
          cat("*** Total Sn:", length(unique(exp.dat$sn)), "\n\n")
          cat("\n\n ******* Summary Stats on Predictor Variables ******* \n")
          print(summary(exp.dat[,c("Sum_Median", "Mean_Median")]))
          cat("\n\nSD Sum: ", sd(exp.dat[,c("Sum_Median")]))
          cat("\nSD Average: ", sd(exp.dat[,c("Mean_Median")]))


          cat("\n\n ******* Default Model ******* \n")
          cat("**** Summary \n")
          print(summary(model.default))
          cat("\n**** ANOVA \n")
          print(model.default.anova)
          cat("\n**** Coefficients \n")
          print(model.default.coef)
          cat("\n**** Marginal and Conditional R2 \n")
          print(model.default.r2)
          cat("\n**** BIC \n")
          print(model.default.BIC)

          cat("\n\n******************* Strong Models *******************\n\n")
          cat("\n ******* SUM Model ******* \n")
          cat("\n Approximate Overall R2: ", model.sum.r2.overall)
          cat("\n Approximate Marginal R2: ", model.sum.r2.marg)
          cat("\n BIC: ", model.sum.strong.BIC)
          cat("\n\n Test whether the Random Effects Intercepts == 0:  T-Test,  Mu = 0 \n")
          cat("If the model fits well, the Random Effects should have a mean around 0 \n")
          print(model.sum.ranInt.tt)
          cat("\n\n Test whether the Random Effects Slopes == 0:  T-Test,  Mu = 0 \n")
          cat("If the model fits well, the Random Effects should have a mean around 0 \n")
          print(model.sum.ranSlope.tt)
          cat("\n\n When we fit the model to each individual participant, do the r2's differ from 0? \n")
          cat("If the model fits well, the r2's should be greater than 0. \n")
          cat("We do a non-parametric Sign Test,  Mu = 0, because there are some outliers \n")
          print(sum.signTest)
          cat("\n\n ******* Average Model ******* \n")
          cat("\n Approximate Overall R2: ", model.ave.r2.overall)
          cat("\n Approximate Marginal R2: ", model.ave.r2.marg)
          cat("\n BIC: ", model.ave.strong.BIC)
          cat("\n\n Test whether the Random Effects Intercepts == 0:  T-Test,  Mu = 0 \n")
          cat("If the model fits well, the Random Effects should have a mean around 0 \n")
          print(model.ave.ranInt.tt)
          cat("\n\n Test whether the Random Effects Slopes == 0:  T-Test,  Mu = 0 \n")
          cat("If the model fits well, the Random Effects should have a mean around 0 \n")
          print(model.ave.ranSlope.tt)
          cat("\n\n When we fit the model to each individual participant, do the r2's differ from 0? \n")
          cat("If the model fits well, the r2's should be greater than 0. \n")
          cat("We do a non-parametric Sign Test,  Mu = 0, because there are some outliers \n")
          print(ave.signTest)
          cat("\n\n ******* Max Ave Average Model ******* \n")
          cat("\n Approximate Overall R2: ", model.maxAve.r2.overall)
          cat("\n Approximate Marginal R2: ", model.maxAve.r2.marg)
          cat("\n BIC: ", model.maxAve.strong.BIC)
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

          cat("\n\n Does the Max Ave model explain significantly more variance than the Average model in individual participants? \n")
          cat("We do a non-parametric Sign Test,  Max Ave Vs Ave Sign Test, because there are some outliers \n")
          print(maxAveVsAve.signTest)


          cat("\n\n******************* Weaker Models *******************\n\n")

          cat("\n\n ******* Biased Average Model ******* \n")
          cat("***** Average Plus Max ***** \n\n")
          cat("**** Summary \n")
          print(summary(model.biasedMean))
          cat("\n**** ANOVA \n")
          print(model.biasedMean.anova)
          cat("\n**** Coefficients \n")
          print(model.biasedMean.coef)
          cat("\n**** Marginal and Conditional R2 \n")
          print(model.biasedMean.r2)
          cat("\n**** Overall R2 \n")
          print(model.biasedMean.r2.overall)
          cat("\n**** BIC \n")
          print(model.biasedMean.BIC)

        sink(NULL)

#### Model Graphs

        strong.exp$biasedMeanFit <- predict(model.biasedMean)

        df.fit <- strong.exp %>% group_by(Mean_Median_Round) %>% summarise(mRespS = mean(respS),n = length(respS), biasedMeanFit = mean(biasedMeanFit), mSum = mean(Sum_Median_Round), mMedian = mean(Mean_Median), mMax = mean(Max_Median), strongFit.maxAve = mean(fit.marg.maxAve), strongFit.sum = mean(fit.marg.sum), strongFit.ave = mean(fit.marg.ave))
        df.group <- strong.exp %>% group_by(groupSize, Mean_Median_Round) %>% summarise(mRespS = mean(respS),n = length(respS), biasedMeanFit = mean(biasedMeanFit), mSum = mean(Sum_Median_Round), mMedian = mean(Mean_Median), mMax = mean(Max_Median), strongFit.maxAve = mean(fit.marg.maxAve), strongFit.sum = mean(fit.marg.sum), strongFit.ave = mean(fit.marg.ave))

        pdf("exploratory Predictor x Response.pdf")
          df.tmp <- df.group[df.group$n>minN,]
          for(i in 6:length(df.group)) {
            maxX <- ifelse(max(df.tmp[[i]]) < 4, 4, 20)
            plot(df.tmp$mRespS ~ df.tmp[[i]], ylim=c(0,5), xlim=c(0,maxX),pch=16,las=1,frame=F,tck = 0,ylab="Group Value Response",xlab=colnames(df.tmp[i]), main = colnames(df.tmp[i]))
            abline(0,1, col=gray(.5))
          }

          df.tmp <- df.group[df.group$n>minN,]
          for(i in 6:length(df.group)) {
            maxX <- ifelse(max(df.tmp[[i]]) < 4, 4, 20)
            plot(NULL, ylim=c(0,5), xlim=c(0,maxX),pch=16,las=1,frame=F,tck = 0,ylab="Group Value Response",xlab=colnames(df.tmp[i]), main = colnames(df.tmp[i]))
            points(df.tmp[df.tmp$groupSize==5,][["mRespS"]] ~ df.tmp[df.tmp$groupSize==5,][[i]],col=gray(0.25), pch=16)
            points(df.tmp[df.tmp$groupSize==4,][["mRespS"]] ~ df.tmp[df.tmp$groupSize==4,][[i]],col=gray(0.40), pch=16)
            points(df.tmp[df.tmp$groupSize==3,][["mRespS"]] ~ df.tmp[df.tmp$groupSize==3,][[i]],col=gray(0.55), pch=16)
            points(df.tmp[df.tmp$groupSize==2,][["mRespS"]] ~ df.tmp[df.tmp$groupSize==2,][[i]],col=gray(0.70), pch=16)
            abline(0,1, col=gray(.5))
          }

          with(df.group[df.group$n>minN,], plot(NULL, ylim=c(0,10), xlim=c(0,4),pch=16,las=1,frame=F,tck = 0,ylab="Prediction",xlab="Average Group Value", main="raw summary stats"))
          with(df.fit[df.fit$n > 10,], points(mRespS ~ Mean_Median_Round, pch=16, col="grey85"))
          with(df.fit[df.fit$n > 10,], lines(mMedian ~ Mean_Median_Round, pch=109, col="grey20", lty = 3))
          with(df.fit[df.fit$n > 10,], lines(mMax ~ Mean_Median_Round, pch=115, col="grey50", lty = 3))
          with(df.fit[df.fit$n > 10,], lines(mSum ~ Mean_Median_Round, pch=109, col="grey75", lty = 1))
          abline(0,1, lty=5)

          with(df.group[df.group$n>minN,], plot(NULL, ylim=c(-5,5), xlim=c(0,4),pch=16,las=1,frame=F,tck = 0,ylab="Prediction",xlab="Average Group Value", main="Stronger Models"))
          with(df.fit[df.fit$n > 10,], points(mRespS ~ Mean_Median_Round, pch=16, col="grey85"))
          with(df.fit[df.fit$n > 10,], lines(strongFit.sum ~ Mean_Median_Round, pch=83, col="grey25", lty = 3))
          with(df.fit[df.fit$n > 10,], lines(strongFit.ave ~ Mean_Median_Round, pch=109, col="grey50", lty = 3))
          with(df.fit[df.fit$n > 10,], lines(strongFit.maxAve ~ Mean_Median_Round, pch=109, col="grey75", lty = 3))
          abline(0,1, lty=5)

          with(df.group[df.group$n>minN,], plot(NULL, ylim=c(0,15), xlim=c(0,4),pch=16,las=1,frame=F,tck = 0,ylab="Prediction",xlab="Average Group Value", main="Strong Model Predictions"))
          with(exp.dat[exp.dat$groupSize == 2, ], lines(Sum_Median ~ Mean_Median, type = "l", col=gray(0.70), lty = 1, lwd=2))
          with(exp.dat[exp.dat$groupSize == 3, ], lines(Sum_Median ~ Mean_Median, col=gray(0.55), lty = 1, lwd=2))
          with(exp.dat[exp.dat$groupSize == 4, ], lines(Sum_Median ~ Mean_Median, col=gray(0.40), lty = 1, lwd=2))
          with(exp.dat[exp.dat$groupSize == 5, ], lines(Sum_Median ~ Mean_Median, col=gray(0.25), lty = 1, lwd=2))
          abline(0,1, lty=3, lwd=2)
          abline(0.5,1, lty=5, lwd=2)

          df.fit <- df.fit[order(df.fit$Mean_Median_Round),]
          with(df.group[df.group$n>minN,], plot(NULL, ylim=c(0,15), xlim=c(0,4),pch=16,las=1,frame=F,tck = 0,ylab="Prediction",xlab="Average Group Value", main="Strong Model Predictions"))
          with(df.group[df.group$groupSize == 2 & df.group$n>minN, ], points(mRespS ~ Mean_Median_Round, pch=16, col=gray(0.70)))
          with(df.group[df.group$groupSize == 3 & df.group$n>minN, ], points(mRespS ~ Mean_Median_Round, pch=16, col=gray(0.55)))
          with(df.group[df.group$groupSize == 4 & df.group$n>minN, ], points(mRespS ~ Mean_Median_Round, pch=16, col=gray(0.40)))
          with(df.group[df.group$groupSize == 5 & df.group$n>minN, ], points(mRespS ~ Mean_Median_Round,, pch=16, col=gray(0.25)))
          with(df.fit[df.fit$n>minN, ], lines(biasedMeanFit ~ Mean_Median_Round, lwd=2, col=gray(0.25)))
          abline(0,1, lty=3, lwd=2)

      dev.off()



  ######## IQR analysis ##############
    exp.dat$mean_25th<- rowMeans(cbind(exp.dat$tValueQ251,exp.dat$tValueQ252,exp.dat$tValueQ253,exp.dat$tValueQ254,exp.dat$tValueQ255),na.rm = T)
    exp.dat$mean_75th<- rowMeans(cbind(exp.dat$tValueQ751,exp.dat$tValueQ752,exp.dat$tValueQ753,exp.dat$tValueQ754,exp.dat$tValueQ755),na.rm = T)
    exp.dat$IQR <- rowDiffs(cbind(exp.dat$mean_25th, exp.dat$mean_75th),na.rm=T)
    exp.dat <- transform(exp.dat, Z = mean_25th + 0.5 * IQR)
    exp.dat <- transform(exp.dat, QR = (respS - Z)/IQR)

    #Histogram collapsed across group sizes
    pdf("QR Histogram for All Probes.pdf")
      b1 <- seq(-39.5, 39.5, 1)
       with(exp.dat, hist(QR,tck=0,las=1,cex.axis=0.9,main="", xlim=c(-10,10),ylim=c(0,0.6),ylab="Percentage of Responses",xlab="Standardized IQR",breaks = b1, freq=F))
      abline(v= -.5,col="red", lwd=2, lty=3)
      abline(v= 0.5,col="red", lwd=2, lty=3)
    dev.off()

      QR_df <- data.frame(exp.dat %>% group_by(groupSize) %>% summarise(
        resp_SD=sd(respS),
        Q25 = (quantile(QR, .25, na.rm=T)),
        Q75 = (quantile(QR, .75, na.rm=T)),
        Q50 = (quantile(QR, .50, na.rm=T)),
        QR_range = (Q75 - Q25)))
    write.table(QR_df,"QR_quantiles.txt",row.names=F,quote=F,sep="\t")


##### Boxplot of estimates by groupsize and the predictions of sum and average ####
    op <- par(mfrow=c(2,2), bty="n", mar=c(6,4,4,1), oma=c(1,1,1,1))
      with(exp.dat, boxplot(respS ~ groupSize, main = "Data", ylim=c(0,15), ylab="Psychological Value", xlab="Group Size"))
      with(exp.dat, boxplot(predict(model.maxAve.a) ~ groupSize, main = "Predicted\nBiased Average", ylim=c(0,15), ylab="Psychological Value", xlab="Group Size"))
      with(exp.dat, boxplot(Sum_Median ~ groupSize, main = "Predicted\nSum", ylim=c(0,15), ylab="Psychological Value", xlab="Group Size"))
      with(exp.dat, boxplot(Mean_Median ~ groupSize, main = "Predicted\nAverage", ylim=c(0,15), ylab="Psychological Value", xlab="Group Size"))
    dev.copy(pdf,"exploratory data groupSize Boxplots.pdf", width=8, height=8)
    dev.off()
    par(op)

    sink("Exploratory Analysis.txt",append=T)
      cat("\n******************* SD Analysis *******************\n")
      print(QR_df)
    sink(NULL)
