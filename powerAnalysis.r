library(chutils)
library(lmerTest)
require(MuMIn)

#### Power Assumptions
#this is the random slope SD
sd1 <- 0.5
#This is the added error to each response SD
sd2 <- 4
#This is the significance alpha
alpha <- 0.05

### Task Parameters
#we will have about 300 trials per participant
nTrials <- 300
#we will have 4 group sizes
groupSize <- c(2,3,4,5)

#### simulation Parameters
#Number of bootstrap loops
loops <- 100
#vary number of participants per group.
nSubs <- seq(1, 11, 2)


df.power <- NULL
#For each number of subjects
for(nSs in nSubs) {
  df.out <- NULL
  #run a bootstrap
  for(run in 1:loops) {
    #generate simulated data
    df.all <- NULL
    #For each group
    subID <- 1
    for(i in groupSize) {
      #generate the data for each subject
      for(j in 1:nSs) {
        #log the group size and subject number
        grp <- rep(i, nTrials)
        sn <- rep(subID, nTrials)

        #get the random slope
        slp <- 1 + rnorm(1, 0, sd1)
        #generate random values for the stimuli
        val <- rnorm(nTrials, 2, 0.6)
        #generate random error
        error <- rnorm(nTrials, 0, sd2)
        #generate response
        resp <- slp*grp * val + error
        #create dataframe
        df.tmp <- data.frame(sn = sn, size = grp, value = val, resp = resp, slope = slp, error = error)
        df.all <- ch.rbind(df.all, df.tmp)
        #increase subject ID
        subID <- subID + 1
      }
    }

    #run mixed model on the data
    df.all$sn <- as.factor(df.all$sn)
    fit <- lmer(resp ~ size + (-1+size|sn), data = df.all)
    #extract marginal r2
    r2marg <- r.squaredGLMM(fit)[1]
    #extract p-value of critical parameter
    pVal <- anova(fit)[["Pr(>F)"]]
    #determine if it is significant
    sig <- ifelse (pVal < alpha, 1, 0)
    #create dataframe of results
    df.tmp1 <- data.frame(runNum = run, r2marginal = r2marg, pValue = pVal, significant = sig)
    df.out <- ch.rbind(df.out, df.tmp1)
  }
  #create a summary of the bootstrap results for each nSubs
  r2margMean <- mean(df.out$r2marginal, na.rm=T)
  r2margSD <- sd(df.out$r2marginal, na.rm=T)
  power <- mean(df.out$significant, na.rm=T)
  df.tmp2 <- data.frame(numSubsPerGroup = nSs, r2marginalMean = r2margMean, r2marginalSD = r2margSD, power = power)
  df.power <- ch.rbind(df.power, df.tmp2)
}

sink("powerAnalysisOutput.txt")
  print(df.power)
sink(NULL)
