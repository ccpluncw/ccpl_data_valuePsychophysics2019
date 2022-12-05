library(chutils)
library(dplyr)
library(matrixStats)

source("combineValueDatasets4.r")

mainDir <- getwd()

rawDataFile <- "data.grouping.txt"
probeCategoryFile <- "probecategories.txt"

#read in rawdata
data.grouping.tmp <- read.table(rawDataFile,sep="\t", header=T, quote="\"")
#read in probe categories
dat.cat<-read.table(probeCategoryFile,sep="\t", header=T, quote="\"")
##Read in individual probe stats from previous experiment
stats.temp<-read.table("ItemStats1.txt",sep="\t", header=T, quote="\"")
stats.temp2<-read.table("ItemStats2.txt",sep="\t", header=T, quote="\"")
stats.temp3<-read.table("ItemStats3.txt",sep="\t", header=T, quote="\"")
stats.temp4<-read.table("ItemStats4.txt",sep="\t", header=T, quote="\"")
stats.temp5<-read.table("ItemStats5.txt",sep="\t", header=T, quote="\"")

subDir <- "analysis"
ch.newDir(mainDir,subDir)
  filterSsBasedOnRT <- FALSE

  source("../powerAnalysis.r")
  source("../grouping_DataClean.r")
  #set seed to the randomization that splits the data into confirmatory and exploratory datasets can be replicated
  seed <- 952003
  source("../splitDat.r")
  source("../exploratory_analysis.r")
  source("../confirmatory_analysis.r")
setwd(mainDir)
