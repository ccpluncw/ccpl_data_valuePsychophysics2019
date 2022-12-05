## Mode function
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

dat.tmp <- read.table("dat.trim.txt",sep="\t", header=T, quote="\"")

# remove single, common items
dat.trim <- dat.tmp[dat.tmp$groupSize > 1, ]

  exp.dat <- NULL
  confirm.dat <- NULL
  for(i in 2:5) {
    #get the data for groupSize == i
    df.tmp<-dat.trim[dat.trim$groupSize == i,]
    #identify the subjects in that group
    tmp2<-unique(df.tmp$sn)
    #count the subjects in that group
    tmp2length<-length(tmp2)
    #adjust seed for each sample
    set.seed(seed + i)
    #randomly round up or down so that we do not bias one group with more or less participants
    if(sample(c(0,1),1) == 0) {
      #sample half the subjects and round up
      sns<-sample(tmp2,ceiling(tmp2length/2))
    } else {
      #sample half the subjects and round down
      sns<-sample(tmp2,floor(tmp2length/2))
    }
    #get the exploratory subjects
    tmp.out.exp<-df.tmp[df.tmp$sn %in% sns,]
    #the remaining subjects are confirmatory
    tmp.out.conf<-df.tmp[!(df.tmp$sn %in% sns),]
    exp.dat <- ch.rbind(exp.dat, tmp.out.exp)
    confirm.dat <- ch.rbind(confirm.dat, tmp.out.conf)
  }

    write.table(exp.dat,"exp.dat.txt",row.names=F,quote=F,sep="\t")
    write.table(confirm.dat,"confirm.dat.txt",row.names=F,quote=F,sep="\t")
