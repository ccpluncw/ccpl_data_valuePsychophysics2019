## Mode function
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

dat.trim <- read.table("dat.trim.txt",sep="\t", header=T, quote="\"")

###Randomly select half of participants from each group size
  #Create dataframes for participants by each group size
  groupSize_sn<-data.frame(dat.trim %>% group_by(sn) %>% summarise(Group_mode=Mode(groupSize)))
  sn_2<-groupSize_sn[groupSize_sn$Group_mode==2,"sn"]
  sn_3<-groupSize_sn[groupSize_sn$Group_mode==3,"sn"]
  sn_4<-groupSize_sn[groupSize_sn$Group_mode==4,"sn"]
  sn_5<-groupSize_sn[groupSize_sn$Group_mode==5,"sn"]
  df_2.tmp<-dat.trim[dat.trim$sn %in% sn_2,]
  df_3.tmp<-dat.trim[dat.trim$sn %in% sn_3,]
  df_4.tmp<-dat.trim[dat.trim$sn %in% sn_4,]
  df_5.tmp<-dat.trim[dat.trim$sn %in% sn_5,]

  #Randomly select half of participants from each group size
    tmp2<-unique(df_2.tmp$sn)
    tmp2length<-length(tmp2)
    df_2<-sample(tmp2,tmp2length/2)

    tmp3<-unique(df_3.tmp$sn)
    tmp3length<-length(tmp3)
    df_3<-sample(tmp3,tmp3length/2)

    tmp4<-unique(df_4.tmp$sn)
    tmp4length<-length(tmp4)
    df_4<-sample(tmp4,tmp4length/2)

    tmp5<-unique(df_5.tmp$sn)
    tmp5length<-length(tmp5)
    df_5<-sample(tmp5,tmp5length/2)

    #Get data from selected participants for exploratory dataset
    tmp_2<-dat.trim[dat.trim$sn %in% df_2,]
    tmp_3<-dat.trim[dat.trim$sn %in% df_3,]
    tmp_4<-dat.trim[dat.trim$sn %in% df_4,]
    tmp_5<-dat.trim[dat.trim$sn %in% df_5,]
    #bind data to create exploratory data frame then write to a file
    exp.dat<-rbind(tmp_2,tmp_3,tmp_4,tmp_5)
    write.table(exp.dat,"exp.dat.txt",row.names=F,quote=F,sep="\t")

    #create confirmatory dataframe with remaining data
    confirm_2<-df_2.tmp[!(df_2.tmp$sn %in% df_2),]
    confirm_3<-df_3.tmp[!(df_3.tmp$sn %in% df_3),]
    confirm_4<-df_4.tmp[!(df_4.tmp$sn %in% df_4),]
    confirm_5<-df_5.tmp[!(df_5.tmp$sn %in% df_5),]
    confirm.dat<-rbind(confirm_2,confirm_3,confirm_4,confirm_5)
    write.table(confirm.dat,"confirm.dat.txt",row.names=F,quote=F,sep="\t")
