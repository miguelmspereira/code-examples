#Get dosage data from a subset of SNPs retrieved from the the UK Biobank imputed data files
#01Aug2017

require(data.table)

setwd('/export101/work/mmd14/shared/imputed/generalscripts')


#Subject ID's
id<-fread('/export101/work/mmd14/shared/imputed/ukb1913_imp_chr1_v2_s487398.sample')
eid<-id$ID_1[-1]


sub<-fread('subset.gen')
info<-subset(sub,select=1:6)
write.csv(info,'snpinfo.csv',quote=F,row.names=F)


b<-subset(sub,select=seq(8,ncol(sub),by=3))
c<-subset(sub,select=seq(9,ncol(sub),by=3))

cols<-sub$V3
rm(sub)



b1<-subset(b,select=1:150000)
b2<-subset(b,select=150001:300000)
b3<-subset(b,select=300001:ncol(b))


c1<-subset(c,select=1:150000)
c2<-subset(c,select=150001:300000)
c3<-subset(c,select=300001:ncol(c))


#(ncol(b1)+ncol(b2)+ncol(b3))==ncol(b) #only a check
#(ncol(c1)+ncol(c2)+ncol(c3))==ncol(c) #only a check
rm(b)
rm(c)

dosages1<-b1+2*c1
rm(b1)
rm(c1)

dosages2<-b2+2*c2
rm(b2)
rm(c2)

dosages3<-b3+2*c3
rm(b3)
rm(c3)




tdosages1<-dosages1[, data.table(t(.SD), keep.rownames=F)]
rm(dosages1)

tdosages2<-dosages2[, data.table(t(.SD), keep.rownames=F)]
rm(dosages2)

tdosages3<-dosages3[, data.table(t(.SD), keep.rownames=F)]
rm(dosages3)



tdosages<-rbind(tdosages1,tdosages2,tdosages3)
rm(tdosages1)
rm(tdosages2)
rm(tdosages3)

colnames(tdosages)<-cols
#head(colnames(tdosages)) #to see if colnames are OK

#Final data matrix
dos<-as.data.table(cbind(eid,tdosages))


#To save the files
save(dos,file='dosages_newrelease.RData')
rm(tdosages)

write.csv(dos,file='dosages_newlrelease.csv',row.names=F,quote=F)






