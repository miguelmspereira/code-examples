#Get dosage data from a subset of SNPs retrieved from the the UK Biobank imputed data files
#01Aug2017

require(data.table)

#Working directory
setwd('/Users/miguelmspereira/Box Sync/Congenica/scripts/')


#Subject ID's
id<-fread('ukb1913_imp_chr1_v2_s487398.sample')
eid<-id$ID_1[-1] #vector with the IDs
rm(id)


#GEN file with the subset of SNPs in 3-column format (this file is for 42 SNPs)
sub<-fread('subset.gen')

#42 rows and 152,249 subjects - corresponds to (456753-6)/3 colums
dim(sub)

#Check data
sub[1:10,1:10]


#----------------------------------------------------------
#File with SNP general info
info<-subset(sub,select=1:6)
write.csv(info,'snpinfo.csv',quote=F,row.names=F)
#----------------------------------------------------------

#Converting to dosage data
b<-subset(sub,select=seq(8,ncol(sub),by=3)) #Takes the 2nd column of each subject
c<-subset(sub,select=seq(9,ncol(sub),by=3)) #Takes the 3rd column of each subject

#dim(b)==dim(c) #Checks if b and c are the same size


#----------------------------------------------------------
#Column names for the final file
#SNP RSID's
cols<-sub$V3
rm(sub) #removed to save memory


#----------------------------------------------------------
#Splits the b and c data matrices in 3 parts - for matrix transposition
b1<-subset(b,select=1:50000)
b2<-subset(b,select=50001:100000)
b3<-subset(b,select=100001:ncol(b))


c1<-subset(c,select=1:50000)
c2<-subset(c,select=50001:100000)
c3<-subset(c,select=100001:ncol(c))


#(ncol(b1)+ncol(b2)+ncol(b3))==ncol(b) #checks if the columns add up to the total
#(ncol(c1)+ncol(c2)+ncol(c3))==ncol(c) #checks if 
rm(b) #for memory
rm(c) #for memory


#Calculating dosage data - done 3 times for each of the 3 matrices
dosages1<-b1+2*c1
rm(b1)
rm(c1)

dosages2<-b2+2*c2
rm(b2)
rm(c2)

dosages3<-b3+2*c3
rm(b3)
rm(c3)



#Matrix transposition
tdosages1<-dosages1[, data.table(t(.SD), keep.rownames=F)]
rm(dosages1)

tdosages2<-dosages2[, data.table(t(.SD), keep.rownames=F)]
rm(dosages2)

tdosages3<-dosages3[, data.table(t(.SD), keep.rownames=F)]
rm(dosages3)



#Joining the 3 transposed matrices to get the final matrix
tdosages<-rbind(tdosages1,tdosages2,tdosages3)
rm(tdosages1)
rm(tdosages2)
rm(tdosages3)

#Adds column names to the final matrix with dosage data
colnames(tdosages)<-cols
head(colnames(tdosages)) #to see if colnames are OK
dim(tdosages)  #checks size
nrow(tdosages)==(ncol(sub)-6)/3 #compares with original size

#Final data matrix
dos<-as.data.table(cbind(eid,tdosages)) #adds the subject IDs from the sample file
rm(tdosages)

#To save the files
save(dos,file='dosages_newrelease.RData')

#Optional - save as CSV - generates big files if for a big subset of SNPs
write.csv(dos,file='dosages_newlrelease.csv',row.names=F,quote=F)






