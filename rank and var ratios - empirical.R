#Data from 6 true SNPs associated with BMI and 24 random SNPs from the ECRHS imputed dataset + all the SNPs in the same LD blocks with them (2,839 unique SNPs from 30 LD blocks)
#This Script joins the data with the phenotype data from ECRHS and does the analysis of the data
#Revised for thesis on 15Dec2017

#Working directory
setwd('/Users/miguelmspereira/Box Sync/Congenica/scripts/')

#Phenotype data
fulldata0<-read.table('fulldata.csv',header=T,stringsAsFactors=T,sep=',')
fulldata<-fulldata0[,-1]
fulldata[1:10,1:10]
dim(fulldata)

#__________________________________________________________________________________________________________________
#Prior SNP knowledge
prior0<-read.table('priorAnswers.csv', header=T, stringsAsFactors=F, sep=',')
prior<-prior0[,-1]
dim(prior) #2,614 SNPs
head(prior)

#__________________________________________________________________________________________________________________
#LD block data
ldblock0<-read.table('snp6_24_ldblock.csv',header=T,stringsAsFactors=F,sep=',')
ldblock<-ldblock0[,-1]
head(ldblock)
dim(ldblock) #2,866 - includes duplicated SNPs

#LD block matrix without duplicate SNPs (it is important to do this because some SNPs will have more than one Alt.Allele)
length(which(duplicated(ldblock$snp)==F)) #number of SNPs that are not duplicated
ldblockDup<-ldblock[-which(duplicated(ldblock$snp)),]
dim(ldblockDup)
head(ldblockDup$coords)


#__________________________________________________________________________________________________________________
#Split genomic coordinates to sort data
coord.list<-strsplit(ldblockDup$coords,split=':')
head(coord.list)

#Splitting SNP coordinates for ordering
posSplit<-unlist(strsplit(ldblockDup[,2],split=':'))
chr<-seq(from=2,to=length(posSplit),by=3) #Gets the chr
snpPos<-seq(from=3,to=length(posSplit),by=3) #Gets the SNP position

#Final LD block matrix with sorted SNPs (by Chromossome and Genomic coordinate)
ldblock.final0<-cbind(ldblockDup$snp,ldblockDup$coords,ldblockDup$ldblock,posSplit[chr],posSplit[snpPos],ldblockDup$refAll,ldblockDup$altAll)
colnames(ldblock.final0)<-c('snp','coords','ldblock','chr','position','Ref','Alt')
head(ldblock.final0)

#ORDERING the SNPs by chr and position
ldblock.final<-ldblock.final0[order(as.numeric(ldblock.final0[,4]),as.numeric(ldblock.final0[,5])),]
head(ldblock.final)
unique(ldblock.final[,4]) #checks that the chromossomes are in the correct order 

#Assigning a number to each LDblock - index facilitate operations
ind<-seq(from=1,to=length(unique(ldblock.final[,3])))
ind.ldblock<-cbind(ind,unique(ldblock.final[,3]))
ind.ldblock


#Matching LD block vector with the assigned numbers
index.ldblock<-rep(99,times=nrow(ldblock.final))
for(i in 1:nrow(ldblock.final)){
  index.ldblock[i]<-as.numeric(ind.ldblock[which(ind.ldblock[,2]==ldblock.final[i,3]),1])
}
index.ldblock

#Joining the index to the matrix
ldblock.final1<-cbind(ldblock.final,index.ldblock)
head(ldblock.final1)


#Subseting the matrix to only the SNPs with prior knowledge n=2,614
ldblock.final2<-ldblock.final1[which(ldblock.final1[,1] %in% prior$X.1),]
dim(ldblock.final2)
head(ldblock.final2)


#__________________________________________________________________________________________________________________
#Reordering columns in the data
orderedSnp<-ldblock.final2[,1]
head(orderedSnp)
which(match(colnames(fulldata),orderedSnp)!='NA')
sortedData<-cbind(fulldata[,1:15],fulldata[orderedSnp])

colnames(fulldata[,1:15])
head(colnames(fulldata[orderedSnp]))


#__________________________________________________________________________________________________________________
#Calculating BMI - outcome measure
bmi<-sortedData$weight/(sortedData$height^2)
hist(bmi)
summary(bmi)

bmidata0<-cbind(bmi,sortedData)
bmidata0[1:10,1:20]



#__________________________________________________________________________________________________________________
#Standardizing the data
bmidata<-as.data.frame(scale(bmidata0))
attach(bmidata)
colnames(bmidata)[16:30] #check to see if RSID's are OK and start where they are supposed to
bmidata[1:10,1:10] #quick look at the data

summary(bmidata$PC4) #To check if it is standardized


#__________________________________________________________________________________________________________________
#Linear regression results (output matrix does not include information regarding centre)
#Standard analysis - Performs linear regression to estimate single SNP effects and adjusts the for the covariates of interest

#Output matrix
output<-data.frame(matrix(nrow=nrow(ldblock.final2),ncol=23))
colnames(output)<-c('snp','group','MAF','ref.allele','alt.allele','snp.beta','snp.ste','snp.pval','age','age.ste','age.pval','sex','sex.ste','sex.pval','pc1','pc1.ste','pc1.pval','pc2','pc2.ste','pc2.pval','pc3','pc3.ste','pc3.pval')
head(output)

output$snp<-orderedSnp
output$group<-ldblock.final2[,8]
output$MAF<-snpinfo$MAF
output$ref.allele<-ldblock.final2[,6]
output$alt.allele<-ldblock.final2[,7]


for(i in 1:nrow(output)){
  #SNP
  output[i,6]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[2,1]
  output[i,7]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[2,2]
  output[i,8]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[2,4]
  
  #age
  output[i,9]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[3,1]
  output[i,10]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[3,2]
  output[i,11]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[3,4]
  
  #Sex
  output[i,12]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[4,1]
  output[i,13]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[4,2]
  output[i,14]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[4,4]
  
  #PC1
  output[i,15]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[5,1]
  output[i,16]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[5,2]
  output[i,17]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[5,4]
  
  #PC2
  output[i,18]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[6,1]
  output[i,19]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[6,2]
  output[i,20]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[6,4]
  
  #PC3
  output[i,21]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[7,1]
  output[i,22]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[7,2]
  output[i,23]<-coef(summary(lm(bmi~bmidata[,16+i]+age+sex.x+PC1+PC2+PC3+as.character(centre))))[7,4]
}

head(output)

#Writes the output to file
write.csv(output,'classical analysis ecrhs (standardized).csv')

#__________________________________________________________________________________________________________________
#Standard analysis without covariates
output.noCovar<-data.frame(matrix(nrow=nrow(ldblock.final2),ncol=8))
colnames(output.noCovar)<-c('snp','group','MAF','ref.allele','alt.allele','snp.beta','snp.ste','snp.pval')
head(output.noCovar)

output.noCovar$snp<-orderedSnp
output.noCovar$group<-ldblock.final2[,8]
output.noCovar$MAF<-snpinfo$MAF
output.noCovar$ref.allele<-ldblock.final2[,6]
output.noCovar$alt.allele<-ldblock.final2[,7]


for(i in 1:nrow(output.noCovar)){
  #SNP
  output.noCovar[i,6]<-coef(summary(lm(bmi~bmidata[,16+i])))[2,1]
  output.noCovar[i,7]<-coef(summary(lm(bmi~bmidata[,16+i])))[2,2]
  output.noCovar[i,8]<-coef(summary(lm(bmi~bmidata[,16+i])))[2,4]
}

head(output)
head(output.noCovar)

#Writes the output to file
write.csv(output.noCovar,'classical analysis ecrhs no covariates (only SNPs and standardized).csv')


#---------------------------------------------------
#Average ranking of the true LD blocks

#__________________
#Support function to get the standard error
std <- function(x) sd(x)/sqrt(length(x))
#__________________

order.output.noCovar<-output.noCovar[order(output.noCovar$snp.pval),]
head(order.output.noCovar)
order.ldblock.classical<-as.numeric(unique(order.output.noCovar$group))
order.ldblock.classical

#Average ranking of the true LD blocks (no covariates)
rank.classical<-mean(which(order.ldblock.classical %in% c(17,19,24,25,28,29)))
rank.classical.sd<-std(which(order.ldblock.classical %in% c(17,19,24,25,28,29)))

#Average ranking of the true LD blocks (WITH covariates)
order.ldblock.classical.withCovar<-as.numeric(unique(output[order(output$snp.pval),2]))
rank.classical.withCovar<-mean(which(order.ldblock.classical.withCovar %in% c(17,19,24,25,28,29)))
rank.classical.withCovar.sd<-std(which(order.ldblock.classical.withCovar %in% c(17,19,24,25,28,29)))



#----------------------------------------------------------------------------------
#Calculating the prior knowledge score - qSum
head(prior)
dim(prior)
qSum<-rowSums(prior[,2:11])
for(i in 1:length(qSum)){
  if(qSum[i]>4) qSum[i]=4
}
table(qSum)

#With weights - optional
#prior2<-data.frame(cbind(1.4*prior$q1,4*prior$q2,7.7*prior$q3,9.6*prior$q4,1.2*prior$q5,2.4*prior$q6,5.7*prior$q7,9.5*prior$q12,3.4*prior$q13,2.5*prior$q15))
#head(prior2)


#prior3<-data.frame(cbind(log(1.4)*prior$q1,log(4)*prior$q2,log(7.7)*prior$q3,log(9.6)*prior$q4,log(1.2)*prior$q5,log(2.4)*prior$q6,log(5.7)*prior$q7,log(9.5)*prior$q12,log(3.4)*prior$q13,log(2.5)*prior$q15))
#qSum<-rowSums(prior3)
#table(qSum)



#----------------------------------------------------------------------------------
#Bayesian modelling using bglm() - Yi et al. 2011
library(BhGLM) #This can be download from my github account

bmidata.bglm<-bmidata[,-(2:16)]

#Genotype data (uses function round() because other functions in the BhGLM package do not accept 'dosage' data)
geno<-round(bmidata0[,17:ncol(bmidata0)])
geno[1:10,1:10]
head(colnames(geno))

#Allele frequencies - check the genotypes frequencies
freq<-geno.freq(geno,verbose=F)


#Design matrix for the SNPs (the function make.main() transforms the matrix appropriatelly for bglm())
x.m<-make.main(geno=geno,model='additive')
x.m[1:10,1:10]

#-----------------
identical(colnames(bmidata.bglm)[-1],prior[,1]) #little check
#-----------------
#Groups
groups<-list()
for (i in 1:max(index.ldblock)){
  groups[[i]]<-ldblock.final2[which(as.numeric(ldblock.final2[,8])==i),1]
}
str(groups)

ldblock.groupCorrespondence<-cbind(ldblock.final1[-which(duplicated(ldblock.final1[,8])),3],paste('G',unique(ldblock.final1[,8]),sep=''))
ldblock.groupCorrespondence




#---------------------------------------------------------------------------------------
#No inclusion of prior knowledge
#Shrinkage values to be tested
sca<-c(100000,10000,5,2.5,1,0.5,0.1,0.01,0.001,0.0001,0.00001,0.000001) #shrinkage parameters to be tested

rank.all<-cbind(sca,rep(99,times=length(sca)))

empirical.effects<-c(99,99)
empirical.effects0<-c(99,99)


ptm <- proc.time() #to get a time estimate
for(i in 1:length(sca)){
  mod1<-bglm(bmi ~ ., data=bmidata.bglm, family = gaussian, prior = "t", prior.mean=0, mean.update=F,prior.scale = sca[i], scale.update=F, group=groups,verbose=T)
  
  prov.output<-matrix(, nrow = nrow(ldblock.final2), ncol = 0)
  prov.output<-cbind(coef(summary.bglm(mod1))[-1,4],as.numeric(ldblock.final2[,8]))
  #head(prov.output)
  
  prov.orderedOutput<-prov.output[order(prov.output[,1]),]
  
  prov.singleOutput<-unique(prov.orderedOutput[,2])
  
  rank<-mean(which(prov.singleOutput %in% c(17,19,24,25,28,29)))
  rank.all[i,2]<-rank
  
  #Output with effect size and block (top 6 blocks)
  new.output<-cbind(coef(summary.bglm(mod1))[-1,c(1,4)],as.numeric(ldblock.final2[,8]),qSum)
  new.orderedOutput<-new.output[order(new.output[,2]),]
  new.singleOutput<-new.orderedOutput[-which(duplicated(new.orderedOutput[,3])==T),]
  final.matrix<-cbind(new.singleOutput[1:6,1],ldblock.groupCorrespondence[new.singleOutput[1:6,3],1],new.singleOutput[1:6,4],new.singleOutput[1:6,3])
  assign((paste("final.matrix.noprior",i,sep='')),final.matrix)
  
  
  print(sca[i])
  
  ##Variance Ratio
  #Real blocks
  prov.output3<-coef(summary.bglm(mod1))[-1,1] #Estimated effects
  
  sd.true.fg<-sd(prov.output3[which(as.numeric(ldblock.final2[,8]) %in% c(17,19,24,25,28,29))])
  sd.false.fg<-sd(prov.output3[-which(as.numeric(ldblock.final2[,8]) %in% c(17,19,24,25,28,29))])
  
  
  #Top 6 blocks
  top.prov<-cbind(coef(summary.bglm(mod1))[-1,4],as.numeric(ldblock.final2[,8]))
  top.prov.orderedOutput<-top.prov[order(top.prov[,1]),]
  top<-unique(top.prov.orderedOutput[,2])
  
  sd.true.fg.top<-sd(prov.output3[which(as.numeric(ldblock.final2[,8]) %in% top[1:6])])
  sd.false.fg.top<-sd(prov.output3[-which(as.numeric(ldblock.final2[,8]) %in% top[1:6])])
  
  
  empirical.effects<-rbind(empirical.effects,c((sd.true.fg^2)/((sd.false.fg^2)),(sd.true.fg.top^2)/(sd.false.fg.top^2)))
  rownames(empirical.effects)[nrow(empirical.effects)]<-sca[i]
  
  #------------------------------------
  #Calculating the variance ratios centered at zero
  #Real blocks
  sd.true.fg0<-mean(prov.output3[which(as.numeric(ldblock.final2[,8]) %in% c(17,19,24,25,28,29))]^2)
  sd.false.fg0<-mean(prov.output3[-which(as.numeric(ldblock.final2[,8]) %in% c(17,19,24,25,28,29))]^2)

  
  #Top blocks
  sd.true.fg.top0<-mean(prov.output3[which(as.numeric(ldblock.final2[,8]) %in% top[1:6])]^2)
  sd.false.fg.top0<-mean(prov.output3[-which(as.numeric(ldblock.final2[,8]) %in% top[1:6])]^2)
  
  empirical.effects0<-rbind(empirical.effects0,c(sd.true.fg0/sd.false.fg0,sd.true.fg.top0/sd.false.fg.top0))
  rownames(empirical.effects0)[nrow(empirical.effects0)]<-sca[i]
}

proc.time() - ptm
empirical.effects
empirical.effects0


var.ratios<-cbind(empirical.effects[,2],empirical.effects[,1],empirical.effects0[,2],empirical.effects0[,1])[-1,]
colnames(var.ratios)<-c("top.blocks","true.blocks","top.blocks at 0","true.blocks at 0")

#Little table with the results
print.xtable(xtable(var.ratios,digits=3), type="html", file="var ratios.html")


#---------------------------------------------------------------------------------
#Adding prior knowledge as differential shrinkage

#shrinkage value to be tested
#S.min corresponds to the absence of prior knowledge
#s.max corresponds to maximum prior knowledge
#Other shrinkage values are interpolated linearly of exponentially in between
s.min<-c(0.01,0.001,0.0001,0.0001,0.0001,0.0001,0.00001,0.00001,0.00001,0.00001,0.00001,0.000001,0.000001,0.000001,0.000001)
s.max<-c(1,1,1,0.1,0.01,0.001,1,0.1,0.01,0.001,0.0001,0.1,0.01,0.001,0.0001)

rank.all.prior<-cbind(s.min,s.max,rep(99,times=length(s.min)))

empirical.effects.prior<-c(99,99)
empirical.effects.prior0<-c(99,99)

ptm <- proc.time()
for(i in 1:length(s.min)){
  #Linear interpolation
  shrinkValues<-seq(from=s.min[i],to=s.max[i],by=(s.max[i]-s.min[i])/(length(unique(qSum))-1))

  #Exponential interpolation
  #shrinkValues<-s.min[i]*(exp(log(s.max[i]/s.min[i])/max(qSum)))^(unique(qSum)[order(unique(qSum))])

  shrinkValues.index<-unique(qSum)[order(unique(qSum))]
  shrinkContinuous.Exp<-rep(99,times=length(qSum))
  for(j in 1:length(qSum)){
    shrinkContinuous.Exp[j]<-shrinkValues[which(shrinkValues.index==qSum[j])]
  }

  mod<-bglm(bmi ~ ., data=bmidata.bglm, family = gaussian, prior = "t", prior.mean=0, mean.update=F,prior.scale = shrinkContinuous.Exp, scale.update=F, group=groups,verbose=T)
  
  prov.output<-matrix(, nrow = nrow(ldblock.final2), ncol = 0)
  prov.output<-cbind(coef(summary.bglm(mod))[-1,4],as.numeric(ldblock.final2[,8]))

  prov.orderedOutput<-prov.output[order(prov.output[,1]),]

  prov.singleOutput<-unique(prov.orderedOutput[,2])
  
  rank<-mean(which(prov.singleOutput %in% c(17,19,24,25,28,29)))
  rank.all.prior[i,3]<-rank

  #Output with effect size and block (top 6 blocks)
  new.output<-cbind(coef(summary.bglm(mod))[-1,c(1,4)],as.numeric(ldblock.final2[,8]),qSum)
  new.orderedOutput<-new.output[order(new.output[,2]),]
  new.singleOutput<-new.orderedOutput[-which(duplicated(new.orderedOutput[,3])==T),]
  final.matrix<-cbind(new.singleOutput[1:6,1],ldblock.groupCorrespondence[new.singleOutput[1:6,3],1],new.singleOutput[1:6,4],new.singleOutput[1:6,3])
  assign((paste("final.matrix.noprior",i,sep='')),final.matrix)

  
  
  ##Variance Ratio
  #Real blocks
  prov.output3<-coef(summary.bglm(mod))[-1,1] #Estimated effects
  
  sd.true.fg<-sd(prov.output3[which(as.numeric(ldblock.final2[,8]) %in% c(17,19,24,25,28,29))])
  sd.false.fg<-sd(prov.output3[-which(as.numeric(ldblock.final2[,8]) %in% c(17,19,24,25,28,29))])
  
  
  #Top 6 blocks
  top.prov<-cbind(coef(summary.bglm(mod))[-1,4],as.numeric(ldblock.final2[,8]))
  top.prov.orderedOutput<-top.prov[order(top.prov[,1]),]
  top<-unique(top.prov.orderedOutput[,2])
  
  sd.true.fg.top<-sd(prov.output3[which(as.numeric(ldblock.final2[,8]) %in% top[1:6])])
  sd.false.fg.top<-sd(prov.output3[-which(as.numeric(ldblock.final2[,8]) %in% top[1:6])])
  
  empirical.effects.prior<-rbind(empirical.effects.prior,c((sd.true.fg^2)/((sd.false.fg^2)),(sd.true.fg.top^2)/(sd.false.fg.top^2)))
  rownames(empirical.effects.prior)[nrow(empirical.effects.prior)]<-paste(s.min[i],s.max[i],sep='-')
  
  #------------------------------------
  #Calculating the variance ratios centered at zero
  #Real blocks
  sd.true.fg0<-mean(prov.output3[which(as.numeric(ldblock.final2[,8]) %in% c(17,19,24,25,28,29))]^2)
  sd.false.fg0<-mean(prov.output3[-which(as.numeric(ldblock.final2[,8]) %in% c(17,19,24,25,28,29))]^2)
  
  
  #Top blocks
  sd.true.fg.top0<-mean(prov.output3[which(as.numeric(ldblock.final2[,8]) %in% top[1:6])]^2)
  sd.false.fg.top0<-mean(prov.output3[-which(as.numeric(ldblock.final2[,8]) %in% top[1:6])]^2)
  
  empirical.effects.prior0<-rbind(empirical.effects.prior0,c(sd.true.fg0/sd.false.fg0,sd.true.fg.top0/sd.false.fg.top0))
  rownames(empirical.effects.prior0)[nrow(empirical.effects.prior0)]<-paste(s.min[i],s.max[i],sep='-')
  
  
  print(c(i,s.min[i],s.max[i]))
}
proc.time() - ptm
rank.all.prior
empirical.effects.prior
empirical.effects.prior0

#Little table with the resutls
var.ratios.prior<-cbind(empirical.effects.prior[,2],empirical.effects.prior[,1],empirical.effects.prior0[,2],empirical.effects.prior0[,1])[-1,]
colnames(var.ratios.prior)<-c("top.blocks","true.blocks","top.blocks at 0","true.blocks at 0")
print.xtable(xtable(var.ratios.prior,digits=3), type="html", file="var ratios prior knowledge.html")


