#Assessing prior knowledge for the SNPs that were replicated
#---The idea is to make a comparison between prior knowledge of the SNPs detected
#---The standard analysis and the baysian method
#---02Nov2017

library(BhGLM)
library(data.table)
library(ggplot2)

#Working directory
setwd('/Users/miguelmspereira/Box Sync/Congenica/scripts/')

#Top 50,000 SNPs from the standard SNP analysis
#FVC
fvc.threshold<-read.table('fvctop500000withLDblocks.tsv',header=T,sep='\t')
head(fvc.threshold)
dim(fvc.threshold)
length(unique(fvc.threshold$snp))
length(unique(fvc.threshold$chr))
length(unique(fvc.threshold$LD.block))


#Top 50,000 SNPs from the standard SNP analysis
#Ratio
ratio.threshold<-read.table('ratiotop500000withLDblocks.tsv',header=T,sep='\t')
head(ratio.threshold)
dim(ratio.threshold)
length(unique(ratio.threshold$snp))
length(unique(ratio.threshold$chr))
length(unique(ratio.threshold$LD.block))



###################################################################################
#Subsetting to the top 20.000 SNPs
#FVC
fvc.threshold25<-fvc.threshold[which(fvc.threshold$snp %in% unique(fvc.threshold$snp)[1:20000]),]
length(unique(fvc.threshold25$snp))
dim(fvc.threshold25)

#Ratio
ratio.threshold25<-ratio.threshold[which(ratio.threshold$snp %in% unique(ratio.threshold$snp)[1:20000]),]
length(unique(ratio.threshold25$snp))
dim(ratio.threshold25)



###################################################################################
#Prior knowledge matrix
q.fvc<-read.csv('priorknowledge_fvc.csv')
q.final.fvc<-q.fvc[which(q.fvc$snp %in% fvc.threshold25$snp),]
head(q.final.fvc)
table(q.final.fvc$qSum) #Prior knowledge distribution

q.ratio<-read.csv('/Users/miguelmspereira/Desktop/lungdevprior/ratio_dosages/priorknowledge/priorknowledge_ratio.csv')
q.final.ratio<-q.ratio[which(q.ratio$snp %in% ratio.threshold25$snp),]
head(q.final.ratio)
table(q.final.ratio$qSum) #Prior knowledge distribution

length(intersect(q.final.fvc$snp,q.final.ratio$snp)) #3700 SNPs intersection (18.5% intersection)


###################################################################################
#Results - to get the replicated SNPs
fvc.standard<-read.csv('fvc standard results with eaf.csv',header=T, stringsAsFactors = F)
fvc.bayes<-read.csv('fvc bayes results with eaf.csv',header=T, stringsAsFactors = F)
ratio.standard<-read.csv('ratio standard results with eaf.csv',header=T, stringsAsFactors = F)
ratio.bayes<-read.csv('ratio bayes results with eaf.csv',header=T, stringsAsFactors = F)


#List of known signals - to remove the LD blocks in common
fvc.known<-read.csv('fvc known signals.txt',header=T,stringsAsFactors = F,sep='\t')
head(fvc.known)
ratio.known<-read.table('ratio known signals.txt',header=T,stringsAsFactors = F,sep='\t')
head(ratio.known)

#Remove known LD blocks - Ratio
intersect(fvc.standard$LD.block,fvc.known$LD.Block) #just a little check
intersect(fvc.bayes$LD.block,fvc.known$LD.Block) #just a little check
intersect(ratio.standard$LD.block,ratio.known$LD.Block)
intersect(ratio.bayes$LD.block,ratio.known$LD.Block)

ratio.standard2<-ratio.standard[-which(ratio.standard$LD.block %in% intersect(ratio.standard$LD.block,ratio.known$LD.Block)),]
ratio.bayes2<-ratio.bayes[-which(ratio.bayes$LD.block %in% intersect(ratio.bayes$LD.block,ratio.known$LD.Block)),]

#Top 100 signlas for FVC and top 400 signals for FEV1/FVC - number of signals taken to replication according to the power calculations
fvc.standard100<-fvc.standard[1:100,]
fvc.bayes100<-fvc.bayes[1:100,]
ratio.standard400<-ratio.standard2[1:393,]
ratio.bayes400<-ratio.bayes2[1:400,]



#Replicated signals matrices
fvc.standard.replicated<-fvc.standard100[which(fvc.standard100$meta.pval<(0.05/100)),]
fvc.bayes.replicated<-fvc.bayes100[which(fvc.bayes100$meta.pval<(0.05/100)),]
ratio.standard.replicated<-ratio.standard400[which(ratio.standard400$meta.pval<(0.05/393)),]
ratio.bayes.replicated<-ratio.bayes400[which(ratio.bayes400$meta.pval<(0.05/400)),]



#Combined matrices
fvc.comb.gene.0<-rbind(fvc.standard[1:100,],fvc.bayes[1:100,-c(10,11,32,33)])
fvc.comb.gene.1<-fvc.comb.gene.0[-which(duplicated(fvc.comb.gene.0$snp)),]
fvc.comb.gene.2<-fvc.comb.gene.1[order(fvc.comb.gene.1$p.val),][1:100,]
fvc.comb.gene.3<-fvc.comb.gene.2[which(fvc.comb.gene.2$meta.pval[1:100]<(0.05/(100))),]
dim(fvc.comb.gene.3) #32 SNPs
length(unique(fvc.comb.gene.3$gene1)) #19 genes

ratio.comb.gene.0<-rbind(ratio.standard2[1:393,],ratio.bayes2[1:400,-c(10,11,32,33)])
ratio.comb.gene.1<-ratio.comb.gene.0[-which(duplicated(ratio.comb.gene.0$snp)),]
ratio.comb.gene.2<-ratio.comb.gene.1[order(ratio.comb.gene.1$p.val),][1:400,]
ratio.comb.gene.3<-ratio.comb.gene.2[which(ratio.comb.gene.2$meta.pval[1:400]<(0.05/(400))),]
dim(ratio.comb.gene.3) #36 SNPs
length(unique(ratio.comb.gene.3$gene1)) #20 genes



#All genes
unique(fvc.standard$gene1[which(fvc.standard$meta.pval[1:100]<(0.05/(100)))])
gene.list<-list()
gene.list$fvc.standard<-unique(fvc.standard.replicated$gene1)
gene.list$fvc.bayes<-unique(fvc.bayes.replicated$gene1)
gene.list$fvc.combined<-unique(fvc.comb.gene.3$gene1)
gene.list$ratio.standard<-unique(ratio.standard.replicated$gene1)
gene.list$ratio.bayes<-unique(ratio.bayes.replicated$gene1)
gene.list$ratio.combined<-unique(ratio.comb.gene.3$gene1)
str(gene.list)



#Prior knowledge of the replicated SNPs
q.fvc.standard<-q.final.fvc[which(q.final.fvc$snp %in% fvc.standard.replicated$snp),]
q.fvc.bayes<-q.final.fvc[which(q.final.fvc$snp %in% fvc.bayes.replicated$snp),]
q.fvc.comb<-q.final.fvc[which(q.final.fvc$snp %in% fvc.comb.gene.3$snp),]
q.ratio.standard<-q.final.ratio[which(q.final.ratio$snp %in% ratio.standard.replicated$snp),]
q.ratio.bayes<-q.final.ratio[which(q.final.ratio$snp %in% ratio.bayes.replicated$snp),]
q.ratio.comb<-q.final.ratio[which(q.final.ratio$snp %in% ratio.comb.gene.3$snp),]



#p-values of replicated SNPs and scores
q.pval.fvc.standard<-merge(fvc.standard.replicated[,c(1,4,12)],q.final.fvc,by='snp')
q.pval.fvc.bayes<-merge(fvc.bayes.replicated[,c(1,4,14)],q.final.fvc,by='snp')
q.pval.fvc.comb<-merge(fvc.comb.gene.3[,c(1,4,12)],q.final.fvc,by='snp')
q.pval.ratio.standard<-merge(ratio.standard.replicated[,c(1,4,12)],q.final.ratio,by='snp')
q.pval.ratio.bayes<-merge(ratio.bayes.replicated[,c(1,4,14)],q.final.ratio,by='snp')
q.pval.ratio.comb<-merge(ratio.comb.gene.3[,c(1,4,12)],q.final.ratio,by='snp')

q.pval.fvc.all<-merge(fvc.threshold25[-which(duplicated(fvc.threshold25$snp)),c(1,16,11)],q.final.fvc,by='snp')
q.pval.ratio.all<-merge(ratio.threshold25[-which(duplicated(ratio.threshold25$snp)),c(1,16,11)],q.final.ratio,by='snp',all.x = F)



#Box plot p-values by score
#FVC
q.pval.data.fvc<-as.data.frame(cbind(
  c(rep('Standard Analysis',times=nrow(q.pval.fvc.standard)),rep('Bayesian Analysis',times=nrow(q.pval.fvc.bayes)),rep('Combined Analysis',times=nrow(q.pval.fvc.comb)),rep('All SNPs',times=20000)),
  c(q.pval.fvc.standard$p.val,q.pval.fvc.bayes$p.val,q.pval.fvc.comb$p.val,q.pval.fvc.all$p.val),
  c(q.pval.fvc.standard$qSum,q.pval.fvc.bayes$qSum,q.pval.fvc.comb$qSum,q.pval.fvc.all$qSum)
))
colnames(q.pval.data.fvc)<-c('analysis','p.value','score')
q.pval.data.fvc$analysis<-factor(q.pval.data.fvc$analysis,levels=c('Standard Analysis','Bayesian Analysis','Combined Analysis','All SNPs'))
q.pval.data.fvc$p.value<--log10(as.numeric(as.character(q.pval.data.fvc$p.value)))
q.pval.data.fvc$score<-as.factor(q.pval.data.fvc$score)



ggplot(q.pval.data.fvc,aes(x=score,y=p.value))+
  geom_boxplot()+
  labs(x='Score', y='-log P',title='FVC - Only replicated SNPs')+
  facet_grid(. ~ analysis)+
  theme_bw()+
  theme(legend.position='none',panel.grid.minor = element_blank())+
  scale_y_continuous(limits=c(0,12),breaks=seq(2,12,2))



#FEV1/FVC
q.pval.data.ratio<-as.data.frame(cbind(
  c(rep('Standard Analysis',times=nrow(q.pval.ratio.standard)),rep('Bayesian Analysis',times=nrow(q.pval.ratio.bayes)),rep('Combined Analysis',times=nrow(q.pval.ratio.comb)),rep('All SNPs',times=20000)),
  c(q.pval.ratio.standard$p.val,q.pval.ratio.bayes$p.val,q.pval.ratio.comb$p.val,q.pval.ratio.all$p.val),
  c(q.pval.ratio.standard$qSum,q.pval.ratio.bayes$qSum,q.pval.ratio.comb$qSum,q.pval.ratio.all$qSum)
))
colnames(q.pval.data.ratio)<-c('analysis','p.value','score')
q.pval.data.ratio$analysis<-factor(q.pval.data.ratio$analysis,levels=c('Standard Analysis','Bayesian Analysis','Combined Analysis','All SNPs'))
q.pval.data.ratio$p.value<--log10(as.numeric(as.character(q.pval.data.ratio$p.value)))
q.pval.data.ratio$score<-as.factor(q.pval.data.ratio$score)


ggplot(q.pval.data.ratio,aes(x=score,y=p.value))+
  geom_boxplot()+
  labs(x='Score', y='-log P',title='FEV1/FVC - Only replicated SNPs')+
  facet_grid(. ~ analysis)+
  theme_bw()+
  theme(legend.position='none',panel.grid.minor = element_blank())+
  scale_y_continuous(limits=c(0,12),breaks=seq(2,12,2))



#Box plot p-values by question
#FVC
colnames(q.pval.fvc.standard)[4:13]<-c('Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8','Q9','Q10')
colnames(q.pval.fvc.bayes)[4:13]<-c('Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8','Q9','Q10')
colnames(q.pval.fvc.comb)[4:13]<-c('Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8','Q9','Q10')
colnames(q.pval.fvc.all)[4:13]<-c('Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8','Q9','Q10')


boxplot.fvc.standard<-melt(q.pval.fvc.standard[,c(3:13)],id.vars='p.val')
boxplot.fvc.bayes<-melt(q.pval.fvc.bayes[,c(3:13)],id.vars='p.val')
boxplot.fvc.comb<-melt(q.pval.fvc.comb[,c(3:13)],id.vars='p.val')
boxplot.fvc.all<-melt(q.pval.fvc.all[,c(3:13)],id.vars='p.val')
######################################
qSum.pval.data.fvc<-as.data.frame(
  cbind(c(rep('Standard Analysis',times=nrow(boxplot.fvc.standard)),rep('Bayesian Analysis',times=nrow(boxplot.fvc.bayes)),rep('Combined Analysis',times=nrow(boxplot.fvc.comb)),rep('All SNPs',times=nrow(boxplot.fvc.all))),
  rbind(
  boxplot.fvc.standard,
  boxplot.fvc.bayes,
  boxplot.fvc.comb,
  boxplot.fvc.all
)))

colnames(qSum.pval.data.fvc)[1]<-'analysis'
qSum.pval.data.fvc$analysis<-factor(qSum.pval.data.fvc$analysis,levels=c('Standard Analysis','Bayesian Analysis','Combined Analysis','All SNPs'))
qSum.pval.data.fvc$p.val.log<--log10(qSum.pval.data.fvc$p.val)
qSum.pval.data.fvc$value2<-factor(qSum.pval.data.fvc$value,levels=c('0','1'))
######################################
#Without the 20.000 SNPs
qSum.pval.data.fvc<-as.data.frame(
  cbind(c(rep('Standard Analysis',times=nrow(boxplot.fvc.standard)),rep('Bayesian Analysis',times=nrow(boxplot.fvc.bayes)),rep('Combined Analysis',times=nrow(boxplot.fvc.comb))),
        rbind(
          boxplot.fvc.standard,
          boxplot.fvc.bayes,
          boxplot.fvc.comb
        )))

colnames(qSum.pval.data.fvc)[1]<-'analysis'
qSum.pval.data.fvc$analysis<-factor(qSum.pval.data.fvc$analysis,levels=c('Standard Analysis','Bayesian Analysis','Combined Analysis'))
qSum.pval.data.fvc$p.val.log<--log10(qSum.pval.data.fvc$p.val)
qSum.pval.data.fvc$value2<-factor(qSum.pval.data.fvc$value,levels=c('0','1'))
class(qSum.pval.data.fvc$variable)
######################################
ggplot(qSum.pval.data.fvc,aes(x=value2,y=p.val.log,fill=analysis))+
  geom_boxplot()+
  facet_wrap(~variable,ncol=5,nrow=2)+
  labs(x='Score', y='-log P',title='FVC - Only replicated SNPs')+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),legend.title = element_blank())+
  scale_y_continuous(limits=c(0,12),breaks=seq(2,12,2))

#ggplot(qSum.pval.data.fvc,aes(x=0,y=p.val.log,fill=analysis))+
#  geom_boxplot()+
#  facet_grid(variable~value2)

#ggplot(qSum.pval.data.fvc,aes(x=0,y=p.val.log,fill=analysis))+
#  geom_boxplot()+
#  facet_wrap(value2~variable)



#FEV1/FVC
colnames(q.pval.ratio.standard)[4:13]<-c('Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8','Q9','Q10')
colnames(q.pval.ratio.bayes)[4:13]<-c('Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8','Q9','Q10')
colnames(q.pval.ratio.comb)[4:13]<-c('Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8','Q9','Q10')
colnames(q.pval.ratio.all)[4:13]<-c('Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8','Q9','Q10')


boxplot.ratio.standard<-melt(q.pval.ratio.standard[,c(3:13)],id.vars='p.val')
boxplot.ratio.bayes<-melt(q.pval.ratio.bayes[,c(3:13)],id.vars='p.val')
boxplot.ratio.comb<-melt(q.pval.ratio.comb[,c(3:13)],id.vars='p.val')
boxplot.ratio.all<-melt(q.pval.ratio.all[,c(3:13)],id.vars='p.val')
######################################
qSum.pval.data.ratio<-as.data.frame(
  cbind(c(rep('Standard Analysis',times=nrow(boxplot.ratio.standard)),rep('Bayesian Analysis',times=nrow(boxplot.ratio.bayes)),rep('Combined Analysis',times=nrow(boxplot.ratio.comb)),rep('All SNPs',times=nrow(boxplot.ratio.all))),
        rbind(
          boxplot.ratio.standard,
          boxplot.ratio.bayes,
          boxplot.ratio.comb,
          boxplot.ratio.all
        )))

colnames(qSum.pval.data.ratio)[1]<-'analysis'
qSum.pval.data.ratio$analysis<-factor(qSum.pval.data.ratio$analysis,levels=c('Standard Analysis','Bayesian Analysis','Combined Analysis','All SNPs'))
qSum.pval.data.ratio$p.val.log<--log10(qSum.pval.data.ratio$p.val)
qSum.pval.data.ratio$value2<-factor(qSum.pval.data.ratio$value,levels=c('0','1'))
######################################
#Without the 20.000 SNPs
qSum.pval.data.ratio<-as.data.frame(
  cbind(c(rep('Standard Analysis',times=nrow(boxplot.ratio.standard)),rep('Bayesian Analysis',times=nrow(boxplot.ratio.bayes)),rep('Combined Analysis',times=nrow(boxplot.ratio.comb))),
        rbind(
          boxplot.ratio.standard,
          boxplot.ratio.bayes,
          boxplot.ratio.comb
        )))

colnames(qSum.pval.data.ratio)[1]<-'analysis'
qSum.pval.data.ratio$analysis<-factor(qSum.pval.data.ratio$analysis,levels=c('Standard Analysis','Bayesian Analysis','Combined Analysis'))
qSum.pval.data.ratio$p.val.log<--log10(qSum.pval.data.ratio$p.val)
qSum.pval.data.ratio$value2<-factor(qSum.pval.data.ratio$value,levels=c('0','1'))
######################################
ggplot(qSum.pval.data.ratio,aes(x=value2,y=p.val.log,fill=analysis))+
  geom_boxplot()+
  facet_wrap(~variable,ncol=5,nrow=2)+
  labs(x='Score', y='-log P',title='FEV1/FVC - Only replicated SNPs')+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),legend.title = element_blank())+
  scale_y_continuous(limits=c(0,12),breaks=seq(2,12,2))

#ggplot(qSum.pval.data.ratio,aes(x=0,y=p.val.log,fill=analysis))+
#  geom_boxplot()+
#  facet_grid(variable~value2)

#ggplot(qSum.pval.data.ratio,aes(x=0,y=p.val.log,fill=analysis))+
#  geom_boxplot()+
#  facet_wrap(value2~variable)


