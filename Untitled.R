# Larks Analysis 
library(ggplot2)
library(RColorBrewer)
library(data.table)

#Read Dengue data.
infile="~/GitHub/LARK_Analysis/Viral_LARK_Analysis/LabStrains/177214_DENV2.csv"

cutoff<-data.frame(STGGYG=4.0,GYNGFG=10.0,SYSGYS=-3.0)

LARKwindowFile<-fread(infile,sep = "\t",header = F,skip = 1)

colnames(LARKwindowFile)<-c("pos","seq","zip","STGGYG","GYNGFG","SYSGYS")
LARK=melt(LARKwindowFile,id.vars = c(1:2),measure.vars = c(4:6))
data.table(t(apply(LARKwindowFile[,4:6,by=V1],MARGIN = 1,function(X){X>cutoff})))

# Read in other data (fitness, structure, regions, proteins)
load("/Users/ptdolan/Google Drive/DenguePanelsAndOutput4-18/fitnesstable_All.Rdata")
cast=dcast(Wtable[!is.na(pos)],pos+wtRes+PONDR+TMstat+Class+reg~.,value.var = "wrel.ciUpper",fun.aggregate = max)
mergedLarks<-merge(LARK,cast,by="pos")

# Do some plotting
LARKS<-ggplot(mergedLarks)+
  geom_line(aes(pos,value,group=variable,col=ifelse(PONDR>0.6,variable,"Structured")))+
  geom_vline(xintercept = mergedLarks[,min(pos),by=reg]$V1)+
  geom_label(data=mergedLarks[,min(pos),by=reg],aes(V1+35,1200,label=reg))+
  scale_color_manual(values = c(brewer.pal(name = "Dark2",n=3),"grey"),"Backbones")

PONDR<-ggplot(mergedLarks)+
  geom_line(aes(pos,PONDR,group=variable,col=ifelse(PONDR>0.6,"Unstructured","Structured")))+
  geom_vline(xintercept = mergedLarks[,min(pos),by=reg]$V1)+
  geom_label(data=mergedLarks[,min(pos),by=reg],aes(V1+35,1,label=reg))+scale_color_brewer(palette = "Set1",direction = -1,"PONDR Pred.")
TMstat<-ggplot(mergedLarks)+geom_point(aes(pos,TMstat,col=TMstat))
cowplot::plot_grid(LARKS,PONDR,TMstat,nrow = 3)




