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
cast=dcast(Wtable[!is.na(pos)],pos+wtRes+TMstat+Class+reg~.,value.var = "wrel.ciUpper",fun.aggregate = max)
mergedLarks<-merge(LARK,cast,by="pos")

disoPred<-fread("~/GitHub/LARK_Analysis/Viral_LARK_Analysis/LabStrains/pondr.txt")
mergedLarks<-merge(mergedLarks,disoPred,by.x="pos",by.y="ResNum")

mergedLarks[,co:=as.numeric(cutoff[variable])]
mergedLarks[,LARK:=value<co]
mergedLarks[,struct:=as.factor(ifelse(((PONDRFIT-PONDERFITERROR)>0.5),"Unstructured","Structured"))]

# Do some plotting
LARKS<-ggplot(mergedLarks)+
  geom_line(aes(pos,value,group=variable,col=variable))+
  geom_point(cex=1,data=mergedLarks[LARK==T],aes(pos,value,group=variable,alpha=struct,col=variable))+
  geom_vline(xintercept = mergedLarks[,min(pos),by=reg]$V1)+
  scale_color_manual(values = c(brewer.pal(name = "Dark2",n=3),"black","grey"),"Backbones")

PONDR<-ggplot(mergedLarks)+
  geom_line(aes(pos,PONDRFIT,group=variable,col=ifelse(PONDRFIT-PONDERFITERROR>0.5,"Unstructured","Structured")))+
  geom_vline(xintercept = mergedLarks[,min(pos),by=reg]$V1)+
  geom_label(data=mergedLarks[,min(pos),by=reg],aes(V1+35,1,label=reg))+scale_color_brewer(palette = "Set1",direction = -1,"PONDR Pred.")

TMstat<-ggplot(mergedLarks)+geom_point(aes(pos,TMstat,col=TMstat))

StrClass<-ggplot(mergedLarks)+geom_point(aes(pos,Class,col=TMstat))

ggsave(filename = paste(sep = "","~/GitHub/LARK_Analysis/Viral_LARK_Analysis/LabStrains/polyprotein.pdf"),cowplot::plot_grid(LARKS,PONDR,TMstat,StrClass,nrow = 4))

for (L in levels(mergedLarks$reg)){
  LARKS<-ggplot(mergedLarks[reg==L])+
    geom_line(aes(pos,value,group=variable,col=variable,alpha=struct))+
    geom_point(cex=1,data=mergedLarks[reg==L&LARK==T],aes(pos,value,group=variable,alpha=struct,col=variable))+
    geom_vline(xintercept = mergedLarks[reg==L,min(pos),by=reg]$V1)+
    scale_color_manual(values = c(brewer.pal(name = "Dark2",n=3),"black","grey"),"Backbones")
    
  PONDR<-ggplot(mergedLarks[reg==L])+
    geom_line(aes(pos,PONDRFIT,group=variable,col=struct))+
    geom_vline(xintercept = mergedLarks[reg==L,min(pos),by=reg]$V1)+
    scale_color_brewer(palette = "Set1",direction = -1,"PONDR Pred.")
  
  TMstat<-ggplot(mergedLarks[reg==L])+geom_point(cex=0.2,aes(pos,TMstat,col=TMstat))
  StrClass<-ggplot(mergedLarks[reg==L])+geom_point(cex=0.2,aes(group=Class,pos,Class,col=TMstat))  
  ggsave(width=5,height=7,filename = paste(sep = "","~/GitHub/LARK_Analysis/Viral_LARK_Analysis/LabStrains/",L,".pdf"),cowplot::plot_grid(LARKS,PONDR,TMstat,StrClass,nrow = 4))
  
}





