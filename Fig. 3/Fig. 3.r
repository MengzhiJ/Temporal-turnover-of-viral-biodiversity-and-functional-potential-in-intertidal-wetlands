library(ggplot2)
library(dplyr)
library(readxl)
library(vegan)
library(ggpmisc)
library(betapart)
library(ggpubr)
library(reshape2)
library(tidyr)
library(nlme)
setwd('D:/Tidal_time/Rdata/Fig. 3')
#alpha diversity (the methods used for mOTUs and vPCs were similar to vOTUs)
otu <- read.delim('vOTU_bray.txt',row.names = 1)
otu1 <- data.frame(t(otu))
richness<- rowSums(otu1 > 0)
shannon_index <- diversity(otu1, index = 'shannon')
richness <- data.frame(richness)
shannon <- data.frame(shannon_index)
write.table (richness,file ="richness_vOTU.csv") 
write.table (shannon,file ="shannon_vOTU.csv")

#Fig.3a
data<-read_excel('Fig. 3.xlsx')
# calculate AIC to determine the best fit model of each variable
vOTU<-data[1:36,]
mOTU<-data[37:72,]
vPC<-data[73:108,]
object2<-gls(Richness~poly(Time,2), data=mOTU)
object1<-gls(Richness~poly(Time), data=mOTU)
summary(object1)
summary(object2)
#richness
ggplot(data,aes(x=Time,y=Richness))+
  geom_point(aes(fill=Type),size = 3.5, alpha = 0.7,shape=21)+
  labs(x="Time (month)",y="Alpha diversity (richness)")+
  geom_smooth(data=data,aes(color=Type),method = 'lm',se=F,size=1.3,fullrange=T,formula=y~poly(x,2))+
  theme_test()+ theme(legend.position ="none")+scale_color_manual(breaks=c("vPCs","vOTUs","Host microbes"),
                                                                  values=c("#2570AE","#9D5D2E","gray"))+
  scale_fill_manual(breaks=c("vPCs","vOTUs","Host microbes"),
                    values=c("#2570AE","#9D5D2E","gray"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+scale_x_continuous(breaks = seq(0,10,2))+ylim(100,600)+
  stat_poly_eq(aes(color=Type,label = paste(..rr.label..,p.value.label,sep = "~`,`~"),size=2),
             formula = y ~ poly(x,2),parse = TRUE)

#shannon
ggplot(data, mapping = aes(x=Time,y=Shannon))+
  geom_point(aes(fill=Type),size = 3.5, alpha = 0.7,shape=21)+
  labs(x="Time (month)",y="Alpha diversity (shannon index)")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(aes(color=Type),method = 'lm',se=F,size=1.3,fullrange=T,formula=y~poly(x,2))+
  scale_color_manual(breaks=c("vPCs","vOTUs","Host microbes"),values=c("#2570AE","#9D5D2E","gray"))+
  scale_fill_manual(breaks=c("vPCs","vOTUs","Host microbes"),
                    values=c("#2570AE","#9D5D2E","gray"))+
  theme_test()+ theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),strip.text = element_text(size = 14),
        legend.title = element_text(size=16))+scale_x_continuous(breaks = seq(0,10,2))+facet_wrap(~Type,scales="free_y")+
stat_poly_eq(aes(color=Type,label = paste(..rr.label..,p.value.label,sep = "~`,`~"),size=2),
             formula = y ~ poly(x,2),parse = TRUE)

#Fig. 3b beta-diversity
votu<-read.delim('votu_sor.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#sorensen
distv<-beta.pair(t(votu),index.family = "sorensen")
sorensen<-as.matrix(distv$beta.sor)
write.csv(sorensen,"sorensen_votu.csv")
#bray
distv<-vegdist(t(votu), method="bray")
bray<-as.matrix(distv)
write.csv(bray,"bray_votu.csv")
#sorensen
data<-read_excel('Fig. 3.xlsx')
data$Time <- factor(data$Time,levels=c('0','2','4','6','8','10'))
ggplot(data, mapping = aes(x=Time,y=Beta))+
  geom_point(aes(fill=Type),size = 3.5, alpha = 0.7,shape=21)+
  labs(x="Time (month)",y="Beta diversity (community distance)")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(aes(color=Type,group=Type),method = 'lm',se=F,size=1.3,formula = y ~ poly(x,2))+
  theme_test()+ theme(legend.position ="none")+scale_color_manual(breaks=c("vPCs","vOTUs","Host microbes"),
                                                                  values=c("#2570AE","#9D5D2E","gray"))+
  scale_fill_manual(breaks=c("vPCs","vOTUs","Host microbes"),
                    values=c("#2570AE","#9D5D2E","gray"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+
stat_poly_eq(aes(color=Type,group=Type,label = paste(..rr.label..,p.value.label,sep = "~`,`~"),size=2),
             formula = y ~ poly(x,2),parse = TRUE)  

#bray
data<-read_excel("Fig. 3.xlsx")
data$Time <- factor(data$Time,levels=c('0','2','4','6','8','10'))
ggplot(data, mapping = aes(x=Time,y=Beta))+
  geom_point(aes(fill=Type),size = 3.5, alpha = 0.7,shape=21)+
  labs(x="Time (month)",y="Beta diversity (Bray-Curtis distance)")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(aes(color=Type),method = 'lm',se=F,size=1.3,formula = y ~ poly(x,2))+
  theme_test()+ theme(legend.position ="none")+scale_color_manual(breaks=c("vPCs","vOTUs","Host microbes"),
                                                                  values=c("#2570AE","#9D5D2E","gray"))+
  scale_fill_manual(breaks=c("vPCs","vOTUs","Host microbes"),
                    values=c("#2570AE","#9D5D2E","gray"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+
  stat_poly_eq(aes(color=Type,label = paste(..rr.label..,p.value.label,sep = "~`,`~"),size=2),
               formula = y ~ poly(x,2),parse = TRUE)  

#Fig. 3c
data<-read_excel("Fig. 3.xlsx")
ggplot(data = data,aes(x=dist_loca_num,y=distv_num))+
  geom_point(aes(fill=Type),size = 3.5, alpha = 0.7,shape=21)+
  geom_smooth(aes(color=Type),method = 'lm',se=F,size=1.3)+
  labs(x="Time (month)",y="Cumulative richness")+
  theme_test()+ theme(legend.position ="none")+scale_color_manual(breaks=c("vPCs","vOTUs","Host microbes"),
                                                                  values=c("#2570AE","#9D5D2E","gray"))+
  scale_fill_manual(breaks=c("vPCs","vOTUs","Host microbes"),
                    values=c("#2570AE","#9D5D2E","gray"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+scale_x_continuous(breaks = seq(0,10,2))+ylim(230,300)+
  stat_poly_eq(aes(color=Type,label = paste(..rr.label..,p.value.label,sep = "~`,`~"),size=2),
               formula = y ~ x,parse = TRUE)  

#Fig. 3d
#sorensen
votu<-read.delim('votu_sor.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
distv<-beta.pair(t(votu),index.family = "sorensen")
distv_num<-1-as.numeric(distv$beta.sor)
distv_num<-log10(distv_num)
distv_num<-as.numeric(distv_num)

loca<-read_excel('time.xlsx')
dist_loca<-vegdist(as.numeric(loca$Time),method = 'euclidean')
dist_loca_log<-log10(dist_loca)
dist_loca_num<-as.numeric(dist_loca_log)

data_new<-data.frame(distv_num,dist_loca_num)
write.csv(data_new,'TDR_vPCs_bray.csv')

data<-read_excel('Fig. 3.xlsx')

ggplot(data = data,aes(x=dist_loca_num,y=distv_num))+
  geom_point(aes(fill=Type),size = 3.5, alpha = 0.5,shape=21)+
  geom_smooth(aes(color=Type),method = 'lm',se=F,size=1.3)+
  labs(x="ln[Time (month)]",y="ln[Community similarity]")+
  theme_test()+ theme(legend.position ="none")+scale_color_manual(breaks=c("vPCs","vOTUs","Host microbes"),
                                                                  values=c("#2570AE","#9D5D2E","gray"))+
  scale_fill_manual(breaks=c("vPCs","vOTUs","Host microbes"),
                    values=c("#2570AE","#9D5D2E","gray"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+scale_x_continuous(breaks = seq(0,1,0.2))+
  stat_poly_eq(aes(color=Type,label = paste(..eq.label..,p.value.label,sep = "~`,`~"),size=2),
               formula = y ~ x,parse = TRUE)  

#Fig. 3e
data<-t(read.delim('votu_sor.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE))
data.core.s<-betapart.core(data) 
data.dist.sor<-beta.pair(data.core.s, index.family="sor")
beta <- as.matrix(data.dist.sor$beta.sor)
diag(beta) <- 0
beta[upper.tri(beta)] <- 0
beta <- reshape2::melt(beta)
beta <- subset(beta, value != 0)

#turnover
turnover <- as.matrix(data.dist.sor$beta.sim)
diag(turnover) <- 0
turnover[upper.tri(turnover)] <- 0  
turnover <- reshape2::melt(turnover)
turnover <- subset(turnover, value != 0)

#nestedness
nestedness <- as.matrix(data.dist.sor$beta.sne)
diag(nestedness) <- 0  
nestedness[upper.tri(nestedness)] <- 0   
nestedness <- reshape2::melt(nestedness)
nestedness <- subset(nestedness, value != 0)

#merge
dff=merge(beta,turnover,by=c("Var1","Var2"))
dff=merge(dff,nestedness,by=c("Var1","Var2"))
colnames(dff)=c("site1","site2","beta","turnover","nestedness")
dff$site=paste(dff$site1,dff$site2,sep="_")

write.csv(dff,"votu_betapart_sor.csv")


#bray
data<-t(read.delim('vOTU_bray.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE))
data.core.s<-betapart.core.abund(data) 
data.dist.bray<-beta.pair.abund(data.core.s, index.family="bray")
beta <- as.matrix(data.dist.bray$beta.bray)
diag(beta) <- 0
beta[upper.tri(beta)] <- 0
beta <- reshape2::melt(beta)
beta <- subset(beta, value != 0)

#turnover
turnover <- as.matrix(data.dist.bray$beta.bray.bal)
diag(turnover) <- 0
turnover[upper.tri(turnover)] <- 0  
turnover <- reshape2::melt(turnover)
turnover <- subset(turnover, value != 0)
#nestedness
nestedness <- as.matrix(data.dist.bray$beta.bray.gra)
diag(nestedness) <- 0  
nestedness[upper.tri(nestedness)] <- 0   
nestedness <- reshape2::melt(nestedness)
nestedness <- subset(nestedness, value != 0)
#merge
dff=merge(beta,turnover,by=c("Var1","Var2"))
dff=merge(dff,nestedness,by=c("Var1","Var2"))
colnames(dff)=c("site1","site2","beta","turnover","nestedness")
dff$site=paste(dff$site1,dff$site2,sep="_")
write.csv(dff,"votu_betapart_bray.csv")

#plot
data<-read.csv("betapart.csv")
ggtern(data,aes(1-beta,turnover,nestedness,fill=Type))+
  geom_point(alpha=0.5,size=2,shape=21)+ scale_fill_manual(breaks=c("vPCs","Host microbes","vOTUs"),
                                                  values=c("#2570AE","gray","#9D5D2E"))+
  theme_rgbw(base_size = 14, base_family = "")+
  theme(panel.border = element_rect(size=0.1))+
  labs(x ="",y = "",z = "",
       title = "")+ theme(legend.position ="none")

#Fig. 3f
data<-read.csv("betapart_dis.csv")
data$Type <- factor(data$Type,levels=c('vPCs','mOTUs','vOTUs'))
ggboxplot(data,notch = F,bxp.errorbar = T, x='component', y='beta', fill = 'Type', palette = c("#2570AE","gray","#9D5D2E"),width = 0.5,size=0.3,outlier.size=0.2)+
  labs(x="",y="Beta diversity (bray-curtis distance)",fill="Type")+
  theme_test()+  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
                       axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
                       axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
                       axis.text.y=element_text(size=16,colour = 'black'),
                       panel.border = element_rect(size=1.2),
                       legend.text = element_text(size=15),
                       legend.title = element_text(size=16))+
  theme(legend.position = "none")
