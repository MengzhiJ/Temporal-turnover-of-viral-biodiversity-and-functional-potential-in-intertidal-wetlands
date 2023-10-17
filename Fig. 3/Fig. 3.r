library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(vegan)
library(ggpmisc)
library(betapart)
library(ggpubr)
library(reshape2)
library(ggtern)
#Fig.3a
data<-read_excel("E:/Tidal_time/R_for_github/Fig. 3/Fig. 3.xlsx")
ggplot(data, mapping = aes(x=Time,y=Richness))+
  geom_point(aes(fill=Type),size = 3.5, alpha = 0.7,shape=21)+
  labs(x="Time (month)",y="Alpha diversity (richness)")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(aes(color=Type),method = 'gam',se=F,size=1.3,fullrange=T,formula=y~s(x,bs = "ts",k = 3))+
  theme_test()+ theme(legend.position ="none")+scale_color_manual(breaks=c("vPCs","vOTUs","mOTUs"),
                                                                  values=c("#2570AE","#9D5D2E","gray"))+
  scale_fill_manual(breaks=c("vPCs","vOTUs","mOTUs"),
                    values=c("#2570AE","#9D5D2E","gray"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+scale_x_continuous(breaks = seq(0,10,2))+
stat_poly_eq(aes(color=Type,label = paste(..rr.label..,p.value.label,sep = "~`,`~"),size=2),
             formula = y ~ poly(x,2),parse = TRUE)
#Fig.3b
data<-read_excel("E:/Tidal_time/R_for_github/Fig. 3/Fig. 3.xlsx")
data$Time <- factor(data$Time,levels=c('0','2','4','6','8','10'))
ggplot(data, mapping = aes(x=Time,y=Beta))+
  geom_point(aes(fill=Type),size = 3.5, alpha = 0.7,shape=21)+
  labs(x="Time (month)",y="Beta diversity (community distance)")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(aes(color=Type,group=Type),method = 'lm',se=F,size=1.3,formula = y ~ poly(x,2))+
  theme_test()+ theme(legend.position ="none")+scale_color_manual(breaks=c("vPCs","vOTUs","mOTUs"),
                                                                  values=c("#2570AE","#9D5D2E","gray"))+
  scale_fill_manual(breaks=c("vPCs","vOTUs","mOTUs"),
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
#Fig.3c
data_new<-read_excel("E:/Tidal_time/R_for_github/Fig. 3/Fig. 3.xlsx")
ggplot(data = data_new,aes(x=dist_loca_num,y=distv_num))+
  geom_point(aes(fill=Type),size = 3.5, alpha = 0.7,shape=21)+
  geom_smooth(aes(color=Type,group=Type),method = 'lm',se=F,size=1.3)+
  labs(x="Time (month)",y="Cumulative OTU richness")+
  theme_test()+ theme(legend.position ="none")+scale_color_manual(breaks=c("vPCs","vOTUs","mOTUs"),
                                                                  values=c("#2570AE","#9D5D2E","gray"))+
  scale_fill_manual(breaks=c("vPCs","vOTUs","mOTUs"),
                    values=c("#2570AE","#9D5D2E","gray"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+ylim(200,300)+scale_x_continuous(breaks = seq(0,10,2))+
  stat_poly_eq(aes(color=Type,group=Type,label = paste(..rr.label..,p.value.label,sep = "~`,`~"),size=2),
               formula = y ~ x,parse = TRUE)  
#Fig. 3d
#calculate Sorensen index and time distance
votu<-read.delim('host_beta.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
distv<-beta.pair(t(votu),index.family = "sorensen")
distv_num<-1-as.numeric(distv$beta.sor)
distv_num<-log10(distv_num)
distv_num<-as.numeric(distv_num)
loca<-read_excel('G:/Tidal_time/time.xlsx')
dist_loca<-vegdist(as.numeric(loca$Time),method = 'euclidean')
dist_loca_log<-log10(dist_loca)
dist_loca_num<-as.numeric(dist_loca_log)
data_new<-data.frame(distv_num,dist_loca_num)
write.csv(data_new,'G:/Tidal_time/TDR_host.csv')

data_new<-read_excel('E:/Tidal_time/R_for_github/Fig. 3/Fig. 3.xlsx')
ggplot(data = data_new,aes(x=dist_loca_num,y=distv_num))+
  geom_point(aes(fill=Type),size = 3.5, alpha = 0.5,shape=21)+
  geom_smooth(aes(color=Type,group=Type),method = 'lm',se=F,size=1.3)+
  labs(x="ln[Time (month)]",y="ln[Community similarity]")+
  theme_test()+ theme(legend.position ="none")+scale_color_manual(breaks=c("vPCs","vOTUs","mOTUs"),
                                                                  values=c("#2570AE","#9D5D2E","gray"))+
  scale_fill_manual(breaks=c("vPCs","vOTUs","mOTUs"),
                    values=c("#2570AE","#9D5D2E","gray"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+scale_x_continuous(breaks = seq(0,1,0.2))+ylim(-0.45,0)+
  stat_poly_eq(aes(color=Type,group=Type,label = paste(..eq.label..,p.value.label,sep = "~`,`~"),size=2),
               formula = y ~ poly(x),parse = TRUE)  
#Fig. 3e
data<-t(read.delim('beta.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE))
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

write.csv(dff,"fig.3e.csv")
#plot
data<-read_excel('E:/Tidal_time/R_for_github/Fig. 3/Fig. 3.xlsx')
ggtern(data,aes(1-beta,turnover,nestedness,fill=Type))+
  geom_point(alpha=0.5,size=2,shape=21)+ scale_fill_manual(breaks=c("vPCs","mOTUs","vOTUs"),
                                                           values=c("#2570AE","gray","#9D5D2E"))+
  theme_rgbw(base_size = 14, base_family = "")+
  theme(panel.border = element_rect(size=0.1))+
  labs(x ="",y = "",z = "",
       title = "")+ theme(legend.position ="none")
#Fig. 3f
data<-read_excel('E:/Tidal_time/R_for_github/Fig. 3/Fig. 3.xlsx')
data$Type <- factor(data$Type,levels=c('vPCs','mOTUs','vOTUs'))
ggboxplot(data,notch = F,bxp.errorbar = T, x='component', y='beta', fill = 'Type', palette = c("#2570AE","gray","#9D5D2E"),width = 0.5,size=0.3,outlier.size=0.2)+
  labs(x="",y="Beta diversity (community distance)",fill="Type")+
  theme_test()+theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
                     axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
                     axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
                     axis.text.y=element_text(size=16,colour = 'black'),
                     panel.border = element_rect(size=1.2),
                     legend.text = element_text(size=15),
                     legend.title = element_text(size=16))+
  theme(legend.position = "none")
