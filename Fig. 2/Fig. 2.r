library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(vegan)
library(ggpmisc)
library(gcookbook)
library(nlme)
setwd('D:/Tidal_time/Rdata/Fig. 2')
#Fig.2b
data<-read_excel("Fig. 2.xlsx")
data$Date<-as.Date(data$Date)
data<-melt(data,id="Date")
data$Type=factor(data$Type, levels=c('Alphaproteobacteria',"Betaproteobacteria",'Gammaproteobacteria',"Deltaproteobacteria","Other Proteobacteria" ,'Actinobacteria','Bacteroidetes',"Chloroflexi",
                                     "Acidobacteria","Firmicutes","Planctomycetes",
                                     "Gemmatimonadetes","Candidatus Latescibacteria","Candidatus Aminicenantes","Thermodesulfobacteria","candidate division WWE3",
                                     "Other Archaea","Other Bacteria"))
ggplot(data, aes(x =Date, y = value,fill=Type))+
  geom_area(position="stack",alpha=0.9)+
  theme_test()+
  scale_x_date(date_labels = "%Y",date_breaks = "1 day")+
  scale_fill_manual(values = c('#2570AE','#2C85CE', '#539DDA',"#8DB2F1","#BBFFFF",'#43CD80', '#9370DB', '#FFD700',
                               "#5F9EA0","#9D5D2E","#CD6839",
                               "#EE5C42","#EE3A8C","#A79F14","#BBC4FF",'#FF9632',
                               "black","gray"))+labs(x="",y="")+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        panel.background = element_blank())+scale_x_date(date_breaks = "1 day")+theme(legend.position = "right")+guides(fill = guide_legend(ncol = 1, byrow = TRUE))

#Fig.2c
data<-read_excel("Fig. 2.xlsx")
ggplot(data, aes(x=Host,y=Virus))+
  geom_point(fill='gray',size = 3, alpha = 0.85,shape=21)+
  labs(x="Normalized host microbial abundance",y="Normalized viral abundance")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(color="#539DDA",method = 'lm',se=F,size=1,fullrange=T,formula = y ~ poly(x,2))+
  geom_smooth(color="#6C6C6C",method = 'lm',se=F,size=1,fullrange=T,formula = y ~ x,linetype=2)+
  theme_bw()+ 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        axis.title.x = element_text(size=16.5),axis.title.y = element_text(size=16.5),
        axis.text.x = element_text(hjust =0.5,size=15,colour = 'black'),
        axis.text.y=element_text(size=15,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+ylim(2500,13000)+xlim(1800,8000)+
  stat_poly_eq(aes(label = paste(..rr.label..,p.value.label,sep = "~`,`~"),size=2),
               formula = y ~ poly(x,2),parse = TRUE)+
  stat_poly_eq(aes(label = paste(..rr.label..,p.value.label,sep = "~`,`~"),size=2),
               formula = y ~ x,parse = TRUE)
# calculate AICC to choose lm model
object2<-gls(Virus~poly(Host,2), data=data)
object1<-gls(Virus~poly(Host), data=data)
summary(object)

#Fig.2d
data<-read_excel("Fig. 2.xlsx")
ggplot(data, mapping = aes(x=Host,y=Virus))+
  geom_point(aes(fill=Type),size = 2.5, alpha = 0.8,shape=21)+
  labs(x="Normalized host microbial abundance",y="Normalized viral abundance")+
  geom_smooth(aes(color=Type),method = 'lm',se=F,size=0.7,fullrange=F,alpha=0.8)+
  theme_bw()+ 
  scale_color_manual(breaks=c('Alphaproteobacteria',"Betaproteobacteria",'Gammaproteobacteria',"Deltaproteobacteria","Other Proteobacteria" ,'Actinobacteria','Bacteroidetes',"Chloroflexi",
                              "Acidobacteria","Firmicutes","Planctomycetes",
                              "Gemmatimonadetes","Candidatus Latescibacteria","Candidatus Aminicenantes","Thermodesulfobacteria","candidate division WWE3",
                              "Other Archaea","Other Bacteria"),
                     values=c('#2570AE','#2C85CE', '#539DDA',"#8DB2F1","#BBFFFF",'#43CD80', '#9370DB', '#FFD700',
                              "#5F9EA0","#9D5D2E","#CD6839",
                              "#EE5C42","#EE3A8C","#A79F14","#BBC4FF",'#FF9632',
                              "black","gray"))+
  scale_fill_manual(breaks=c('Alphaproteobacteria',"Betaproteobacteria",'Gammaproteobacteria',"Deltaproteobacteria","Other Proteobacteria" ,'Actinobacteria','Bacteroidetes',"Chloroflexi",
                             "Acidobacteria","Firmicutes","Planctomycetes",
                             "Gemmatimonadetes","Candidatus Latescibacteria","Candidatus Aminicenantes","Thermodesulfobacteria","candidate division WWE3",
                             "Other Archaea","Other Bacteria"),
                    values=c('#2570AE','#2C85CE', '#539DDA',"#8DB2F1","#BBFFFF",'#43CD80', '#9370DB', '#FFD700',
                             "#5F9EA0","#9D5D2E","#CD6839",
                             "#EE5C42","#EE3A8C","#A79F14","#BBC4FF",'#FF9632',
                             "black","gray"))+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        axis.title.x = element_text(size=19),axis.title.y = element_text(size=19),
        axis.text.x = element_text(hjust =0.5,size=16.5,colour = 'black'),
        axis.text.y=element_text(size=16.5,colour = 'black'),
        panel.border = element_rect(size=1),
        legend.text = element_text(size=12),
        legend.title = element_blank())+
  stat_poly_eq(aes(label = paste(..rr.label..,p.value.label,sep = "~`,`~"),size=2),
               formula = y ~ x,parse = TRUE)
# calculate AICC to choose lm model
object2<-gls(Virus~poly(Host,2), data=data)
object1<-gls(Virus~poly(Host), data=data)
summary(object)

#Fig.2e
data<-read_excel("Fig. 2.xlsx")
ggplot(data, mapping = aes(x=Host,y=Viruses))+
  geom_point(aes(fill=lifestyle),size = 3, alpha = 0.85,shape=21)+
  labs(x="Normalized host microbial abundance",y="Relative abundance of viruses (%)")+
  geom_smooth(aes(color=lifestyle),method = 'lm',se=F,size=1.1,fullrange=F)+
  theme_test()+ theme(legend.position ="none")+scale_fill_manual(breaks=c("Lytic","Lysogenic"),
                                                                 values=c("#2570AE","#9D5D2E"))+
  scale_color_manual(breaks=c("Lytic","Lysogenic"),values=c("#2570AE","#9D5D2E"))+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        axis.title.x = element_text(size=16.5),axis.title.y = element_text(size=16.5),
        axis.text.x = element_text(hjust =0.5,size=15,colour = 'black'),
        axis.text.y=element_text(size=15,colour = 'black'),
        panel.border = element_rect(size=1),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))+ylim(30,80)+xlim(5000,13000)+
  stat_poly_eq(aes(color=lifestyle,label = paste(..rr.label..,p.value.label,sep = "~`,`~"),size=2),
               parse = TRUE)
