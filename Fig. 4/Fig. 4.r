library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(vegan)
library(ggpubr)
#Fig.4a
data<-read_excel("E:/Tidal_time/R_for_github/Fig. 4/Fig. 4.xlsx")
ggplot(data, mapping = aes(x=Time,y=Pi))+
  geom_point(size = 3, alpha = 0.85,shape=21,fill="gray")+
  labs(x="Time (month)",y="Microdiversity (p value)")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(method = 'lm',se=T,size=1,fullrange=T,color='#539DDA',fill="#539DDA",alpha=0.2)+
  theme_test()+ theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+scale_x_continuous(breaks = seq(0,10,2))+ylim(0.0003,0.0018)+
stat_cor(method = "pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
         label.x = )
#Fig. 4b
data<-read_excel("E:/Tidal_time/R_for_github/Fig. 4/Fig. 4.xlsx")
ggplot(data, mapping = aes(x=Time,y=gene_pi))+
  geom_point(size = 3, alpha = 0.85,shape=21,fill="gray")+
  labs(x="Time (month)",y="Microdiversity (p value)")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(method = 'lm',se=T,size=1,fullrange=T,color='#539DDA',fill="#539DDA",alpha=0.2)+
  theme_test()+ theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+scale_x_continuous(breaks = seq(0,10,2))+ylim(0.0003,0.0018)+
stat_cor(method = "pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
         label.x = )  
#Fig. 4c
data<-read_excel("E:/Tidal_time/R_for_github/Fig. 4/Fig. 4.xlsx")
ggplot(data, mapping = aes(x=Pi,y=diversity))+
  geom_point(size = 3, alpha = 0.85,shape=21,fill="gray")+
  labs(x="Microdiversity (p value)",y="Macrodiversity (shannon index)")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(method = 'lm',se=T,size=1,fullrange=T,color='#539DDA',fill="#539DDA",alpha=0.2)+
  theme_test()+ 
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+ylim(6,7)+xlim(0.0003,0.0017)
stat_cor(method = "pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
         label.x = )+ylim(6,7)

#Fig. 4d
data<-read_excel("E:/Tidal_time/R_for_github/Fig. 4/Fig. 4.xlsx")
ggplot(data, mapping = aes(x=Time,y=pNpS))+
  geom_point(size = 3, alpha = 0.85,shape=21,fill="gray")+
  labs(x="Time (month)",y="Gene selection pressure (pN/pS)")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(method = 'lm',se=T,size=1,fullrange=T,color='#539DDA',fill="#539DDA",alpha=0.2)+
  theme_test()+ 
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+ylim(0.02,0.18)+scale_x_continuous(breaks = seq(0,10,2))
stat_cor(method = "pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
         label.x = )

#Fig.4e
data<-read_excel("E:/Tidal_time/R_for_github/Fig. 4/Fig. 4.xlsx")
ggplot(data, mapping = aes(x=Time,y=TjD))+
  geom_point(size = 3, alpha = 0.85,shape=21,fill="gray")+
  labs(x="Time (month)",y="|Tajima's D|")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(method = 'lm',se=T,size=1,fullrange=T,color='#539DDA',fill="#539DDA",alpha=0.2)+
  theme_test()+ 
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+scale_x_continuous(breaks = seq(0,10,2))+ylim(0.5,2)+
  stat_cor(method = "pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = )
