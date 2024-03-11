library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(vegan)
library(ggpubr)
library(nlme)
library(ggforce)
library(ggsci)
setwd('D:/Tidal_time/Rdata/Fig. 4/')
#Fig.4a
data<-read_excel("Fig. 4.xlsx")
ggplot(data, mapping = aes(x=Time,y=Pi))+
  geom_point(size = 3, alpha = 0.85,shape=21,fill="gray")+
  labs(x="Time (month)",y="Microdiversity (p value)")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(method = 'lm',se=T,size=1,fullrange=T,color='#539DDA',fill="#539DDA",alpha=0.2,formula = y ~ poly(x,2))+
  theme_test()+ theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+scale_x_continuous(breaks = seq(0,10,2))+ylim(0.0003,0.0018)+
stat_poly_eq(aes(label = paste(..rr.label..,p.value.label,sep = "~`,`~"),size=2),
             formula = y ~ poly(x,2),parse = TRUE)

#Fig. 4b
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
  stat_poly_eq(aes(label = paste(..rr.label..,p.value.label,sep = "~`,`~"),size=2),
               formula = y ~ x,parse = TRUE)
#Fig. 4c
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
        legend.title = element_text(size=16))+ylim(6,7)+xlim(0.0003,0.0017)+
stat_cor(method = "pearson",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
         label.x = )

#Fig. 4d
ggplot(data, mapping = aes(x=Time,y=pNpS))+
  geom_boxplot(aes(group=Time),color="black",width=1,outlier.color="white",position = position_dodge(0.8))+
  geom_jitter(aes(x=Time,y=pNpS),fill="gray",width=0.4,shape=21,size=2,alpha=0.7)+
  labs(x="Time (month)",y="Pnps")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(method = 'loess',se=T,size=1,fullrange=T,color='#539DDA',fill="#539DDA",alpha=0.2,aes(group=1))+
  theme_test()+ theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+scale_x_continuous(breaks = seq(0,10,2))+ylim(0.02,0.18)

#Fig.4e
ggplot(data, mapping = aes(x=Time,y=TjD))+
  geom_boxplot(aes(group=Time),color="black",width=1,outlier.color="white",position = position_dodge(0.8))+
  geom_jitter(aes(x=Time,y=TjD),fill="gray",width=0.4,shape=21,size=2,alpha=0.7)+
  labs(x="Time (month)",y="|Tajima's D|")+
  geom_smooth(method = 'loess',se=T,size=1,fullrange=T,color='#539DDA',fill="#539DDA",alpha=0.2,aes(group=1))+
  theme(axis.line = element_line(color="black"))+
  theme_test()+ 
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1.2),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+scale_x_continuous(breaks = seq(0,10,2))+ylim(0.5,2)
#Fig. 4f
sales <- c(70,51)
names<-c("a","b")
share<-sales/sum(sales)*100
data <- data.frame(
  sales,share,names)
ggplot()+
  geom_arc_bar(data=data,aes(x0 = 0, y0 = 0, r0 = 0, r = 1,amount=sales,explode=c(0,0.03),fill=names,color=names),stat="pie")+
  coord_fixed()+theme_void()+theme(legend.position = "none")+scale_color_manual(breaks=c("a","b"),values=c("lightgray","#66B3FF"))+
  scale_fill_manual(breaks=c("a","b"),values=c("gray","#66B3FF"))

data<-data.frame(variable=c(23,
                            14,
                            6,
                            4,
                            3), group = paste0("a", 1:5))
ggplot(data, aes(x = 3, y = variable, fill = group))+ geom_col() +
  coord_polar(theta = "y") +xlim(c(1, 4.5))+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("a1","a2","a3","a4","a5","a6"),values=c("Lightgray","#6FB7B7",'#84C1FF', '#FFD306', '#d3a4ff'))
