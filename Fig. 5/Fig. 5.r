library(vegan)
library(htmltools)
library(linkET)
library(Hmisc)
library(RColorBrewer)
library(ggplot2)
library(ggpmisc)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(randomForest)
library(rfUtilities)
library(rfPermute)
library(tidyverse)
library(patchwork)
library(reshape2)
library(tidyr)
library(psych)
#Fig.5a mantel test
env<-read.csv('Environment.csv',row.names = 1) # sample as row, env as column
vOTUs<-t(read.delim('vOTU.txt',row.names = 1))
mOTUs<-t(read.delim('host.txt',row.names = 1))
vPCs<-t(read.delim('PCs.txt',row.names = 1))
env<-apply(env[,1:10],2,as.numeric) #transfer numeric
#calculate distance
vOTUs<-vegdist(vOTUs,method = "bray")
vPCs<-vegdist(vPCs,method = "bray")
mOTUs<-vegdist(mOTUs,method = "bray")
#mantel test between env,geo,OTUs
#calculate distance
env<-vegdist(env,method = 'euclidean',na.rm = T) 
#partial mantel test
mantel(vOTUs,env,permutations = 9999,method="pearson",na.rm=T)

#mantel test for each env with OTUs
mantel<-list()
for (i in 1:10){
  mantel[[i]]<-vegdist(env[,i],method='euclidean',upper = FALSE,na.rm = T)
}
names(mantel)<-colnames(env)[1:10]

#mantel test for env and vOTU
r<-c()
p<-c()
for (i in 1:10){
  r[i]<-mantel(vOTUs,mantel[[i]],permutations = 999,method="spearman",na.rm=T)$statistic
  p[i]<-mantel(vOTUs,mantel[[i]],permutations = 999,method="spearman",na.rm=T)$signif
}
p<-as.data.frame(p)
r<-as.data.frame(r)
rownames(p)<-colnames(env)
rownames(r)<-colnames(env)

#mantel test for env and vPCs
r1<-c()
p1<-c()
for (i in 1:10){
  r1[i]<-mantel(vPCs,mantel[[i]],permutations = 999,method="spearman",na.rm=T)$statistic
  p1[i]<-mantel(vPCs,mantel[[i]],permutations = 999,method="spearman",na.rm=T)$signif
}
p1<-as.data.frame(p1)
r1<-as.data.frame(r1)
rownames(p1)<-colnames(env)
rownames(r1)<-colnames(env)

#mantel test for env and mOTU
r2<-c()
p2<-c()
for (i in 1:10){
  r2[i]<-mantel(mOTUs,mantel[[i]],permutations = 999,method="spearman",na.rm=T)$statistic
  p2[i]<-mantel(mOTUs,mantel[[i]],permutations = 999,method="spearman",na.rm=T)$signif
}
p2<-as.data.frame(p2)
r2<-as.data.frame(r2)
rownames(p2)<-colnames(env)
rownames(r2)<-colnames(env)

write.csv(p,file ="p.csv", row.names = T, quote =FALSE)    
write.csv(r,file ="r.csv", row.names = T, quote =FALSE) 
write.csv(p1,file ="p1.csv", row.names = T, quote =FALSE)    
write.csv(r1,file ="r1.csv", row.names = T, quote =FALSE) 
write.csv(p2,file ="p2.csv", row.names = T, quote =FALSE)    
write.csv(r2,file ="r2.csv", row.names = T, quote =FALSE)

mantel<-read.delim('mantel.txt')
qcorrplot(rcorr(env), type = "upper", diag = F, grid_size = 0.4,grid_col = "gray") +
  geom_square(linetype=0) +
  geom_couple(aes(colour = p_value, size = r_value), data = mantel, curvature = 0.1) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),
                       limits = c(-1, 1),
                       breaks = seq(-1,1,0.5))+
  scale_size_manual(values = c(0.2, 0.4, 1))+ 
  scale_colour_manual(values=c("#87CEEB","gray","#9ACD32")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 1.5), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))+
  geom_mark(
    only_mark = T,
    size = 5, 
    sig_level = c(0.05, 0.01, 0.001), 
    sig_thres = 0.05
  ) 

#Fig. 5b
env<-read.csv('Environment.csv',row.names = 1)# sample as row, env as column
env<-apply(env[,1:10],2,as.numeric)
env<-vegdist(env,method = 'euclidean',na.rm =T) 
distv_env<-as.numeric(env)

vPCs<-read.delim('vPCs.txt',row.names = 1)
vOTUs<-read.delim('vOTU.txt',row.names = 1)
mOTUs<-read.delim('host.txt',row.names = 1)
distv<-vegdist(t(mOTUs),method = 'bray')
distv_num<-1-as.numeric(distv)

data_new<-data.frame(distv_env,distv_num)
write.csv(data_new,file='env-mOTUs.csv',row.names = T)
data1<-read.csv('env-com.csv')# sample as row, env as column

ggplot(data1, mapping = aes(x=distv_env,y=distv_num))+
  geom_point(aes(color=Type),size = 2, alpha = 0.7,shape=16)+
  labs(x="Environmental heterogeneity",y="Community similarity")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(aes(color=Type),method = 'lm',se=T,size=1.3,fullrange=F,fill='lightgray')+
  theme_test()+ theme(legend.position ="none")+scale_color_manual(breaks=c("vPCs","vOTUs","mOTUs"),
                                                                  values=c("#2570AE","#9D5D2E","gray"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+ylim(0.15,1.19)+
stat_poly_eq(aes(color=Type,label = paste(..eq.label..,p.value.label, sep = "~`,`~")),
             formula = y ~ x,parse = TRUE)
 
#VPA 
votu_matrix<-t(read.delim('vOTU.txt',row.names = 1))
votu_matrix<-t(read.delim('vPCs.txt',row.names = 1))
votu_matrix<- decostand(votu_matrix, method = 'hellinger') 

env<-read.csv('Environment.csv',row.names = 1) 
print(row.names(votu_matrix)==row.names(env))
res_all<- cca(votu_matrix~.,data=env)
res_null<-cca(votu_matrix~1,data=env)
env_frwd<-ordistep(res_null,scope=formula(res_all),direction = 'forward',permutations = how(nperm = 999)) #fit model<0.05
env_frwd #no-linear envs

envdat_raw<-read.csv('Environment.csv',row.names = 1) #factors in no-linear geo/env/VMR
print(row.names(votu_matrix)==row.names(envdat_raw))
env<-envdat_raw[,1:3]
cca_all<-cca(votu_matrix,env) #obtain one or two factors common proportion
RsquareAdj(cca_all)
decorana(votu_matrix)

#Fig.5c
myro <- read.table("map.txt",header = T,row.names = 1,sep="\t",check.names = F)
spearman <- corr.test(myro[,1:10], myro[,11:44], method = 'spearman', adjust = 'none')
r <- data.frame(spearman$r) 
r$myro <- rownames(r) 
r <- melt(r, id = 'myro')
spearman <- cbind(r)
spearman
ggplot() +
  geom_tile(data = spearman, aes(x = myro, y = variable, fill = value),color="#9D9D9D") + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),
                       limits = c(-1, 1),
                       breaks = seq(-1,1,0.5))+
  theme_classic()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =1,angle = 45,size=10,colour = 'black'),
        axis.text.y=element_text(size=10,colour = 'black'),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16)) +
  labs(y = '', x = '', fill = 'Correlation')
#importance for all factors
set.seed(123)
Sat_forest <- randomForest(Sat~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
CysC_forest <- randomForest(CysC~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
CysH_forest <- randomForest(CysH~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
DsrC_forest <- randomForest(DsrC~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
PL1_forest <- randomForest(PL1~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
GH1_forest <- randomForest(GH1~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
GH16_forest <- randomForest(GH16~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
GH18_forest <- randomForest(GH18~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
GH26_forest <- randomForest(GH26~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
GH5_forest <- randomForest(GH5~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
GH51_forest <- randomForest(GH51~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
GT1_forest <- randomForest(GT1~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
GT17_forest <- randomForest(GT17~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
GT2_forest <- randomForest(GT2~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
GT4_forest <- randomForest(GT4~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
GT6_forest <- randomForest(GT6~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
PhoA_forest <- randomForest(PhoA~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
PhoD_forest <- randomForest(PhoD~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Acidobacteria_forest <- randomForest(Acidobacteria~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Actinobacteria_forest <- randomForest(Actinobacteria~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Alphaproteobacteria_forest <- randomForest(Alphaproteobacteria~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Bacteroidetes_forest <- randomForest(Bacteroidetes~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Betaproteobacteria_forest <- randomForest(Betaproteobacteria~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
candidate_division_WWE3_forest <- randomForest(candidate_division_WWE3~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Candidatus_Aminicenantes_forest <- randomForest(Candidatus_Aminicenantes~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Candidatus_Latescibacteria_forest <- randomForest(Candidatus_Latescibacteria~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Chloroflexi_forest <- randomForest(Chloroflexi~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Deltaproteobacteria_forest <- randomForest(Deltaproteobacteria~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Firmicutes_forest <- randomForest(Firmicutes~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Gammaproteobacteria_forest <- randomForest(Gammaproteobacteria~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Gemmatimonadetes_forest <- randomForest(Gemmatimonadetes~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Other_Proteobacteria_forest <- randomForest(Other_Proteobacteria~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Planctomycetes_forest <- randomForest(Planctomycetes~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Thermodesulfobacteria_forest <- randomForest(Thermodesulfobacteria~Temperature+Salinity+pH+NH4+NO3+NO2+TN+TP+SO4+TOC, data =myro, importance = TRUE, ntree = 500)
Sat<- data.frame(importance(Sat_forest), check.names = FALSE)
CysC<- data.frame(importance(CysC_forest), check.names = FALSE)
CysH<- data.frame(importance(CysH_forest), check.names = FALSE)
DsrC<- data.frame(importance(DsrC_forest), check.names = FALSE)
PL1<- data.frame(importance(PL1_forest), check.names = FALSE)
GH1<- data.frame(importance(GH1_forest), check.names = FALSE)
GH16<- data.frame(importance(GH16_forest), check.names = FALSE)
GH18<- data.frame(importance(GH18_forest), check.names = FALSE)
GH26<- data.frame(importance(GH26_forest), check.names = FALSE)
GH5<- data.frame(importance(GH5_forest), check.names = FALSE)
GH51<- data.frame(importance(GH51_forest), check.names = FALSE)
GT1<- data.frame(importance(GT1_forest), check.names = FALSE)
GT17<- data.frame(importance(GT17_forest), check.names = FALSE)
GT2<- data.frame(importance(GT2_forest), check.names = FALSE)
GT4<- data.frame(importance(GT4_forest), check.names = FALSE)
GT6<- data.frame(importance(GT6_forest), check.names = FALSE)
PhoA<- data.frame(importance(PhoA_forest), check.names = FALSE)
PhoD<- data.frame(importance(PhoD_forest), check.names = FALSE)
Acidobacteria<- data.frame(importance(Acidobacteria_forest), check.names = FALSE)
Actinobacteria<- data.frame(importance(Actinobacteria_forest), check.names = FALSE)
Alphaproteobacteria<- data.frame(importance(Alphaproteobacteria_forest), check.names = FALSE)
Bacteroidetes<- data.frame(importance(Bacteroidetes_forest), check.names = FALSE)
Betaproteobacteria<- data.frame(importance(Betaproteobacteria_forest), check.names = FALSE)
candidate_division_WWE3<- data.frame(importance(candidate_division_WWE3_forest), check.names = FALSE)
Candidatus_Aminicenantes<- data.frame(importance(Candidatus_Aminicenantes_forest), check.names = FALSE)
Candidatus_Latescibacteria<- data.frame(importance(Candidatus_Latescibacteria_forest), check.names = FALSE)
Chloroflexi<- data.frame(importance(Chloroflexi_forest), check.names = FALSE)
Deltaproteobacteria<- data.frame(importance(Deltaproteobacteria_forest), check.names = FALSE)
Firmicutes<- data.frame(importance(Firmicutes_forest), check.names = FALSE)
Gammaproteobacteria<- data.frame(importance(Gammaproteobacteria_forest), check.names = FALSE)
Gemmatimonadetes<- data.frame(importance(Gemmatimonadetes_forest), check.names = FALSE)
Other_Proteobacteria<- data.frame(importance(Other_Proteobacteria_forest), check.names = FALSE)
Planctomycetes<- data.frame(importance(Planctomycetes_forest), check.names = FALSE)
Thermodesulfobacteria<- data.frame(importance(Thermodesulfobacteria_forest), check.names = FALSE)

Importance<- data.frame(cbind(Sat$`%IncMSE`, CysC$`%IncMSE`,CysH$`%IncMSE`,DsrC$`%IncMSE`,PL1$`%IncMSE`,
                              GH1$`%IncMSE`,GH16$`%IncMSE`,GH18$`%IncMSE`, GH26$`%IncMSE`,GH5$`%IncMSE`,GH51$`%IncMSE`,GT1$`%IncMSE`,
                              GT17$`%IncMSE`,GT2$`%IncMSE`,GT4$`%IncMSE`, GT6$`%IncMSE`,PhoA$`%IncMSE`,PhoD$`%IncMSE`,
                              Acidobacteria$`%IncMSE`, Actinobacteria$`%IncMSE`,Alphaproteobacteria$`%IncMSE`,Bacteroidetes$`%IncMSE`,Betaproteobacteria$`%IncMSE`,
                              candidate_division_WWE3$`%IncMSE`,Candidatus_Aminicenantes$`%IncMSE`,Candidatus_Latescibacteria$`%IncMSE`, Chloroflexi$`%IncMSE`,Deltaproteobacteria$`%IncMSE`,Firmicutes$`%IncMSE`,
                              Gammaproteobacteria$`%IncMSE`,Gemmatimonadetes$`%IncMSE`,Other_Proteobacteria$`%IncMSE`, Planctomycetes$`%IncMSE`,Thermodesulfobacteria$`%IncMSE`))
Importance[Importance<1] <- 0
Importance[is.na(Importance)] <- 0
write.csv(Importance,file = "Importance.csv")

impdata <- read.csv("Importance.CSV",header = TRUE)

data1=melt(impdata ,
           id.vars='Items', 
           measure.vars=measure_name,
           variable.name = "sample", 
           value.name = "expr")
data1$sample <- factor(data1$sample,levels=c('Temperature','Salinity','pH','NH4','NO3','NO2','TN','TP','SO4','TOC'))
spearman$myro<-factor(spearman$myro,levels=c('Temperature','Salinity','pH','NH4','NO3','NO2','TN','TP','SO4','TOC'))
ggplot() +
  geom_tile(data = spearman, aes(x = myro, y = variable, fill = value),color="#ADADAD") + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),
                       limits = c(-1, 1),
                       breaks = seq(-1,1,0.5))+theme_bw()+
  theme(panel.border=element_blank(),panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =1,angle = 45,size=10,colour = 'black'),
        axis.text.y=element_text(size=10,colour = 'black'),
        legend.text = element_text(size=8.5),
        legend.title = element_text(size=10),axis.ticks = element_blank()) +
  labs(y = '', x = '', fill = 'Pearson’s r')+
  geom_point(data = data1, aes(x = sample, y = Items, size = expr), shape = 21) +
  scale_size_continuous(range = c(1,6)) +
  labs(size = 'Importance')
