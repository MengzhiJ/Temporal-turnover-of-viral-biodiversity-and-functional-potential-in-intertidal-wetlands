#Fig. 5a NST ratio
wd="D:/Tidal_time/Rdata/Fig. 5"
save.wd="D:/Tidal_time/Rdata/Fig. 5/NST-vOTUs" 
com.file="D:/Tidal_time/Rdata/Fig. 5/vOTU.txt" 
group.file="D:/Tidal_time/Rdata/Fig. 5/group.txt"
library(vegan)
library(ape)
library(iCAMP)
library(parallel)
library(bigmemory)
library(NST) 
library(permute)
library(DirichletReg)
setwd(wd)
comm=t(read.delim(com.file, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE))
group=read.delim(group.file, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

samp.ck=NST::match.name(rn.list=list(comm=comm,group=group))
comm=samp.ck$comm
comm=comm[,colSums(comm)>0,drop=FALSE]
group=samp.ck$group

groupi=group[,1,drop=FALSE]
prefixi=paste0(prefix,"")
meta.groupi=groupi
dist.method="bray"

setwd(save.wd)
tnst=tNST(comm=comm, group=groupi, meta.group=meta.groupi, meta.com=NULL,
          dist.method=dist.method, abundance.weighted=TRUE, rand=1000,
          output.rand=TRUE, nworker=28, LB=FALSE, null.model="PF",dirichlet=F,
          between.group=T, SES=TRUE, RC=TRUE)

write.table(tnst$index.grp,file = paste0(prefixi,".tNST.summary.csv"), quote = FALSE,sep = ",")
write.table(tnst$index.pair.grp,file = paste0(prefixi,".tNST.pairwise.csv"),quote = FALSE,sep = ",")
write.table(tnst$index.pair,file = paste0(prefixi,".tNST.pair.csv"),quote = FALSE,sep = ",")

data<-read_excel("Fig. 5.xlsx")
data$Time <- factor(data$Time,levels=c('Aug','Oct','Dec','Feb','April','June'))
data$Type <- factor(data$Type,levels=c('vOTUs','Host microbes'))
ggplot(data,aes(x=Time,y=NST))+
  stat_summary(aes(fill=Type), alpha=0.65,fun = mean,geom="bar",color="black",width=0.55,size=0.3)+
  stat_summary(fun.data = mean_se,geom="errorbar",width=.08)+
  labs(x="",y="Normalized stochastic ratio (%)",fill="Type")+
  scale_fill_manual(breaks=c("vOTUs","Host microbes"),
                    values=c("#9D5D2E","gray"))+
  theme_test()+theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
                     axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
                     axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
                     axis.text.y=element_text(size=12,colour = 'black'),
                     panel.border = element_rect(size=1),
                     legend.text = element_text(size=10),
                     legend.title = element_text(size=10))+ylim(0,100)
scale_y_continuous(breaks = seq(0,80,10))+
  theme(legend.position = "none")+facet_wrap(~Type)


#Fig.5b mantel test
library(vegan)
setwd('D:/Tidal_time/Rdata/Fig. 5')
env<-read.csv('Environment.csv',row.names = 1) # sample as row, env as column
vOTUs<-t(read.delim('vOTU.txt',row.names = 1))
vPCs<-t(read.delim('vPCs.txt',row.names = 1))
host<-t(read.delim('host.txt',row.names = 1))
env<-apply(env[,1:10],2,as.numeric) #transfer numeric
#calculate distance
vOTUs<-vegdist(vOTUs,method = "bray")
vPCs<-vegdist(vPCs,method = "bray")
host<-vegdist(HOST,method = "bray")
#mantel test between env,geo,OTUs
#calculate distance
env<-vegdist(env,method = 'euclidean',na.rm = T) 

#partial mantel test
mantel(vOTUs,env,permutations = 9999,method="spearman",na.rm=T)
mantel(VPCs,env,permutations = 9999,method="spearman",na.rm=T)

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


#mantel test for env and host
r2<-c()
p2<-c()
for (i in 1:10){
  r2[i]<-mantel(host,mantel[[i]],permutations = 999,method="spearman",na.rm=T)$statistic
  p2[i]<-mantel(host,mantel[[i]],permutations = 999,method="spearman",na.rm=T)$signif
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

#Fig. 5b plot
library(htmltools)
library(linkET)
library(Hmisc)
library(RColorBrewer)
library(ggplot2)
mantel<-read.delim('mantel.txt')
qcorrplot(rcorr(env,type=c("spearman")), type = "upper", diag = F, grid_size = 0.4,grid_col = "gray") +
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
         fill = guide_colorbar(title = "Spearman's rho", order = 3))

#Fig. S9
library(vegan)
library(ggpmisc)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
env<-read.csv('Environment.csv',row.names = 1)# sample as row, env as column
env<-apply(env[,1:10],2,as.numeric)
env<-vegdist(env,method = 'euclidean',na.rm =T) 
distv_env<-as.numeric(env)

vOTUs<-read.delim('vOTU.txt',row.names = 1)
host<-read.delim('host.txt',row.names = 1)
distv<-vegdist(t(host),method = 'bray')
distv_num<-1-as.numeric(distv)

data_new<-data.frame(distv_env,distv_num)
write.csv(data_new,file='env-host.csv',row.names = T)
data<-read.csv('env-com.CSV')# sample as row, env as column

ggplot(data, mapping = aes(x=distv_env,y=distv_num))+
  geom_point(aes(color=Type),size = 2, alpha = 0.7,shape=16)+
  labs(x="Environmental heterogeneity",y="Community similarity")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(aes(color=Type),method = 'lm',se=T,size=1.3,fullrange=F,fill='lightgray')+
  theme_test()+ theme(legend.position ="none")+scale_color_manual(breaks=c("vPCs","vOTUs","Host microbes"),
                                                                  values=c("#2570AE","#9D5D2E","gray"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=18.5),axis.title.y = element_text(size=18.5),
        axis.text.x = element_text(hjust =0.5,size=16,colour = 'black'),
        axis.text.y=element_text(size=16,colour = 'black'),
        panel.border = element_rect(size=1),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16))+
stat_poly_eq(aes(color=Type,label = paste(..eq.label..,p.value.label, sep = "~`,`~")),
             formula = y ~ x,parse = TRUE)
 
#VPA 
library(vegan)
library(ggplot2)
votu_matrix<-t(read.delim('vOTU.txt',row.names = 1))
host_matrix<-t(read.delim('host.txt',row.names = 1))
votu_matrix<- decostand(votu_matrix, method = 'hellinger') 
host_matrix<- decostand(host_matrix, method = 'hellinger') 
env<-read.csv('Environment.csv',row.names = 1) 

print(row.names(host_matrix)==row.names(env))
res_all<- cca(host_matrix~.,data=env)
res_null<-cca(host_matrix~1,data=env)
env_frwd<-ordistep(res_null,scope=formula(res_all),direction = 'forward',permutations = how(nperm = 999)) #fit model<0.05
env_frwd #no-linear envs
envdat_raw<-read.csv('Environment.csv',row.names = 1) #factors in no-linear geo/env/VMR
print(row.names(host_matrix)==row.names(envdat_raw))
env<-envdat_raw[,1:3]
cca_all<-cca(host_matrix,env) #obtain one or two factors common proportion
RsquareAdj(cca_all)

#Fig.5c
library(vegan)
library(randomForest)
library(rfUtilities)
library(rfPermute)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(reshape2)
library(tidyr)
library(ggplot2)
library(psych)
setwd('D:/Tidal_time/Rdata/Fig. 5/random')
myro <- read.table("map.txt",header = T,row.names = 1,sep="\t",check.names = F)
#check model significance (screened)
otu_forest.pval <- a3(Sat~Temperature+Salinity+NO2+TOC, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(CysC~Temperature+NO2+SO4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(CysH~Temperature+Salinity+TN+TP+NH4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(DsrC~pH+TN, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(PL1~NH4+NO2+TN, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(GH1~Temperature+NH4+TN+pH, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(GH16~Temperature+Salinity+pH+SO4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(GH18~Temperature+Salinity+NO2, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(GH26~Temperature+Salinity+NH4+NO3+NO2+SO4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(GH5~NO3+NO2+TN+Salinity+NH4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(GH51~Temperature+Salinity+NO3+SO4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(GT1~NO3+TP, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(GT17~Temperature+NH4+NO3+NO2+SO4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(GT2~Temperature+Salinity+pH+NO3+NO2+TP+SO4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(GT6~Temperature+Salinity+pH+NH4+TN, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(PhoA~Temperature+NH4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(PhoD~pH+NO3+TN, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Acidobacteria~SO4+TOC+TN, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Actinobacteria~Temperature+Salinity+TOC, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Alphaproteobacteria~Temperature+pH+TP+SO4+TOC, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Bacteroidetes~Temperature+Salinity+pH+NH4+NO3+TN+TP, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Betaproteobacteria~Temperature+TP+SO4+TOC, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(candidate_division_WWE3~Temperature+SO4+TP, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Candidatus_Aminicenantes~Temperature+pH+TN+TP, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Candidatus_Latescibacteria~Temperature+pH+NH4+TN+SO4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Chloroflexi~Temperature+Salinity+NO2+SO4+NH4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Deltaproteobacteria~Salinity+pH+SO4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Firmicutes~Temperature+NH4+NO2+SO4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Gammaproteobacteria~Temperature+Salinity+TN, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Gemmatimonadetes~Temperature+Salinity+pH+NO2+TN+SO4, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
Otu_forest.pval <- a3(Other_Proteobacteria~Temperature+pH+NO3+SO4+TOC, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Planctomycetes~Temperature+NO3+NO2+SO4+TOC, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))
otu_forest.pval <- a3(Thermodesulfobacteria~Temperature+Salinity+NO3+SO4+TOC, data =myro, model.fn = randomForest, p.acc = 0.05, model.args = list(importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10))

#coefficient
spearman <- corr.test(myro[,1:10], myro[,11:43], method = 'spearman', adjust = 'none')
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
Sat_forest <- randomForest(Sat~Temperature+Salinity+NO2+TOC, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
CysC_forest <- randomForest(CysC~Temperature+NO2+SO4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
CysH_forest <- randomForest(CysH~Temperature+Salinity+TN+TP+NH4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
DsrC_forest <- randomForest(DsrC~pH+TN, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
PL1_forest <- randomForest(PL1~NH4+NO2+TN, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
GH1_forest <- randomForest(GH1~Temperature+NH4+TN+pH, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
GH16_forest <- randomForest(GH16~Temperature+Salinity+pH+SO4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
GH18_forest <- randomForest(GH18~Temperature+Salinity+NO2, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
GH26_forest <- randomForest(GH26~Temperature+Salinity+NH4+NO3+NO2+SO4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
GH5_forest <- randomForest(GH5~NO3+NO2+TN+Salinity+NH4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
GH51_forest <- randomForest(GH51~Temperature+Salinity+NO3+SO4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
GT1_forest <- randomForest(GT1~NO3+TP, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
GT17_forest <- randomForest(GT17~Temperature+NH4+NO3+NO2+SO4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
GT2_forest <- randomForest(GT2~Temperature+Salinity+pH+NO3+NO2+TP+SO4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
GT6_forest <- randomForest(GT6~Temperature+Salinity+pH+NH4+TN, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
PhoA_forest <- randomForest(PhoA~Temperature+NH4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
PhoD_forest <- randomForest(PhoD~pH+NO3+TN, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Acidobacteria_forest <- randomForest(Acidobacteria~SO4+TOC+TN, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Actinobacteria_forest <- randomForest(Actinobacteria~Temperature+Salinity+TOC, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Alphaproteobacteria_forest <- randomForest(Alphaproteobacteria~Temperature+pH+TP+SO4+TOC, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Bacteroidetes_forest <- randomForest(Bacteroidetes~Temperature+Salinity+pH+NH4+NO3+TN+TP, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Betaproteobacteria_forest <- randomForest(Betaproteobacteria~Temperature+TP+SO4+TOC, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
candidate_division_WWE3_forest <- randomForest(candidate_division_WWE3~Temperature+SO4+TP, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Candidatus_Aminicenantes_forest <- randomForest(Candidatus_Aminicenantes~Temperature+pH+TN+TP, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Candidatus_Latescibacteria_forest <- randomForest(Candidatus_Latescibacteria~Temperature+pH+NH4+TN+SO4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Chloroflexi_forest <- randomForest(Chloroflexi~Temperature+Salinity+NO2+SO4+NH4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Deltaproteobacteria_forest <- randomForest(Deltaproteobacteria~Salinity+pH+SO4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Firmicutes_forest <- randomForest(Firmicutes~Temperature+NH4+NO2+SO4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Gammaproteobacteria_forest <- randomForest(Gammaproteobacteria~Temperature+Salinity+TN, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Gemmatimonadetes_forest <- randomForest(Gemmatimonadetes~Temperature+Salinity+pH+NO2+TN+SO4, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Other_Proteobacteria_forest <- randomForest(Other_Proteobacteria~Temperature+pH+NO3+SO4+TOC, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Planctomycetes_forest <- randomForest(Planctomycetes~Temperature+NO3+NO2+SO4+TOC, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)
Thermodesulfobacteria_forest <- randomForest(Thermodesulfobacteria~Temperature+Salinity+NO3+SO4+TOC, data =myro, importance = TRUE, tree = 1000, nrep = 1000, num.cores = 10)

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
write.csv(Importance,file = "Importance.csv")

impdata <- read.csv("Importance.CSV",header = TRUE)
data1=melt(impdata,
           id.vars='Items',
           variable.name = "sample", 
           value.name = "expr")
data1$sample <- factor(data1$sample,levels=c('Temperature','Salinity','pH','NH4','NO3','NO2','TN','TP','SO4','TOC'))
spearman$myro<-factor(spearman$myro,levels=c('Temperature','Salinity','pH','NH4','NO3','NO2','TN','TP','SO4','TOC'))
p1<-ggplot() +
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
  labs(y = '', x = '', fill = 'Spearman rho')+
  geom_point(data = data1, aes(x = sample, y = Items, size = expr), shape = 21) +
  scale_size_continuous(range = c(0.5,6.5)) +
  labs(size = 'Importance')+theme(legend.position="none")
data2<- read.csv("explanation.CSV",header = TRUE)
data2$Items<-factor(data2$Items,levels=impdata$Items)
p2<-ggplot()+geom_bar(data=data2,aes(x=Explanation,y=Items),stat="identity",width=0.75, position="dodge",fill="#BEBEBE",alpha=0.7)+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(color="black")) + 
  theme_test()+theme(legend.position="none")+labs(x="",y="")+xlim(0,50)+
  theme(
    legend.text=element_text(size=12,face="plain",color="black"),
    axis.title=element_text(size=12.5,face="plain",color="black"),axis.text.x = element_text(angle=30,hjust =1,size=10,colour = 'black'),
    axis.text.y = element_blank(),axis.ticks.y=element_blank(),legend.title = element_text(size=13,face="plain",color="black"),panel.border = element_rect(size=1))
p<-p1 +p2 +  plot_layout(nrow= 1, widths = c(2.3, 1))