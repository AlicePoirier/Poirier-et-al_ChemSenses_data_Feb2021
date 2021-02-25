

#### R code for Poirier et al., Chemical Senses, 2021
#### Scent marks signal species, sex and reproductive status in tamarins (Saguinus spp., Neotropical primates)



#### 1)  BEHAVIOUR   ----

### Summary statistics  ----
library(dplyr)

# Dataset: all recorded scent-marking events: Scent-Marking_behaviour_dataset_all_Feb2021.csv
SMdata<- read.csv(file.choose(), header=T, na.strings = "NA", sep=",")
SMdata[,-10]<- lapply(SMdata[,-10], factor)
head(SMdata)


# Summary individual
range(SMdata$Value)
mean(SMdata$Value)
sd(SMdata$Value)

coplot(Value~Gland | Sex + Repro, data=SMdata)


# Summary species, sex and repro
SMdata.species<- SMdata %>%
  group_by(Day, Bout, Species, Sex, Repro) %>%
  summarise(Value=sum(Value))
head(SMdata.species)

SMdata.species %>% group_by(Species) %>% 
  summarise(n=sum(Value), mean=mean(Value), med=median(Value), sd=sd(Value), IQR=IQR(Value),
            se=round(sd(Value)/sqrt(sum(n())), digits=3), min=min(Value), max=max(Value))

SMdata.species %>% group_by(Species, Sex) %>% 
  summarise(n=sum(Value), mean=mean(Value), med=median(Value), sd=sd(Value), IQR=IQR(Value),
            se=round(sd(Value)/sqrt(sum(n())), digits=3), min=min(Value), max=max(Value))

SMdata.species %>% group_by(Species, Repro) %>% 
  summarise(n=sum(Value), mean=mean(Value), med=median(Value), sd=sd(Value), IQR=IQR(Value),
            se=round(sd(Value)/sqrt(sum(n())), digits=3), min=min(Value), max=max(Value))


# Summary group
SMdata.group<- SMdata %>%
  group_by(Day, Bout, Group) %>%
  summarise(Value=sum(Value))
head(SMdata.group)

SMdata.group %>% group_by(Group) %>% 
  summarise(n=sum(Value), mean=mean(Value), med=median(Value), sd=sd(Value), IQR=IQR(Value),
            se=round(sd(Value)/sqrt(sum(n())), digits=3), min=min(Value), max=max(Value))

 
# Summary Gland
SMdata.gland<- SMdata %>%
  group_by(Day, Bout, Species, Sex, Repro, Gland) %>%
  summarise(Value=sum(Value))
head(SMdata.gland)

SMdata.gland %>% group_by(Species, Gland) %>% 
  summarise(n=sum(Value), mean=mean(Value), med=median(Value), sd=sd(Value), IQR=IQR(Value),
            se=round(sd(Value)/sqrt(sum(n())), digits=3), min=min(Value), max=max(Value))



### Behaviour linear model (glmmTMB) ----
library(glmmTMB)

# Dataset: sternal gland (n = 12) removed for GLMM: Scent-Marking_behaviour_dataset_noSTG_Feb2021.csv
SMdata2<- read.csv(file.choose(), header=T, na.strings = "NA", sep=",")
SMdata2[,-10]<- lapply(SMdata2[,-10], factor)
head(SMdata2)


# Null model
m0.null<- glmmTMB(Value~1+(1|Group:Indiv)+(1|Group:Day)+(1|Group:Day:Bout), data=SMdata2, ziformula=~1, family=poisson(link="log"))
AIC(m0.null) #2427
summary(m0.null)
# DHARMa OK!


# Species only model ----
m1<- glmmTMB(Value~Species+(1|Group:Indiv)+(1|Group:Day)+(1|Group:Day:Bout), data=SMdata2, ziformula=~1, family=poisson(link="log"))
AIC(m1) #2429
summary(m1)
# DHARMa OK!


# Full model following reviewers comments ----
m2.full<- glmmTMB(Value~Species+Sex+Repro+Gland+(1|Group:Indiv)+(1|Group:Day)+(1|Group:Day:Bout), data=SMdata2, ziformula=~1, family=poisson(link="log"))
AIC(m2.full) #2193
summary(m2.full)
# DHARMa OK!


# Check residuals (DHARMa)
library(DHARMa)

mymodel=m2.full
simulationOutput<- simulateResiduals(fittedModel=mymodel, n=1000, refit=F, use.u=F)

plot(simulationOutput)
testResiduals(simulationOutput)



####  2)  SEMIOCHEMISTRY  ----

### Summary statistics ----

# All anogenital gland samples: Scent-Marking_semiochem_all_Feb2021.csv
Comp.all<-read.csv(file.choose(), header=T, sep=",")
Comp.all[,-c(10,11)]<- lapply(Comp.all[,-c(10,11)], factor)
head(Comp.all)

Comp.all %>% 
  group_by(Sample) %>% tally() %>%
  summarise(min=min(n), max=max(n), mean=mean(n), sd=sd(n))


# unique compounds and compounds found in half the samples
library(reshape2)

Comp.list<- melt(xtabs(~Comp, data=Comp.all, addNA=T))
head(Comp.list)

Comp.unique<- length(which(Comp.list$value==1))
Comp.unique
Comp.unique * 100 / nrow(Comp.list)

Comp.half<- length(which(Comp.list$value>length(unique(Comp.all$Sample))/2))
Comp.half
Comp.half * 100 / nrow(Comp.list)



### Calculate log(x+1) relative peak area ----
# All samples
Comp.all.transfo<- data.frame(Comp.all %>% group_by(Sample) %>%
          summarise(Area=Area, SArea=sum(Area), logRPA=round(log(((Area/sum(Area))*100)+1), digits=3)) %>%
          add_tally(name="Ncomp") )
head(Comp.all.transfo)

Comp.all.transfo<- Comp.all.transfo %>% arrange(Sample, Area)
Comp.all<- Comp.all %>% arrange(Sample, Area)

Comp.all2<- cbind(Comp.all, Comp.all.transfo[,-(1:2)])
head(Comp.all2)



### Bray-Curtis dissimilarity index and Non-metric Multidimensional Scaling  ----
library(vegan)

head(Comp.all2)
All.cla<-unique(Comp.all2[,c(1:8)])
All.cla2<-arrange(All.cla, Sample)
All.mat<- xtabs(logRPA~Sample+Comp, Comp.all2)
All.vegdist<-vegdist(All.mat, method="bray")
All.nmds<- metaMDS(All.vegdist, k=2, autotransform=F, noshare=F, wascores=F, expand=T, trymax=200, trace=F)
All.nmds.sco<- scores(All.nmds)
All.dissim.dataset<- as.data.frame(cbind(All.cla2, All.nmds.sco))
head(All.dissim.dataset)
All.nmds$stress



### Figure 1: NMDS ----
library(ggplot2)

fig1a<- 
  ggplot(All.dissim.dataset, mapping=aes(x=NMDS1, y=NMDS2, shape=Species, linetype=Species))+
  geom_point(size=2.2, colour="black")+
  scale_shape_manual(values=c(2,17))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2), level=0.75, colour="black")+ 
  scale_linetype_manual(values=c("dashed","solid"))+
  xlab("NMDS Axis 1") + ylab("NMDS Axis 2")+ 
  labs(title=expression(paste("a) Species")))+
  guides(colour=F, shape=F, linetype=F)+
  scale_x_continuous(breaks=seq(-0.6, 0.6, 0.3))+
  scale_y_continuous(breaks=seq(-0.6, 0.6, 0.3))+
  theme_classic(base_size=11)+ 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))

fig1b<- 
  ggplot(All.dissim.dataset, mapping=aes(x=NMDS1, y=NMDS2, shape=Sex, linetype=Sex))+
  geom_point(size=2.2, colour="black")+
  scale_shape_manual(values=c(16,1))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2), level=0.75, colour="black")+ 
  scale_linetype_manual(values=c("solid","dashed"))+
  xlab("NMDS Axis 1") + ylab("NMDS Axis 2")+ 
  labs(title=expression(paste("b) Sex")))+
  guides(colour=F, shape=F, linetype=F)+
    scale_x_continuous(breaks=seq(-0.6, 0.6, 0.3))+
    scale_y_continuous(breaks=seq(-0.6, 0.6, 0.3))+
  theme_classic(base_size=11)+ 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))

fig1c<- 
  ggplot(All.dissim.dataset, mapping=aes(x=NMDS1, y=NMDS2, shape=Repro, linetype=Repro))+
  geom_point(size=2.2, colour="black")+
  scale_shape_manual(values=c(0,15))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2), level=0.75, colour="black")+ 
  scale_linetype_manual(values=c("dashed","solid"))+
  xlab("NMDS Axis 1") + ylab("NMDS Axis 2")+ 
  labs(title=expression(paste("c) Reproductive status")))+
  guides(colour=F, shape=F, linetype=F)+
  scale_x_continuous(breaks=seq(-0.6, 0.6, 0.3))+
  scale_y_continuous(breaks=seq(-0.6, 0.6, 0.3))+
  theme_classic(base_size=11)+ 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))

#dev.new()
library(cowplot)

fig1<- plot_grid(fig1a, fig1b, fig1c, rel_heights=c(1,1,1), ncol=3, align="h")

ggsave(plot=fig1, filename="Fig1_NMDS_Feb2021.tiff", device="tiff", scale=1.2, dpi=600, units="cm", width=18.5, height=6)




### PerMANOVA  ----

# Permutests
permutest(betadisper(All.vegdist, All.dissim.dataset$Species)) #ns
permutest(betadisper(All.vegdist, All.dissim.dataset$Sex)) #ns
permutest(betadisper(All.vegdist, All.dissim.dataset$Repro)) #ns


# colinearity tests
adonis2(All.vegdist~Species+Sex, data=All.dissim.dataset, permutations=999, method="bray", arg="margins")
adonis2(All.vegdist~Sex+Species, data=All.dissim.dataset, permutations=999, method="bray", arg="margins")

adonis2(All.vegdist~Species+Repro, data=All.dissim.dataset, permutations=999, method="bray", arg="margins")
adonis2(All.vegdist~Repro+Species, data=All.dissim.dataset, permutations=999, method="bray", arg="margins")

adonis2(All.vegdist~Repro+Sex, data=All.dissim.dataset, permutations=999, method="bray", arg="margins")
adonis2(All.vegdist~Sex+Repro, data=All.dissim.dataset, permutations=999, method="bray", arg="margins")
# no colinearity


# Global model with interactions
adonis2(All.vegdist~Species*Sex*Repro, data=All.dissim.dataset, permutations=999, method="bray", arg="margins") # interactions ns


# Global model without interaction (final)
adonis2(All.vegdist~Species+Sex+Repro, data=All.dissim.dataset, permutations=999, method="bray", arg="margins") 



### SIMPER ----
summary(simper(All.mat, All.dissim.dataset$Species), ordered=T)
summary(simper(All.mat, All.dissim.dataset$Sex), ordered=T)
summary(simper(All.mat, All.dissim.dataset$Repro), ordered=T)
# manually save these results as dataset


# Calculate mean contribution score
# Open SIMPER results: Scent-Marking_semiochem_SIMPER_allComp_Feb2021.csv
simper.results<-read.csv(file.choose(), header=T, sep=",", na.strings="BL")
head(simper.results)

simper.results.sum<- data.frame(simper.results %>% group_by(comp) %>%
                     summarise(Mean.ave=mean(average), SD.ave=sd(average), Mean.rank=mean(rank)))
simper.results.sum2<- simper.results.sum %>% arrange(desc(Mean.ave))
simper.results.sum2<- mutate(simper.results.sum2, General.rank=1:n())
head(simper.results.sum2)

simper.results.interest<- simper.results.sum2[which(simper.results.sum2$Mean.ave>0.0044),]
# Save compound contribution results and manually add compound names to dataset
write.table(simper.results.interest,"Scent-Marking_semiochem_SIMPER_Feb2021.csv", sep=";")



# Figure 2: 41 comp of interest ----
# Open compound contribution dataset: Scent-Marking_semiochem_SIMPER_Feb2021.csv
Comp.interest<-read.csv(file.choose(), header=T, sep=";", na.strings="BL")
Comp.interest[,c(1,3,4)]<- lapply(Comp.interest[,c(1,3,4)], factor)
head(Comp.interest)


#dev.new()
fig2<-
ggplot(Comp.interest, aes(x=name, y=average))+
  geom_point(aes(shape=factor), size=2)+
  scale_shape_manual(name="Category", values=c(2,1,0), breaks=c("Species","Sex","Repro. status"))+
  scale_x_discrete(limits=rev(levels(Comp.interest$name)))+
  xlab("")+
  ylab("Relative contribution to chemical dissimilarity")+
  coord_flip()+
  theme_bw(base_size=10)+
  theme(legend.position="top")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))


ggsave(plot=fig2, filename="Fig2_contrib_600dpi_Feb2021.tiff", device="tiff", scale=1, dpi=600, units="cm", width=18.3, height=15.2)


