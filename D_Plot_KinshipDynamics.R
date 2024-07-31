rm(list=ls())
library(data.table)
library(rethinking)
library(ggplot2)
library(ggpubr)


#### Community ####

##### Predicted Kinship Dynamics #####
trot=fread("./Data/PredictedKinship_Community.csv")

p1=ggplot(trot[(N*0.1)>+0.45], aes(x=N*0.1, y=rff)) + geom_line(col="darkorchid1", lwd=2) + geom_line(aes(x=N*0.1, y=rmf), lwd=2,col="orange1") +
  ylab("Average relatedness") + #xlab("Age relative to mean generation time") +
  xlab(" ") +
  #ylim(c(0,0.04)) +
  scale_x_continuous(breaks=c(0.5,1.6,2.6), labels=c("13", "40", "65"), limits=c(0.4,2.6)) +
  scale_y_continuous(breaks=c(0.0, 0.02,0.04,0.06,0.08, 0.1), limits=c(0.0,0.1)) +
    theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        panel.grid.major.x = element_line(colour="black", linetype=2),
        axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))


##### Observed Kinship Dynamics #####

Ngogo <- fread(file = "./Data/AdultFemale_MeanR_OverYears.csv")
Ngogo=copy(Ngogo[year>=2008 & year<=2015])

Ngogo$male=ifelse(Ngogo$sex=="M", 1, 0)
meanAge=mean(Ngogo$age) 
sdAge=sd(Ngogo$age)

dat1=list(
  avg_r=Ngogo$Mean.GenR.Adult,
  age=Ngogo$age,
  age_s=(Ngogo$age-meanAge)/sdAge,
  male=Ngogo$male,
  ID=as.integer(as.factor(Ngogo$IndF))
)

load(file="./Data/FemaleRelatedness_Average_Genetic_Adults_Beta.rdata")

precis(m1)
post <- extract.samples(m1)
p.females.Slope <- (post$b_age)
p.males.Slope <- (post$b_age + post$b_age_male)
precis(data.frame(p.females.Slope, p.males.Slope))

#Plot
########### Female
#zeros for varying effects to just plot main effect - the average male.female chimp
a_id_z = matrix(0, 1000, length(unique(dat1$ID)))

########### Female
dat1=as.data.table(dat1)
age_s.F <- seq(from = min(dat1[male==0]$age_s), to = max(dat1[male==0]$age_s), length.out = 1000)
Ngogo_female <- as.data.frame(age_s.F)
colnames(Ngogo_female)[1] ="age_s"
Ngogo_female$male <- 0
Ngogo_female$ID <- 1
Ngogo_female_dis <- link (m1, n=1000, data = Ngogo_female, replace=list(a_id=a_id_z, b_age_id=a_id_z, b_male_id=a_id_z, b_age_male_id=a_id_z))
Ngogo_female_mean <- apply( Ngogo_female_dis , 2 , mean )
Ngogo_female_PI <- apply( Ngogo_female_dis , 2 , PI, prob=0.89)

########### Male 
age_s.M <- seq(from = min(dat1[male==1]$age_s), to = max(dat1[male==1]$age_s), length.out = 1000)
Ngogo_male <- as.data.frame(age_s.M )
colnames(Ngogo_male)[1] ="age_s"
Ngogo_male$male <- 1
Ngogo_male$ID <- 1 
Ngogo_male_dis <- link (m1, n=1000, data = Ngogo_male, replace=list(a_id=a_id_z, b_age_id=a_id_z, b_male_id=a_id_z, b_age_male_id=a_id_z))
Ngogo_male_mean <- apply( Ngogo_male_dis , 2 , mean )
Ngogo_male_PI <- apply( Ngogo_male_dis , 2 , PI , prob=0.89 )

femaleline=data.table(age_s=age_s.F, Mean.GenR.Adult=Ngogo_female_mean, sex="F", PI_L=Ngogo_female_PI["5%",], PI_U=Ngogo_female_PI["94%",])
maleline=data.table(age_s=age_s.M, Mean.GenR.Adult=Ngogo_male_mean, sex="M", PI_L=Ngogo_male_PI["5%",], PI_U=Ngogo_male_PI["94%",])
femaleline$age=meanAge+femaleline$age_s*sdAge
maleline$age=meanAge+maleline$age_s*sdAge


#plot observed for community
p2=ggplot(Ngogo, aes(x=age/25, y=Mean.GenR.Adult, col=sex)) + geom_point(alpha=0.1) +
  #ylab("Average relatedness") + xlab("Age relative to mean generation time") +
  xlab(" ") + ylab(" ") +
  labs(col="Sex") +
  #ylim(c(0,0.04)) +
  scale_x_continuous(breaks=c(0.5,1.6,2.6), labels=c("13", "40", "65"), limits=c(0.4,2.6)) +
  scale_color_manual(values=c("darkorchid1", "orange1"), name="Sex", labels=c("Female", "Male")) +
  scale_y_continuous(breaks=c(0.0, 0.02,0.04,0.06,0.08, 0.1), limits=c(0.0,0.1)) +
    #Posterior mean lines for each sex
  geom_line(data=maleline, aes(x=age/25, y=Mean.GenR.Adult), lwd=2) +
  geom_line(data=femaleline, aes(x=age/25, y=Mean.GenR.Adult), lwd=2) +
  
  #Posterior 89% compatability interval
  geom_ribbon(data=femaleline, aes(x=age/25, ymin=PI_L, ymax=PI_U), alpha=0.2, col=NA, fill="darkorchid1") + 
  geom_ribbon(data=maleline, aes(x=age/25, ymin=PI_L, ymax=PI_U), alpha=0.2, col=NA, fill="orange1") +
  
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        panel.grid.major.x = element_line(colour="black", linetype=2),
        axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = "none",
        legend.title = element_text(size=12)) #+ facet_grid(.~sex)


#### Neighbourhood ####

##### Predicted Kinship Dynamics #####
trotN = fread("./Data/PredictedKinship_Neighbourhood.csv")

p3=ggplot(trotN[(N*0.1)>=0.5], aes(x=N*0.1, y=rff)) + geom_line(col="darkorchid1", lwd=2) + geom_line(aes(x=N*0.1, y=rmf), lwd=2,col="orange1") +
  ylab("Average relatedness") + #xlab("Age (Relative to mean generation time)") +
  xlab(" ") +
  #ylim(c(0,0.1)) +
  scale_x_continuous(breaks=c(0.5,1.6,2.6), labels=c("13", "40", "65"), limits=c(0.4,2.6)) +
  scale_y_continuous(breaks=c(0.0, 0.02,0.04,0.06,0.08, 0.1), limits=c(0.0,0.1)) +
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        panel.grid.major.x = element_line(colour="black", linetype=2),
        axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))

##### Observed Kinship Dynamics #####

NgogoN <- fread("./Data/AdultFemale_MeanR_OverYears_Neighbourhood.csv")

NgogoN$male=ifelse(NgogoN$sex=="M", 1, 0)
meanAgeN=mean(NgogoN$age) 
sdAgeN=sd(NgogoN$age)

dat1N=list(
  avg_r=NgogoN$Mean.GenR.Adult,
  age=NgogoN$age,
  age_s=(NgogoN$age-meanAgeN)/sdAgeN,
  male=NgogoN$male,
  ID=as.integer(as.factor(NgogoN$IndF))
)

load(file="./Data/FemaleRelatedness_Average_Genetic_Adults_Neighbourhoods_Beta.rdata")

precis(m3)
postN <- extract.samples(m3)
p.females.slope.N <- (postN$b_age)
p.males.slope.N <- (postN$b_age + postN$b_age_male)
precis(data.frame(p.females.slope.N, p.males.slope.N))


#Plot
########### Female
#zeros for varying effects to just plot main effect - the average male.female chimp
a_id_z_N = matrix(0, 1000, length(unique(dat1N$ID)))

########### Female
dat1N=as.data.table(dat1N)
age_s.F_N <- seq(from = min(dat1N[male==0]$age_s), to = max(dat1N[male==0]$age_s), length.out = 1000)
Ngogo_female_N <- as.data.frame(age_s.F_N)
colnames(Ngogo_female_N)[1] ="age_s"
Ngogo_female_N$male <- 0
Ngogo_female_N$ID <- 1
Ngogo_female_dis_N <- link (m3, n=1000, data = Ngogo_female_N, replace=list(a_id=a_id_z_N, b_age_id=a_id_z_N, b_male_id=a_id_z_N, b_age_male_id=a_id_z_N))
Ngogo_female_mean_N <- apply( Ngogo_female_dis_N , 2 , mean )
Ngogo_female_PI_N <- apply( Ngogo_female_dis_N , 2 , PI, prob=0.89)

########### Male 
age_s.M_N <- seq(from = min(dat1N[male==1]$age_s), to = max(dat1N[male==1]$age_s), length.out = 1000)
Ngogo_male_N <- as.data.frame(age_s.M_N )
colnames(Ngogo_male_N)[1] ="age_s"
Ngogo_male_N$male <- 1
Ngogo_male_N$ID <- 1 
Ngogo_male_dis_N <- link (m3, n=1000, data = Ngogo_male_N, replace=list(a_id=a_id_z_N, b_age_id=a_id_z_N, b_male_id=a_id_z_N, b_age_male_id=a_id_z_N))
Ngogo_male_mean_N <- apply( Ngogo_male_dis_N , 2 , mean )
Ngogo_male_PI_N <- apply( Ngogo_male_dis_N , 2 , PI , prob=0.89 )

femaleline_N=data.table(age_s=age_s.F_N, Mean.GenR.Adult=Ngogo_female_mean_N, sex="F", PI_L=Ngogo_female_PI_N["5%",], PI_U=Ngogo_female_PI_N["94%",])
maleline_N=data.table(age_s=age_s.M_N, Mean.GenR.Adult=Ngogo_male_mean_N, sex="M", PI_L=Ngogo_male_PI_N["5%",], PI_U=Ngogo_male_PI_N["94%",])
femaleline_N$age=meanAgeN+femaleline_N$age_s*sdAgeN
maleline_N$age=meanAgeN+maleline_N$age_s*sdAgeN

#plot observed for neighbourhood
p4=ggplot(NgogoN, aes(x=age/25, y=Mean.GenR.Adult, col=sex)) + geom_point(alpha=0.2) +
  #ylab("Average relatedness") + xlab("Age (Relative to mean generation time)") +
  xlab(" ") + ylab(" ") +
  labs(col="Sex") +
  #xlim(c(0.5,3)) +
  #ylim(c(0,0.1)) +
  scale_x_continuous(breaks=c(0.5,1.6,2.6), labels=c("13", "40", "65"), limits=c(0.4,2.6)) +
  scale_y_continuous(breaks=c(0.0, 0.02,0.04,0.06,0.08, 0.1), limits=c(0.0,0.1)) +
  scale_color_manual(values=c("darkorchid1", "orange1"), name="Sex", labels=c("Female", "Male")) +
  
  #Posterior mean lines for each sex
  geom_line(data=maleline_N, aes(x=age/25, y=Mean.GenR.Adult), lwd=2) +
  geom_line(data=femaleline_N, aes(x=age/25, y=Mean.GenR.Adult), lwd=2) +
  
  #Posterior 89% compatability interval
  geom_ribbon(data=femaleline_N, aes(x=age/25, ymin=PI_L, ymax=PI_U), alpha=0.2, col=NA, fill="darkorchid1") + 
  geom_ribbon(data=maleline_N, aes(x=age/25, ymin=PI_L, ymax=PI_U), alpha=0.2, col=NA, fill="orange1") +
  
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        panel.grid.major.x = element_line(colour="black", linetype=2),
        axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = "none",
        legend.title = element_text(size=12)) #+ facet_grid(.~sex)


#### Make Figure 1 ####
require(grid) 
figure=ggarrange(p1 + rremove("ylab") + rremove("xlab"), p2 + rremove("ylab") + rremove("xlab"), 
                 p3 + rremove("ylab") + rremove("xlab"), p4 + rremove("ylab") + rremove("xlab"),
                 labels="AUTO", label.x=c(0.12, 0.12, 0.12, 0.12))
#Warnings due to restricting plots to age at sexual maturity (13) and oldest age in the dataset (65)

annotate_figure(figure, left = textGrob("Local Relatedness\n ", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Age (Years)", gp = gpar(cex = 1.3)))



#### Supplementary Figure 1: Lower local mating rate
trotNS = fread("./Data/PredictedKinship_Neighbourhood_LowerLocalMating.csv")

ggplot(trotNS[(N*0.1)>=0.5], aes(x=N*0.1, y=rff)) + geom_line(col="darkorchid1", lwd=2) + geom_line(aes(x=N*0.1, y=rmf), lwd=2,col="orange1") +
  ylab("Average relatedness") + #xlab("Age (Relative to mean generation time)") +
  xlab(" ") +
  #ylim(c(0,0.1)) +
  scale_x_continuous(breaks=c(0.5,1.6,2.6), labels=c("13", "40", "65"), limits=c(0.4,2.6)) +
  scale_y_continuous(breaks=c(0.0, 0.02,0.04,0.06,0.08, 0.1), limits=c(0.0,0.1)) +
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        panel.grid.major.x = element_line(colour="black", linetype=2),
        axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))



#### Average relatedness values at age 13 vs age 65 for females: ####

##### Community #####

###### Predicted ######
######## Females #######
min(trot[(N*0.1)>=0.5]$rff)
max(trot[(N*0.1)<=2.6]$rff)
####### Males #######
min(trot[(N*0.1)>=0.5]$rmf)
max(trot[(N*0.1)<=2.6]$rmf)

###### Observed ######
####### Males #######
maleline=as.data.table(maleline)
maleline[which.min(abs(age/25-0.5))]
maleline[which.min(abs(age/25-2.6))]
######## Females #######
femaleline=as.data.table(femaleline)
femaleline[which.min(abs(age/25-0.5))]
femaleline[which.min(abs(age/25-2.6))]


##### Neighbourhood #####
###### Predicted ######
######## Females #######
min(trotN[(N*0.1)>=0.5]$rff)
max(trotN[(N*0.1)<=2.6]$rff)
####### Males #######
min(trotN[(N*0.1)>=0.5]$rmf)
max(trotN[(N*0.1)<=2.6]$rmf)

###### Observed ######
####### Males #######
maleline_N=as.data.table(maleline_N)
maleline_N[which.min(abs(age/25-0.5))]
maleline_N[which.min(abs(age/25-2.6))]
######## Females #######
femaleline_N=as.data.table(femaleline_N)
femaleline_N[which.min(abs(age/25-0.5))]
femaleline_N[which.min(abs(age/25-2.6))]



#lower local mating rate
#females
min(trotNS[(N*0.1)>=0.5]$rff)
max(trotNS[(N*0.1)<=2.6]$rff)

#males
min(trotNS[(N*0.1)>=0.5]$rmf)
max(trotNS[(N*0.1)<=2.6]$rmf)
                                      