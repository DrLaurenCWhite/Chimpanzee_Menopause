rm(list=ls())
library(data.table)
library(rethinking)
library(ggplot2)
library(ggpubr)
library(grid)

#### Community ####

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



##Cost/benefit for observed
M=1 
u=0.076
df=0.5
dm=0
nf=66
nm=44

femaleline2=data.table(age_s=age_s.F, avg_r_F=Ngogo_female_mean)
maleline2=data.table(age_s=age_s.M, avg_r_M=Ngogo_male_mean)
line=merge(femaleline2, maleline2, by="age_s")
line$age=meanAge+line$age_s*sdAge
line$Gen=line$age/25
hum=data.table()
for (i in line$Gen){
  RM=line[Gen==i]$avg_r_M
  RF=line[Gen==i]$avg_r_F
  hhat=(1/2)*((1-df)^2+(1-dm)^2)
  brkt=(1/nf)*(1+RM)+((nf-1)/nf)*(RF+RM)
  Res=((RF+RM)-hhat*brkt)/((1+RM)-hhat*brkt)
  new=data.table(Ns=i, res=Res)
  hum=rbind(hum,new)
}
hum<-hum[15:1000,]

cols = c("Help"="pink", "Harm"="blue")


p1=ggplot(hum, aes(x=Ns, y=res)) + geom_line(lwd=2) + #xlab('Age relative to mean generation time') + 
  #ylab("Individual-cost / Group-benefit") +
  xlab(" ") + ylab(" ") +
  geom_hline(yintercept=0, lwd=1, lty=2) + 
  scale_x_continuous(breaks=c(0.55,1.6,2.6), labels=c("13", "40", "65"), limits=c(0.4,2.6)) +
  scale_y_continuous(limits = c(-0.005,0.01), breaks=c(-0.004,-0.002, 0, 0.002,0.004,0.006,0.008)) +
  #geom_area(data=hum[res<0], aes(fill="Harm"), alpha=0.5) +
  geom_area(data=hum[res>0], aes(fill="Help"), alpha=0.5) +
  scale_fill_manual(values=cols) +
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        panel.grid.major.x = element_line(colour="black", linetype=2),
        axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position="none") +
  NULL

#### Subgroup ####

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


##Cost/benefit for observed
M=1
u=0.076
df=0.9
dm=0.1
nf=18
nm=12

femaleline2_N=data.table(age_s=age_s.F_N, avg_r_F=Ngogo_female_mean_N)
maleline2_N=data.table(age_s=age_s.M_N, avg_r_M=Ngogo_male_mean_N)
line_N=merge(femaleline2_N, maleline2_N, by="age_s")
line_N$age=meanAgeN+line_N$age_s*sdAgeN
line_N$Gen=line_N$age/25
hum_N=data.table()
for (i in line_N$Gen){
  RM=line_N[Gen==i]$avg_r_M
  RF=line_N[Gen==i]$avg_r_F
  hhat=(1/2)*((1-df)^2+(1-dm)^2)
  brkt=(1/nf)*(1+RM)+((nf-1)/nf)*(RF+RM)
  Res=((RF+RM)-hhat*brkt)/((1+RM)-hhat*brkt)
  new=data.table(Ns=i, res=Res)
  hum_N=rbind(hum_N,new)
}


p2=ggplot(hum_N, aes(x=Ns, y=res)) + geom_line(lwd=2) + #xlab('Age (Relative to mean generation time)') + 
  #ylab("Individual-cost / Group-benefit") +
  xlab(" ") + ylab(" ") +
  geom_hline(yintercept=0, lwd=1, lty=2) + 
  scale_x_continuous(breaks=c(0.55,1.6,2.6), labels=c("13", "40", "65"), limits=c(0.4,2.6)) +
  scale_y_continuous(limits = c(-0.005,0.01), breaks=c(-0.004,-0.002, 0, 0.002,0.004,0.006,0.008)) +
  geom_area(data=hum_N[res<0], aes(fill="Harm"), alpha=0.5) +
  geom_area(data=hum_N[res>0], aes(fill="Help"), alpha=0.5) +
  scale_fill_manual(values=cols) +
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        panel.grid.major.x = element_line(colour="black", linetype=2),
                axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position="none") +
  NULL

fig=ggarrange(p1, p2, labels="AUTO", label.x=c(0.17, 0.18))
#Warnings due to restricting plots to age at sexual maturity (13) and oldest age in the dataset (65)

annotate_figure(fig, left = textGrob("       Individual-cost / Group-benefit", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Age (Years)", gp = gpar(cex = 1.3)))

