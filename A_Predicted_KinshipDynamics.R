rm(list=ls())
library(data.table)
library(ggplot2)
library(ggpubr)



#### Predicted Kinship Dynamics at the Community Scale ####
M=1 #local breeding rate
u=0.076 #turn over rate
df=0.5 #female dispersal rate
dm=0 #male dispersal rate
nf=66 #number of adult females
nm=44 #number of adult males

#To start off, we set G(ff), probability of two individuals sharing allele to 0.05
Gff=0.05
Gmm=0.05
Gfm=0.05
Gf=(1+M*Gfm)/2
Gm=Gf

#Now we iterate this a few times to get the stable values
fox=data.table(GFM=Gfm, GFF=Gff, GMM=Gmm, GS=Gm)
for (i in 1:150) {
  
  Gff<-(u^2)*((1-df)^2)*(1/4)*(( (1/nf)*Gf+((nf-1)/nf)*Gff) +2*M*Gfm + (M^2)*( (1/nm)*Gm+((nm-1)/nm)*Gmm))+
    2*u*(1-u)*(1-df)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff) + M*Gfm)+
    ((1-u)^2)*Gff
  
  Gmm<-(u^2)*((1-dm)^2)*(1/4)*((M^2)*( (1/nm)*Gm+((nm-1)/nm)*Gmm) +2*M*Gfm + (M^2)*( (1/nf)*Gf+((nf-1)/nf)*Gff))+
    2*u*(1-u)*(1-dm)*(1/2)*(M*((1/nm)*Gm +((nm-1)/nm)*Gmm) + Gfm)+
    ((1-u)^2)*Gmm
  
  Gfm<-(u^2)*(1-df)*(1-dm)*(1/4)*(( (1/nf)*Gf+((nf-1)/nf)*Gff)+2*M*Gfm + (M^2)*( (1/nm)*Gm+((nm-1)/nm)*Gmm)) +
    u*(1-u)*(1-df)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm) + Gfm) +
    u*(1-u)*(1-dm)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff) + M*Gfm) +
    (1-u)*(1-u)*Gfm
  
  Gf=(1+(M*Gfm))/2
  Gm=Gf
  
  NEW=data.table(GFM=Gfm, GFF=Gff, GMM=Gmm, GS=Gm)
  fox=rbind(NEW, fox)
  
}


#relatedness among youngest individuals
g0ff<-u*((1-df)^2)*(1/4)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+2*M*Gfm+(M^2)*((1/nm)*Gm+((nm-1)/nm)*Gmm))+
  (1-u)*(1-df)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+M*Gfm)
g0mm<-u*((1-dm)^2)*(1/4)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+2*M*Gfm+(M^2)*((1/nm)*Gm+((nm-1)/nm)*Gmm))+
  (1-u)*(1-dm)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm)+Gfm)
g0mf<-u*(1-df)*(1-dm)*(1/4)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+2*M*Gfm+(M^2)*((1/nm)*Gm+((nm-1)/nm)*Gmm))+
  (1-u)*(1-df)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm)+Gfm)
g0fm<-u*(1-df)*(1-dm)*(1/4)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+2*Gfm+(M^2)*((1/nm)*Gm+((nm-1)/nm)*Gmm))+
  (1-u)*(1-dm)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+M*Gfm)

Gff_N=g0ff
Gmm_N=g0mm
Gmf_N=g0mf
Gfm_N=g0fm
trot=data.table(N=0, gff=Gff_N, gmm=Gmm_N, gfm=Gfm_N, gmf=Gmf_N)

for (i in 1:30){
  Gff_N1<-u*(1-df)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff_N)+M*Gmf_N)+(1-u)*Gff_N
  Gmm_N1<-u*(1-dm)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm_N)+Gfm_N)+(1-u)*Gmm_N
  Gmf_N1<-u*(1-dm)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff_N)+M*Gmf_N)+(1-u)*Gmf_N
  Gfm_N1<-u*(1-df)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm_N)+Gfm_N)+(1-u)*Gfm_N
  new=data.table(N=i, gff=Gff_N1, gmm=Gmm_N1, gfm=Gfm_N1, gmf=Gmf_N1)
  trot=rbind(trot, new)
  Gff_N=Gff_N1
  Gmm_N=Gmm_N1
  Gmf_N=Gmf_N1
  Gfm_N=Gfm_N1
}

trot$rff=trot$gff/Gf
trot$rmf=trot$gmf/Gf

write.csv(trot, file="./Data/PredictedKinship_Community.csv")



#### Predicted Kinship Dynamics at the Neighbourhood Scale ####
M=1 #local breeding rate
u=0.076 #turn over rate
df=0.9 #female dispersal rate
dm=0.1 #male dispersal rate
nf=18 #number of females
nm=12 #number of males

#To start off, we set G(ff), probability of two individuals sharing allele to 0.05
Gff=0.05
Gmm=0.05
Gfm=0.05
Gf=(1+M*Gfm)/2
Gm=Gf

#Now we iterate this a few times to get the stable values
fox=data.table(GFM=Gfm, GFF=Gff, GMM=Gmm, GS=Gm)
for (i in 1:150) {
  
  Gff<-(u^2)*((1-df)^2)*(1/4)*(( (1/nf)*Gf+((nf-1)/nf)*Gff) +2*M*Gfm + (M^2)*( (1/nm)*Gm+((nm-1)/nm)*Gmm))+
    2*u*(1-u)*(1-df)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff) + M*Gfm)+
    ((1-u)^2)*Gff
  
  Gmm<-(u^2)*((1-dm)^2)*(1/4)*((M^2)*( (1/nm)*Gm+((nm-1)/nm)*Gmm) +2*M*Gfm + (M^2)*( (1/nf)*Gf+((nf-1)/nf)*Gff))+
    2*u*(1-u)*(1-dm)*(1/2)*(M*((1/nm)*Gm +((nm-1)/nm)*Gmm) + Gfm)+
    ((1-u)^2)*Gmm
  
  Gfm<-(u^2)*(1-df)*(1-dm)*(1/4)*(( (1/nf)*Gf+((nf-1)/nf)*Gff)+2*M*Gfm + (M^2)*( (1/nm)*Gm+((nm-1)/nm)*Gmm)) +
    u*(1-u)*(1-df)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm) + Gfm) +
    u*(1-u)*(1-dm)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff) + M*Gfm) +
    (1-u)*(1-u)*Gfm
  
  Gf=(1+(M*Gfm))/2
  Gm=Gf
  
  NEW=data.table(GFM=Gfm, GFF=Gff, GMM=Gmm, GS=Gm)
  fox=rbind(NEW, fox)
  
}


#relatedness among youngest individuals
g0ff<-u*((1-df)^2)*(1/4)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+2*M*Gfm+(M^2)*((1/nm)*Gm+((nm-1)/nm)*Gmm))+
  (1-u)*(1-df)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+M*Gfm)
g0mm<-u*((1-dm)^2)*(1/4)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+2*M*Gfm+(M^2)*((1/nm)*Gm+((nm-1)/nm)*Gmm))+
  (1-u)*(1-dm)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm)+Gfm)
g0mf<-u*(1-df)*(1-dm)*(1/4)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+2*M*Gfm+(M^2)*((1/nm)*Gm+((nm-1)/nm)*Gmm))+
  (1-u)*(1-df)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm)+Gfm)
g0fm<-u*(1-df)*(1-dm)*(1/4)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+2*Gfm+(M^2)*((1/nm)*Gm+((nm-1)/nm)*Gmm))+
  (1-u)*(1-dm)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+M*Gfm)

Gff_N=g0ff
Gmm_N=g0mm
Gmf_N=g0mf
Gfm_N=g0fm
trotN=data.table(N=0, gff=Gff_N, gmm=Gmm_N, gfm=Gfm_N, gmf=Gmf_N)
for (i in 1:30){
  Gff_N1<-u*(1-df)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff_N)+M*Gmf_N)+(1-u)*Gff_N
  Gmm_N1<-u*(1-dm)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm_N)+Gfm_N)+(1-u)*Gmm_N
  Gmf_N1<-u*(1-dm)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff_N)+M*Gmf_N)+(1-u)*Gmf_N
  Gfm_N1<-u*(1-df)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm_N)+Gfm_N)+(1-u)*Gfm_N
  new=data.table(N=i, gff=Gff_N1, gmm=Gmm_N1, gfm=Gfm_N1, gmf=Gmf_N1)
  trotN=rbind(trotN, new)
  Gff_N=Gff_N1
  Gmm_N=Gmm_N1
  Gmf_N=Gmf_N1
  Gfm_N=Gfm_N1
}

trotN$rff=trotN$gff/Gf
trotN$rmf=trotN$gmf/Gf

write.csv(trotN, file="./Data/PredictedKinship_Neighbourhood.csv")


#### Predicted Kinship Dynamics at the Neighbourhood Scale With Reduced Local Breeding Rate ####
M=0.5 #local breeding rate
u=0.076 #turnover rate
df=0.9 #female dispersal rate
dm=0.1 #male dispersal rate
nf=18 #number of females
nm=12 #number of males

#To start off, we set G(ff), probability of two individuals sharing allele to 0.05
Gff=0.05
Gmm=0.05
Gfm=0.05
Gf=(1+M*Gfm)/2
Gm=Gf

#Now we iterate this a few times to get the stable values
fox=data.table(GFM=Gfm, GFF=Gff, GMM=Gmm, GS=Gm)
for (i in 1:150) {
  
  Gff<-(u^2)*((1-df)^2)*(1/4)*(( (1/nf)*Gf+((nf-1)/nf)*Gff) +2*M*Gfm + (M^2)*( (1/nm)*Gm+((nm-1)/nm)*Gmm))+
    2*u*(1-u)*(1-df)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff) + M*Gfm)+
    ((1-u)^2)*Gff
  
  Gmm<-(u^2)*((1-dm)^2)*(1/4)*((M^2)*( (1/nm)*Gm+((nm-1)/nm)*Gmm) +2*M*Gfm + (M^2)*( (1/nf)*Gf+((nf-1)/nf)*Gff))+
    2*u*(1-u)*(1-dm)*(1/2)*(M*((1/nm)*Gm +((nm-1)/nm)*Gmm) + Gfm)+
    ((1-u)^2)*Gmm
  
  Gfm<-(u^2)*(1-df)*(1-dm)*(1/4)*(( (1/nf)*Gf+((nf-1)/nf)*Gff)+2*M*Gfm + (M^2)*( (1/nm)*Gm+((nm-1)/nm)*Gmm)) +
    u*(1-u)*(1-df)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm) + Gfm) +
    u*(1-u)*(1-dm)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff) + M*Gfm) +
    (1-u)*(1-u)*Gfm
  
  Gf=(1+(M*Gfm))/2
  Gm=Gf
  
  NEW=data.table(GFM=Gfm, GFF=Gff, GMM=Gmm, GS=Gm)
  fox=rbind(NEW, fox)
  
}


#relatedness among youngest individuals
g0ff<-u*((1-df)^2)*(1/4)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+2*M*Gfm+(M^2)*((1/nm)*Gm+((nm-1)/nm)*Gmm))+
  (1-u)*(1-df)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+M*Gfm)
g0mm<-u*((1-dm)^2)*(1/4)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+2*M*Gfm+(M^2)*((1/nm)*Gm+((nm-1)/nm)*Gmm))+
  (1-u)*(1-dm)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm)+Gfm)
g0mf<-u*(1-df)*(1-dm)*(1/4)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+2*M*Gfm+(M^2)*((1/nm)*Gm+((nm-1)/nm)*Gmm))+
  (1-u)*(1-df)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm)+Gfm)
g0fm<-u*(1-df)*(1-dm)*(1/4)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+2*Gfm+(M^2)*((1/nm)*Gm+((nm-1)/nm)*Gmm))+
  (1-u)*(1-dm)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff)+M*Gfm)

Gff_N=g0ff
Gmm_N=g0mm
Gmf_N=g0mf
Gfm_N=g0fm
trotN=data.table(N=0, gff=Gff_N, gmm=Gmm_N, gfm=Gfm_N, gmf=Gmf_N)
for (i in 1:30){
  Gff_N1<-u*(1-df)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff_N)+M*Gmf_N)+(1-u)*Gff_N
  Gmm_N1<-u*(1-dm)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm_N)+Gfm_N)+(1-u)*Gmm_N
  Gmf_N1<-u*(1-dm)*(1/2)*(((1/nf)*Gf+((nf-1)/nf)*Gff_N)+M*Gmf_N)+(1-u)*Gmf_N
  Gfm_N1<-u*(1-df)*(1/2)*(M*((1/nm)*Gm+((nm-1)/nm)*Gmm_N)+Gfm_N)+(1-u)*Gfm_N
  new=data.table(N=i, gff=Gff_N1, gmm=Gmm_N1, gfm=Gfm_N1, gmf=Gmf_N1)
  trotN=rbind(trotN, new)
  Gff_N=Gff_N1
  Gmm_N=Gmm_N1
  Gmf_N=Gmf_N1
  Gfm_N=Gfm_N1
}

trotN$rff=trotN$gff/Gf
trotN$rmf=trotN$gmf/Gf

max(trotN$rmf) # 0.03140293
max(trotN$rmf)/max(trotN$rff) # 8.720855


write.csv(trotN, file="./Data/PredictedKinship_Neighbourhood_LowerLocalMating.csv")

