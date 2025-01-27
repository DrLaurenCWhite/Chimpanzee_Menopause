rm(list=ls())
library(data.table)
library(rethinking)
library(ggplot2)


#Modeling female's average relatedness to all other adult males and females in her group as a function of that females age


#### Community Scale ####
Ngogo <- fread(file = "./Data/AdultFemale_MeanR_OverYears.csv")
Ngogo=copy(Ngogo[year>=2008 & year<=2015])

Ngogo$male=ifelse(Ngogo$sex=="M", 1, 0)
meanAge=mean(Ngogo$age) 
sdAge=sd(Ngogo$age)

length(unique(Ngogo$IndF))
unique(Ngogo[, c("IndF", "ImmiOrNatal")])[, .N, by=ImmiOrNatal]
Ngogo[, .N/2, by=IndF]

dat1=list(
  avg_r=Ngogo$Mean.GenR.Adult + 0.000001,
  age=Ngogo$age,
  age_s=(Ngogo$age-meanAge)/sdAge,
  male=Ngogo$male,
  ID=as.integer(as.factor(Ngogo$IndF))
)

# Beta regression parameters (using dbeta2)
# p = average probability
# theta = shape parameter

#Can load model results rather than re-run 
load("./Data/FemaleRelatedness_Average_Genetic_Adults_Beta.rdata")

#run model
# m1=map2stan(
#   alist(
#     avg_r ~ dbeta2( p, theta ),
#     logit(p) <- a + b_age*age_s + b_male*male + b_age_male*age_s*male +
#       
#       a_id[ID] + b_age_id[ID]*age_s + b_male_id[ID]*male + b_age_male_id[ID]*age_s*male,
# 
#     c(a_id, b_age_id, b_male_id, b_age_male_id)[ID] ~ dmvnormNC(sigma_s, Rho),
# 
#     a ~ dnorm ( -2, 1.5 ),
#     c(b_age, b_male, b_age_male) ~ dnorm ( 0, 1 ) ,
#     sigma_s ~ dexp(1),
#     Rho ~ dlkjcorr(2),
#     theta ~ dcauchy ( 2, 1 ) #2=above zero
#   ),
#   data=dat1,
#   iter = 2000,
#   chains = 2,
#   cores = 4,
#   control = list(adapt_delta = 0.95, max_treedepth = 15),
# )
# save(m1, file="./Data/FemaleRelatedness_Average_Genetic_Adults_Beta.rdata")

precis(m1, depth=3, pars=c("a", "b_age", "b_male", "b_age_male", "theta", "Rho"))

#CHECK PRIORS
source("./Z_ExtractPrior_LCWedit.R")
prior1=extract.prior(m1, n=1000)

dat1=as.data.table(dat1)
a_id_z = matrix(0, 1000, length(unique(dat1$ID)))
prior_pred=data.frame(age_s=seq(from=min(dat1$age_s), to=max(dat1$age_s), length.out=1000), male=0, ID=1)
prior_pred$age=meanAge+prior_pred$age_s*sdAge

prior_1=link(m1, post=prior1, data=prior_pred, 
             replace=list(a_id=a_id_z, b_age_id=a_id_z, b_male_id=a_id_z, b_age_male_id=a_id_z))

plot(NULL, xlim=c(-1.5,3), ylim=c(0,1))
for (i in 101:200){
  lines(prior_pred$age_s, prior_1[i,], col=col.alpha("black", alpha=0.5))
}

