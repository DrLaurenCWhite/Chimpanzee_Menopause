GetRel=function(df, dm, nf, nm, M, u){
  Gff=0.05
  Gmm=0.05
  Gfm=0.05
  Gf=(1+M*Gfm)/2
  Gm=Gf
  
  #Now we iterate this a few times to get the stable values
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
  
  for (i in 1:40){
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
  trot$rmm=trot$gmm/Gm
  trot$rfm=trot$gfm/Gm
  trot$Ns=trot$N*0.1
  return(trot)
}


GetCB <- function(df, dm, nf, nm, m, r){
  if (m==1) {
    fox=data.table()
    for (i in r$Ns){
      RM=r[Ns==i]$rmf
      RF=r[Ns==i]$rff
      hhat=(1/2)*((1-df)^2+(1-dm)^2)
      brkt=(1/nf)*(1+RM)+((nf-1)/nf)*(RF+RM)
      Res=((RF+RM)-hhat*brkt)/((1+RM)-hhat*brkt)
      new=data.table(Ns=i, res=Res)
      fox=rbind(fox,new)
    }
    return(fox)
  } else {
    fox=data.table()
    for (i in r$Ns){
      RM=r[Ns==i]$rmf
      RF=r[Ns==i]$rff
      hhat=(1/2)*((1-df)^2+(1-dm)^2)
      
      brkt=(1/nf)+((nf-1)/nf)*RF
      Res=((nf-1)*RF + nm*RM - (nf-1)*hhat*brkt)/((nf+nm-1)*(1-hhat*brkt))
      new=data.table(Ns=i, res=Res)
      fox=rbind(fox,new)
    }
    return(fox)
  }

}

