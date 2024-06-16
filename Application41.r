# rm
rm(list=ls())
wd = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

## Function PGEV
pgev <- function(x,eta) {
  out <- NA
  w <- rep(NA,length(x))
  for (i in 1:length(x))
    w[i] <- max(1+eta*x[i],0)
  out <- exp(-w^(-1/eta))
  return(out)
}

## Data simulation
n=500
set.seed(10000003)
x = rnorm(n)
alfasim = 0
betasim = 1
etasim = alfasim+betasim*x

# RGEV link
prob<-pgev(etasim,5)
y<-rbinom(n,1,prob)
mean(y)
plot(etasim,prob)

# GEV link
#prob<-1-pgev(-etasim,5)
#y<-rbinom(n,1,prob)
#mean(y)
#plot(etasim,prob)

## Initial values
clas=glm(y ~ x, family=binomial(link="probit"))$coefficients
alpha0=clas[1]
beta0=clas[2]
alpha0.star= alpha0-beta0*mean(x)

## For WigBugs
datas <- list("n","x","y")
library(R2WinBUGS)
library(coda)
bd="C:/WinBUGS14/"
wd1 = "outwinbugs/GEV"
wd2 = "outwinbugs/SN"
wd3 = "outwinbugs/ST"
wd4 = "outwinbugs/SSL"
wd5 = "outwinbugs/SCN"
wd6 = "outwinbugs/N"
wd7 = "outwinbugs/T"
wd8 = "outwinbugs/SL"
wd9 = "outwinbugs/CN"

## MCMC values
niter=3000
nburnin=1000
nthin=2

## MCMC values on the paper
#niter=30000
#nburnin=10000
#nthin=20

## RGEV model or GEV model
###########################################################
prog.bug<-file.path(wd, "bugsmodels/gev.txt")
parameters<-c("alpha","beta","eta")
numpar<-length(parameters)
inits <- function(){list(alpha.star=alpha0.star, beta=beta0, eta=0.5)}

timesgev<-system.time(
    resu.simgev<- bugs(datas, inits, parameters, prog.bug,
    n.chains=2, n.iter=niter, n.burnin=nburnin, n.thin=nthin, debug=F, codaPkg=F,
    bugs.directory=bd,bugs.seed=12,
    working.directory=wd1, digits = 5, clearWD=F))
tim_gev = timesgev
sum_gev = resu.simgev$summary

# Criteria GEV
t = bugs.log(paste0(wd1,"/log.txt"))
alfamean = t$stat[1,1]
betamean = t$stat[2,1]
etamean = alfamean+betamean*x
etagev = t$stat[4,1]
ppgev = 1-pgev(-etamean,etagev)
sesgopgev = sum(prob-ppgev)
sesgopgev_mean = sesgopgev/n
Dbarmgev = t$DIC[2,1]
Dicmgev = t$DIC[2,4]
likauxmedia = ((ppgev)^y)*(1-ppgev)^(1-y)
Dhatgev = -2*sum(log(likauxmedia))
DICmgevo = 2*Dbarmgev-Dhatgev
EAICmgev = Dbarmgev+2*numpar
EBICmgev = Dbarmgev+ log(n)*numpar

## auxiliar functions
source("auxiliar/CumulativeFunctions.r")
source("auxiliar/MeanVarCorrecFunctions.r")
source("auxiliar/MeasureFunctions.r")

## SN model
###########################################################
prog.bug = file.path(wd, "bugsmodels/sn.txt")
parameters<-c("alpha","beta","delta","sse")
numpar<-length(parameters)-1
inits <- function(){list(alpha.star=alpha0.star, beta=beta0, delta=0,v=rep(0.5,n))}
timessn<-system.time(
  resu.simsn<- bugs(datas, inits, parameters, prog.bug,
                    n.chains=2, n.iter=niter, n.burnin=nburnin, n.thin=nthin, debug=F,
                    codaPkg=F, bugs.directory=bd,bugs.seed=12,
                    working.directory=wd2, digits = 5, clearWD=F))
tim_sn = timessn
sum_sn = resu.simsn$summary

# Criteria SN
t = bugs.log(paste0(wd2,"/log.txt"))
alfamean=t$stat[1,1]
betamean=t$stat[2,1]
deltamean=t$stat[3,1]
etamean=alfamean+betamean*x
aux=MedvarSNI(deltamean,nu,type="SN")
psn<-cdfSNI(etamean,aux$med,aux$var,deltamean,nu,type="SN")
sesgopsn=sum(prob-psn)
sesgopsn_mean = sesgopsn/n
Dbarmsn=t$DIC[2,1]
Dicmsn=t$DIC[2,4]
likauxmedia<-((psn)^y)*(1-psn)^(1-y)
Dhatsn=-2*sum(log(likauxmedia))
DICmsno=2*Dbarmsn-Dhatsn
EAICmsn<-Dbarmsn+2*numpar
EBICmsn<-Dbarmsn+ log(n)*numpar
# aprox
measn=MeasuresSNI(resu.simsn,y,x,numpar,n,type="SN")
DICsna=2*measn$Dbara-Dhatsn

# figures
obser<-seq(1,n,1)
library(calibrate)
postscript('output/Figures/cposn.eps', width = 12.0, height = 6.0, horizontal=F)
par(mfcol=c(1,1), xaxs='r', yaxs='r', bty='l' , cex=0.8)
par(mfrow=c(1,2))
plot(log(measn$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(measn$CPOa),obser)
plot(measn$KLa,xlab="Observation",ylab="KL",type="l")
textxy(obser,measn$KLa,obser)
graphics.off()
jpeg(filename = "output/Figures/cpoklsn.jpg", width = 480, height = 480, 
     units = "px", pointsize = 12, quality = 100,
     bg = "white",  res = NA, restoreConsole = TRUE)
par(mfrow=c(1,2))
plot(log(measn$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(measn$CPOa),obser)
plot(measn$KLa,xlab="Observation",ylab="KL",,type="l")
textxy(obser,measn$KLa,obser)
dev.off()

## ST model
###########################################################
prog.bug<-file.path(wd, "bugsmodels/st.txt")
parameters<-c("alpha","beta","delta","nu","sse")
numpar<-length(parameters)-1
inits <- function(){list(alpha.star=alpha0.star, beta=beta0,
                         delta=0,nu=4,u=rep(0.5,n),v=rep(0.5,n))}
timesst<-system.time(
  resu.simst<- bugs(datas, inits, parameters, prog.bug,
                    n.chains=2, n.iter=niter, n.burnin=nburnin, n.thin=nthin, debug=F,
                    codaPkg=F, bugs.directory=bd,bugs.seed=12,
                    working.directory=wd3, digits = 5, clearWD=F))
tim_st = timesst
sum_st = resu.simst$summary

# Criteria ST
t=bugs.log(paste0(wd3,"/log.txt"))
alfamean=t$stat[1,1]
betamean=t$stat[2,1]
deltamean=t$stat[3,1]
numean=t$stat[5,1]
etamean=alfamean+betamean*x
aux=MedvarSNI(deltamean,numean,type="ST")
pst<-cdfSNI(etamean,aux$med,aux$var,deltamean,numean,type="ST")
sesgopst=sum(prob-pst)
Dbarmst=t$DIC[2,1]
Dicmst=t$DIC[2,4]
likauxmedia<-((pst)^y)*(1-pst)^(1-y)
Dhatst=-2*sum(log(likauxmedia))
DICmsto=2*Dbarmst-Dhatst
EAICmst<-Dbarmst+2*numpar
EBICmst<-Dbarmst+ log(n)*numpar
# aprox
meast=MeasuresSNI(resu.simst,y,x,numpar,n,type="ST")
DICsta=2*meast$Dbara-Dhatst
library(calibrate)
postscript('output/Figures/cpost.eps', width = 12.0, height = 6.0, horizontal=F)
par(mfcol=c(1,1), xaxs='r', yaxs='r', bty='l' , cex=0.8)
par(mfrow=c(1,2))
plot(log(meast$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(meast$CPOa),obser)
plot(meast$KLa,xlab="Observation",ylab="KL",,type="l")
textxy(obser,meast$KLa,obser)
graphics.off()
jpeg(filename = "output/Figures/cpoklst.jpg", width = 480, height = 480, 
     units = "px", pointsize = 12, quality = 100,
     bg = "white",  res = NA, restoreConsole = TRUE)
par(mfrow=c(1,2))
plot(log(meast$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(meast$CPOa),obser)
plot(meast$KLa,xlab="Observation",ylab="KL",type="l")
textxy(obser,meast$KLa,obser)
dev.off()


## SSL model
###########################################################
prog.bug<-file.path(wd, "bugsmodels/ssl.txt")
parameters<-c("alpha","beta","delta","nu","sse")
numpar<-length(parameters)-1
inits <- function(){list(alpha.star=alpha0.star, beta=beta0,
                         delta=0,nu=4,v=rep(0.5,n),u=rep(0.5,n))}
timesssl<-system.time(
  resu.simssl<- bugs(datas, inits, parameters, prog.bug,
                     n.chains=2, n.iter=niter, n.burnin=nburnin, n.thin=nthin,
                     debug=F, codaPkg=F, bugs.directory=bd,bugs.seed=12,
                     working.directory=wd4, digits = 5, clearWD=F))
tim_ssl = timesssl
sum_ssl = resu.simssl$summary

# Criteria SSL
t=bugs.log(paste0(wd4,"/log.txt"))
alfamean=t$stat[1,1]
betamean=t$stat[2,1]
deltamean=t$stat[3,1]
numean=t$stat[5,1]
etamean=alfamean+betamean*x
aux=MedvarSNI(deltamean,numean,type="SSL")
system.time(
  pssl<-cdfSNI(etamean,aux$med,aux$var,deltamean,numean,type="SSL")
)
sesgopssl=sum(prob-pssl)
Dbarmssl=t$DIC[2,1]
Dicmssl=t$DIC[2,4]
likauxmedia<-((pssl)^y)*(1-pssl)^(1-y)
Dhatssl=-2*sum(log(likauxmedia))
DICmsslo=2*Dbarmssl-Dhatssl
EAICmssl<-Dbarmssl+2*numpar
EBICmssl<-Dbarmssl+ log(n)*numpar
# aprox
meassl=MeasuresSNI(resu.simssl,y,x,numpar,n,type="SSL")
DICssla=2*meassl$Dbara-Dhatssl
library(calibrate)
postscript('output/Figures/cpossl.eps', width = 12.0, height = 6.0, horizontal=F)
par(mfcol=c(1,1), xaxs='r', yaxs='r', bty='l' , cex=0.8)
par(mfrow=c(1,2))
plot(log(meassl$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(meassl$CPOa),obser)
plot(meassl$KLa,xlab="Observation",ylab="KL",,type="l")
textxy(obser,meassl$KLa,obser)
graphics.off()
jpeg(filename = "output/Figures/cpoklssl.jpg", width = 480, height = 480, 
     units = "px", pointsize = 12, quality = 100,
     bg = "white",  res = NA, restoreConsole = TRUE)
par(mfrow=c(1,2))
plot(log(meassl$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(meassl$CPOa),obser)
plot(meassl$KLa,xlab="Observation",ylab="KL",type="l")
textxy(obser,meassl$KLa,obser)
dev.off()


## SCN model
###########################################################
prog.bug<-file.path(wd, "bugsmodels/scn.txt")
parameters<-c("alpha","beta","delta","nu","gamma","sse")
numpar<-length(parameters)-1
inits <- function(){list(alpha.star=alpha0.star, beta=beta0,
                         delta=0,nu=0.5,gamma=0.5,v=rep(0.5,n))}
timesscn<-system.time(
  resu.simscn<- bugs(datas, inits, parameters, prog.bug,
                     n.chains=2, n.iter=niter, n.burnin=nburnin, n.thin=nthin, debug=F, codaPkg=F,
                     bugs.directory=bd,bugs.seed=12,
                     working.directory=wd5, digits = 5, clearWD=F))
tim_scn = timesscn
sum_scn = resu.simscn$summary

# Criteria SCN
t=bugs.log(paste0(wd5,"/log.txt"))
alfamean=t$stat[1,1]
betamean=t$stat[2,1]
deltamean=t$stat[3,1]
gammamean=t$stat[5,1]
numean=t$stat[6,1]
nuvmean=c(numean,gammamean)
etamean=alfamean+betamean*x
aux=MedvarSNI(deltamean,nuvmean,type="SCN")
system.time(
  pscn<-cdfSNI(etamean,aux$med,aux$var,deltamean,nuvmean,type="SCN")
)
sesgopscn=sum(prob-pscn)
likauxmedia<-((pscn)^y)*(1-pscn)^(1-y)
Dhatscn=-2*sum(log(likauxmedia))
Dbarmscn=NA
Dicmscn=NA
DICmscno=NA
EAICmscn=NA
EBICmscn=NA
# aprox
meascn=MeasuresSNI(resu.simscn,y,x,numpar,n,type="SCN")
DICscna=2*meascn$Dbara-Dhatscn
library(calibrate)
postscript('output/Figures/cposcn.eps', width = 12.0, height = 6.0, horizontal=F)
par(mfcol=c(1,1), xaxs='r', yaxs='r', bty='l' , cex=0.8)
par(mfrow=c(1,2))
plot(log(meascn$CPO),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(meascn$CPO),obser)
plot(meascn$KL,xlab="Observation",ylab="KL",,type="l")
textxy(obser,meascn$KL,obser)
graphics.off()

jpeg(filename = "output/Figures/cpoklscn.jpg", width = 480, height = 480, 
     units = "px", pointsize = 12, quality = 100,
     bg = "white",  res = NA, restoreConsole = TRUE)
par(mfrow=c(1,2))
plot(log(meascn$CPO),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(meascn$CPO),obser)
plot(meascn$KL,xlab="Observation",ylab="KL",type="l")
textxy(obser,meascn$KL,obser)
dev.off()


## N model
###########################################################
prog.bug<-file.path(wd, "bugsmodels/n.txt")
parameters<-c("alpha","beta","sse")
numpar<-length(parameters)-1
inits <- function(){list(alpha.star=alpha0.star, beta=beta0)}
timesn<-system.time(
  resu.simn<- bugs(datas, inits, parameters, prog.bug,
                   n.chains=2, n.iter=niter, n.burnin=nburnin, n.thin=nthin, debug=F, codaPkg=F,
                   bugs.directory=bd,bugs.seed=12,
                   working.directory=wd6, digits = 5, clearWD=F))
tim_n = timesn
sum_n = resu.simn$summary

# Criteria N
t=bugs.log(paste0(wd6,"/log.txt"))
alfamean=t$stat[1,1]
betamean=t$stat[2,1]
etamean=alfamean+betamean*x
numean=1
aux=MedvarNI(numean,type="N")
system.time(
  pn<-cdfNI(etamean,aux$med,aux$var,numean,type="N")
)
sesgopn=sum(prob-pn)
Dbarmn=t$DIC[2,1]
Dicmn=t$DIC[2,4]
likauxmedia<-((pn)^y)*(1-pn)^(1-y)
Dhatn=-2*sum(log(likauxmedia))
DICmno=2*Dbarmn-Dhatn
EAICmn<-Dbarmn+2*numpar
EBICmn<-Dbarmn+ log(n)*numpar
#aprox
mean=MeasuresNI(resu.simn,y,x,numpar,n,type="N")
DICta=2*mean$Dbara-Dhatn

library(calibrate)
postscript('output/Figures/cpon.eps', width = 12.0, height = 6.0, horizontal=F)
par(mfcol=c(1,1), xaxs='r', yaxs='r', bty='l' , cex=0.8)
par(mfrow=c(1,2))
plot(log(mean$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(mean$CPOa),obser)
plot(mean$KLa,xlab="Observation",ylab="KL",,type="l")
textxy(obser,mean$KLa,obser)
graphics.off()
jpeg(filename = "output/Figures/cpokln.jpg", width = 480, height = 480, 
     units = "px", pointsize = 12, quality = 100,
     bg = "white",  res = NA, restoreConsole = TRUE)
par(mfrow=c(1,2))
plot(log(mean$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(mean$CPOa),obser)
plot(mean$KLa,xlab="Observation",ylab="KL",type="l")
textxy(obser,mean$KLa,obser)
dev.off()

## T model
###########################################################
prog.bug<-file.path(wd, "bugsmodels/t.txt")
parameters<-c("alpha","beta","nu","sse")
numpar<-length(parameters)-1
inits <- function(){list(alpha.star=alpha0.star, beta=beta0, nu=4,u=rep(0.5,n))}
timest<-system.time(
  resu.simt<- bugs(datas, inits, parameters, prog.bug,
                   n.chains=2, n.iter=niter, n.burnin=nburnin, n.thin=nthin, debug=F, codaPkg=F,
                   bugs.directory=bd,bugs.seed=12,
                   working.directory=wd7, digits = 5, clearWD=F))
tim_t = timest
sum_t = resu.simt$summary

# Criteria T
t=bugs.log(paste0(wd7,"/log.txt"))
alfamean=t$stat[1,1]
betamean=t$stat[2,1]
etamean=alfamean+betamean*x
numean=t$stat[4,1]
aux=MedvarNI(numean,type="N")
system.time(
  pt<-cdfNI(etamean,aux$med,aux$var,numean,type="T")
)
sesgopt=sum(prob-pt)
Dbarmt=t$DIC[2,1]
Dicmt=t$DIC[2,4]
likauxmedia<-((pt)^y)*(1-pt)^(1-y)
Dhatt=-2*sum(log(likauxmedia))
DICmto=2*Dbarmt-Dhatt
EAICmt<-Dbarmt+2*numpar
EBICmt<-Dbarmt+ log(n)*numpar
#aprox
meat=MeasuresNI(resu.simt,y,x,numpar,n,type="T")
DICta=2*meat$Dbara-Dhatt

library(calibrate)
postscript('output/Figures/cpot.eps', width = 12.0, height = 6.0, horizontal=F)
par(mfcol=c(1,1), xaxs='r', yaxs='r', bty='l' , cex=0.8)
par(mfrow=c(1,2))
plot(log(meat$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(meat$CPOa),obser)
plot(meat$KLa,xlab="Observation",ylab="KL",,type="l")
textxy(obser,meat$KLa,obser)
graphics.off()
jpeg(filename = "output/Figures/cpoklt.jpg", width = 480, height = 480, 
     units = "px", pointsize = 12, quality = 100,
     bg = "white",  res = NA, restoreConsole = TRUE)
par(mfrow=c(1,2))
plot(log(meat$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(meat$CPOa),obser)
plot(meat$KLa,xlab="Observation",ylab="KL",,type="l")
textxy(obser,meat$KLa,obser)
dev.off()

## SL model
###########################################################
prog.bug<-file.path(wd, "bugsmodels/sl.txt")
parameters<-c("alpha","beta","nu","sse")
numpar<-length(parameters)-1
inits <- function(){list(alpha.star=alpha0.star, beta=beta0,nu=4,u=rep(0.5,n))}
timessl<-system.time(
  resu.simsl<- bugs(datas, inits, parameters, prog.bug,
                    n.chains=2, n.iter=niter, n.burnin=nburnin, n.thin=nthin,
                    debug=F, codaPkg=F, bugs.directory=bd,bugs.seed=12,
                    working.directory=wd8, digits = 5, clearWD=F))
tim_sl = timessl
sum_sl = resu.simsl$summary

# Criteria SL
t=bugs.log(paste0(wd8,"/log.txt"))
alfamean=t$stat[1,1]
betamean=t$stat[2,1]
etamean=alfamean+betamean*x
numean=t$stat[4,1]
aux=MedvarNI(numean,type="SL")
system.time(
  psl<-cdfNI(etamean,aux$med,aux$var,numean,type="SL")
)
sesgopsl=sum(prob-psl)
Dbarmsl=t$DIC[2,1]
Dicmsl=t$DIC[2,4]
likauxmedia<-((psl)^y)*(1-psl)^(1-y)
Dhatsl=-2*sum(log(likauxmedia))
DICmslo=2*Dbarmsl-Dhatsl
EAICmsl<-Dbarmsl+2*numpar
EBICmsl<-Dbarmsl+ log(n)*numpar
#aprox
measl=MeasuresNI(resu.simsl,y,x,numpar,n,type="SL")
DICsla=2*measl$Dbara-Dhatsl

library(calibrate)
postscript('output/Figures/cposl.eps', width = 12.0, height = 6.0, horizontal=F)
par(mfcol=c(1,1), xaxs='r', yaxs='r', bty='l' , cex=0.8)
par(mfrow=c(1,2))
plot(log(measl$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(measl$CPOa),obser)
plot(measl$KLa,xlab="Observation",ylab="KL",type="l")
textxy(obser,measl$KLa,obser)
graphics.off()
jpeg(filename = "output/Figures/cpoklsl.jpg", width = 480, height = 480, 
     units = "px", pointsize = 12, quality = 100,
     bg = "white",  res = NA, restoreConsole = TRUE)
par(mfrow=c(1,2))
plot(log(measl$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(measl$CPOa),obser)
plot(measl$KLa,xlab="Observation",ylab="KL",type="l")
textxy(obser,measl$KLa,obser)
dev.off()

# CN model
###########################################################
prog.bug<-file.path(wd, "bugsmodels/cn.txt")
parameters<-c("alpha","beta","nu","gamma","sse")
numpar<-length(parameters)-1
inits <- function(){list(alpha.star=alpha0.star, beta=beta0,nu=0.5,gamma=0.5)}
timescn<-system.time(
  resu.simcn<- bugs(datas, inits, parameters, prog.bug,
                    n.chains=2, n.iter=niter, n.burnin=nburnin, n.thin=nthin, debug=F, codaPkg=F,
                    bugs.directory=bd,bugs.seed=12,
                    working.directory=wd9, digits = 5, clearWD=F))
tim_cn = timescn
sum_cn = resu.simcn$summary

# Criteria CN
t=bugs.log(paste0(wd9,"/log.txt"))
alfamean=t$stat[1,1]
betamean=t$stat[2,1]
etamean=alfamean+betamean*x
numean=t$stat[5,1]
gammamean=t$stat[4,1]
nuvmean=c(numean,gammamean)
aux=MedvarNI(nuvmean,type="CN")
system.time(
  pcn<-cdfNI(etamean,aux$med,aux$var,nuvmean,type="CN")
)
sesgopcn=sum(prob-pcn)
likauxmedia<-((pcn)^y)*(1-pcn)^(1-y)
Dhatcn=-2*sum(log(likauxmedia))
Dbarmcn=NA
Dicmcn=NA
DICmcno=NA
EAICmcn=NA
EBICmcn=NA
#aprox
meacn=MeasuresNI(resu.simscn,y,x,numpar,n,type="CN")
DICcna=2*meacn$Dbara-Dhatcn

library(calibrate)
postscript('output/Figures/cpocn.eps', width = 12.0, height = 6.0, horizontal=F)
par(mfcol=c(1,1), xaxs='r', yaxs='r', bty='l' , cex=0.8)
par(mfrow=c(1,2))
plot(log(meacn$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(meacn$CPOa),obser)
plot(meacn$KLa,xlab="Observation",ylab="KL",type="l")
textxy(obser,meacn$KLa,obser)
graphics.off()
jpeg(filename = "output/Figures/cpoklcn.jpg", width = 480, height = 480, 
     units = "px", pointsize = 12, quality = 100,
     bg = "white",  res = NA, restoreConsole = TRUE)
par(mfrow=c(1,2))
plot(log(meacn$CPOa),xlab="Observation",ylab="log(CPO)",type="l")
textxy(obser,log(meacn$CPOa),obser)
plot(meacn$KLa,xlab="Observation",ylab="KL",type="l")
textxy(obser,meacn$KLa,obser)
dev.off()


## Complet results
sesgofin=c(sesgopgev,sesgopsn,sesgopst,sesgopssl,sesgopscn,sesgopn,sesgopt,sesgopsl,sesgopcn)
mediasesgo=sesgofin/n
Dbarap=c(Dbarmgev,measn$Dbara,meast$Dbara,meassl$Dbara,meascn$Dbara,mean$Dbara,meat$Dbara,measl$Dbara,meacn$Dbara)
DICap=c(Dicmgev,measn$DICa,meast$DICa,meassl$DICa,meascn$DICa,mean$DICa,meat$DICa,measl$DICa,meacn$DICa)
EAICap=c(EAICmgev,measn$EAICa,meast$EAICa,meassl$EAICa,meascn$EAICa,mean$EAICa,meat$EAICa,measl$EAICa,meacn$EAICa)
EBICap=c(EBICmgev,measn$EBICa,meast$EBICa,meassl$EBICa,meascn$EBICa,mean$EBICa,meat$EBICa,measl$EBICa,meacn$EBICa)
Bap=c(0,measn$Ba,meast$B,meassl$B,meascn$B,mean$B,meat$B,measl$B,meacn$B)
Dbarmfin=c(Dbarmgev,Dbarmsn,Dbarmst,Dbarmssl,Dbarmscn,Dbarmn,Dbarmt,Dbarmsl,Dbarmcn)
Dhatfin=c(Dhatgev,Dhatsn,Dhatst,Dhatssl,Dhatscn,Dhatn,Dhatt,Dhatsl,Dhatcn)
DICmfin=c(Dicmgev,Dicmsn,Dicmst,Dicmssl,Dicmscn,Dicmn,Dicmt,Dicmsl,Dicmcn)
DICmfino=c(DICmgevo,DICmsno,DICmsto,DICmsslo,DICmscno,DICmno,DICmto,DICmslo,DICmcno)
EAICmfin=c(EAICmgev,EAICmsn,EAICmst,EAICmssl,EAICmscn,EAICmn,EAICmt,EAICmsl,EAICmcn)
EBICmfin=c(EBICmgev,EBICmsn,EBICmst,EBICmssl,EBICmscn,EBICmn,EBICmt,EBICmsl,EBICmcn)
DICapfin=2*Dbarap-Dhatfin

resultsb=rbind(sesgofin,mediasesgo,Dbarap,DICap,EAICap,EBICap,Bap,Dbarmfin,Dhatfin,DICmfin,DICmfino,EAICmfin,EBICmfin,DICapfin)
colnames(resultsb) <- c("GEV","SN","ST","SSL","SCN","N","T","SL","CN")

## Reducid results
resultsf=rbind(sesgofin,mediasesgo,Dbarap,EAICap,EBICap,Bap,Dhatfin,DICapfin)
colnames(resultsf) <- c("GEV","SN","ST","SSL","SCN","N","T","SL","CN")

## SAVE RDA
save.image(file="ProgSimulation.rda")






