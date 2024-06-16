
## NI: CPO KL, DIC EAIC and EBIC
###########################################################

source("auxiliar/VerossFunctions.r")

MeasuresNI<-function(sim,y,x,numpar,n,type="N")
{

dd<-attach.bugs(sim,overwrite = T)
y=y
x=x
numpar=numpar
n=n
espac<-10

iter<-length(dd$beta)/espac
likaux<-matrix(0,n,iter)
CPOaux<-matrix(0,n,iter)

for(k in 1:iter){
i<-espac*k
betass<-dd$beta[i]
alfass<-dd$alpha[i]
nuss<-dd$nu[i]
likaux[,k]<- VerossSMSN1(y,x,alfass,betass,nuss,type="N")
CPOaux[,k]<-1/likaux[,k]
}

CPO<-1/apply(CPOaux,1,mean)
KL<--log(CPO)+apply(log(likaux),1,mean)


Dbar=mean(apply(-2*log(likaux),2,sum))
B<-sum(log(1/(apply(CPOaux,1,mean))))


DIC<-mean(apply(-2*log(likaux),2,sum))+sim$pD


EAIC<-mean(apply(-2*log(likaux),2,sum))+2*numpar
EBIC<-mean(apply(-2*log(likaux),2,sum))+ log(n)*numpar


return(list(CPOa=CPO,KLa=KL,Ba=B,Dbara=Dbar,DICa=DIC,EAICa=EAIC,EBICa=EBIC))
}



## SNI: CPO KL, DIC EAIC and EBIC
###########################################################

#source("auxiliar/VerosimilitudSNI.r")

MeasuresSNI<-function(sim,y,x,numpar,n,type="N")
{
  dd<-attach.bugs(sim,overwrite = T)
  
  y=y
  x=x
  numpar=numpar
  n=n
  espac<-10
  
  iter<-length(dd$beta)/espac
  likaux<-matrix(0,n,iter)
  CPOaux<-matrix(0,n,iter)
  
  for(k in 1:iter){
    i<-espac*k
    betass<-dd$beta[i]
    alfass<-dd$alpha[i]
    deltass<-dd$delta[i]
    nuss<-dd$nu[i]
    likaux[,k]<- VerossSMSN2(y,x,alfass,betass,deltass,nuss,type="SN")
    CPOaux[,k]<-1/likaux[,k]
  }
  
  CPO<-1/apply(CPOaux,1,mean)
  KL<--log(CPO)+apply(log(likaux),1,mean)
  
  Dbar<-mean(apply(-2*log(likaux),2,sum))
  
  B<-sum(log(1/(apply(CPOaux,1,mean))))
  DIC<-mean(apply(-2*log(likaux),2,sum))+sim$pD
  EAIC<-mean(apply(-2*log(likaux),2,sum))+2*numpar
  EBIC<-mean(apply(-2*log(likaux),2,sum))+ log(n)*numpar
  
  return(list(CPOa=CPO,KLa=KL,Ba=B,Dbara=Dbar,DICa=DIC,EAICa=EAIC,EBICa=EBIC))
}