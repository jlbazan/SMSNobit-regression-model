
## likelihood of NI
########################################################### 

source("auxiliar/CumulativeFunctions.r")
source("auxiliar/MeanVarCorrecFunctions.r")
  
VerossSMSN1<-function(y,x,alfa,beta,nu,type = "N"){
n<-length(y)
resp<-matrix(0,n,1)
aux=MedvarNI(nu,type)
eta<-alfa+beta*x
prob<-cdfNI(eta,aux$med,aux$var,nu,type)
resp <-((prob/(1-prob))^y)*(1-prob)
return(resp)
}   


## likelihood of SNI
########################################################### 

#source("auxiliar/AcumuladasdasSNI.r")
#source("auxiliar/MediasyvarianzascorregidasSNI.r")

VerossSMSN2<-function(y,x,alfa,beta,delta,nu,type = "SN"){
  n<-length(y)
  resp<-matrix(0,n,1)
  aux=MedvarSNI(delta,nu,type)
  eta<-alfa+beta*x
  prob<-cdfSNI(eta,aux$med,aux$var,delta,nu,type)
  resp <-((prob/(1-prob))^y)*(1-prob)
  return(resp)
}   