
## NI: Calculating corrected means and variances
###########################################################

MedvarNI<-function(nu,type="N"){
if(type=="N"){
k1<-1
k2<-1
}

if(type=="T"){
k1<-(nu/2)^(1/2)*gamma((nu-1)/2)/gamma(nu/2)
k2<-(nu/2)*gamma((nu-2)/2)/gamma(nu/2)
}

if(type=="SL"){
k1<-2*nu/(2*nu-1)
k2<-2*nu/(2*nu-2)
}

if(type=="CN"){
k1<-nu[1]/nu[2]^(1/2)+1-nu[1]
k2<-nu[1]/nu[2]+1-nu[1]
}

Auxmed<-0
Auxvar<-k2
med<- Auxmed
var<- 1/Auxvar
return(list(med=med,var=var))
}


## SNI: Calculating corrected means and variances
###########################################################

MedvarSNI<-function(delta,nu,type="SN"){
  if(type=="SN"){
    k1<-1
    k2<-1
  }
  
  if(type=="ST"){
    k1<-(nu/2)^(1/2)*gamma((nu-1)/2)/gamma(nu/2)
    k2<-(nu/2)*gamma((nu-2)/2)/gamma(nu/2)
  }
  
  if(type=="SSL"){
    k1<-2*nu/(2*nu-1)
    k2<-2*nu/(2*nu-2)
  }
  
  if(type=="SCN"){
    k1<-nu[1]/nu[2]^(1/2)+1-nu[1]
    k2<-nu[1]/nu[2]+1-nu[1]
  }
  Auxmed<-sqrt(2)*k1*delta
  Auxvar<-pi*k2-2*k1^2*delta^2
  med<- -Auxmed/sqrt(Auxvar)
  var<-pi/Auxvar
  return(list(med=med,var=var))
}
