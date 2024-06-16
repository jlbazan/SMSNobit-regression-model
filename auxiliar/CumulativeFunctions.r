
## NI 
## Cumulative distribution function
###########################################################

cdfNI<-function(x,mu,sigma2,nu,type="N"){
resp<-matrix(0,length(x),1)
if(type=="N"){

cdf<- function(x){dnorm((x-mu)/sqrt(sigma2))/sqrt(sigma2)}
          }

if(type=="T"){
cdf<- function(y){
z=(y-mu)/sqrt(sigma2)
cdf=dt(z,nu)/sqrt(sigma2)
                 }
          }

if(type=="SL"){

cdf<- function(x){
f <- function(u){ nu*u^(nu - 1)*dnorm(x,mu,sqrt(u^(-1)*sigma2))}
cdf <- integrate(Vectorize(f),0,1)$value
                 }
          }


if(type=="CN"){
cdf<- function(x){
z=(x-mu)/sqrt(sigma2)
z2=z*sqrt(nu[2])
cdf=nu[1]*dnorm(z2)/sqrt(sigma2/nu[2])+(1-nu[1])*dnorm(z)/sqrt(sigma2)
                 }
           }

for(i in 1:length(x)) {
resp[i]<-integrate(Vectorize(cdf),-Inf,x[i])$value
}
return(resp)
}



## SNI 
## Cumulative distribution function
###########################################################

cdfSNI<-function(x,mu,sigma2,delta,nu,type="SN"){
  lambda<-delta/sqrt(1-delta^2)
  resp<-matrix(0,length(x),1)
  if(type=="SN"){
    cdf<- function(x){2*dnorm((x-mu)/sqrt(sigma2))*pnorm(lambda*(x-mu)/sqrt(sigma2))/sqrt(sigma2)}
    
  }
  
  if(type=="ST"){
    cdf<- function(x){
      z=(x-mu)/sqrt(sigma2)
      cdf=2*dt(z,nu)*pt(sqrt(nu+1)*lambda*z/sqrt(nu+z^2),nu+1)/sqrt(sigma2)
    }
  }
  
  if(type=="SSL"){
    require(sn)
    cdf<- function(x){
      f <- function(u){ 2*nu*u^(nu - 1)*dnorm(x,mu,sqrt(u^(-1)*sigma2))*pnorm(u^(1/2)*lambda*(x-mu)/sqrt(sigma2))}
      cdf <- integrate(Vectorize(f),0,1)$value
    }
  }
  
  
f <- function(x) 1

integrate(Vectorize(f),0,1)$value
  
  if(type=="SCN"){
    cdf<- function(x){
      z=(x-mu)/sqrt(sigma2)
      z2=z*sqrt(nu[2])
      cdf=2*nu[1]*dnorm(z2)*pnorm(lambda*z2)/sqrt(sigma2/nu[2])+2*(1-nu[1])*dnorm(z)*pnorm(lambda*z)/sqrt(sigma2)
    }
  }
  
  
  for(i in 1:length(x)) {
    resp[i]<-integrate(Vectorize(cdf),-Inf,x[i])$value
  }
  return(resp)
}