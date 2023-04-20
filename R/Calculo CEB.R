###############################################################################
###
###              Calculo CEB (Acompo) (Acompo)              
###              
###                 Julio 2017    
### Autor: Agustin Perez
###############################################################################

#Modificaciones:
#Diciembre 2017
#Lola y Agus
#Generacion de normal multivariante por Cholesky, correcion de errores y mejora de funcionamiento

ceb <- function(D,X,y,Vudc,beta.est,Ved,L) {
  
  DosL <- L*2
  udc <- array(rep(0,3*D*DosL),dim=c(3,D,DosL))
  
  
  
  Ad <- mucdk.gorro <- array(rep(0,3*D*DosL), dim=c(3,D,DosL))
  Bd <- array(rep(0,D*DosL),dim=c(D,DosL))
  SumAd <- array(rep(0,3*D),dim=c(3,D))
  SumBd <- vector()
  pebdk.gorro <- array(rep(0,3*D),dim=c(3,D))
  
  Xd <- X
  yd <- y
  
  
  medias <- c(0,0,0)     
  
  # Simulacion de las ud,muc, Ad y Bd
  for(d in 1:D){ 
    for(l in 1:L){
      udc[,d,l] <- var.norm.mult(1,medias,Vudc)
      udc[,d,l+L] <- -udc[,d,l]
    }
  }
  
  
  for(d in 1:D){ 
    AUX <- solve(Ved[[d]])
    for(l in 1:DosL){
      mucdk.gorro[,d,l] <- X[[d]]%*%beta.est + udc[,d,l] 
      aux <- t(yd[[d]]-X[[d]]%*%beta.est-udc[,d,l])%*%AUX
      Bd[d,l] <- exp(-0.5*(aux%*%(yd[[d]]-X[[d]]%*%beta.est- udc[,d,l])))
      Ad[,d,l] <- exp(mucdk.gorro[,d,l])*Bd[d,l]/(1+sum(exp(mucdk.gorro[,d,l])))
    }
  }
  
  # Calculo de las sumas
  for(d in 1:D){ 
    for(l in 1:DosL){
      SumAd[,d] <- SumAd[,d]+Ad[,d,l]
    }
  }
  
  SumBd <- addmargins(Bd,2)[,DosL+1]
  
  for(d in 1:D)
    pebdk.gorro[,d] <- SumAd[,d]/SumBd[d]
  
  return(pebdk.gorro)
}
