###############################################################################
###
###                 Transformacion Vu a Thetas y Thetas a Vu     
###                  Derivadas
###                   Mayo 2017
                               
###  Autor: Lola y Agus

UveU <- function(thetas) {
  Vuu<-matrix(0,nrow=3,ncol=3)
  sup<-c(thetas[4]*sqrt(thetas[1])*sqrt(thetas[2]),thetas[5]*sqrt(thetas[1])*sqrt(thetas[3]),thetas[6]*sqrt(thetas[2])*sqrt(thetas[3]))
  Vuu[upper.tri(Vuu)]<-sup
  Vuu<-Vuu+t(Vuu)
  diag(Vuu)<-c(thetas[1],thetas[2],thetas[3])
  
  return(Vuu)
}


Thetas<- function(Vus) {
Zethas<-vector()
Zethas<-c(Vus[1,1],Vus[2,2],Vus[3,3],Vus[1,2]/(sqrt(Vus[1,1]*Vus[2,2])),Vus[1,3]/(sqrt(Vus[1,1]*Vus[3,3])),Vus[2,3]/(sqrt(Vus[2,2]*Vus[3,3])))
  
return(Zethas)
}



FirstDer<-function(thetas){
    Vud1<-matrix(0,nrow=3,ncol=3)
    Vud1[1,1]<-1
    Vud1[1,2]<-Vud1[2,1]<- thetas[4]*sqrt(thetas[2])/(2*sqrt(thetas[1]))
    Vud1[1,3]<-Vud1[3,1]<- thetas[5]*sqrt(thetas[3])/(2*sqrt(thetas[1]))
    
    Vud2<-matrix(0,nrow=3,ncol=3)
    Vud2[2,2]<-1
    Vud2[1,2]<-Vud2[2,1]<- thetas[4]*sqrt(thetas[1])/(2*sqrt(thetas[2]))
    Vud2[2,3]<-Vud2[3,2]<- thetas[6]*sqrt(thetas[3])/(2*sqrt(thetas[2]))
    	
    Vud3<-matrix(0,nrow=3,ncol=3)
    Vud3[3,3]<-1
    Vud3[1,3]<-Vud3[3,1]<- thetas[5]*sqrt(thetas[1])/(2*sqrt(thetas[3]))
    Vud3[2,3]<-Vud3[3,2]<- thetas[6]*sqrt(thetas[2])/(2*sqrt(thetas[3]))
    
    Vud4<-matrix(0,nrow=3,ncol=3)
    Vud4[1,2]<-Vud4[2,1]<-sqrt(thetas[1])*sqrt(thetas[2])
   
    Vud5<-matrix(0,nrow=3,ncol=3)
    Vud5[1,3]<-Vud5[3,1]<-sqrt(thetas[1])*sqrt(thetas[3])
    
    Vud6<-matrix(0,nrow=3,ncol=3)
    Vud6[2,3]<-Vud6[3,2]<-sqrt(thetas[2])*sqrt(thetas[3])	
	
	return(list(Vud1,Vud2,Vud3,Vud4,Vud5,Vud6))
	
	
	
	
}