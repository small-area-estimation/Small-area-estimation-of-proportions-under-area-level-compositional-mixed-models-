###############################################################################
###############################################################################
###
###                                      
###                                                         
###


### Editado por:
### con fecha:
### cambios realizados:



BETA.U.compo <- function(X, y, D,Ved, theta.hat) {

	
    Vd.inv<-Vd.gorro<-list()
    Xd<-yd<-list()
    p <- ncol(X[[1]])
   
    Xd<-X
    yd<-y
  
   Vudbeta<-UveU(theta.hat)
   for(d in 1:D){
   Vd.gorro[[d]] <- Vudbeta+Ved[[d]] 
   Vd.inv[[d]] <- solve(Vd.gorro[[d]]) 
   }
  # Calculo beta gorro
   
    Q.inv <- matrix(0, nrow=p, ncol=p)
    XVy <- rep(0,p)
    for(d in 1:D) {
       Q.inv <- Q.inv + t(Xd[[d]])%*%Vd.inv[[d]]%*%Xd[[d]]         # Inversa de Q. Posteriormente se calcula Q
       XVy <- XVy + t(Xd[[d]])%*%Vd.inv[[d]]%*%yd[[d]]             # Producto X^t_d  por  V^-1_d  por  y_d para las d submatrices
     }
    Q <- solve(Q.inv)
    
    betax <- Q%*%XVy
    
  # Calculo u gorro
    u <- list()
    
    for(d in 1:D){
        u[[d]] <-Vudbeta%*%Vd.inv[[d]]%*%(yd[[d]]-X[[d]]%*%betax) 
        }
   # u <- as.matrix(unlist(u))
    
    
    return(list(betax,u))
}

