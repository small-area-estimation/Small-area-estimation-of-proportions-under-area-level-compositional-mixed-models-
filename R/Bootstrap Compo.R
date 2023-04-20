###############################################################################
###############################################################################
###
###                   Bootstrap modelo area AR(1)
##                  Modelos área con efectos de tiempo independientes                    
###                                Febrero 2012                             
###

### AUTOR: Agustin Perez

### Editado por:        Agustin Perez Martin
### con fecha:          Octubre 2017
### cambios realizados: modificacion del beta en la generacion de yb[[d]]. Pequenos cambios de formato



BOOT.compo <-function(X, D, Vedb, thetasb,betasb,i,B=50){
  source("REML.R")
  source("Estimacion BETA.R")
  
  
  Xd <- X
  Vudb <- UveU(thetasb)
  
  edb <- udb <- array(rep(0,D*3), dim=c(D,3))
  theta.gorro.ast <- vector()
  beta.gorro.ast <- list()
  pidk.gorro.ast <- mudk.gorro.ast <- array(rep(0,3*D*B), dim=c(3,D,B))
  
  b <- BadTot2 <- difmu <- difpi <- 0
  medias <- c(0, 0, 0)
  mudb <- pidkb <- array(rep(0,3*D*B), dim=c(3,D,B))
  yb <- list()
  excepcion <- vector()

    while(b<B){
    b <- b+1
    cat("Comienza la muestra bootstrap", b, "\n")
    
    # Simulacion de las ud y ed
    for(d in 1:D){ 
      edb[d,] <- mvrnorm(1, medias, Vedb[[d]])
      udb[d,] <- mvrnorm(1, medias, Vudb)
    }
    
    # Simulacion de la variable objetivo yd y de mud 
    for(d in 1:D){   
      mudb[,d,b] <- X[[d]]%*%betasb + udb[d,] 
      yb[[d]] <- X[[d]]%*%betasb + udb[d,] + edb[d,]  ####OJO, he cambiado beta por betasb
      
    }
    for(d in 1:D){   
      pidkb[,d,b] <- exp(mudb[,d,b])/(1+sum(exp(mudb[,d,b])))
    }     
    
    #Ajuste modelo bootstrap
    fit2 <- try(REMLcompo(X, yb, D, Vedb, Vudb), TRUE)
    
    # Si hay algun error en el ajuste REML tiramos la simulacuon y usamos otra nueva
    if(class(fit2)=="try-error"){
      excepcion <- c(excepcion, b)
      write.table(data.frame(class(fit2),D,i,b), file="WARNING.txt", append=TRUE, col.names=FALSE)
      b <- b-1
      BadTot2 <- BadTot2+1
    }
    else{
      cat("D=",D,"Muestra Bootstrap" , b,"\n")
      if(fit2[[3]]<100) {  	# si en el REML no se ha llegado a 100 iteraciones nos vale la estiomacion
        theta.gorro.ast[b] <- fit2[1]
        theta.ast <- theta.gorro.ast[[b]]
        
        # Calculo de los beta´s y las u´s mediante la funcion BETA.U.compo definida en Estimacion BETA.R)    
        beta.gorro.ast[[b]] <- BETA.U.compo(X, yb, D, Vedb ,theta.ast)
        betas.gorro.ast <- as.array(beta.gorro.ast[[b]][[1]], dim=c(p*1)) 
        u <- beta.gorro.ast[[b]][[2]]   
        
        # Estimador EBLUP del parametro poblacional
        for(d in 1:D) { 
          mudk.gorro.ast[,d,b] <- X[[d]]%*%betas.gorro.ast+ u[[d]]
        }
        # Obtencion de laos valors de los pidk a partir de los valores anteriores
        for(d in 1:D){   
          pidk.gorro.ast[,d,b] <- exp(mudk.gorro.ast[,d,b])/(1+sum(exp(mudk.gorro.ast[,d,b])))
        }     
        
        #### FUNCION CEB
        # pebdk.gorro.ast[,,b]<-t(ceb(D,X,y,theta.hat,betas.gorro,Ved))
        
        difmu <- difmu + (mudk.gorro.ast[,,b]-mudb[,,b])^2
        difpi <- difpi + (pidk.gorro.ast[,,b]-pidkb[,,b])^2
      }
      else {
        b <- b-1
        BadTot2<-BadTot2+fit2[[4]]
      }
      cat( "\n")
    }
  }
  
  mse_mudk_ast <- difmu/B
  mse_pidk_ast <- difpi/B
  
  return(list(mse_mudk_ast,mse_pidk_ast))
}




###