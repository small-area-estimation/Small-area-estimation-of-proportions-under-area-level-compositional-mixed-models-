###############################################################################
###############################################################################
###
###          Experimento de simulacion 2 estudio de variacion de L               
###                                                            
###                            Diciembre 2017                            
###


### Autor:  Lola y Agus   


# Llamada a las funciones Generacion X, REMLcompo y BETA.U.compo
source("Generacion X.R") # Funcion para generar los valores de las variables explicativas
source("REML.R") # Funcion ajuste REML
source("Estimacion BETA.R") # Calculo de beta´s y u´s estimados
source("Transformaciones Vu_Theta y Theta_Vu.R")  
source("Calculo CEB.R")  
source("rnorM.R")

library(Matrix)
library(MASS)
#library(mvtnorm)


set.seed(19032005)
I <- 500   # Numero de vueltas para las simulaciones
L <- 10000  # Se hace para c(10000, 1000, 100, 10)
Dominios <- 100

# Declaraciones de variables

udk <- list()
xd1 <- xd2 <- xd3 <- X <- list()
Badg <- vector()
Ved <- list()

for(g in 1:length(L)) { # Bucle para los distintos tamaños L
  D <- Dominios
  cat("\n Comienza la simulacion", g, " con D =", D,  "\n")
  
  # Calculo de las variables explicativas (X) mediante la funcion equis definida en Generacion X.R
  Xp <- equis(D)  
  
  for(d in 1:D){
    X[[d]] <- Xp[[1]][[d]]
  }   
  p <- Xp[[2]]
  
  # Beta, Vud, Ved: (Valores iniciales)    
  beta <- as.vector(rep(1,p))
  
  for(d in 1:D){ 
    Ved[[d]] <- matrix(-0.25,nrow=3,ncol=3)
    diag(Ved[[d]]) <- 1
  }
  
  thetas <- c(1,3/2,2,-1/2,-1/2,0)
  
  # Matriz Vud
  Vud <- UveU(thetas)
  mediasu <- c(0,0,0)
  
  
  # Declaraciones de variables
  theta.gorro <- vector()
  beta.gorro <- betas.gorro <- list()
  y <- list()
  
  ud <- ed <- array(rep(0,D*3), dim=c(D,3))
  
  mud <- pidk.gorro <- pi <- mudi.gorro.comp <- pebdk <- array(rep(0,3*D*I),dim=c(3,D,I))
  sydk <- pdk <- array(rep(0,3*D*I),dim=c(3,D,I))
  Bpebdk <- Epebdk <- Bpidk <- Epidk <- Bmudk <- Emudk <- pidk <- mudk <- array(rep(0,D*3),dim=c(D,3))
  BMdk <- EMdk <- array(rep(0,D*3),dim=c(D,3))
  Bpdk <- Epdk <- array(rep(0,D*3),dim=c(D,3))
  
  
  i <- BadTot <- 0
  
  # Empieza el bucle de las simulaciones:
  while(i<I){
    i <- i+1
    cat("Comienza la muestra", i, "\n")
    
    # Simulacion de las ud y ed
    for(d in 1:D){ 
      ed[d,] <- var.norm.mult(1,mediasu,Ved[[d]])
      ud[d,] <- var.norm.mult(1,mediasu,Vud)
    }
    
    # Simulacion de la variable objetivo yd y de mud 
    for(d in 1:D){   
      y[[d]] <- X[[d]]%*%beta + ud[d,] + ed[d,]  
      mud[,d,i] <- X[[d]]%*%beta + ud[d,] 
    }
    
    # Calculo de las variables pid a partir de la mud
    for(d in 1:D){   
      pi[,d,i] <- exp(mud[,d,i])/(1+sum(exp(mud[,d,i])))
    }
    
    
    # AJUSTE DEL MODELO DE AREA MEDIANTE REML (funcion REML.compo definida en REML.R)
    fit <- try(REMLcompo(X, y, D, Ved, Vud, MAXITER=100), TRUE)
    if(class(fit)=="try-error"){
      write.table(data.frame(class(fit),L[g],D,i), file="WARNING indep.txt", append=TRUE, row.names=FALSE)
      cat("\t Muestra", i, " rechazada por try-error\n")
      i <- i-1
      BadTot <- BadTot+1
    }
    else{
      cat("Grupo L= ", L[g], "Muestra" , i,"\n")
      write.table(data.frame(fit[[3]],D,i), file="ITER.txt", append=TRUE, row.names=FALSE, col.names=FALSE)
      
      if(fit[[3]]<100) {  	# si en el REML no se ha llegado a 100 iteraciones nos vale la estiomacion
        theta.gorro[i] <-fit[1]
        theta.hat<-theta.gorro[[i]]
        
        # Calculo de los beta´s y las u´s mediante la funcion BETA.U.compo definida en Estimacion BETA.R)    
        beta.gorro[[i]] <- BETA.U.compo(X, y, D, Ved ,theta.hat)
        betas.gorro <- as.array(beta.gorro[[i]][[1]],dim=c(p*1)) 
        u<- beta.gorro[[i]][[2]]   
        
        # Estimador EBLUP del parametro poblacional
        for(d in 1:D)
          mudi.gorro.comp[,d,i] <- X[[d]]%*%betas.gorro+ u[[d]]
        
        # Obtencion de laos valors de los pidk a partir de los valores anteriores
        for(d in 1:D)
          pidk.gorro[,d,i] <- exp(mudi.gorro.comp[,d,i])/(1+sum(exp(mudi.gorro.comp[,d,i])))
        
        #### FUNCION CEB
        t1 <- proc.time()
        Vudc <- UveU(theta.hat)
        pebdk[,,i] <- ceb(D,X,y,Vudc,betas.gorro,Ved,L[g])
        t2 <- proc.time(); cat("\n Tiempo de ceb", (t2-t1)[3], "\n")
        write.table(data.frame(L[g], i, (t2-t1)[3]), file="tiempo CEB.txt", append=TRUE, row.names=FALSE, col.names=FALSE)
        
      }     #fin bucle if(fit[[3]]) del limite de 100 iteraciones
      
    }# fin bucle else
    
  }# Fin bucle I (simulaciones)
  
  # Calculo de errores y sesgos absolutos y relativos
  
  for(d in 1:D){ 
    for(k in 1:3){
      mudk[d,k]<-mean(mud[k,d,],na.rm=TRUE)
      pidk[d,k]<-mean(pi[k,d,],na.rm=TRUE)
      Emudk[d,k]<-(mean((mudi.gorro.comp[k,d,]-mud[k,d,])^2,na.rm=TRUE))^0.5
      Bmudk[d,k]<-mean((mudi.gorro.comp[k,d,]-mud[k,d,]),na.rm=TRUE)
      Epidk[d,k]<-(mean((pidk.gorro[k,d,]-pi[k,d,])^2,na.rm=TRUE))^0.5
      Bpidk[d,k]<-mean(pidk.gorro[k,d,]-pi[k,d,],na.rm=TRUE)
      Epebdk[d,k]<-(mean((pebdk[k,d,]-pi[k,d,])^2,na.rm=TRUE))^0.5
      Bpebdk[d,k]<-mean(pebdk[k,d,]-pi[k,d,],na.rm=TRUE)
    }
  }
  
  
  Emuk<-REmuk<-Bmuk<-RBmuk<-EMk<-REMk<-BMk<-RBMk<-Epik<-REpik<-Bpik<-RBpik<-Epk<-REpk<-Bpk<-RBpk<-Bpebk<-RBpebk<-Epebk<-REpebk<-vector()
  
  for(k in 1:3){
    Emuk[k] <- mean(Emudk[,k],na.rm=TRUE)
    Bmuk[k] <- mean(Bmudk[,k],na.rm=TRUE)
    BMk[k] <- mean(BMdk[,k],na.rm=TRUE)
    Epik[k] <- mean(Epidk[,k],na.rm=TRUE)
    Bpik[k] <- mean(Bpidk[,k],na.rm=TRUE)
    Epebk[k] <- mean(Epebdk[,k],na.rm=TRUE)
    Bpebk[k] <- mean(Bpebdk[,k],na.rm=TRUE)
  }
  
  escritura <- paste("L", L[g], ".txt", sep="")
  # Escritutra de resultados en ficheros
  write.table(Emudk, file=paste("Emudk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
  write.table(Epidk, file=paste("Epidk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
  write.table(Bmudk, file=paste("Bmudk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
  write.table(Bpidk, file=paste("Bpidk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
  write.table(Epebdk, file=paste("Epebdk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
  write.table(Bpebdk, file=paste("Bpebdk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
  
}# Fin bucle g





###