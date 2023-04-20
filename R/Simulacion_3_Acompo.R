###############################################################################
###############################################################################
###
###        Experimento de simulacion 3  Compositional parametrizacion con rho               
###                            Bootstrap                                
###                            Marzo 2017                            
###


### Autor:  Agustin Perez 


# Llamada a las funciones Generacion X, REMLcompo y BETA.U.compo
source("Generacion X.R") # Funcion para generar los valores de las variables explicativas
source("REML.R") # Funcion ajuste REML
source("Estimacion BETA.R") # Calculo de beta´s y u´s estimados
source("Bootstrap Compo.R")
source("Transformaciones Vu_Theta y Theta_Vu.R")  
#source("Calculo CEB.R")  

#library(mme)
library(Matrix)
library(MASS)
library(mvtnorm)
 

set.seed(19032005)
I <-100 # Numero de vueltas para las simulaciones
B<-400
Dominios <-c(100)#,75,100,150,200,300) # c(50,75,100,150,200,300)Tamaños de los dominios D

# Declaraciones de variables

  udk<-list()
  xd1<-xd2<-xd3<-X<-list()
  Badg <- vector()
  Ved<-list()

   for(g in 1:length(Dominios)) { # Bucle para los distintos tamaños D
   D <- Dominios[g]
  
   cat("\n Comienza la simulacion", g, " con D =", D,  "\n")
   
   
   # Calculo de las variables explicativas (X) mediante la funcion equis definida en Generacion X.R
   
   Xp<-equis(D)  
   
    for(d in 1:D){
   	X[[d]]<-Xp[[1]][[d]]
   }   
   p<-Xp[[2]]
   
   # Beta, Vud, Ved: (Valores iniciales)    
   
   beta<-as.vector(rep(1,p))

   for(d in 1:D){ 
   Ved[[d]]<-matrix(-0.25,nrow=3,ncol=3)
   diag(Ved[[d]])<-1
   }
   
   thetas<-c(1,3/2,2,-1/2,-1/2,0)
   Vud<-UveU(thetas)
  
   mediasu<-c(0,0,0)
   
   
   # Declaraciones de variables
   
   theta.gorro<-vector()
   beta.gorro<-betas.gorro<-list()
   

 	 y<-list()
   
   ud<-ed<-array(rep(0,D*3),dim=c(D,3))
  
   mud<-pidk.gorro<-pi<-mudi.gorro.comp<-array(rep(0,3*D*I),dim=c(3,D,I))
   #sydk<-pdk<-array(rep(0,3*D*I),dim=c(3,D,I))
  # Bpebdk<-Bpidk<-Epebdk<-Epidk<-RBmudk<-REmudk<-Bmudk<-Emudk<-pidk<-mudk<-array(rep(0,D*3),dim=c(D,3))
    Bpidk<-Epidk<-Bmudk<-Emudk<-pidk<-mudk<-array(rep(0,D*3),dim=c(D,3))
   #BMdk<-EMdk<-REMdk<-RBMdk<-array(rep(0,D*3),dim=c(D,3))
   #Bpdk<-Epdk<-REpdk<-RBpdk<-array(rep(0,D*3),dim=c(D,3))
   
   i<-BadTot<-0
    mse_ast_mudk<- mse_ast_pidk<-E_ast_mudk<-E_ast_pidk<-B_ast_mudk<-B_ast_pidk<-0
    
    
    # Lectura
    if(g==99) {
        mse_mudk <- read.table(file="Emudk_D50.txt")
        mse_pidk <- read.table(file="Epidk_D50.txt")
       # mse_pebdk<- read.table(file="Epebdk_D50.txt")
    }
    if(g==2) {
       mse_mudk <- read.table(file="Emudk_D75.txt")
       mse_pidk<- read.table(file="Epidk_D75.txt")
      # mse_pebdk<- read.table(file="Epebdk_D75.txt")
    }
    if(g==1) {
       mse_mudk <- read.table(file="Emudk_D100.txt")
       mse_pidk <- read.table(file="Epidk_D100.txt")
      # mse_pebdk<- read.table(file="Epebdk_D100.txt")
    }
    if(g==4) {
       mse_mudk <- read.table(file="Emudk_D150.txt")
       mse_pidk <- read.table(file="Epidk_D150.txt")
      # mse_pebdk<- read.table(file="Epebdk_D150.txt")
    }
    if(g==5) {
      mse_mudk <- read.table(file="Emudk_D200.txt")
      mse_pidk <- read.table(file="Epidk_D200.txt")
      #mse_pebdk<- read.table(file="Epebdk_D200.txt")
    }
    if(g==6) {
       mse_mudk <- read.table(file="Emudk_D300.txt")
       mse_pidk <- read.table(file="Epidk_D300.txt")
      # mse_pebdk<- read.table(file="Epebdk_D300.txt")
    }
   
   
    mse_mudk<-mse_mudk^2
    mse_pidk<- mse_pidk^2
   # mse_pebdk<- mse_pebdk^2
    
    #mse_ast_mudk<- mse_ast_pidk<- mse_ast_pebdk<-E_ast_mudk<-	B_ast_mudk<-E_ast_pidk<-B_ast_pidk<-E_ast_pebdk<-B_ast_pebdk<-array(rep(0,D*3),dim=c(D,3))
    
     mse_ast_mudk<- mse_ast_pidk<- E_ast_mudk<-	B_ast_mudk<-E_ast_pidk<-B_ast_pidk<-B_ast_pebdk<-array(rep(0,D*3),dim=c(D,3))
    
   # Empieza el bucle de las simulaciones:
   
    while(i<I){
   	i<-i+1
   
   		cat("Comienza la muestra", i, "\n")
   	
   	# Simulacion de las ud y ed
   	
   	for(d in 1:D){ 
   		ed[d,] <- mvrnorm(1,mediasu,Ved[[d]])
   		ud[d,] <- mvrnorm(1,mediasu,Vud)
   	}
   	
   	
   	# Simulacion de la variable objetivo yd y de mud 
   	
  for(d in 1:D){   
   		y[[d]] <- X[[d]]%*%beta + ud[d,] + ed[d,]  
   		#mud[,d,i] <- X[[d]]%*%beta + ud[d,] 
   	}
   	
   	# Calculo de las variables pid a partir de la mud
   		
  # 		for(d in 1:D){   
  # 		pi[,d,i] <- exp(mud[,d,i])/(1+sum(exp(mud[,d,i])))
   #		}
   
   	
   			
   		# AJUSTE DEL MODELO DE AREA MEDIANTE REML (funcion REML.compo definida en REML.R)
   	
   	  	fit <- try(REMLcompo(X, y, D, Ved, Vud), TRUE)
   	
   	# Si hay algun error en el ajuste REML tiramos la simulacuon y usamos otra nueva
   	if(class(fit)=="try-error"){
   		write.table(data.frame(class(fit),Dominios[g],D,i), file="WARNING indep.txt", append=TRUE, col.names=FALSE)
   		cat("\t Muestra", i, " rechazada por try-error\n")
   		i <- i-1
   		BadTot <- BadTot+1
   		}
   		else{
   		cat("Dominio", D, "Muestra" , i,"\n")
   		
            		
            		
            	  if(fit[[3]]<100) {  	# si en el REML no se ha llegado a 100 iteraciones nos vale la estiomacion
            	
            	 	theta.gorro[i] <-fit[1]
            	
            		theta.hat<-theta.gorro[[i]]
            
      # Calculo de los beta´s y las u´s mediante la funcion BETA.U.compo definida en Estimacion BETA.R)    
            		
        beta.gorro[[i]] <- BETA.U.compo(X, y, D, Ved ,theta.hat)
        
        
        betas.gorro<-as.array(beta.gorro[[i]][[1]],dim=c(p*1)) 
         
         u<- beta.gorro[[i]][[2]]   
        
          	  
       # Estimador EBLUP del parametro poblacional
            
         for(d in 1:D) { 
               mudi.gorro.comp[,d,i] <- X[[d]]%*%betas.gorro+ u[[d]]
         }
                
         # Obtencion de laos valors de los pidk a partir de los valores anteriores
                
         for(d in 1:D){   
   	        	pidk.gorro[,d,i] <- exp(mudi.gorro.comp[,d,i])/(1+sum(exp(mudi.gorro.comp[,d,i])))
   	        	}     
   	     
       
                
   	    bootcompo <- BOOT.compo(X, D, Ved, theta.hat, betas.gorro,i,B)  
   	    
   	    mse_ast_mudk<- mse_ast_mudk+ t(bootcompo[[1]])
   	    mse_ast_pidk<- mse_ast_pidk+ t(bootcompo[[2]])
   	    #mse_ast_pebdk<- mse_ast_pebdk+ t(bootcompo[[3]])
   	    
   	    E_ast_mudk<-E_ast_mudk+(t(bootcompo[[1]])-mse_mudk)^2
   	    B_ast_mudk<-B_ast_mudk+(t(bootcompo[[1]])-mse_mudk)
   	    
   	    E_ast_pidk<-E_ast_pidk+( t(bootcompo[[2]])-mse_pidk)^2
   	    B_ast_pidk<-B_ast_pidk+(t(bootcompo[[2]])-mse_pidk)
   	    
   	   # E_ast_pebdk<-E_ast_pebdk+( t(bootcompo[[3]])-mse_pebdk)^2
   	   # B_ast_pebdk<-B_ast_pebdk+(t(bootcompo[[3]])-mse_pebdk)
	        
            	  }# fin bucle if(fit[[3]])
   		
   			
   				}# fin bucle else
   	  	
   	  	
   	   }# Fin bucle I (simulaciones)
   
   
    mse_ast_mudk<-mse_ast_mudk/I
    mse_ast_pidk<-mse_ast_pidk/I
   # mse_ast_pebdk<-mse_ast_pebdk/I
   
    E_ast_mudk<-(E_ast_mudk/I)^0.5
   	B_ast_mudk<-B_ast_mudk/I
   	    
    E_ast_pidk<-(E_ast_pidk/I)^0.5
    B_ast_pidk<-B_ast_pidk/I
    
   # E_ast_pebdk<-(E_ast_pebdk/I)^0.5
   # B_ast_pebdk<-B_ast_pebdk/I
	        
   
   # RE_ast_mudk<-10^2*E_ast_mudk/mse_mudk
   # RB_ast_mudk<-10^2*B_ast_mudk/mse_mudk
     
   # RE_ast_pidk<-10^2*E_ast_pidk/mse_pidk
   # RB_ast_pidk<-10^2*B_ast_pidk/mse_pidk
    
   # RE_ast_pebdk<-10^2*E_ast_pebdk/mse_pebdk
   # RB_ast_pebdk<-10^2*B_ast_pebdk/mse_pebdk
   
    # E_ast_muk<-RE_ast_muk<-B_ast_muk<-RB_ast_muk<-E_ast_pik<-RE_ast_pik<-B_ast_pik<-RB_ast_pik<-E_ast_pebk<-RE_ast_pebk<-B_ast_pebk<-RB_ast_pebk<-vector()
      E_ast_muk<-B_ast_muk<-E_ast_pik<-B_ast_pik<-vector()
   for(k in 1:3){
   	E_ast_muk[k]<-mean(E_ast_mudk[,k],na.rm=TRUE)
  # 	RE_ast_muk[k]<-mean(RE_ast_mudk[,k],na.rm=TRUE)
   	B_ast_muk[k]<-mean(B_ast_mudk[,k],na.rm=TRUE)
  # 	RB_ast_muk[k]<-mean(RB_ast_mudk[,k],na.rm=TRUE)
   	E_ast_pik[k]<-mean(E_ast_pidk[,k],na.rm=TRUE)
  # 	RE_ast_pik[k]<-mean(RE_ast_pidk[,k],na.rm=TRUE)
   	B_ast_pik[k]<-mean(B_ast_pidk[,k],na.rm=TRUE)
  # 	RB_ast_pik[k]<-mean(RB_ast_pidk[,k],na.rm=TRUE)
   #	E_ast_pebk[k]<-mean(E_ast_pebdk[,k],na.rm=TRUE)
  # 	RE_ast_pebk[k]<-mean(RE_ast_pebdk[,k],na.rm=TRUE)
   #	B_ast_pebk[k]<-mean(B_ast_pebdk[,k],na.rm=TRUE)
   #	RB_ast_pebk[k]<-mean(RB_ast_pebdk[,k],na.rm=TRUE)
   }
   
   
   
   
   
    # Escritutra de resultados en ficheros
    
    if(g==99) {
        write.table(E_ast_mudk, file="E_ast_mudk_D50.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_pidk, file="E_ast_pidk_D50.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_pidk, file="E_ast_pebdk_D50.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_mudk, file="B_ast_mudk_D50.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_pidk, file="B_ast_pidk_D50.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_pebdk, file="B_ast_pebdk_D50.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RE_ast_mudk, file="RE_ast_mudk_D50.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RE_ast_pidk, file="RE_ast_pidk_D50.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RE_ast_pebdk, file="RE_ast_pebdk_D50.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RB_ast_mudk, file="RB_ast_mudk_D50.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RB_ast_pidk, file="RB_ast_pidk_D50.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RB_ast_pebdk, file="RB_ast_pebdk_D50.txt", row.names=FALSE, col.names=FALSE)
    	  resultados <- data.frame(E_ast_muk, E_ast_pik, E_ast_pebk,B_ast_muk,B_ast_pik,B_ast_pebk, RE_ast_muk,RE_ast_pik,RE_ast_pebk,RB_ast_muk,RB_ast_pik,RB_ast_pebk)
        resultados <- as.data.frame(t(resultados))
        write.table(resultados, file="Resultados_Simulacion3_Compo_D50.txt", sep="\t")
      

    }
    if(g==2) {
        write.table(E_ast_mudk, file="E_ast_mudk_D75.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_pidk, file="E_ast_pidk_D75.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_pebdk, file="E_ast_pebdk_D75.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_mudk, file="B_ast_mudk_D75.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_pidk, file="B_ast_pidk_D75.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_pebdk, file="B_ast_pebdk_D75.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RE_ast_mudk, file="RE_ast_mudk_D75.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RE_ast_pidk, file="RE_ast_pidk_D75.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RE_ast_pebdk, file="RE_ast_pebdk_D75.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RB_ast_mudk, file="RB_ast_mudk_D75.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RB_ast_pidk, file="RB_ast_pidk_D75.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RB_ast_pebdk, file="RB_ast_pebdk_D75.txt", row.names=FALSE, col.names=FALSE)
    	  resultados <- data.frame(E_ast_muk, E_ast_pik, E_ast_pebk,B_ast_muk,B_ast_pik,B_ast_pebk, RE_ast_muk,RE_ast_pik,RE_ast_pebk,RB_ast_muk,RB_ast_pik,RB_ast_pebk)
        resultados <- as.data.frame(t(resultados))
        write.table(resultados, file="Resultados_Simulacion3_Compo_D75.txt", sep="\t")
    }
    if(g==1) {
        write.table(E_ast_mudk, file="E_ast_mudk_D100.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_pidk, file="E_ast_pidk_D100.txt", row.names=FALSE, col.names=FALSE)
    	 # write.table(E_ast_pebdk, file="E_ast_pebdk_D100.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_mudk, file="B_ast_mudk_D100.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_pidk, file="B_ast_pidk_D100.txt", row.names=FALSE, col.names=FALSE)
    	 # write.table(B_ast_pebdk, file="B_ast_pebdk_D100.txt", row.names=FALSE, col.names=FALSE)
    	 # write.table(RE_ast_mudk, file="RE_ast_mudk_D100.txt", row.names=FALSE, col.names=FALSE)
    	 # write.table(RE_ast_pidk, file="RE_ast_pidk_D100.txt", row.names=FALSE, col.names=FALSE)
    	 # write.table(RE_ast_pebdk, file="RE_ast_pebdk_D100.txt", row.names=FALSE, col.names=FALSE)
    	 # write.table(RB_ast_mudk, file="RB_ast_mudk_D100.txt", row.names=FALSE, col.names=FALSE)
    	 # write.table(RB_ast_pidk, file="RB_ast_pidk_D100.txt", row.names=FALSE, col.names=FALSE)
    	 # write.table(RB_ast_pebdk, file="RB_ast_pebdk_D100.txt", row.names=FALSE, col.names=FALSE)
    	  resultados <- data.frame(E_ast_muk, E_ast_pik,B_ast_muk,B_ast_pik)
        resultados <- as.data.frame(t(resultados))
        write.table(resultados, file="Resultados_Simulacion3_Compo_D100.txt", sep="\t")
    }
    
     if(g==4) {
        write.table(E_ast_mudk, file="E_ast_mudk_D150.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_pidk, file="E_ast_pidk_D150.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_pebdk, file="E_ast_pebdk_D150.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_mudk, file="B_ast_mudk_D150.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_pidk, file="B_ast_pidk_D150.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_pebdk, file="B_ast_pebdk_D150.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_mudk, file="RE_ast_mudk_D150.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_pidk, file="RE_ast_pidk_D150.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_pebdk, file="RE_ast_pebdk_D150.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_mudk, file="RB_ast_mudk_D150.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_pidk, file="RB_ast_pidk_D150.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_pebdk, file="RB_ast_pebdk_D150.txt", row.names=FALSE, col.names=FALSE)
    	  resultados <- data.frame(E_ast_muk, E_ast_pik, E_ast_pebk,B_ast_muk,B_ast_pik,B_ast_pebk, RE_ast_muk,RE_ast_pik,RE_ast_pebk,RB_ast_muk,RB_ast_pik,RB_ast_pebk)
        write.table(resultados, file="Resultados_Simulacion3_Compo_D150.txt", sep="\t")
    }
    if(g==5) {
        write.table(E_ast_mudk, file="E_ast_mudk_D200.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_pidk, file="E_ast_pidk_D200.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_pebdk, file="E_ast_pebdk_D200.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_mudk, file="B_ast_mudk_D200.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_pidk, file="B_ast_pidk_D200.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_pebdk, file="B_ast_pebdk_D200.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RE_ast_mudk, file="RE_ast_mudk_D200.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RE_ast_pidk, file="RE_ast_pidk_D200.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RE_ast_pebdk, file="RE_ast_pebdk_D200.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RB_ast_mudk, file="RB_ast_mudk_D200.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RB_ast_pidk, file="RB_ast_pidk_D200.txt", row.names=FALSE, col.names=FALSE)
    	  rite.table(RB_ast_pebdk, file="RB_ast_pebdk_D200.txt", row.names=FALSE, col.names=FALSE)
    	  resultados <- data.frame(E_ast_muk, E_ast_pik, E_ast_pebk,B_ast_muk,B_ast_pik,B_ast_pebk, RE_ast_muk,RE_ast_pik,RE_ast_pebk,RB_ast_muk,RB_ast_pik,RB_ast_pebk)
        resultados <- as.data.frame(t(resultados))
        write.table(resultados, file="Resultados_Simulacion3_Compo_D200.txt", sep="\t")
    }
    if(g==6) {
        write.table(E_ast_mudk, file="E_ast_mudk_D300.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_pidk, file="E_ast_pidk_D300.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(E_ast_pebdk, file="E_ast_pebdk_D300.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_mudk, file="B_ast_mudk_D300.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(B_ast_pidk, file="B_ast_pidk_D300.txt", row.names=FALSE, col.names=FALSE)
    	  rite.table(B_ast_pebdk, file="B_ast_pebdk_D300.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RE_ast_mudk, file="E_ast_mudk_D300.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RE_ast_pidk, file="E_ast_pidk_D300.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RE_ast_pebdk, file="E_ast_pebdk_D300.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RB_ast_mudk, file="B_ast_mudk_D300.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RB_ast_pidk, file="B_ast_pidk_D300.txt", row.names=FALSE, col.names=FALSE)
    	  write.table(RB_ast_pebdk, file="B_ast_pebdk_D300.txt", row.names=FALSE, col.names=FALSE)
    	  resultados <- data.frame(E_ast_muk, E_ast_pik, E_ast_pebk,B_ast_muk,B_ast_pik,B_ast_pebk, RE_ast_muk,RE_ast_pik,RE_ast_pebk,RB_ast_muk,RB_ast_pik,RB_ast_pebk)
        resultados <- as.data.frame(t(resultados))
        write.table(resultados, file="Resultados_Simulacion3_Compo_D300.txt", sep="\t")
    }
    


   
   }# Fin bucle g
   