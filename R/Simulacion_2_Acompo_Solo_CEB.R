###############################################################################
###############################################################################
###
###        Experimento de simulacion 2  Compositional parametrizacion con rho               
###                                                            
###                            Marzo 2017                            
###


### Autor:  Lola    


# Llamada a las funciones Generacion X, REMLcompo y BETA.U.compo
source("Generacion X.R") # Funcion para generar los valores de las variables explicativas
source("REML.R") # Funcion ajuste REML
source("Estimacion BETA.R") # Calculo de beta´s y u´s estimados
source("Transformaciones Vu_Theta y Theta_Vu.R")  
#source("Multmme.R")  
source("Calculo CEB.R")  

library(Matrix)
library(MASS)
library(mvtnorm)
 

set.seed(19032005)
I <- 1000 # Numero de vueltas para las simulaciones
L <- 1000
Dominios <- c(300)#, 75,150,200,300)#,75,100,150,200,300)  #Tamaños de los dominios D

# Declaraciones de variables

udk <- list()
xd1 <- xd2 <- xd3 <- X <- list()
Badg <- vector()
Ved <- list()

for (g in 1:length(Dominios)) {
	
	# Bucle para los distintos tamaños D
	
	D <- Dominios[g]
	
	cat("\n Comienza la simulacion", g, " con D =", D,  "\n")
	
	
	# Calculo de las variables explicativas (X) mediante la funcion equis definida en Generacion X.R
	
	Xp <- equis(D)
	
	for (d in 1:D) {
		X[[d]] <- Xp[[1]][[d]]
	}
	p <- Xp[[2]]
	
   # Beta, Vud, Ved: (Valores iniciales)    
	
	beta <- as.vector(rep(1, p))

	for(d in 1:D) {
		Ved[[d]] <- matrix(-0.25, nrow = 3, ncol = 3)
		diag(Ved[[d]]) <- 1
	}
   
	thetas <- c(1, 3/2, 2, -1/2, -1/2, 0)
   
   # Matriz Vud
   
	Vud <- UveU(thetas)
   
   mediasu <- as.vector(rep(0,3))
   
   
   # Declaraciones de variables
   
   theta.gorro <- vector()
   beta.gorro <- betas.gorro <- list()
   y <- list()
   
   ud <- ed <-array(rep(0,D*3),dim=c(D,3))
  
   mud <- pidk.gorro <- pi <- mudi.gorro.comp <- pebdk <-array(rep(0,3*D*I),dim=c(3,D,I))
   
 
   pidk<-mudk<- array(rep(0,D*3),dim=c(D,3))
 
   Bpebdk<-Epebdk<- array(rep(0,D*3),dim=c(D,3))
   Bmudk<-Emudk<- array(rep(0,D*3),dim=c(D,3))
   
   i<-BadTot<-0
   
   # Empieza el bucle de las simulaciones:
   
    while(i<I){
   	i<-i+1
   
   		cat("Comienza la muestra", i, "\n")
   	
   	# Simulacion de las ud y ed
   	
   	for(d in 1:D){ 
   		ed[d,] <- rmvnorm(1,mediasu,Ved[[d]])
   		ud[d,] <- rmvnorm(1,mediasu,Vud)
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
   	
   	  	fit <- try(REMLcompo(X, y, D, Ved, Vud), TRUE)
   	
   
      # Si hay algun error en el ajuste REML tiramos la simulacion y usamos otra nueva
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
        
        
        betas.gorro <- as.array(beta.gorro[[i]][[1]], dim = c(p * 1)) 
         
         u<- beta.gorro[[i]][[2]]   
        

       # Estimador EBLUP del parametro poblacional
            
         for(d in 1:D) { 
               mudi.gorro.comp[,d,i] <- X[[d]]%*%betas.gorro+ u[[d]]
         }
                
         # Obtencion de laos valors de los pidk a partir de los valores anteriores
                
         for(d in 1:D){   
   	        	pidk.gorro[,d,i] <- exp(mudi.gorro.comp[,d,i])/(1+sum(exp(mudi.gorro.comp[,d,i])))
   	        	}     
   	   
   	   #### FUNCION CEB
         
         Vudc <- UveU(theta.hat)
        
         
         pebdk[,,i]<-ceb(D,X,y,Vudc,betas.gorro,Ved,L)
         
     
        }# fin bucle if(fit[[3]])
        
      }# fin bucle else segundo
    #}# fin bucle else primero (por error espacio parametrico)
    
  }# Fin bucle I (simulaciones)
            
   
   # Calculo de errores y sesgos absolutos y relativos
   
  	 for(d in 1:D){ 
  	 	for(k in 1:3){
  	 		Epebdk[d,k]<-(mean((pebdk[k,d,]-pi[k,d,])^2,na.rm=TRUE))^0.5
  	 	  Bpebdk[d,k]<-mean(pebdk[k,d,]-pi[k,d,],na.rm=TRUE)
  	 	  Emudk[d,k]<-(mean((mudi.gorro.comp[k,d,]-mud[k,d,])^2,na.rm=TRUE))^0.5
  	 	  Bmudk[d,k]<-mean((mudi.gorro.comp[k,d,]-mud[k,d,]),na.rm=TRUE)
  	 	}}
   
   
   
   Epebk <- Bpebk <- Emud <- Bmud <- vector()
    
   for(k in 1:3){
    Epebk[k]<-mean(Epebdk[,k],na.rm=TRUE)
    Bpebk[k]<-mean(Bpebdk[,k],na.rm=TRUE)
    Emud[k]<-mean(Emudk[,k],na.rm=TRUE)
    Bmud[k]<-mean(Bmudk[,k],na.rm=TRUE)
   }
    
    
    # Escritutra de resultados en ficheros
   
   
 escritura <- paste("D", Dominios[g], ".txt", sep="")
 # Escritutra de resultados en ficheros
 # write.table(Emudk, file=paste("Emudk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
 # write.table(Epidk, file=paste("Epidk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
 # write.table(Bmudk, file=paste("Bmudk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
 # write.table(Bpidk, file=paste("Bpidk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
 # write.table(EMdk, file=paste("EMdk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
 # write.table(Epdk, file=paste("Epdk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
 # #write.table(BMdk, file=paste("BMdk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
 # write.table(Bpdk, file=paste("Bpdk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
 write.table(Epebdk, file=paste("Epebdk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
 write.table(Bpebdk, file=paste("Bpebdk", escritura,sep="_"), row.names=FALSE, col.names=FALSE)
 resultados <- data.frame(Bmud,Bpebk,Emud,Epebk)
 resultados <- as.data.frame(t(resultados))
 write.table(resultados, file=paste("Resultados_Simulacion2_Compo_CEB", escritura,sep="_"), sep="\t")
    
    # if(g==99) {
    #    
    # 	 write.table(Epebdk, file="Epebdk_D50.txt", row.names=FALSE, col.names=FALSE)
    # 	 write.table(Bpebdk, file="Bpebdk_D50.txt", row.names=FALSE, col.names=FALSE)
    # 	 resultados <- data.frame(Epebk,Bpebk)
    #     #resultados <- as.data.frame(t(resultados))
    #     write.table(resultados, file="Resultados_Simulacion2_Compo_D50_CEB.txt", sep="\t")
    # 
    # }
    # if(g==2) {
    #     write.table(Epebdk, file="Epebdk_D75.txt", row.names=FALSE, col.names=FALSE)
    # 	 write.table(Bpebdk, file="Bpebdk_D75.txt", row.names=FALSE, col.names=FALSE)
    # 	 resultados <- data.frame(Epebk,Bpebk)
    #     #resultados <- as.data.frame(t(resultados))
    #     write.table(resultados, file="Resultados_Simulacion2_Compo_D75_CEB.txt", sep="\t")
    # 
    # }
    # if(g==1) {
    #     write.table(Epebdk, file="Epebdk_D100.txt", row.names=FALSE, col.names=FALSE)
    # 	 write.table(Bpebdk, file="Bpebdk_D100.txt", row.names=FALSE, col.names=FALSE)
    # 	 resultados <- data.frame(Epebk,Bpebk)
    #     #resultados <- as.data.frame(t(resultados))
    #     write.table(resultados, file="Resultados_Simulacion2_Compo_D100_CEB.txt", sep="\t")
    # }
    # 
    #  if(g==4) {
    #     write.table(Epebdk, file="Epebdk_D150.txt", row.names=FALSE, col.names=FALSE)
    # 	 write.table(Bpebdk, file="Bpebdk_D150.txt", row.names=FALSE, col.names=FALSE)
    # 	 resultados <- data.frame(Epebk,Bpebk)
    #     #resultados <- as.data.frame(t(resultados))
    #     write.table(resultados, file="Resultados_Simulacion2_Compo_D150_CEB.txt", sep="\t")
    # 
    # }
    # if(g==5) {
    #    write.table(Epebdk, file="Epebdk_D200.txt", row.names=FALSE, col.names=FALSE)
    # 	 write.table(Bpebdk, file="Bpebdk_D200.txt", row.names=FALSE, col.names=FALSE)
    # 	 resultados <- data.frame(Epebk,Bpebk)
    #     #resultados <- as.data.frame(t(resultados))
    #     write.table(resultados, file="Resultados_Simulacion2_Compo_D200_CEB.txt", sep="\t")
    # }
    # if(g==6) {
    #      write.table(Epebdk, file="Epebdk_D300.txt", row.names=FALSE, col.names=FALSE)
    # 	 write.table(Bpebdk, file="Bpebdk_D300.txt", row.names=FALSE, col.names=FALSE)
    # 	 resultados <- data.frame(Epebk,Bpebk)
    #     #resultados <- as.data.frame(t(resultados))
    #     write.table(resultados, file="Resultados_Simulacion2_Compo_D300_CEB.txt", sep="\t")
    # }
    # 
    # 

   
   }# Fin bucle g
   