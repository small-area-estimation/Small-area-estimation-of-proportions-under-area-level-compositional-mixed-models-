 sink("estado.txt")
###############################################################################
###############################################################################
###
###        Experimento de simulacion 1  Compositional Data (parametrizacion con rho)               
###                                                            
###                            Marzo 2017                            
###


### Autor:  Lola    



# Llamada a las funciones Generacion X, REMLcompo y BETA.U.compo
source("Generacion X.R")  # Funcion para generar los valores de las variables explicativas
source("REML.R") # Funcion ajuste REML
source("Estimacion BETA.R")  # Calculo de beta´s y u´s estimados
source("Transformaciones Vu_Theta y Theta_Vu.R")  

library(Matrix)
library(MASS)
library(mvtnorm)
 

set.seed(19032005)
# Numero de simulaciones
K <-1000                                                                                                                                                                                                                                              
# Tamaños Dominios D
Dominios <-c(50,75,100,150,200,300)

# Inicializacion de variables 

  BIAS.theta1 <- BIAS.theta2 <- BIAS.theta3 <- BIAS.theta4 <- BIAS.theta5 <- BIAS.theta6 <- 
	BIAS <- vector()
  
  BIAS.betaest<-vector()
  RBIAS.betaest <-vector()
  MSE.betaest<-RMSE.betaest<-vector()
  RBIAS.theta1 <- RBIAS.theta2 <- RBIAS.theta3 <- RBIAS.theta4 <- RBIAS.theta5 <- RBIAS.theta6 <- 
	RBIAS.beta33 <- RBIAS <- vector()
  MSE.theta1 <- MSE.theta2 <- MSE.theta3 <- MSE.theta4 <- MSE.theta5 <- MSE.theta6 <- 
	MSE.beta11 <- MSE.beta12 <- MSE.beta21 <- MSE.beta22 <- MSE.beta31 <- MSE.beta32 <- 
	MSE.beta33 <- MSE <- vector()
  RMSE.theta1 <- RMSE.theta2 <- RMSE.theta3 <- RMSE.theta4 <- RMSE.theta5 <- RMSE.theta6 <- 
  MSE <- vector()
  
  udk<-list()
  xd1<-xd2<-xd3<-X<-list()
  Badg <- vector()
  Ved<-list()
 
 
 # Inicio del bucle de dominios

   for(g in 1:length(Dominios)) {
   D <- Dominios[g]
   
   cat("\n Comienza la simulacion", g, " con D =", D,  "\n")
   
   
   # Simulacion de las variables explicativas (funcion equis del fixhero Generacion X.R)
   
   Xp<-equis(D)  
   
   for(d in 1:D){
   	X[[d]]<-Xp[[1]][[d]]
   }   
   
   p<-Xp[[2]]
   
   # Valores iniciales parametros
   
   # Beta, thetai,Vud, Ved    
   
   beta<-as.vector(rep(1,p))
   
   for(d in 1:D){ 
   Ved[[d]]<-matrix(-0.25,nrow=3,ncol=3)
   diag(Ved[[d]])<-1
   }
   

   thetas<-c(1,3/2,2,-1/2,-1/2,0)
   Vud<-UveU(thetas)
   
   mediasu<-rep(0,3)
   
   
   # Inicializacion variables (parametros)
   theta.gorro<-vector()
   beta.gorro<-list()
   betas.gorro<-matrix(0,nrow=p,ncol=K)
  
   ed<-ud<-y<-list()
   
   mud <- mud.gorro <- list()
   
   k<-BadTot<-0
   
   # Inicio bucle simulaciones (k)
   
   while(k<K){
   	k<-k+1
   	
   	cat("Comienza la muestra", k, "\n")
   	
   	# Simulacion de los ed y ud (normales multivariantes)
   	
   	for(d in 1:D){ 
   		ed[[d]] <- rmvnorm(1,mediasu,Ved[[d]])
   		ud[[d]] <- rmvnorm(1,mediasu,Vud)
   	}
   	
   	
   	# Simulacion de la variable objetivo (y)
   	
   	for(d in 1:D){   
   		y[[d]] <- X[[d]]%*%beta + t(ud[[d]]) + t(ed[[d]])                      
   	}
   
   		
   # Ajuste mediante REML. Obtencion de parametros (funcion REML.compo emn el fichero REML.R)
   	
   	
   	fit <- try(REMLcompo(X, y, D, Ved, Vud), TRUE)
   	
   	
   	if(class(fit)=="try-error"){
   		write.table(data.frame(class(fit),Dominios[g],D,k), file="WARNING indep.txt", append=TRUE, col.names=FALSE)
   		cat("\t Muestra", k, " rechazada por try-error\n")
   		k <- k-1
   		BadTot <- BadTot+1
   		}
   		else{
   		cat("Dominio", D, "Muestra" , k,"\n")
   		
            		
            		
            	  if(fit[[3]]<100) {  	
            		theta.gorro[k] <-fit[1]
            	
            		theta.hat<-theta.gorro[[k]]  # Valores de los paramtros obtenidos mediante REML
            		
      # Calculo de los valores estimados de los beta  (funcion BBETA.U.compo en el fichero Estimacion BETA.R)   
            			
        beta.gorro[[k]] <- BETA.U.compo(X, y, D, Ved ,theta.hat)
            }
            else {
                k <- k-1
                BadTot <- BadTot+fit[[4]]
            }  
       
         
   		}
   
        cat("\n")
    } # Fin bucle simulaciones
        
    Badg[g] <- BadTot
    

    # Calculo de medidas absolutas y relativas para los beta
   
    for(k in 1:K){
     for(l in 1:p){
		 	           betas.gorro[l,k]<-beta.gorro[[k]][[1]][l]
            	  }}
    
    
    
    
    for(l in 1:p){
    BIAS.betaest[l] <- mean(betas.gorro[l,]-beta[l])
    RBIAS.betaest[l] <- 100*BIAS.betaest[l]/abs(beta[l])
    MSE.betaest[l] <- sqrt(mean((betas.gorro[l,]-beta[l])^2))
    RMSE.betaest[l] <- 100*sqrt(MSE.betaest[l])/abs(beta[l])
    }
    
    # Calculo de medidas absolutas y relativas para los theta
    
  theta1.gorro<-theta2.gorro<-theta3.gorro<-theta4.gorro<-theta5.gorro<-theta6.gorro<-vector()
    	
   for (j in 1:K){
    	theta1.gorro[j]<-theta.gorro[[j]][1]
    	theta2.gorro[j]<-theta.gorro[[j]][2]
    	theta3.gorro[j]<-theta.gorro[[j]][3]
    	theta4.gorro[j]<-theta.gorro[[j]][4]
    	theta5.gorro[j]<-theta.gorro[[j]][5]
    	theta6.gorro[j]<-theta.gorro[[j]][6]
    }
    
    
    BIAS.theta1[g] <- mean(theta1.gorro-thetas[1])
    BIAS.theta2[g] <- mean(theta2.gorro-thetas[2])
    BIAS.theta3[g] <- mean(theta3.gorro-thetas[3])
    BIAS.theta4[g] <- mean(theta4.gorro-thetas[4])
    BIAS.theta5[g] <- mean(theta5.gorro-thetas[5])
    BIAS.theta6[g] <- mean(theta6.gorro-thetas[6])
   
    MSE.theta1[g] <- sqrt(mean((theta1.gorro-thetas[1])^2))
    MSE.theta2[g] <- sqrt(mean((theta2.gorro-thetas[2])^2))
    MSE.theta3[g] <- sqrt(mean((theta3.gorro-thetas[3])^2))
    MSE.theta4[g] <- sqrt(mean((theta4.gorro-thetas[4])^2))
    MSE.theta5[g] <- sqrt(mean((theta5.gorro-thetas[5])^2))
    MSE.theta6[g] <- sqrt(mean((theta6.gorro-thetas[6])^2))
    
    RMSE.theta1[g] <- 100*sqrt(MSE.theta1[g])/abs(thetas[1])
    RMSE.theta2[g] <- 100*sqrt(MSE.theta2[g])/abs(thetas[2])
    RMSE.theta3[g] <- 100*sqrt(MSE.theta3[g])/abs(thetas[3])
    RMSE.theta4[g] <- 100*sqrt(MSE.theta4[g])/abs(thetas[4])
    RMSE.theta5[g] <- 100*sqrt(MSE.theta5[g])/abs(thetas[5])
    RMSE.theta6[g] <- 100*sqrt(MSE.theta6[g])/abs(thetas[6])
    
    RBIAS.theta1[g] <- 100*BIAS.theta1[g]/abs(thetas[1])
    RBIAS.theta2[g] <- 100*BIAS.theta2[g]/abs(thetas[2])
    RBIAS.theta3[g] <- 100*BIAS.theta3[g]/abs(thetas[3])
    RBIAS.theta4[g] <- 100*BIAS.theta4[g]/abs(thetas[4])
    RBIAS.theta5[g] <- 100*BIAS.theta5[g]/abs(thetas[5])
    RBIAS.theta6[g] <- 100*BIAS.theta6[g]/abs(thetas[6])
    
# Escritura de resultados en ficheros para los distintos dominios (D (g))

if(g==1) {
resultadosthetas <- data.frame(BIAS.theta1, MSE.theta1,BIAS.theta2, MSE.theta2,BIAS.theta3, MSE.theta3, BIAS.theta4, MSE.theta4,BIAS.theta5, MSE.theta5,BIAS.theta6, MSE.theta6)
resultadosbetas <- data.frame(BIAS.betaest, MSE.betaest)
resultadosthetasR <- data.frame(RBIAS.theta1, RMSE.theta1,RBIAS.theta2, RMSE.theta2,RBIAS.theta3, RMSE.theta3, RBIAS.theta4, RMSE.theta4,RBIAS.theta5, RMSE.theta5,RBIAS.theta6, RMSE.theta6)
resultadosbetasR <- data.frame(RBIAS.betaest, RMSE.betaest)
        resultadosthetas <- as.data.frame(t(resultadosthetas))
        resultadosbetas <- as.data.frame(t(resultadosbetas))
        resultadosthetasR <- as.data.frame(t(resultadosthetasR))
        resultadosbetasR <- as.data.frame(t(resultadosbetasR))
        write.table(resultadosthetas, file="Resultados_Simulacion1_Compo_50thetas.txt", sep="\t")
        write.table(resultadosbetas, file="Resultados_Simulacion1_Compo_50betas.txt", sep="\t")
				write.table(resultadosthetasR, file="Resultados_Simulacion1_Compo_50thetas_R.txt", sep="\t")
	      write.table(resultadosbetasR, file="Resultados_Simulacion1_Compo_50betas_R.txt", sep="\t")
    }
    if(g==2) {
resultadosthetas <- data.frame(BIAS.theta1, MSE.theta1,BIAS.theta2, MSE.theta2,BIAS.theta3, MSE.theta3, BIAS.theta4, MSE.theta4,BIAS.theta5, MSE.theta5,BIAS.theta6, MSE.theta6)
resultadosbetas <- data.frame(BIAS.betaest, MSE.betaest)
resultadosthetasR <- data.frame(RBIAS.theta1, RMSE.theta1,RBIAS.theta2, RMSE.theta2,RBIAS.theta3, RMSE.theta3, RBIAS.theta4, RMSE.theta4,RBIAS.theta5, RMSE.theta5,RBIAS.theta6, RMSE.theta6)
resultadosbetasR <- data.frame(RBIAS.betaest, RMSE.betaest)
        resultadosthetas <- as.data.frame(t(resultadosthetas))
        resultadosbetas <- as.data.frame(t(resultadosbetas))
        resultadosthetasR <- as.data.frame(t(resultadosthetasR))
        resultadosbetasR <- as.data.frame(t(resultadosbetasR))
        write.table(resultadosthetas, file="Resultados_Simulacion1_Compo_75thetas.txt", sep="\t")
        write.table(resultadosbetas, file="Resultados_Simulacion1_Compo_75betas.txt", sep="\t")
				write.table(resultadosthetasR, file="Resultados_Simulacion1_Compo_75thetas_R.txt", sep="\t")
	      write.table(resultadosbetasR, file="Resultados_Simulacion1_Compo_75betas_R.txt", sep="\t")
    }
    if(g==3) {
resultadosthetas <- data.frame(BIAS.theta1, MSE.theta1,BIAS.theta2, MSE.theta2,BIAS.theta3, MSE.theta3, BIAS.theta4, MSE.theta4,BIAS.theta5, MSE.theta5,BIAS.theta6, MSE.theta6)
resultadosbetas <- data.frame(BIAS.betaest, MSE.betaest)
resultadosthetasR <- data.frame(RBIAS.theta1, RMSE.theta1,RBIAS.theta2, RMSE.theta2,RBIAS.theta3, RMSE.theta3, RBIAS.theta4, RMSE.theta4,RBIAS.theta5, RMSE.theta5,RBIAS.theta6, RMSE.theta6)
resultadosbetasR <- data.frame(RBIAS.betaest, RMSE.betaest)
        resultadosthetas <- as.data.frame(t(resultadosthetas))
        resultadosbetas <- as.data.frame(t(resultadosbetas))
        resultadosthetasR <- as.data.frame(t(resultadosthetasR))
        resultadosbetasR <- as.data.frame(t(resultadosbetasR))
        write.table(resultadosthetas, file="Resultados_Simulacion1_Compo_100thetas.txt", sep="\t")
        write.table(resultadosbetas, file="Resultados_Simulacion1_Compo_100betas.txt", sep="\t")
				write.table(resultadosthetasR, file="Resultados_Simulacion1_Compo_100thetas_R.txt", sep="\t")
	      write.table(resultadosbetasR, file="Resultados_Simulacion1_Compo_100betas_R.txt", sep="\t")
    }
    
    if(g==4) {
 resultadosthetas <- data.frame(BIAS.theta1, MSE.theta1,BIAS.theta2, MSE.theta2,BIAS.theta3, MSE.theta3, BIAS.theta4, MSE.theta4,BIAS.theta5, MSE.theta5,BIAS.theta6, MSE.theta6)
resultadosbetas <- data.frame(BIAS.betaest, MSE.betaest)
resultadosthetasR <- data.frame(RBIAS.theta1, RMSE.theta1,RBIAS.theta2, RMSE.theta2,RBIAS.theta3, RMSE.theta3, RBIAS.theta4, RMSE.theta4,RBIAS.theta5, RMSE.theta5,RBIAS.theta6, RMSE.theta6)
resultadosbetasR <- data.frame(RBIAS.betaest, RMSE.betaest)
        resultadosthetas <- as.data.frame(t(resultadosthetas))
        resultadosbetas <- as.data.frame(t(resultadosbetas))
        resultadosthetasR <- as.data.frame(t(resultadosthetasR))
        resultadosbetasR <- as.data.frame(t(resultadosbetasR))
        write.table(resultadosthetas, file="Resultados_Simulacion1_Compo_150thetas.txt", sep="\t")
        write.table(resultadosbetas, file="Resultados_Simulacion1_Compo_150betas.txt", sep="\t")
				write.table(resultadosthetasR, file="Resultados_Simulacion1_Compo_150thetas_R.txt", sep="\t")
	      write.table(resultadosbetasR, file="Resultados_Simulacion1_Compo_150betas_R.txt", sep="\t")
    }
    if(g==5) {
resultadosthetas <- data.frame(BIAS.theta1, MSE.theta1,BIAS.theta2, MSE.theta2,BIAS.theta3, MSE.theta3, BIAS.theta4, MSE.theta4,BIAS.theta5, MSE.theta5,BIAS.theta6, MSE.theta6)
resultadosbetas <- data.frame(BIAS.betaest, MSE.betaest)
resultadosthetasR <- data.frame(RBIAS.theta1, RMSE.theta1,RBIAS.theta2, RMSE.theta2,RBIAS.theta3, RMSE.theta3, RBIAS.theta4, RMSE.theta4,RBIAS.theta5, RMSE.theta5,RBIAS.theta6, RMSE.theta6)
resultadosbetasR <- data.frame(RBIAS.betaest, RMSE.betaest)
        resultadosthetas <- as.data.frame(t(resultadosthetas))
        resultadosbetas <- as.data.frame(t(resultadosbetas))
        resultadosthetasR <- as.data.frame(t(resultadosthetasR))
        resultadosbetasR <- as.data.frame(t(resultadosbetasR))
        write.table(resultadosthetas, file="Resultados_Simulacion1_Compo_200thetas.txt", sep="\t")
        write.table(resultadosbetas, file="Resultados_Simulacion1_Compo_200betas.txt", sep="\t")
				write.table(resultadosthetasR, file="Resultados_Simulacion1_Compo_200thetas_R.txt", sep="\t")
	      write.table(resultadosbetasR, file="Resultados_Simulacion1_Compo_200betas_R.txt", sep="\t")
    }
    if(g==6) {
resultadosthetas <- data.frame(BIAS.theta1, MSE.theta1,BIAS.theta2, MSE.theta2,BIAS.theta3, MSE.theta3, BIAS.theta4, MSE.theta4,BIAS.theta5, MSE.theta5,BIAS.theta6, MSE.theta6)
resultadosbetas <- data.frame(BIAS.betaest, MSE.betaest)
resultadosthetasR <- data.frame(RBIAS.theta1, RMSE.theta1,RBIAS.theta2, RMSE.theta2,RBIAS.theta3, RMSE.theta3, RBIAS.theta4, RMSE.theta4,RBIAS.theta5, RMSE.theta5,RBIAS.theta6, RMSE.theta6)
resultadosbetasR <- data.frame(RBIAS.betaest, RMSE.betaest)
        resultadosthetas <- as.data.frame(t(resultadosthetas))
        resultadosbetas <- as.data.frame(t(resultadosbetas))
        resultadosthetasR <- as.data.frame(t(resultadosthetasR))
        resultadosbetasR <- as.data.frame(t(resultadosbetasR))
        write.table(resultadosthetas, file="Resultados_Simulacion1_Compo_300thetas.txt", sep="\t")
        write.table(resultadosbetas, file="Resultados_Simulacion1_Compo_300betas.txt", sep="\t")
				write.table(resultadosthetasR, file="Resultados_Simulacion1_Compo_300thetas_R.txt", sep="\t")
	      write.table(resultadosbetasR, file="Resultados_Simulacion1_Compo_300betas_R.txt", sep="\t")
    }
    
mal <- data.frame(Badg)
write.table(mal, file="Fallos.txt", sep="\t")
}



 