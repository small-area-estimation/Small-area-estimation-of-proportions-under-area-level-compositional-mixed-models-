sink("Simul4.txt")
###############################################################################
###############################################################################
###
###        Experimento de simulacion 4a  Compositional Data (parametrizacion con rho)               
###                                                            
###                            Septiembre 2017                            
###         Estudiamos el efecto de distintas multinomiales


### Autor:  MJ Lombardia    

#getwd()
#setwd(/Users/mjose/Documents/Files/FilesI/Investigacion/GrupoElche/Datos composicionales/R)


# Llamada a las funciones Generacion X, REMLcompo y BETA.U.compo
source("Generacion X.R") # Funcion para generar los valores de las variables explicativas
source("REML.R") # Funcion ajuste REML
source("Estimacion BETA.R") # Calculo de beta´s y u´s estimados
source("Transformaciones Vu_Theta y Theta_Vu.R")  
source("Calculo CEB.R")  


library(mme)
library(Matrix)
library(MASS)
library(mvtnorm)

set.seed(19032005)
I <- 100# Numero de vueltas para las simulaciones
Dominios <-100 # Tamaños de los dominios D
#Dominios <-c(50,75,100) # Tamaños de los dominios D
#Dominios <-c(50,75,100,150,200,300) # Tamaños de los dominios D


# Inicializacion de variables 


udk<-list()
xd1<-xd2<-xd3<-X<-list()

H0<-Ved<-list()


# Inicio del bucle de dominios

for(g in 1:length(Dominios)) {
  D <- Dominios[g]
  
  cat("\n Comienza la simulacion", g, " con D =", D,  "\n")
  
  
  # Simulacion de las variables explicativas (funcion equis del fichero Generacion X.R)

  ##Opcion 1: Construimos las x's reales
  Xp<-equis(D)   
  
  for(d in 1:D){
    X[[d]]<-Xp[[1]][[d]]
  }   
  
  p<-Xp[[2]]
  
 
  
  # Valores iniciales parametros
  
  # Beta, thetai,Vud, Ved    
  
  betaM<-as.vector(c(1,-0.12,1,-0.13,1,-0.14)) #se fijan de acuerdo a la Sec.7 del documento Acompo
  theta<-rep(0.1,3) #sigma_u
  mediasu<-rep(0,3)
  nud<-100 #Sec7.Line 7
  q<-4 #numero de categorias
  
  
  # Inicializacion variables (parametros) para el ajuste mme

  ud<-ed<-etad<-array(rep(0,D*3),dim=c(D,3))  
  muMd<-array(rep(0,3*D*I),dim=c(3,D,I))
  zd.bar<-zd<-y<-Y<-list()
  pMi<-pMi.est<-array(rep(0,3*D*I),dim=c(3,D,I))
  covMd<-array(rep(0,3*3),dim=c(3,3))
  Vud<-Ved<-list()  
 
 ## Inicializamos el vector de errores y sesgos
 pMdk<-pMdk.est<-array(rep(0,D*3),dim=c(D,3))
 BpMdk<-EpMdk<-RBpMdk<-REpMdk<-array(rep(0,D*3),dim=c(D,3))
 
 # Inicializacion variables (parametros) para el ajuste TFH
 theta.gorro<-vector()
 beta.gorro<-betas.gorro<-list()
 pebdk<-pidk.gorro<- mudi.gorro.comp<-array(rep(0,3*D*I),dim=c(3,D,I))
 pdk.gorro<-array(rep(0,D*3),dim=c(D,3))
 ## Inicializamos el vector de errores y sesgos
 Epidk<-Bpidk<-Epebdk<-Bpebdk<-array(rep(0,D*3),dim=c(D,3))
 
 ## Calculo de varianzas para este ajuste TFH

 #Calculo de H0 para la varianza Ved del ajuste TFH
 k<-q-1
 unos<-as.vector(rep(1,k))
 H0<-array(rep(0,k*k),dim=c(k,k))
 diag(H0)<-1
 H0<-q*(H0+unos%*%t(unos))
 
 
 # Matriz Vud
 Vud<-array(rep(0,(q-1)*(q-1)),dim=c((q-1),(q-1))) #Vud para el modelo TFH
# diag(Vud)<-1 
 
 
  # Inicio bucle simulaciones (k)
i<-BadTot <-0  
  while(i<I){
    i<-i+1
    
    cat("Comienza la muestra", i, "\n")
    
    # Simulacion de los ud (normales univariantes), zd, pMd para el modelo MLM (k=1,2,3)
  
    for(d in 1:D) {

      y[[d]]<-zd.bar[[d]]<-zd[[d]]<-Ved[[d]]<-vector()
      
      ud[d,]<-rnorm(3,mediasu,theta)
      etad[d,]<-X[[d]]%*%betaM+ ud[d,]
      pMi[,d,i] <-exp(etad[d,])/(1+sum(exp(etad[d,]))) #Probabilidades de la multinomial teoricas
      prob<-c(pMi[1,d,i],pMi[2,d,i],pMi[3,d,i],1-pMi[1,d,i]-pMi[2,d,i]-pMi[3,d,i])
      zd[[d]]<-rmultinom(1,nud,prob) 
      
      zd.bar[[d]]<-zd[[d]]/nud
      y[[d]]<-log(zd.bar[[d]]/zd.bar[[d]][4]) #Revisar calculo de los datos en logaritmo:aparecen 0's e INF?????
      muMd[,d,i] <- nud*pMi[,d,i] #calculo de las verdaderas medias
      covMd<-nud*(diag(pMi[,d,i])-pMi[,d,i]%*%t(pMi[,d,i])) #calculo de las cov. para Ved
      Ved[[d]]<-nud^(-2)*H0%*%covMd%*%t(H0) #Ved para el modelo TFH
     
      if(zd[[d]][1,1]==0||zd[[d]][2,1]==0||zd[[d]][3,1]==0||zd[[d]][4,1]==0){
        #        zd[[d]]<-rmultinom(1,nud,prob)
        cat("\t iter y area:",i, d, " rechazada por ceros")
        d<-d-1
        BadTot <- BadTot+1
        next
      }
      
      
      
      } #Fin del bucle D:recorre las areas


# Utilizacion del paquete de R mme para el ajuste    
   #Needed matrix and initial values
    Y[[d]]<-X1<-X2<-X3<-Z1<-Z2<-Z3<-Z4<-Area<-Time<-Sample<-Population<-vector()
     szdk<-array(rep(0,3*D*I),dim=c(3,D,I))
#    De1<-De2<-De3<-vector() #Opcion1 para calculo Vu
     Y1<-Y2<-Y3<-vector()    #Opcion2 para calculo Vu
     
     
     
   for(d in 1:D)	{
     Area[d]=d
     Time[d]=1
     Sample[d]=nud 
     Population[d]=1000 ## no se utiliza en 'mme' a menos que se haga un estudio poblacional
     
     X1[d]<-X[[d]][1,2] #se construye el fichero de datos desde las X's y los conteos de la Multinomial
     X2[d]<-X[[d]][2,4]
     X3[d]<-X[[d]][3,6]
     Z1[d]<-zd[[d]][1,1]
     Z2[d]<-zd[[d]][2,1]
     Z3[d]<-zd[[d]][3,1]
     Z4[d]<-zd[[d]][4,1]
     
     Y[[d]]<-cbind(y[[d]][-4]) #vector de dimension 3 necesario para el ajuste TFH
     
     ##Opcion 1 para el calculo de Vu: Vu=dig(Mediana Ved1, Ved2,Ved3), ajuste TFH
#     De1[d]<-Ved[[d]][1,1] #Guardamos los elementos de la diagonal de Ved 
#     De2[d]<-Ved[[d]][2,2] #construimos la diagonal de Vu a partir de  
#     De3[d]<-Ved[[d]][3,3] #la mediana de los valores diagonales de Ved
     
     ##Opcion 2 para el calculo de Vu: Vu=dig(var(Y1), var(Y2),var(Y3)), ajuste TFH
     Y1[d]<-Y[[d]][1,1] #Guardamos los elementos de la diagonal de Ved 
     Y2[d]<-Y[[d]][2,1] #construimos la diagonal de Vu a partir de  
     Y3[d]<-Y[[d]][3,1] #la mediana de los valores diagonales de Ved
     
   }
                                    
    resmme <- data.frame(Area,Time,Sample, Population,Z1,Z2,Z3,Z4,X1,X2,X3 ) 
    
#    q=4 # number of categories of the response variable (definido arriba)
    pp=c(1,1,1) # vector with the number of auxiliary variables in each category
    mod=1 # Model 1
    # Needed matrix and initial values
    datar=data.mme(resmme,q,pp,mod)
    # Model fit    ###PUEDE FALLAR EN ALGUNAS ITERACIONES###
    result=model(datar$d,datar$t,pp,datar$Xk,datar$X,datar$Z,datar$initial,
                 + datar$y[,1:(q-1)],datar$n,datar$N,mod)
    
    #result
    
    # Predictors of sample totals
    sZ1=result$mean[,1]*datar$n/datar$N 
    sZ2=result$mean[,2]*datar$n/datar$N
    sZ3=result$mean[,3]*datar$n/datar$N
    resultadomme<- data.frame(Z1,Z2,Z3,sZ1,sZ2,sZ3)
    # Predictors of domain*category probabilities
    p1<-result$mean[,1]/datar$N #que calcula????? son las pMd.est??? pero no se parecen...
    p2<-result$mean[,2]/datar$N
    p3<-result$mean[,3]/datar$N

    for(d in 1:D){
      szdk[1,d,i]<-sZ1[d]
      szdk[2,d,i]<-sZ2[d]
      szdk[3,d,i]<-sZ3[d]
      pMi.est[1,d,i]<-p1[d]
      pMi.est[2,d,i]<-p2[d]
      pMi.est[3,d,i]<-p3[d]
    }	
    
    # AJUSTE DEL MODELO DE AREA MEDIANTE TFH   
#    diag(Vud)<-c(median(De1),median(De2),median(De3))  ## Opcion 1 para calcular la Vu
    diag(Vud)<-c(var(Y1)/2,var(Y2)/2,var(Y3)/2)     ## Opcion 2 para calcular la Vu
    # Inicializamos Vud antes del bucle I   

   # AJUSTE DEL MODELO DE AREA MEDIANTE REML (funcion REML.compo definida en REML.R)

      fit <- try(REMLcompo(X, Y, D, Ved, Vud), TRUE) 
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
          
          beta.gorro[[i]] <- BETA.U.compo(X, Y, D, Ved ,theta.hat)
          
          
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
          #### FUNCION CEB
          
#          pebdk[,,i]<-t(ceb(D,X,Y,theta.hat,betas.gorro,Ved))
          
        }# fin bucle if(fit[[3]])
      }# fin bucle else       

    
    }#Fin del bucle I:Iteraciones del Monte Carlo

# Calculo de errores y sesgos absolutos y relativos

for(d in 1:D){ 
  for(k in 1:3){
    pMdk[d,k]<-mean(pMi[k,d,],na.rm=TRUE) #promedio de probabilidades multinomiales
    pMdk.est[d,k]<-mean(pMi.est[k,d,],na.rm=TRUE) #promedio de probabilidades estimadas multinomiales
    pdk.gorro[d,k]<-mean(pidk.gorro[k,d,],na.rm=TRUE)#promedio de probabilidades estimadas TFH
    EpMdk[d,k]<-(mean((pMi.est[k,d,]-pMi[k,d,])^2,na.rm=TRUE))^0.5
    BpMdk[d,k]<-mean(pMi.est[k,d,]-pMi[k,d,],na.rm=TRUE)
    Epidk[d,k]<-(mean((pidk.gorro[k,d,]-pMi[k,d,])^2,na.rm=TRUE))^0.5
    Bpidk[d,k]<-mean(pidk.gorro[k,d,]-pMi[k,d,],na.rm=TRUE)
#    Epebdk[d,k]<-(mean((pebdk[k,d,]-pMi[k,d,])^2,na.rm=TRUE))^0.5
#    Bpebdk[d,k]<-mean(pebdk[k,d,]-pMi[k,d,],na.rm=TRUE)
    
  }
}
REpMdk<-10^2*(EpMdk/pMdk)
RBpMdk<-10^2*(BpMdk/pMdk)
REpidk<-10^2*(Epidk/pMdk)
RBpidk<-10^2*(Bpidk/pMdk)
#REpebdk<-10^2*(Epebdk/pMdk)
#RBpebdk<-10^2*(Bpebdk/pMdk)

EpMk<-REpMk<-BpMk<-RBpMk<-Epik<-REpik<-Bpik<-RBpik<-Epebk<-REpebk<-Bpebk<-RBpebk<-vector()

for(k in 1:3){
   EpMk[k]<-mean(EpMdk[,k],na.rm=TRUE)
   REpMk[k]<-mean(REpMdk[,k],na.rm=TRUE)
   BpMk[k]<-mean(abs(BpMdk[,k]),na.rm=TRUE)
   RBpMk[k]<-mean(abs(RBpMdk[,k]),na.rm=TRUE)
   Epik[k]<-mean(Epidk[,k],na.rm=TRUE)
   REpik[k]<-mean(REpidk[,k],na.rm=TRUE)
   Bpik[k]<-mean(abs(Bpidk[,k]),na.rm=TRUE)
   RBpik[k]<-mean(abs(RBpidk[,k]),na.rm=TRUE)
#   Epebk[k]<-mean(Epebdk[,k],na.rm=TRUE)
#   REpebk[k]<-mean(REpebdk[,k],na.rm=TRUE)
#   Bpebk[k]<-mean(abs(Bpebdk[,k]),na.rm=TRUE)
#   RBpebk[k]<-mean(abs(RBpebdk[,k]),na.rm=TRUE)
}

# Escritutra de resultados en ficheros

if(g==1) {
  write.table(pMdk, file="pMdk_D100.txt", row.names=FALSE, col.names=FALSE)
  resultadosP <- data.frame(pMdk,pMdk.est,pdk.gorro)
  resultadosP <- as.data.frame(resultadosP)
  write.table(resultadosP, file="Probabilidadesdk_D100.txt", row.names=F, col.names=F)
  resultadosE <- data.frame(EpMdk,Epidk)
  resultadosE <- as.data.frame(resultadosE)
  write.table(resultadosE, file="Erroresdk_D100.txt", row.names=F, col.names=F)
  resultadosB <- data.frame(BpMdk,Bpidk)
  resultadosB <- as.data.frame((resultadosB))
  write.table(resultadosB, file="Biasdk_D100.txt", row.names=FALSE, col.names=FALSE)
  resultadosRE <- data.frame(REpMdk,REpidk)
  resultadosRE <- as.data.frame((resultadosRE))
  write.table(resultadosRE, file="RelativeError_100.txt", row.names=FALSE, col.names=FALSE)
  resultadosRB <- data.frame(RBpMdk,RBpidk)
  resultadosRB <- as.data.frame((resultadosRB))
  write.table(resultadosRB, file="RelativeBias_D100.txt", row.names=FALSE, col.names=FALSE)
  
#  write.table(EpMk, file="EpMk_D50.txt", row.names=FALSE, col.names=FALSE)
#  write.table(REpMk, file="REpMk_D50.txt", row.names=FALSE, col.names=FALSE)
#  write.table(BpMk, file="BpMk_D50.txt", row.names=FALSE, col.names=FALSE)
#  write.table(RBpMk, file="RBpMk_D50.txt", row.names=FALSE, col.names=FALSE)
  resultados <- data.frame(EpMk,Epik,REpMk,REpik,BpMk,Bpik, RBpMk,RBpik)
#  resultados <- data.frame(EpMk,Epik,Epebk,REpMk,REpik,REpebk,BpMk,Bpik,Bpebk, RBpMk,RBpik,RBpebk)
  resultados <- as.data.frame(t(resultados))
  write.table(resultados, file="Resultados_Simulacion4_Compo_D100.txt", sep="\t")
  
}
if(g==2) {
  write.table(pMdk, file="pMdk_D75.txt", row.names=FALSE, col.names=FALSE)
  resultadosP <- data.frame(pMdk,pMdk.est,pdk.gorro)
  resultadosP <- as.data.frame(resultadosP)
  write.table(resultadosP, file="Probabilidadesdk_D75.txt", row.names=F, col.names=F)
  #  write.table(resultadosP, file="pest_D50.txt", row.names=FALSE, col.names=FALSE)
  resultadosE <- data.frame(EpMdk,Epidk)
  resultadosE <- as.data.frame(resultadosE)
  write.table(resultadosE, file="Erroresdk_D75.txt", row.names=F, col.names=F)
  resultadosB <- data.frame(BpMdk,Bpidk)
  resultadosB <- as.data.frame((resultadosB))
  write.table(resultadosB, file="Biasdk_D75.txt", row.names=FALSE, col.names=FALSE)
  resultadosRE <- data.frame(REpMdk,REpidk)
  resultadosRE <- as.data.frame((resultadosRE))
  write.table(resultadosRE, file="RelativeError_75.txt", row.names=FALSE, col.names=FALSE)
  resultadosRB <- data.frame(RBpMdk,RBpidk)
  resultadosRB <- as.data.frame((resultadosRB))
  write.table(resultadosRB, file="RelativeBias_D75.txt", row.names=FALSE, col.names=FALSE)
  
  #  write.table(EpMk, file="EpMk_D75.txt", row.names=FALSE, col.names=FALSE)
  #  write.table(REpMk, file="REpMk_D75.txt", row.names=FALSE, col.names=FALSE)
  #  write.table(BpMk, file="BpMk_D75.txt", row.names=FALSE, col.names=FALSE)
  #  write.table(RBpMk, file="RBpMk_D75.txt", row.names=FALSE, col.names=FALSE)
  resultados <- data.frame(EpMk,Epik,REpMk,REpik,BpMk,Bpik, RBpMk,RBpik)
  #  resultados <- data.frame(EpMk,Epik,Epebk,REpMk,REpik,REpebk,BpMk,Bpik,Bpebk, RBpMk,RBpik,RBpebk)
  resultados <- as.data.frame(t(resultados))
  write.table(resultados, file="Resultados_Simulacion4_Compo_D75.txt", sep="\t")
  
}
if(g==3) {
  write.table(pMdk, file="pMdk_D100.txt", row.names=FALSE, col.names=FALSE)
  resultadosP <- data.frame(pMdk,pMdk.est,pdk.gorro)
  resultadosP <- as.data.frame(resultadosP)
  write.table(resultadosP, file="Probabilidadesdk_D100.txt", row.names=F, col.names=F)
  resultadosE <- data.frame(EpMdk,Epidk)
  resultadosE <- as.data.frame(resultadosE)
  write.table(resultadosE, file="Erroresdk_D100.txt", row.names=F, col.names=F)
  resultadosB <- data.frame(BpMdk,Bpidk)
  resultadosB <- as.data.frame((resultadosB))
  write.table(resultadosB, file="Biasdk_D100.txt", row.names=FALSE, col.names=FALSE)
  resultadosRE <- data.frame(REpMdk,REpidk)
  resultadosRE <- as.data.frame((resultadosRE))
  write.table(resultadosRE, file="RelativeError_100.txt", row.names=FALSE, col.names=FALSE)
  resultadosRB <- data.frame(RBpMdk,RBpidk)
  resultadosRB <- as.data.frame((resultadosRB))
  write.table(resultadosRB, file="RelativeBias_D100.txt", row.names=FALSE, col.names=FALSE)
  
  #  write.table(EpMk, file="EpMk_D100.txt", row.names=FALSE, col.names=FALSE)
  #  write.table(REpMk, file="REpMk_D100.txt", row.names=FALSE, col.names=FALSE)
  #  write.table(BpMk, file="BpMk_D100.txt", row.names=FALSE, col.names=FALSE)
  #  write.table(RBpMk, file="RBpMk_D100.txt", row.names=FALSE, col.names=FALSE)
  resultados <- data.frame(EpMk,Epik,REpMk,REpik,BpMk,Bpik, RBpMk,RBpik)
  #  resultados <- data.frame(EpMk,Epik,Epebk,REpMk,REpik,REpebk,BpMk,Bpik,Bpebk, RBpMk,RBpik,RBpebk)
  resultados <- as.data.frame(t(resultados))
  write.table(resultados, file="Resultados_Simulacion4_Compo_D100.txt", sep="\t")
  
}

if(g==4) {
  write.table(pMdk, file="pMdk_D150.txt", row.names=FALSE, col.names=FALSE)
  write.table(EpMdk, file="EpMdk_D150.txt", row.names=FALSE, col.names=FALSE)
  write.table(BpMdk, file="BpMdk_D150.txt", row.names=FALSE, col.names=FALSE)
  write.table(REpMdk, file="REpMdk_D150.txt", row.names=FALSE, col.names=FALSE)
  write.table(RBpMdk, file="RBpMdk_D150.txt", row.names=FALSE, col.names=FALSE)
  write.table(EpMk, file="EpMk_D150.txt", row.names=FALSE, col.names=FALSE)
  write.table(REpMk, file="REpMk_D150.txt", row.names=FALSE, col.names=FALSE)
  write.table(BpMk, file="BpMk_D150.txt", row.names=FALSE, col.names=FALSE)
  write.table(RBpMk, file="RBpMk_D100.txt", row.names=FALSE, col.names=FALSE)
  resultados <- data.frame(EpMk,REpMk,BpMk,RBpMk)
  resultados <- as.data.frame(t(resultados))
  write.table(resultados, file="Resultados_Simulacion4_Compo_D150.txt", sep="\t")
  
}
if(g==5) {
  write.table(pMdk, file="pMdk_D200.txt", row.names=FALSE, col.names=FALSE)
  write.table(EpMdk, file="EpMdk_D200.txt", row.names=FALSE, col.names=FALSE)
  write.table(BpMdk, file="BpMdk_D200.txt", row.names=FALSE, col.names=FALSE)
  write.table(REpMdk, file="REpMdk_D200.txt", row.names=FALSE, col.names=FALSE)
  write.table(RBpMdk, file="RBpMdk_D200.txt", row.names=FALSE, col.names=FALSE)
  write.table(EpMk, file="EpMk_D200.txt", row.names=FALSE, col.names=FALSE)
  write.table(REpMk, file="REpMk_D200.txt", row.names=FALSE, col.names=FALSE)
  write.table(BpMk, file="BpMk_D200.txt", row.names=FALSE, col.names=FALSE)
  write.table(RBpMk, file="RBpMk_D200.txt", row.names=FALSE, col.names=FALSE)
  resultados <- data.frame(EpMk,REpMk,BpMk,RBpMk)
  resultados <- as.data.frame(t(resultados))
  write.table(resultados, file="Resultados_Simulacion4_Compo_D200.txt", sep="\t")
}
if(g==6) {
  write.table(pMdk, file="pMdk_D300.txt", row.names=FALSE, col.names=FALSE)
  write.table(EpMdk, file="EpMdk_D300.txt", row.names=FALSE, col.names=FALSE)
  write.table(BpMdk, file="BpMdk_D300.txt", row.names=FALSE, col.names=FALSE)
  write.table(REpMdk, file="REpMdk_D300.txt", row.names=FALSE, col.names=FALSE)
  write.table(RBpMdk, file="RBpMdk_D300.txt", row.names=FALSE, col.names=FALSE)
  write.table(EpMk, file="EpMk_D300.txt", row.names=FALSE, col.names=FALSE)
  write.table(REpMk, file="REpMk_D300.txt", row.names=FALSE, col.names=FALSE)
  write.table(BpMk, file="BpMk_D300.txt", row.names=FALSE, col.names=FALSE)
  write.table(RBpMk, file="RBpMk_D300.txt", row.names=FALSE, col.names=FALSE)
  resultados <- data.frame(EpMk,REpMk,BpMk,RBpMk)
  resultados <- as.data.frame(t(resultados))
  write.table(resultados, file="Resultados_Simulacion4_Compo_D300.txt", sep="\t")
}




} #Fin del bucle G:dominios distintos







  
  