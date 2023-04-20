###############################################################################
###############################################################################
###
###                 Algoritmo Fisher Scoring REML Compositional data           
###              
###                               febrero 2017                                
###

### AUTOR:  Lola Esteban Lefler
### Modelo:  Ydt = Xd Beta + Ud + ed


### Editado por:
### con fecha: 

REMLcompo <- function(X, y, D, Ved, Vud, MAXITER = 100) {
    
	
    Vud.f  <- Vud
   
    
    Xd<-X
    yd<-y
    
    # Obtencion de los parametros theta a partir de la mariz Vud
    
    theta.f <-Thetas(Vud)
    

    
    # Calculo de las derivadas de la matriz Vud respecto de los parametros (theta)
       
    Vuds<- FirstDer(theta.f)
    
    Bad <- 0
    FLAG <- 0
    
    p <- ncol(X[[1]])
   
    # Inicio algoritmo Fisher Scoring 
    
    for(ITER in 1:MAXITER){
    	
    	#Iinicializacion de las variables necesarias
    	
        Vd.inv <- Vinvyd <- VinvXd <- list()
        
        #1
        VinvVud1 <- XtVinvVud1VinvX <- VinvVud1VinvVud1 <- XtVinvVud1VinvVud1VinvX <- list()
        tr.VinvVud1 <- ytVinvX <- ytVinvVud1Vinvy <- ytVinvVud1VinvX <- SumXtVinvVud1VinvX <- tr.VinvVud1VinvVud1 <- 0
        #2
        VinvVud2 <- XtVinvVud2VinvX <- VinvVud2VinvVud2 <- XtVinvVud2VinvVud2VinvX <- list()
        tr.VinvVud2 <- ytVinvVud2Vinvy <- ytVinvVud2VinvX <- SumXtVinvVud2VinvX <- tr.VinvVud2VinvVud2 <- 0
        #3
        VinvVud3 <- XtVinvVud3VinvX <- VinvVud3VinvVud3 <- XtVinvVud3VinvVud3VinvX <- list()
        tr.VinvVud3 <- ytVinvVud3Vinvy <- ytVinvVud3VinvX <- SumXtVinvVud3VinvX <- tr.VinvVud3VinvVud3 <- 0
        
        #4
        VinvVud4 <- XtVinvVud4VinvX <- VinvVud4VinvVud4 <- XtVinvVud4VinvVud4VinvX <- list()
        tr.VinvVud4 <- ytVinvVud4Vinvy <- ytVinvVud4VinvX <- SumXtVinvVud4VinvX <- tr.VinvVud4VinvVud4 <- 0
        
        #5
        VinvVud5 <- XtVinvVud5VinvX <- VinvVud5VinvVud5 <- XtVinvVud5VinvVud5VinvX <- list()
        tr.VinvVud5 <- ytVinvVud5Vinvy <- ytVinvVud5VinvX <- SumXtVinvVud5VinvX <- tr.VinvVud5VinvVud5 <- 0
        
        #6
        VinvVud6 <- XtVinvVud6VinvX <- VinvVud6VinvVud6 <- XtVinvVud6VinvVud6VinvX <- list()
        tr.VinvVud6 <- ytVinvVud6Vinvy <- ytVinvVud6VinvX <- SumXtVinvVud6VinvX <- tr.VinvVud6VinvVud6 <- 0
        
        # 12
        VinvVud1VinvVud2 <- XtVinvVud1VinvVud2VinvX <- list()
        tr.VinvVud1VinvVud2 <- 0
        
        # 13
        VinvVud1VinvVud3 <- XtVinvVud1VinvVud3VinvX <- list()
        tr.VinvVud1VinvVud3 <- 0
        
        # 14
        VinvVud1VinvVud4 <- XtVinvVud1VinvVud4VinvX <- list()
        tr.VinvVud1VinvVud4 <- 0
        
        # 15
        VinvVud1VinvVud5 <- XtVinvVud1VinvVud5VinvX <- list()
        tr.VinvVud1VinvVud5 <- 0
        
        # 16
        VinvVud1VinvVud6 <- XtVinvVud1VinvVud6VinvX <- list()
        tr.VinvVud1VinvVud6 <- 0
        
        # 23
        VinvVud2VinvVud3 <- XtVinvVud2VinvVud3VinvX <- list()
        tr.VinvVud2VinvVud3 <- 0
        
        # 24
        VinvVud2VinvVud4 <- XtVinvVud2VinvVud4VinvX <- list()
        tr.VinvVud2VinvVud4 <- 0
        
        # 25
        VinvVud2VinvVud5 <- XtVinvVud2VinvVud5VinvX <- list()
        tr.VinvVud2VinvVud5 <- 0
        
        # 26
        VinvVud2VinvVud6 <- XtVinvVud2VinvVud6VinvX <- list()
        tr.VinvVud2VinvVud6 <- 0
        
        # 34
        VinvVud3VinvVud4 <- XtVinvVud3VinvVud4VinvX <- list()
        tr.VinvVud3VinvVud4 <- 0
        
        # 35
        VinvVud3VinvVud5 <- XtVinvVud3VinvVud5VinvX <- list()
        tr.VinvVud3VinvVud5 <- 0
        
        # 36
        VinvVud3VinvVud6 <- XtVinvVud3VinvVud6VinvX <- list()
        tr.VinvVud3VinvVud6 <- 0
        
        # 45
        VinvVud4VinvVud5 <- XtVinvVud4VinvVud5VinvX <- list()
        tr.VinvVud4VinvVud5 <- 0
        
        # 46
        VinvVud4VinvVud6 <- XtVinvVud4VinvVud6VinvX <- list()
        tr.VinvVud4VinvVud6 <- 0
        
        # 56
        VinvVud5VinvVud6 <- XtVinvVud5VinvVud6VinvX <- list()
        tr.VinvVud5VinvVud6 <- 0
       
       
        Q.inv <- matrix(0, nrow=p, ncol=p)
        
          
            
        
        for(d in 1:D) {
    			
          Vd <-Vud.f+Ved[[d]]        # Elementos de la matriz de varianza
        	
            
            if(abs(det(Vd))<0.000000001 || abs(det(Vd))>100000000000) {
                FLAG <- 1
                Bad <- Bad+1
                break
            }
            
            Vd.inv[[d]] <- solve(Vd)                                    # Inversa de la matriz de varianza en d submatrices
            Vinvyd[[d]] <- Vd.inv[[d]]%*%yd[[d]]                        # Producto V^-1_d por y_d para las d submatrices
            VinvXd[[d]] <- Vd.inv[[d]]%*%Xd[[d]]                        # Producto V^-1_d por X_d para las d submatrices
            Q.inv <- Q.inv + t(Xd[[d]])%*%VinvXd[[d]]                   # Inversa de Q. Posteriormente se calcula Q
            
# Calculos necesarios para obtener los scores (s1 a s6)
            
            #S1
            VinvVud1[[d]] <- Vd.inv[[d]]%*%Vuds[[1]]
            tr.VinvVud1 <- tr.VinvVud1 + sum(diag(VinvVud1[[d]]))
            XtVinvVud1VinvX[[d]] <- t(VinvXd[[d]])%*%Vuds[[1]]%*%VinvXd[[d]]
            ytVinvX <- ytVinvX + t(yd[[d]])%*%VinvXd[[d]]
            ytVinvVud1Vinvy <- ytVinvVud1Vinvy + t(Vinvyd[[d]])%*%Vuds[[1]]%*%Vinvyd[[d]]
            ytVinvVud1VinvX <- ytVinvVud1VinvX + t(Vinvyd[[d]])%*%Vuds[[1]]%*%VinvXd[[d]]
            SumXtVinvVud1VinvX <- SumXtVinvVud1VinvX + XtVinvVud1VinvX[[d]]
            
            #S2
            VinvVud2[[d]] <- Vd.inv[[d]]%*%Vuds[[2]]
            tr.VinvVud2 <-  tr.VinvVud2 + sum(diag(VinvVud2[[d]]))
            XtVinvVud2VinvX[[d]] <- t(VinvXd[[d]])%*%Vuds[[2]]%*%VinvXd[[d]]
            ytVinvVud2Vinvy <- ytVinvVud2Vinvy + t(Vinvyd[[d]])%*%Vuds[[2]]%*%Vinvyd[[d]]
            ytVinvVud2VinvX <- ytVinvVud2VinvX + t(Vinvyd[[d]])%*%Vuds[[2]]%*%VinvXd[[d]]
            SumXtVinvVud2VinvX <- SumXtVinvVud2VinvX + XtVinvVud2VinvX[[d]]
            
            #S3
            VinvVud3[[d]] <- Vd.inv[[d]]%*%Vuds[[3]]
            tr.VinvVud3 <-  tr.VinvVud3 + sum(diag(VinvVud3[[d]]))
            XtVinvVud3VinvX[[d]] <- t(VinvXd[[d]])%*%Vuds[[3]]%*%VinvXd[[d]]
            ytVinvVud3Vinvy <- ytVinvVud3Vinvy + t(Vinvyd[[d]])%*%Vuds[[3]]%*%Vinvyd[[d]]
            ytVinvVud3VinvX <- ytVinvVud3VinvX + t(Vinvyd[[d]])%*%Vuds[[3]]%*%VinvXd[[d]]
            SumXtVinvVud3VinvX <- SumXtVinvVud3VinvX + XtVinvVud3VinvX[[d]]
            
            #S4
            VinvVud4[[d]] <- Vd.inv[[d]]%*%Vuds[[4]]
            tr.VinvVud4 <-  tr.VinvVud4 + sum(diag(VinvVud4[[d]]))
            XtVinvVud4VinvX[[d]] <- t(VinvXd[[d]])%*%Vuds[[4]]%*%VinvXd[[d]]
            ytVinvVud4Vinvy <- ytVinvVud4Vinvy + t(Vinvyd[[d]])%*%Vuds[[4]]%*%Vinvyd[[d]]
            ytVinvVud4VinvX <- ytVinvVud4VinvX + t(Vinvyd[[d]])%*%Vuds[[4]]%*%VinvXd[[d]]
            SumXtVinvVud4VinvX <- SumXtVinvVud4VinvX + XtVinvVud4VinvX[[d]]
            
            #S5
            VinvVud5[[d]] <- Vd.inv[[d]]%*%Vuds[[5]]
            tr.VinvVud5 <-  tr.VinvVud5 + sum(diag(VinvVud5[[d]]))
            XtVinvVud5VinvX[[d]] <- t(VinvXd[[d]])%*%Vuds[[5]]%*%VinvXd[[d]]
            ytVinvVud5Vinvy <- ytVinvVud5Vinvy + t(Vinvyd[[d]])%*%Vuds[[5]]%*%Vinvyd[[d]]
            ytVinvVud5VinvX <- ytVinvVud5VinvX + t(Vinvyd[[d]])%*%Vuds[[5]]%*%VinvXd[[d]]
            SumXtVinvVud5VinvX <- SumXtVinvVud5VinvX + XtVinvVud5VinvX[[d]]
            
            #S6
            VinvVud6[[d]] <- Vd.inv[[d]]%*%Vuds[[6]]
            tr.VinvVud6 <-  tr.VinvVud6 + sum(diag(VinvVud6[[d]]))
            XtVinvVud6VinvX[[d]] <- t(VinvXd[[d]])%*%Vuds[[6]]%*%VinvXd[[d]]
            ytVinvVud6Vinvy <- ytVinvVud6Vinvy + t(Vinvyd[[d]])%*%Vuds[[6]]%*%Vinvyd[[d]]
            ytVinvVud6VinvX <- ytVinvVud6VinvX + t(Vinvyd[[d]])%*%Vuds[[6]]%*%VinvXd[[d]]
            SumXtVinvVud6VinvX <- SumXtVinvVud6VinvX + XtVinvVud6VinvX[[d]]
            
# Calculos necesarios para obtener los elementos de la matriz de informacion de Fisher            
            
            #F11
            VinvVud1VinvVud1[[d]] <- VinvVud1[[d]]%*%VinvVud1[[d]]
            tr.VinvVud1VinvVud1 <- tr.VinvVud1VinvVud1 + sum(diag(VinvVud1VinvVud1[[d]]))
            XtVinvVud1VinvVud1VinvX[[d]]  <- t(Xd[[d]])%*%VinvVud1VinvVud1[[d]]%*%VinvXd[[d]]
            
            #F22
            VinvVud2VinvVud2[[d]] <- VinvVud2[[d]]%*%VinvVud2[[d]]
            tr.VinvVud2VinvVud2 <- tr.VinvVud2VinvVud2 + sum(diag(VinvVud2VinvVud2[[d]]))
            XtVinvVud2VinvVud2VinvX[[d]] <- t(Xd[[d]])%*%VinvVud2VinvVud2[[d]]%*%VinvXd[[d]]
            
            #F33
            VinvVud3VinvVud3[[d]] <- VinvVud3[[d]]%*%VinvVud3[[d]]
            tr.VinvVud3VinvVud3 <- tr.VinvVud3VinvVud3 + sum(diag(VinvVud3VinvVud3[[d]]))
            XtVinvVud3VinvVud3VinvX[[d]] <- t(Xd[[d]])%*%VinvVud3VinvVud3[[d]]%*%VinvXd[[d]]
            
            #F44
            VinvVud4VinvVud4[[d]] <- VinvVud4[[d]]%*%VinvVud4[[d]]
            tr.VinvVud4VinvVud4 <- tr.VinvVud4VinvVud4 + sum(diag(VinvVud4VinvVud4[[d]]))
            XtVinvVud4VinvVud4VinvX[[d]] <- t(Xd[[d]])%*%VinvVud4VinvVud4[[d]]%*%VinvXd[[d]]
            
            #F55
            VinvVud5VinvVud5[[d]] <- VinvVud5[[d]]%*%VinvVud5[[d]]
            tr.VinvVud5VinvVud5 <- tr.VinvVud5VinvVud5 + sum(diag(VinvVud5VinvVud5[[d]]))
            XtVinvVud5VinvVud5VinvX[[d]] <- t(Xd[[d]])%*%VinvVud5VinvVud5[[d]]%*%VinvXd[[d]]
            
            #F66
            VinvVud6VinvVud6[[d]] <- VinvVud6[[d]]%*%VinvVud6[[d]]
            tr.VinvVud6VinvVud6 <- tr.VinvVud6VinvVud6 + sum(diag(VinvVud6VinvVud6[[d]]))
            XtVinvVud6VinvVud6VinvX[[d]] <- t(Xd[[d]])%*%VinvVud6VinvVud6[[d]]%*%VinvXd[[d]]
            
            #F12
            VinvVud1VinvVud2[[d]] <- VinvVud1[[d]]%*%VinvVud2[[d]]
            tr.VinvVud1VinvVud2 <- tr.VinvVud1VinvVud2 + sum(diag(VinvVud1VinvVud2[[d]]))
            XtVinvVud1VinvVud2VinvX[[d]] <- t(Xd[[d]])%*%VinvVud1VinvVud2[[d]]%*%VinvXd[[d]]
            
            #F13
            VinvVud1VinvVud3[[d]] <- VinvVud1[[d]]%*%VinvVud3[[d]]
            tr.VinvVud1VinvVud3 <- tr.VinvVud1VinvVud3 + sum(diag(VinvVud1VinvVud3[[d]]))
            XtVinvVud1VinvVud3VinvX[[d]] <- t(Xd[[d]])%*%VinvVud1VinvVud3[[d]]%*%VinvXd[[d]]
            
            #F14
            VinvVud1VinvVud4[[d]] <- VinvVud1[[d]]%*%VinvVud4[[d]]
            tr.VinvVud1VinvVud4 <- tr.VinvVud1VinvVud4 + sum(diag(VinvVud1VinvVud4[[d]]))
            XtVinvVud1VinvVud4VinvX[[d]] <- t(Xd[[d]])%*%VinvVud1VinvVud4[[d]]%*%VinvXd[[d]]
            
            #F15
            VinvVud1VinvVud5[[d]] <- VinvVud1[[d]]%*%VinvVud5[[d]]
            tr.VinvVud1VinvVud5 <- tr.VinvVud1VinvVud5 + sum(diag(VinvVud1VinvVud5[[d]]))
            XtVinvVud1VinvVud5VinvX[[d]] <- t(Xd[[d]])%*%VinvVud1VinvVud5[[d]]%*%VinvXd[[d]]
            
            #F16
            VinvVud1VinvVud6[[d]] <- VinvVud1[[d]]%*%VinvVud6[[d]]
            tr.VinvVud1VinvVud6 <- tr.VinvVud1VinvVud6 + sum(diag(VinvVud1VinvVud6[[d]]))
            XtVinvVud1VinvVud6VinvX[[d]] <- t(Xd[[d]])%*%VinvVud1VinvVud6[[d]]%*%VinvXd[[d]]
            
            #F23
            VinvVud2VinvVud3[[d]] <- VinvVud2[[d]]%*%VinvVud3[[d]]
            tr.VinvVud2VinvVud3 <- tr.VinvVud2VinvVud3 + sum(diag(VinvVud2VinvVud3[[d]]))
            XtVinvVud2VinvVud3VinvX[[d]] <- t(Xd[[d]])%*%VinvVud2VinvVud3[[d]]%*%VinvXd[[d]]
            
            #F24
            VinvVud2VinvVud4[[d]] <- VinvVud2[[d]]%*%VinvVud4[[d]]
            tr.VinvVud2VinvVud4 <- tr.VinvVud2VinvVud4 + sum(diag(VinvVud2VinvVud4[[d]]))
            XtVinvVud2VinvVud4VinvX[[d]] <- t(Xd[[d]])%*%VinvVud2VinvVud4[[d]]%*%VinvXd[[d]]
            
            #F25
            VinvVud2VinvVud5[[d]] <- VinvVud2[[d]]%*%VinvVud5[[d]]
            tr.VinvVud2VinvVud5 <- tr.VinvVud2VinvVud5 + sum(diag(VinvVud2VinvVud5[[d]]))
            XtVinvVud2VinvVud5VinvX[[d]] <- t(Xd[[d]])%*%VinvVud2VinvVud5[[d]]%*%VinvXd[[d]]
            
            #F26
            VinvVud2VinvVud6[[d]] <- VinvVud2[[d]]%*%VinvVud6[[d]]
            tr.VinvVud2VinvVud6 <- tr.VinvVud2VinvVud6 + sum(diag(VinvVud2VinvVud6[[d]]))
            XtVinvVud2VinvVud6VinvX[[d]] <- t(Xd[[d]])%*%VinvVud2VinvVud6[[d]]%*%VinvXd[[d]] 
            
            #F34
            VinvVud3VinvVud4[[d]] <- VinvVud3[[d]]%*%VinvVud4[[d]]
            tr.VinvVud3VinvVud4 <- tr.VinvVud3VinvVud4 + sum(diag(VinvVud3VinvVud4[[d]]))
            XtVinvVud3VinvVud4VinvX[[d]] <- t(Xd[[d]])%*%VinvVud3VinvVud4[[d]]%*%VinvXd[[d]]
            
            #F35
            VinvVud3VinvVud5[[d]] <- VinvVud3[[d]]%*%VinvVud5[[d]]
            tr.VinvVud3VinvVud5 <- tr.VinvVud3VinvVud5 + sum(diag(VinvVud3VinvVud5[[d]]))
            XtVinvVud3VinvVud5VinvX[[d]] <- t(Xd[[d]])%*%VinvVud3VinvVud5[[d]]%*%VinvXd[[d]]
            
            #F36
            VinvVud3VinvVud6[[d]] <- VinvVud3[[d]]%*%VinvVud6[[d]]
            tr.VinvVud3VinvVud6 <- tr.VinvVud3VinvVud6 + sum(diag(VinvVud3VinvVud6[[d]]))
            XtVinvVud3VinvVud6VinvX[[d]] <- t(Xd[[d]])%*%VinvVud3VinvVud6[[d]]%*%VinvXd[[d]] 
            
            #F45
            VinvVud4VinvVud5[[d]] <- VinvVud4[[d]]%*%VinvVud5[[d]]
            tr.VinvVud4VinvVud5 <- tr.VinvVud4VinvVud5 + sum(diag(VinvVud4VinvVud5[[d]]))
            XtVinvVud4VinvVud5VinvX[[d]] <- t(Xd[[d]])%*%VinvVud4VinvVud5[[d]]%*%VinvXd[[d]]
            
            #F46
            VinvVud4VinvVud6[[d]] <- VinvVud4[[d]]%*%VinvVud6[[d]]
            tr.VinvVud4VinvVud6 <- tr.VinvVud4VinvVud6 + sum(diag(VinvVud4VinvVud6[[d]]))
            XtVinvVud4VinvVud6VinvX[[d]] <- t(Xd[[d]])%*%VinvVud4VinvVud6[[d]]%*%VinvXd[[d]]          
            
            #F56
            VinvVud5VinvVud6[[d]] <- VinvVud5[[d]]%*%VinvVud6[[d]]
            tr.VinvVud5VinvVud6 <- tr.VinvVud5VinvVud6 + sum(diag(VinvVud5VinvVud6[[d]]))
            XtVinvVud5VinvVud6VinvX[[d]] <- t(Xd[[d]])%*%VinvVud5VinvVud6[[d]]%*%VinvXd[[d]]  
            
            if(FLAG==1){
                FLAG <- 0
                ITER <- MAXITER 
                break
            }
        }
        
        Q <- solve(Q.inv)
        
        tr.XtVinvVud1VinvXQ <- tr.XtVinvVud2VinvXQ <-tr.XtVinvVud3VinvXQ <- tr.XtVinvVud4VinvXQ<-tr.XtVinvVud5VinvXQ <- tr.XtVinvVud6VinvXQ <- 0
        
        
        tr.XtVinvVud1VinvVud1VinvXQ <- tr.XtVinvVud2VinvVud2VinvXQ <- tr.XtVinvVud3VinvVud3VinvXQ <- tr.XtVinvVud4VinvVud4VinvXQ <- tr.XtVinvVud5VinvVud5VinvXQ <- tr.XtVinvVud6VinvVud6VinvXQ <-0
        tr.XtVinvVud1VinvVud2VinvXQ <- tr.XtVinvVud1VinvVud3VinvXQ <- tr.XtVinvVud1VinvVud4VinvXQ <- tr.XtVinvVud1VinvVud5VinvXQ <- tr.XtVinvVud1VinvVud6VinvXQ <-0
        tr.XtVinvVud2VinvVud3VinvXQ <- tr.XtVinvVud2VinvVud4VinvXQ <- tr.XtVinvVud2VinvVud5VinvXQ <- tr.XtVinvVud2VinvVud6VinvXQ <-0
        tr.XtVinvVud3VinvVud4VinvXQ <- tr.XtVinvVud3VinvVud5VinvXQ <- tr.XtVinvVud3VinvVud6VinvXQ <- 0
        tr.XtVinvVud4VinvVud5VinvXQ <- tr.XtVinvVud4VinvVud6VinvXQ <- 0
        tr.XtVinvVud5VinvVud6VinvXQ <- 0
        
        
         for(d in 1:D){
            tr.XtVinvVud1VinvXQ <- tr.XtVinvVud1VinvXQ + sum(diag(XtVinvVud1VinvX[[d]]%*%Q))
            tr.XtVinvVud2VinvXQ <- tr.XtVinvVud2VinvXQ + sum(diag(XtVinvVud2VinvX[[d]]%*%Q))
            tr.XtVinvVud3VinvXQ <- tr.XtVinvVud3VinvXQ + sum(diag(XtVinvVud3VinvX[[d]]%*%Q))
            tr.XtVinvVud4VinvXQ <- tr.XtVinvVud4VinvXQ + sum(diag(XtVinvVud4VinvX[[d]]%*%Q))
            tr.XtVinvVud5VinvXQ <- tr.XtVinvVud5VinvXQ + sum(diag(XtVinvVud5VinvX[[d]]%*%Q))
            tr.XtVinvVud6VinvXQ <- tr.XtVinvVud6VinvXQ + sum(diag(XtVinvVud6VinvX[[d]]%*%Q))
            
            tr.XtVinvVud1VinvVud1VinvXQ <- tr.XtVinvVud1VinvVud1VinvXQ + sum(diag(XtVinvVud1VinvVud1VinvX[[d]]%*%Q))
            tr.XtVinvVud2VinvVud2VinvXQ <- tr.XtVinvVud2VinvVud2VinvXQ + sum(diag(XtVinvVud2VinvVud2VinvX[[d]]%*%Q))
            tr.XtVinvVud3VinvVud3VinvXQ <- tr.XtVinvVud3VinvVud3VinvXQ + sum(diag(XtVinvVud3VinvVud3VinvX[[d]]%*%Q))
            tr.XtVinvVud4VinvVud4VinvXQ <- tr.XtVinvVud4VinvVud4VinvXQ + sum(diag(XtVinvVud4VinvVud4VinvX[[d]]%*%Q))
            tr.XtVinvVud5VinvVud5VinvXQ <- tr.XtVinvVud5VinvVud5VinvXQ + sum(diag(XtVinvVud5VinvVud5VinvX[[d]]%*%Q))
            tr.XtVinvVud6VinvVud6VinvXQ <- tr.XtVinvVud6VinvVud6VinvXQ + sum(diag(XtVinvVud6VinvVud6VinvX[[d]]%*%Q))
            
            tr.XtVinvVud1VinvVud2VinvXQ <- tr.XtVinvVud1VinvVud2VinvXQ + sum(diag(XtVinvVud1VinvVud2VinvX[[d]]%*%Q))
            tr.XtVinvVud1VinvVud3VinvXQ <- tr.XtVinvVud1VinvVud3VinvXQ + sum(diag(XtVinvVud1VinvVud3VinvX[[d]]%*%Q))
            tr.XtVinvVud1VinvVud4VinvXQ <- tr.XtVinvVud1VinvVud4VinvXQ + sum(diag(XtVinvVud1VinvVud4VinvX[[d]]%*%Q))
            tr.XtVinvVud1VinvVud5VinvXQ <- tr.XtVinvVud1VinvVud5VinvXQ + sum(diag(XtVinvVud1VinvVud5VinvX[[d]]%*%Q))
            tr.XtVinvVud1VinvVud6VinvXQ <- tr.XtVinvVud1VinvVud6VinvXQ + sum(diag(XtVinvVud1VinvVud6VinvX[[d]]%*%Q))
            
            tr.XtVinvVud2VinvVud3VinvXQ <- tr.XtVinvVud2VinvVud3VinvXQ + sum(diag(XtVinvVud2VinvVud3VinvX[[d]]%*%Q))
            tr.XtVinvVud2VinvVud4VinvXQ <- tr.XtVinvVud2VinvVud4VinvXQ + sum(diag(XtVinvVud2VinvVud4VinvX[[d]]%*%Q))
            tr.XtVinvVud2VinvVud5VinvXQ <- tr.XtVinvVud2VinvVud5VinvXQ + sum(diag(XtVinvVud2VinvVud5VinvX[[d]]%*%Q))
            tr.XtVinvVud2VinvVud6VinvXQ <- tr.XtVinvVud2VinvVud6VinvXQ + sum(diag(XtVinvVud2VinvVud6VinvX[[d]]%*%Q))
           
            tr.XtVinvVud3VinvVud4VinvXQ <- tr.XtVinvVud3VinvVud4VinvXQ + sum(diag(XtVinvVud3VinvVud4VinvX[[d]]%*%Q))
            tr.XtVinvVud3VinvVud5VinvXQ <- tr.XtVinvVud3VinvVud5VinvXQ + sum(diag(XtVinvVud3VinvVud5VinvX[[d]]%*%Q))
            tr.XtVinvVud3VinvVud6VinvXQ <- tr.XtVinvVud3VinvVud6VinvXQ + sum(diag(XtVinvVud3VinvVud6VinvX[[d]]%*%Q))
            
            tr.XtVinvVud4VinvVud5VinvXQ <- tr.XtVinvVud4VinvVud5VinvXQ + sum(diag(XtVinvVud4VinvVud5VinvX[[d]]%*%Q))
            tr.XtVinvVud4VinvVud6VinvXQ <- tr.XtVinvVud4VinvVud6VinvXQ + sum(diag(XtVinvVud4VinvVud6VinvX[[d]]%*%Q))
            
            tr.XtVinvVud5VinvVud6VinvXQ <- tr.XtVinvVud5VinvVud6VinvXQ + sum(diag(XtVinvVud5VinvVud6VinvX[[d]]%*%Q))
            
        }
        
        tr.PV1 <- tr.VinvVud1 - tr.XtVinvVud1VinvXQ 
        tr.PV2 <- tr.VinvVud2 - tr.XtVinvVud2VinvXQ
        tr.PV3 <- tr.VinvVud3 - tr.XtVinvVud3VinvXQ 
        tr.PV4 <- tr.VinvVud4 - tr.XtVinvVud4VinvXQ
        tr.PV5 <- tr.VinvVud5 - tr.XtVinvVud5VinvXQ 
        tr.PV6 <- tr.VinvVud6 - tr.XtVinvVud6VinvXQ
        
        SumXtVinvVud1VinvXQ <- SumXtVinvVud1VinvX%*%Q
        SumXtVinvVud2VinvXQ <- SumXtVinvVud2VinvX%*%Q
        SumXtVinvVud3VinvXQ <- SumXtVinvVud3VinvX%*%Q
        SumXtVinvVud4VinvXQ <- SumXtVinvVud4VinvX%*%Q
        SumXtVinvVud5VinvXQ <- SumXtVinvVud5VinvX%*%Q
        SumXtVinvVud6VinvXQ <- SumXtVinvVud6VinvX%*%Q
       
        tr.PV1PV1  <- tr.VinvVud1VinvVud1 - 2*tr.XtVinvVud1VinvVud1VinvXQ + sum(diag(SumXtVinvVud1VinvXQ%*%SumXtVinvVud1VinvXQ))
        tr.PV1PV2  <- tr.VinvVud1VinvVud2 - 2*tr.XtVinvVud1VinvVud2VinvXQ + sum(diag(SumXtVinvVud1VinvXQ%*%SumXtVinvVud2VinvXQ))
        tr.PV1PV3  <- tr.VinvVud1VinvVud3 - 2*tr.XtVinvVud1VinvVud3VinvXQ + sum(diag(SumXtVinvVud1VinvXQ%*%SumXtVinvVud3VinvXQ))
        tr.PV1PV4  <- tr.VinvVud1VinvVud4 - 2*tr.XtVinvVud1VinvVud4VinvXQ + sum(diag(SumXtVinvVud1VinvXQ%*%SumXtVinvVud4VinvXQ))
        tr.PV1PV5  <- tr.VinvVud1VinvVud5 - 2*tr.XtVinvVud1VinvVud5VinvXQ + sum(diag(SumXtVinvVud1VinvXQ%*%SumXtVinvVud5VinvXQ))
        tr.PV1PV6  <- tr.VinvVud1VinvVud6 - 2*tr.XtVinvVud1VinvVud6VinvXQ + sum(diag(SumXtVinvVud1VinvXQ%*%SumXtVinvVud6VinvXQ))
       
        tr.PV2PV2  <- tr.VinvVud2VinvVud2 - 2*tr.XtVinvVud2VinvVud2VinvXQ + sum(diag(SumXtVinvVud2VinvXQ%*%SumXtVinvVud2VinvXQ))
        tr.PV2PV3  <- tr.VinvVud2VinvVud3 - 2*tr.XtVinvVud2VinvVud3VinvXQ + sum(diag(SumXtVinvVud2VinvXQ%*%SumXtVinvVud3VinvXQ))
        tr.PV2PV4  <- tr.VinvVud2VinvVud4 - 2*tr.XtVinvVud2VinvVud4VinvXQ + sum(diag(SumXtVinvVud2VinvXQ%*%SumXtVinvVud4VinvXQ))
        tr.PV2PV5  <- tr.VinvVud2VinvVud5 - 2*tr.XtVinvVud2VinvVud5VinvXQ + sum(diag(SumXtVinvVud2VinvXQ%*%SumXtVinvVud5VinvXQ))
        tr.PV2PV6  <- tr.VinvVud2VinvVud6 - 2*tr.XtVinvVud2VinvVud6VinvXQ + sum(diag(SumXtVinvVud2VinvXQ%*%SumXtVinvVud6VinvXQ))
        
        tr.PV3PV3  <- tr.VinvVud3VinvVud3 - 2*tr.XtVinvVud3VinvVud3VinvXQ + sum(diag(SumXtVinvVud3VinvXQ%*%SumXtVinvVud3VinvXQ))
        tr.PV3PV4  <- tr.VinvVud3VinvVud4 - 2*tr.XtVinvVud3VinvVud4VinvXQ + sum(diag(SumXtVinvVud3VinvXQ%*%SumXtVinvVud4VinvXQ))
        tr.PV3PV5  <- tr.VinvVud3VinvVud5 - 2*tr.XtVinvVud3VinvVud5VinvXQ + sum(diag(SumXtVinvVud3VinvXQ%*%SumXtVinvVud5VinvXQ))
        tr.PV3PV6  <- tr.VinvVud3VinvVud6 - 2*tr.XtVinvVud3VinvVud6VinvXQ + sum(diag(SumXtVinvVud3VinvXQ%*%SumXtVinvVud6VinvXQ))
        
        tr.PV4PV4  <- tr.VinvVud4VinvVud4 - 2*tr.XtVinvVud4VinvVud4VinvXQ + sum(diag(SumXtVinvVud4VinvXQ%*%SumXtVinvVud4VinvXQ))
        tr.PV4PV5  <- tr.VinvVud4VinvVud5 - 2*tr.XtVinvVud4VinvVud5VinvXQ + sum(diag(SumXtVinvVud4VinvXQ%*%SumXtVinvVud5VinvXQ))
        tr.PV4PV6  <- tr.VinvVud4VinvVud6 - 2*tr.XtVinvVud4VinvVud6VinvXQ + sum(diag(SumXtVinvVud4VinvXQ%*%SumXtVinvVud6VinvXQ))
        
        tr.PV5PV5  <- tr.VinvVud5VinvVud5 - 2*tr.XtVinvVud5VinvVud5VinvXQ + sum(diag(SumXtVinvVud5VinvXQ%*%SumXtVinvVud5VinvXQ))
        tr.PV5PV6  <- tr.VinvVud5VinvVud6 - 2*tr.XtVinvVud5VinvVud6VinvXQ + sum(diag(SumXtVinvVud5VinvXQ%*%SumXtVinvVud6VinvXQ))
        
        tr.PV6PV6  <- tr.VinvVud6VinvVud6 - 2*tr.XtVinvVud6VinvVud6VinvXQ + sum(diag(SumXtVinvVud6VinvXQ%*%SumXtVinvVud6VinvXQ))
        
        ytVinvXQ <- ytVinvX%*%Q
        ytPV1Py  <- ytVinvVud1Vinvy - 2*ytVinvVud1VinvX%*%t(ytVinvXQ) + ytVinvXQ%*%SumXtVinvVud1VinvX%*%t(ytVinvXQ)
        ytPV2Py  <- ytVinvVud2Vinvy - 2*ytVinvVud2VinvX%*%t(ytVinvXQ) + ytVinvXQ%*%SumXtVinvVud2VinvX%*%t(ytVinvXQ)
        ytPV3Py  <- ytVinvVud3Vinvy - 2*ytVinvVud3VinvX%*%t(ytVinvXQ) + ytVinvXQ%*%SumXtVinvVud3VinvX%*%t(ytVinvXQ)
        ytPV4Py  <- ytVinvVud4Vinvy - 2*ytVinvVud4VinvX%*%t(ytVinvXQ) + ytVinvXQ%*%SumXtVinvVud4VinvX%*%t(ytVinvXQ)
        ytPV5Py  <- ytVinvVud5Vinvy - 2*ytVinvVud5VinvX%*%t(ytVinvXQ) + ytVinvXQ%*%SumXtVinvVud5VinvX%*%t(ytVinvXQ)
        ytPV6Py  <- ytVinvVud6Vinvy - 2*ytVinvVud6VinvX%*%t(ytVinvXQ) + ytVinvXQ%*%SumXtVinvVud6VinvX%*%t(ytVinvXQ)
        
        
        
        
        
        # Scores 
        
        S1 <- -0.5*tr.PV1 + 0.5*ytPV1Py
        S2 <- -0.5*tr.PV2 + 0.5*ytPV2Py
        S3 <- -0.5*tr.PV3 + 0.5*ytPV3Py
        S4 <- -0.5*tr.PV4 + 0.5*ytPV4Py
        S5 <- -0.5*tr.PV5 + 0.5*ytPV5Py
        S6 <- -0.5*tr.PV6 + 0.5*ytPV6Py
        
       # Matriz de informacion de Fisher
        
        F11  <- 0.5*tr.PV1PV1
        F12  <- 0.5*tr.PV1PV2
        F13  <- 0.5*tr.PV1PV3
        F14  <- 0.5*tr.PV1PV4
        F15  <- 0.5*tr.PV1PV5
        F16  <- 0.5*tr.PV1PV6
        
        F22  <- 0.5*tr.PV2PV2
        F23  <- 0.5*tr.PV2PV3
        F24  <- 0.5*tr.PV2PV4
        F25  <- 0.5*tr.PV2PV5
        F26  <- 0.5*tr.PV2PV6
        
        F33  <- 0.5*tr.PV3PV3
        F34  <- 0.5*tr.PV3PV4
        F35  <- 0.5*tr.PV3PV5
        F36  <- 0.5*tr.PV3PV6
        
        F44  <- 0.5*tr.PV4PV4
        F45  <- 0.5*tr.PV4PV5
        F46  <- 0.5*tr.PV4PV6
        
        F55  <- 0.5*tr.PV5PV5
        F56  <- 0.5*tr.PV5PV6
        
        F66  <- 0.5*tr.PV6PV6
        
        
        if(ITER>1){
         F.sig.prev<-Fsig
        	Q.prev<-Q
        	theta.f.prev<-theta.f
        	}
        
        
        Ssig <- c(S1,S2,S3,S4,S5,S6)
        Fsig <- matrix(c(F11,F12,F13,F14,F15,F16,F12,F22,F23,F24,F25,F26,F13,F23,F33,F34,F35,F36,F14,F24,F34,F44,F45,F46,F15,F25,F35,F45,F55,F56,F16,F26,F36,F46,F56,F66),nrow=6,byrow=TRUE)
    
         
        
        
        
        if(abs(det(Fsig))< 0.000000001 || sum(abs(Fsig))>10000000) {
            #print(SumaF)
            ITER <- MAXITER
            Bad <- Bad+1
            break
        }
        
        Fsig.inv <- solve(Fsig)
        dif <- Fsig.inv%*%Ssig  # formula de actualizacion del algoritmo
        
       
    	
          theta.rml<-theta.f
          	
          	
        # Actualizacion del parametro  	
        
         theta.f <- theta.rml + dif         
          
         Vud.f<-UveU(theta.f)
         
         Vuds<-FirstDer(theta.f)
       
        # Condiciones para la matriz de var-cor. Si no se cumplen, nos quedamos con el valor de la iteracion anterior
        
        if( theta.f[1]<=0 || theta.f[2]<=0 ||  theta.f[3]<=0 ||  theta.f[4]>1 ||  theta.f[5]>1 || theta.f[6]>1 || det(Vud.f)<=0){
        	theta.f<-theta.f.prev
        	Fsig<-F.sig.prev
        	Q<-Q.prev
        	return(list(theta.f,Fsig, ITER, Bad, Q))
        	}
        
			    
			    # Criterio de parada
         identical(as.numeric(abs(dif)<0.000001),rep(1,6))
        if(identical(as.numeric(abs(dif)<0.000001),rep(1,6))){   # Criterio de parada
          
            break
        }
    }
    
   
    return(list(theta.f,Fsig, ITER, Bad, Q))
}

