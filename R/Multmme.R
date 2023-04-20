###############################################################################
###
###                 Aplicacion paquete mme (Acompo)              
###              
###                          Febrero 2017                               
###                         Agustin Perez


MMEcompo <- function(D,X,y,i) {
	
	library(mme)
	
	       # Utilizacion del paquete de R mme para comparacion de resultados
	
        sydk<-pdek<-array(rep(0,3*D),dim=c(3,D))      
        X1<-X2<-X3<-Y1<-Y2<-Y3<-Y4<-Area<-Time<-Sample<-Sampleround <-Population<-vector()
   	    for(d in 1:D)	{
   	  	Area[d]=d
   	  	Time[d]=1
   	  	Sample[d]=100
   	  	Population[d]=1000
   	    X1[d]<-X[[d]][1,2]
   	    X2[d]<-X[[d]][2,4]
   	    X3[d]<-X[[d]][3,6]
   	   
   	    Y1[d]<- Sample[d]*(exp(y[[d]][1])/(1+ sum(exp(y[[d]]))))
   	    Y2[d]<- Sample[d]*(exp(y[[d]][2])/(1+ sum(exp(y[[d]]))))
   	    Y3[d]<- Sample[d]*(exp(y[[d]][3])/(1+ sum(exp(y[[d]]))))
   	    Y4[d]<-Sample[d]-sum(Y1[d],Y2[d],Y3[d])+1
   	
   	    
   	  	Y1[d]<-round(Y1[d],0)
   	  	Y2[d]<-round(Y2[d],0)
   	  	Y3[d]<-round(Y3[d],0)
   	  	Y4[d]<-round(Y4[d],0)
   	  	Sampleround[d]<-sum(Y1[d],Y2[d],Y3[d],Y4[d])
   	  	
   	  	}
   	    
   	  
   	  	resmme <- data.frame(Area,Time,Sampleround, Population,	Y1,Y2,Y3,Y4,X1,X2,X3 )
   	  	
   	  	kk <- 4 # number of categories of the response variable
				pp <- c(1,1,1) # vector with the number of auxiliary variables in each category
				mod <- 1 # Model 1
				# Needed matrix and initial values
				datar <- data.mme(resmme,kk,pp,mod)
				# Model fit
				result <- model(datar$d,datar$t,pp,datar$Xk,datar$X,datar$Z,datar$initial,
                  + datar$y[,1:(kk-1)],datar$n,datar$N,mod)
				#result
     	
				# Predictors of sample totals
					sY1 <- result$mean[,1]*datar$n/datar$N # MuMdk
					sY2 <- result$mean[,2]*datar$n/datar$N
					sY3 <- result$mean[,3]*datar$n/datar$N
					resultadomme<- data.frame(Y1,Y2,Y3,sY1,sY2,sY3)
				# Predictors of domain*category probabilities
					p1<-result$mean[,1]/datar$N #pis
					p2<-result$mean[,2]/datar$N
					p3<-result$mean[,3]/datar$N
					
					
				 for(d in 1:D){
				 	sydk[1,d]<-sY1[d]
				 	sydk[2,d]<-sY2[d]
				 	sydk[3,d]<-sY3[d]
				 	pdek[1,d]<-p1[d]
				 	pdek[2,d]<-p2[d]
				 	pdek[3,d]<-p3[d]
				 }	
return(list(pdek))				
}
	