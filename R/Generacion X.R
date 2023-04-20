###############################################################################
###
###                 Generacion Variables explicativas X (Acompo)              
###              
###                                 Febrero 2017                               
###                                Agustin Perez


# Generacóon de las variables exolicativas (segun documento Acompo)

equis <- function(D) {

set.seed(19032005)	
mux1<-1
mux2<-1
mux3<-1
sigmax11<-1
sigmax22<-3/2
sigmax33<-2

     for(d in 1:D) {
    	udk[[d]]<-vector()
    	for(k in 1:3){
    		udk[[d]][k]<-runif(1,0,1)
    	}}
    
    
    for(d in 1:D) {
    	xd1[[d]]<-xd2[[d]]<-xd3[[d]]<-vector()
    
    	xd1[[d]][1]<-1
     	xd1[[d]][2]<-mux1+sigmax11^0.5*udk[[d]][1]
      xd2[[d]][1]<-1
    	xd2[[d]][2]<-mux2+sigmax22^0.5*udk[[d]][2]
      xd3[[d]][1]<-1
    	xd3[[d]][2]<-mux3+sigmax33^0.5*udk[[d]][3]
     	}
    

for(d in 1:D){
X[[d]]<-as.matrix(bdiag(t(xd1[[d]]),t(xd2[[d]]),t(xd3[[d]])))
}   
p<-ncol(X[[1]])    
    
return(list(X,p))
    }