#####################################################################
#####################################################################
###
###                   CALCULO P-VALORES EN MODELO
###                           TRIVARIANTE                            
###                                                                  
#                            Octubre 2017                          

### AUTOR: Agustin Perez Martin

### Editado por:        Agustin Perez Martin
### con fecha:          31 Octubre 2017
### cambios realizados: introducida funcion de calculo de los IC



pvalue <- function(beta.fitted, fit) {
  if(!inherits(fit,"try-error")){
    
    z <- abs(beta.fitted)/sqrt(as.vector(diag(fit[[5]])))
    p <- pnorm(z, lower.tail=F)
    
    return( 2*p )
  }
  else{
    warning("Only a converging algoritm is allowed", call. = FALSE)
  }
}


IC <- function(beta.fitted, fit, conf.level = 0.95) {
  if(!inherits(fit,"try-error")){
    alpha <- 1-conf.level
    k <- 1-alpha/2
    z <- qnorm(k)
    Finv <- solve(fit[[2]])
    Qinv <- fit[[5]]
    
    lower.beta <- beta.fitted-z*sqrt(diag(Qinv))
    lower.theta <- fit[[1]]-z*sqrt(diag(Finv))
    upper.beta <- beta.fitted+z*sqrt(diag(Qinv))
    upper.theta <- fit[[1]]+z*sqrt(diag(Finv))
    
    return(list(data.frame(l.inf.beta=lower.beta,l.sup.beta=upper.beta),
                data.frame(l.inf.theta=lower.theta,l.sup.theta=upper.theta)))
  }
  else{
    warning("Only a converging algoritm is allowed", call. = FALSE)
  }
}



####