
###############################################################################
###############################################################################
###
###        Metodo de generacion de Normal Multivariante basado en el
###                       algoritmo de Cholesky
###                                                            
###                           Diciembre 2017                           
###

### Autor: Agustin Perez Martin  



### Algoritmo de Cholesky
cholesky <- function(m.cov){
  # Paso 1
  m.A <- matrix(0, nrow = nrow(m.cov), ncol = ncol(m.cov)); m.A
  # Paso 2
  m.A[1,1] <- sqrt(m.cov[1,1])
  for(i in 2:nrow(m.A)) {
    m.A[i,1] <- m.cov[i,1]/sqrt(m.cov[1,1])
  }
  # Paso 3
  for(j in 2:ncol(m.A)) {
    m.A[j,j] <- sqrt(m.cov[j,j] - cumsum(m.A[j,]^2)[j-1])
    for(i in j+1:ncol(m.A)) {
      if(i < ncol(m.A)+1) {     
        m.A[i,j] <- (m.cov[i,j] - cumsum(m.A[i,]*m.A[j,])[j-1])/m.A[j,j]
        #cat("La i es ", i, " y la j es ", j,"\n")
      }
    }    
  }
  return(m.A)
}

# cho <- t(chol(m.cov))
# round(m.A,4)==round(cho,4)


var.norm.mult <- function(n=1, mu, m.cov){
  ###Algoritmo de Generacion
  # Paso 1
  m.A <- cholesky(m.cov)
  # Paso 2
  m.Z <- matrix(0, nrow = n, ncol = ncol(m.A))
  for(j in 1:ncol(m.A))
    m.Z[,j] <- rnorm(n)
  # Paso 3
  return( t(mu + t( m.Z%*%t(m.A)) ) )     #X <- t(mu + t( m.Z%*%t(m.A)) ); X
  
}




# V <- matrix(-0.25, nrow=3, ncol=3)
# diag(V) <- 1
# X <- var.norm.mult(n = 500, mu=c(0,0), m.cov=V)
# 
#####################################################################
###     Comprobaciones y gr?ficos

# det(t(X)%*%X)
# round(solve(t(X)%*%X),2)
# cov(X)
# cor(X)
# 
# par(mfrow = c(2,3))
# for(j in 1:ncol(X)) {
#   hist(X[,j], prob = TRUE, col = j+1, main = paste("Histograma de X", j), xlab = paste("X", j))
#   lines(density(X[,j]), lwd = 2)
#   print(ks.test(X[,j], "pnorm", 0, 1))
#   # print(shapiro.test(X[,j]))
#   cat("Desviacion ", sd(X[,j]))
# }
# 
# summary(X)
# sd(X)






###