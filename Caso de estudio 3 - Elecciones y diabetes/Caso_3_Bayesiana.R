## Parcial 3 Bayesiana ##
setwd("C:/Users/EQUIPO/OneDrive/Documentos/Documentos Ander/UNAL/Estadística Bayesiana/Caso de estudio 3")
getwd()
library(readr)
library(dplyr)
library(readxl)
library(foreach)
library(doParallel)
library(xtable)
library(LaplacesDemon)
library(xtable)
library(ggplot2)
library(gridExtra)

# ALCALDÍA DE BOGOTÁ ####

## Datos ####

Candidato <- c("C. F. Galán", "G. Bolivar", "J. D. Oviedo", "D. Molano", "R. Lara", "J. L. Vargas", "J. E. Robledo", "N. Ramos", "R. A. Quintero", "Voto en Blanco")
Cantidad <- c(493, 257, 227, 48, 41, 38, 28, 11, 3, 54)
Proporcion <- c(0.411, 0.214, 0.189, 0.040, 0.034, 0.032, 0.023, 0.009, 0.003, 0.045)

tabla1 <- data.frame(Candidato, Cantidad, Proporcion)
rm(Candidato, Cantidad, Proporcion)

## Muestreador de Gibbs ####

### Hiperparámetros ####
a <- 1
b <- 1
K <- length(tabla1$Proporcion)
nj <- tabla1$Cantidad
delta <-1

### Número de iteraciones ####
S <- 101000
B <- 10000
ac <- 0
cont <- 0
### Cada 10% anuncie
ncat <- floor(S/10)

### Matriz MCMC
THETA <- matrix(data = NA, nrow = B, ncol = K+1)
LLMULT    <- matrix(data = NA, nrow = B, ncol = 1)


### Inicializar parámetros  ####
set.seed(1010)
alp <- log(rgamma(1, shape = a, rate = b))

### Cadena  ####
ñ <- proc.time() 
set.seed(1010)
for(s in 1:S) {
  # Actualizar theta
  theta <- rdirichlet(n = 1, alpha = exp(alp)+nj)
  
  # Actualizar alpha con Metropolis
  
  # Propuesta
  alp.star   <- rnorm(n = 1, mean = alp, sd = sqrt(delta))
  
  # Tasa de aceptación
  log.r <- (sum(ddirichlet(x = theta, alpha = rep(exp(alp.star),K), log = TRUE) - ddirichlet(x = theta, alpha = rep(exp(alp),K), log = TRUE) 
                + dgamma(x = exp(alp.star), shape = a, rate = b, log = TRUE) - dgamma(x = exp(alp), shape = a, rate = b, log = TRUE)
                + alp.star - alp)) 

  # Actualizar
  if (runif(1) < exp(log.r)){
    alp <- alp.star
    ac <- ac + 1
  }
  # Almacenar
  if(s > 1000 && s%%10 == 0){
    cont <- cont + 1 
    THETA[cont,] <- c(theta, exp(alp))
    LLMULT[cont,] <- dmultinom(x = nj, prob = theta, log = T)
  }
  
  if (s == S) cont <- 0  
    
  # Progreso
  if (s%%ncat == 0) 
    cat(100*round(s/S, 1), "% completado ... ", "Tasa de aceptación del ",  100*ac/S, "% \n",  sep = "")
}
proc.time()-ñ 

### Análisis de convergencia ####
colnames(THETA) <- c(paste0("theta",1:K), "alpha")
colnames(LLMULT) <- c("ll")


neff_THETA <- coda::effectiveSize(THETA)
EMCMULT <- apply(X = THETA, MARGIN = 2, FUN = sd)/sqrt(neff_THETA)
100*abs(EMCMULT/colMeans(THETA))

### Gráfico de logverosimilitud ####
vero<-data.frame(LLMULT)
ggplot(data = vero, aes(x = 1:length(LLMULT), y = LLMULT)) +
  geom_point(shape = 16, size = 1, color = "violetred") +
  ylim(-50, -20) +
  labs(x = "Iteración", y = "Log-verosimilitud", title = "Modelo Elecciones") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.title = element_text(size = 8),
        axis.text = element_text(size=7))

rm(vero)

### Tabla con los resultados ####
THETA <- data.frame(THETA)
medias <- c(mean(THETA$theta1), mean(THETA$theta2),mean(THETA$theta3))*100
IC <- rbind(quantile(THETA$theta1, probs = c(0.025,0.975)),quantile(THETA$theta2, probs = c(0.025,0.975)),quantile(THETA$theta3, probs = c(0.025,0.975)))*100
resultados<-c(0.4902,0.1871,0.2010)*100
out<-cbind(medias,IC, resultados)
colnames(out)<-c("Estimación", "Q2.5%", "Q97.5%","Resultado ofical" )
xtable(out)

### Gráfico con los resultados ####

galan <- data.frame(
  candidato = rep(c("C. F. Galán"), each = 2),
  tipo = rep(c("Estimación"), 1),
  Resultados = rep(c("Estimación", "Resultado Oficial"), 1),
  valor = c(40.89475, 49.02),
  inferior = c(38.19, NA),
  superior = c(43.65, NA)
)
bolivar <- data.frame(
  candidato = rep(c("G. Bolivar"), each = 2),
  tipo = rep(c("Estimación"), 1),
  Resultados = rep(c("Estimación", "Resultado Oficial"), 1),
  valor = c(21.34634, 18.71),
  inferior = c(19.08776 , NA),
  superior = c(23.67598, NA)
)
oviedo <- data.frame(
  candidato = rep(c("J. D. Oviedo"), each = 2),
  tipo = rep(c("Estimación"), 1),
  Resultados = rep(c("Estimación", "Resultado Oficial"), 1),
  valor = c(18.88061, 20.10),
  inferior = c(16.70354 , NA),
  superior = c(21.13293, NA)
)
data <- rbind(galan, bolivar, oviedo)
rm(galan, bolivar, oviedo)
ggplot(data, aes(x = tipo, y = valor, color = Resultados, group = tipo)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = inferior, ymax = superior, color="Intervalo de confianza"), width = 0.3, size = 1) +
  facet_wrap(~candidato, scales = "fixed") +
  labs(
    title = "Resultados de los Candidatos a la Alcaldía",
    y = "Porcentaje de Votos",
    x = "",
  ) +
  theme_minimal() +
  ylim(c(16, 50)) +
  theme(plot.title = element_text(hjust = 0.6,size=17))

# DIABETES ####

## Datos ####
yX.train <- source("http://www2.stat.duke.edu/~pdh10/FCBS/Inline/yX.diabetes.train")[[1]]
yX.test <- source("http://www2.stat.duke.edu/~pdh10/FCBS/Inline/yX.diabetes.test")[[1]]
xtrain <- yX.train[,-1]
ytrain <- yX.train[,1]

xtest <- yX.test[,-1]
ytest <- yX.test[,1]

## Número de parámetros y tamaño de muestra  ####
p <- length(xtrain[1,])
n <- length(xtrain[,1])

## OLS  ####
beta_ols <- solve(t(xtrain)%*%xtrain)%*%t(xtrain)%*%ytrain
sig2_ols <- sum((ytrain - xtrain%*%beta_ols)^2)/(n-p)

## Número de iteraciones  ####
S <- 101000
B <- 10000
cont <- 0

## Previa unitaria ####

### Hiperparámetros ####
beta0 <- beta_ols
sigma0 <- n*sig2_ols*solve(t(xtrain)%*%xtrain)
nu0 <- 1
sig20 <- sig2_ols

### Matriz MCMC 
BETA1 <- matrix(data = NA, nrow = B, ncol = p+1)
LL1    <- matrix(data = NA, nrow = B, ncol = 1)

### Inicializar parámetros y optimizar código ####
beta <- rep(0,p)
sig2 <- 1
isigma0 <- solve(sigma0)
xtx <- t(xtrain)%*%xtrain
xty <- t(xtrain)%*%ytrain
i_n <- diag(nrow = n)

### Cadena ####
ñ <- proc.time() 
set.seed(1010)
for(s in 1:S) {
  # Actualizar beta
  sigman <- solve(isigma0+1/sig2*xtx)
  beta <- mvtnorm::rmvnorm(n = 1, mean = sigman%*%(isigma0%*%beta0+1/sig2*xty), sigma = sigman)
  
  # Actualizar sig2
  yhat <- xtrain%*%t(beta)
  sig2 <- 1/rgamma(n = 1, shape = 0.5*(nu0+n), rate = 0.5*(nu0*sig20+sum((ytrain-yhat)^2)))
  
  # Almacenar
  if(s > 1000 && s%%10 == 0){
    cont <- cont + 1 
    BETA1[cont,] <- c(beta,sig2)
    LL1[cont] <- mvtnorm::dmvnorm(x = ytrain, mean = yhat, sigma = sig2*i_n, log = T)
  }
  
  if (s == S) cont <- 0  
  
  # Progreso
  if (s%%ncat == 0) 
    cat(100*round(s/S, 1), "% completado ... \n",  sep = "")
}
colnames(BETA1) <- c(paste0("beta",1:p),"sig2")
colnames(LL1) <- c("ll")
proc.time()-ñ 

## Previa g ####

### Hiperparámetros  ####
beta0 <- rep(0,p)
g <- n
k <- g*sig2_ols
sigma0 <- k*solve(t(xtrain)%*%xtrain)
nu0 <- 1
sig20 <- sig2_ols

### Matriz MC
BETA2 <- matrix(data = NA, nrow = B, ncol = p+1)
LL2    <- matrix(data = NA, nrow = B, ncol = 1)

### Cantidades importantes ####
Hg    <- (g/(g+1))*(xtrain%*%solve(t(xtrain)%*%xtrain)%*%t(xtrain))
SSRg  <- as.numeric(t(ytrain)%*%(diag(1,n) - Hg)%*%ytrain)
Vbeta <- (g/(g+1))*solve(t(xtrain)%*%xtrain)
Ebeta <- Vbeta%*%t(xtrain)%*%ytrain
i_n <- diag(nrow = n)
ncat2 <- floor(B/10)


### Monte Carlo ####

ñ <- proc.time() 
set.seed(1010)
for(s in 1:B) {
  BETA2[s,p+1] <- 1/rgamma(1, (nu0 + n)/2, (nu0*sig20 + SSRg)/2)
  BETA2[s,c(1:p)] <- c(mvtnorm::rmvnorm(1, Ebeta, BETA2[s,p+1]*Vbeta))
  LL2[s] <- mvtnorm::dmvnorm(x = ytrain, mean = xtrain%*%(BETA2[s,c(1:p)]), sigma = BETA2[s,p+1]*i_n, log = T)
  
  # Progreso
  if (s%%ncat2 == 0) 
    cat(100*round(s/B, 1), "% completado ... \n",  sep = "")
}
colnames(BETA2) <- c(paste0("beta",1:p),"sig2")
colnames(LL2) <- c("ll")
proc.time()-ñ 

## Regresión rígida ####

### Hiperparámetros  ####
beta0 <- rep(0,p)
nu0 <- 1
sig20 <- sig2_ols
alambda <- 1
blambda <- 2
cont <- 0

### Matriz MCMC
BETA3 <- matrix(data = NA, nrow = B, ncol = p+2)
LL3   <- matrix(data = NA, nrow = B, ncol = 1)

### Inicializar parámetros y optimizar código ####
beta <- rep(0,p)
sig2 <- 1
isigma0 <- solve(sigma0)
xtx <- t(xtrain)%*%xtrain
xty <- t(xtrain)%*%ytrain
ip <- diag(nrow = p)
i_n <- diag(nrow = n)
### Necesario para la primera iteración ##
btb <- beta%*%t(beta)
### Cadena ####
ñ <- proc.time() 
set.seed(1010)
for(s in 1:S) {
  # Actualizar lambda
  lambda <- rgamma(n = 1 , shape = 0.5*p+alambda, rate = 0.5*(1/sig2)*btb+blambda)
  
  # Actualizar beta
  sigman <- solve(1/sig2*(lambda*ip+xtx)) 
  beta <- mvtnorm::rmvnorm(n = 1, mean = sigman%*%(1/sig2*xty), sigma = sigman)
  btb <- beta%*%t(beta)
  # Actualizar sig2
  yhat <- xtrain%*%t(beta)
  sig2 <- 1/rgamma(n = 1, shape = 0.5*(nu0+n+p), rate = 0.5*(nu0*sig20+lambda*btb+sum((ytrain-yhat)^2)))
  
  # Almacenar
  if(s > 1000 && s%%10 == 0){
    cont <- cont + 1 
    BETA3[cont,] <- c(beta,sig2,lambda)
    LL3[cont] <- mvtnorm::dmvnorm(x = ytrain, mean = yhat, sigma = sig2*i_n, log = T)
  }
  
  if (s == S) cont <- 0  
  
  # Progreso
  if (s%%ncat == 0) 
    cat(100*round(s/S, 1), "% completado ... \n",  sep = "")
}
colnames(BETA3) <- c(paste0("beta",1:p),"sig2","lambda")
colnames(LL3) <- c("ll")
proc.time()-ñ 

## Regresión con errores correlacionados ####

### Hiperparámetros ####
tau20 <- 50
nu0 <- 1
sig20 <- sig2_ols
arho <- 0
brho <- 1


DY <-abs(outer( (1:n),(1:n) ,"-")) # para construir la matriz de correlacion

### Matriz MCMC
BETA4 <- matrix(data = NA, nrow = B, ncol = p+2)
LL4   <- matrix(data = NA, nrow = B, ncol = 1)

### Número de iteraciones y cantidades importantes ####
S <-51000
B <- 5000
ac    <-0     
ncat3 <- floor(S/10)
delta <- 0.3
cont <- 0

### Inicializar parámetros y optimizar código ####
beta <- rep(0,p)
sig2 <- 1
itau20 <- diag(1/tau20,nrow=p)
rho <- 0.5
xtx <- t(xtrain)%*%xtrain
xty <- t(xtrain)%*%ytrain
ip <- diag(nrow = p)
i_n <- diag(nrow = n)
### Código útil para la primera iteración ####
Cor    <- rho^DY
iCor   <- solve(Cor)

### Funciones auxiliares ####
tr <- function(X) sum(diag(X))
rmvnorm <- function(n,mu,Sigma) 
{
  p<-length(mu)
  res<-matrix(0,nrow=n,ncol=p)
  if( n>0 & p>0 ) {
    E<-matrix(rnorm(n*p),n,p)
    res<-t(  t(E%*%chol(Sigma)) +c(mu))
  }
  res
}

### Cadena ####
ñ <- proc.time()
set.seed(1010)
for (s in 1:S) {
  # Actualizar beta
  V.beta <- solve( t(xtrain)%*%iCor%*%xtrain/sig2 + itau20)
  E.beta <- V.beta%*%( t(xtrain)%*%iCor%*%ytrain/sig2)
  beta   <- t(rmvnorm(1,E.beta,V.beta))
  
  # Actualizar sigma^2
  yhat <- xtrain%*%beta
  sig2 <- 1/rgamma(1,(nu0+n)/2,(nu0*sig20+t(ytrain-yhat)%*%iCor%*%(ytrain-yhat)) /2 )
  
  # Actualizar rho con Metropolis
  
  # 1. Propuesta
  rho.p <- abs(runif(1,rho-delta,rho+delta))
  rho.p <- min(rho.p, 2-rho.p)
  # 2. Tasa de aceptacion
  Cor.p <- rho.p^DY
  iCor.p <- solve(Cor.p)
  lr <- -0.5*( determinant(Cor.p,log=TRUE)$mod - determinant(Cor,log=TRUE)$mod + 
                tr( (ytrain-yhat)%*%t(ytrain-yhat)%*%(iCor.p - iCor) )/sig2 )
  # 3. Actualizar valor
  if( log(runif(1)) < lr ) { 
    rho <-rho.p
    ac <- ac+1 
    Cor    <- Cor.p
    iCor   <- iCor.p
  }
  # Almacenar
  if(s > 1000 && s%%10 == 0) {
    cont <- cont + 1
    BETA4[cont,] <- c(beta,sig2,rho)
    LL4[cont]  <- mvtnorm::dmvnorm(x = ytrain, mean = yhat, sigma = sig2*Cor, log = T)
  # Progreso
  if (s%%ncat3 == 0) 
    cat(100*round(s/S, 1), "% completado ... ", "Tasa de aceptación del ",  100*ac/S, "% \n",  sep = "")
  }
}
colnames(BETA4) <- c(paste0("beta",1:p),"sig2","rho")
colnames(LL4) <- c("ll")
proc.time()-ñ 

## Análisis de convergencia ####

### Modelo 1 ####
neff_BETA1 <- coda::effectiveSize(BETA1)
summary(neff_BETA1[1:64])
neff_BETA1[65]
EMC1 <- apply(X = BETA1, MARGIN = 2, FUN = sd)/sqrt(neff_BETA1)
summary(100*abs(EMC1/colMeans(BETA1))[1:64])
100*abs(EMC1/colMeans(BETA1))[65]

### Modelo 2 ####
neff_BETA2 <- coda::effectiveSize(BETA2)
summary(neff_BETA2[1:64])
neff_BETA2[65]
EMC2 <- apply(X = BETA2, MARGIN = 2, FUN = sd)/sqrt(neff_BETA2)
summary(100*abs(EMC2/colMeans(BETA2))[1:64])
100*abs(EMC2/colMeans(BETA2))[65]

### Modelo 3 ####
neff_BETA3 <- coda::effectiveSize(BETA3)
summary(neff_BETA3[1:64])
neff_BETA3[65:66]
EMC3 <- apply(X = BETA3, MARGIN = 2, FUN = sd)/sqrt(neff_BETA3)
summary(100*abs(EMC3/colMeans(BETA3))[1:64])
100*abs(EMC3/colMeans(BETA3))[65:66]

### Modelo 4 ####
neff_BETA4 <- coda::effectiveSize(BETA4)
summary(neff_BETA4[1:64])
neff_BETA4[65:66]
EMC4 <- apply(X = BETA4, MARGIN = 2, FUN = sd)/sqrt(neff_BETA4)
summary(100*abs(EMC4/colMeans(BETA4))[1:64])
100*abs(EMC4/colMeans(BETA4))[65:66]

### Gráfico de logversimilitud ####
vero1 <- data.frame(LL1)
vero2 <- data.frame(LL2)
vero3 <- data.frame(LL3)
vero4 <- data.frame(LL4)

plot1 <- ggplot(data = vero1, aes(x = 1:length(LL1), y = LL1)) +
  geom_point(shape = 16, size = 0.8, color = "lightcoral") +
  ylim(-400, -320) +
  labs(x = "Iteración", y = "Log-verosimilitud", title = "Modelo 1") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))

plot2 <- ggplot(data = vero2, aes(x = 1:length(LL2), y = LL2)) +
  geom_point(shape = 16, size = 0.8, color = "lightgreen") +
  ylim(-400, -320) +
  labs(x = "Iteración", y = "Log-verosimilitud", title = "Modelo 2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))

plot3 <- ggplot(data = vero3, aes(x = 1:length(LL3), y = LL3)) +
  geom_point(shape = 16, size = 0.8, color = "lightseagreen") +
  ylim(-400, -320) +
  labs(x = "Iteración", y = "Log-verosimilitud", title = "Modelo 3") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))

plot4 <- ggplot(data = vero4, aes(x = 1:length(LL4), y = LL4)) +
  geom_point(shape = 16, size = 0.8, color = "violet") +
  ylim(-400, -320) +
  labs(x = "Iteración", y = "Log-verosimilitud", title = "Modelo 4") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))

# Organizar los gráficos en un panel 2x2
grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)
rm(vero1,vero2,vero3,vero4)

## Pregunta 1 ####

### Generar ytest gorro ####
ntest <- 100

#### Modelo 1 ####
beta1hat <- colMeans(BETA1[,1:p])
ythat1 <- xtest%*%beta1hat
error1 <- 1/ntest*sum(abs(ytest - ythat1))

#### Modelo 2 ####
beta2hat <- colMeans(BETA2[,1:p])
ythat2 <- xtest%*%beta2hat
error2 <- 1/ntest*sum(abs(ytest - ythat2))

#### Modelo 3 ####
beta3hat <- colMeans(BETA3[,1:p])
ythat3 <- xtest%*%beta3hat
error3 <- 1/ntest*sum(abs(ytest - ythat3))

#### Modelo 4 ####
beta4hat <- colMeans(BETA4[,1:p])
ythat4 <- xtest%*%beta4hat
error4 <- 1/ntest*sum(abs(ytest - ythat4))

#### Gráfico ytest gorro vs ytest #### 
layout(matrix(c(1,2,3,4,5,6,7,8), ncol=2, byrow=TRUE), heights = c(0.09,0.9,0.09,0.9,0.09,0.9,0.09,0.9))
par(mar=c(0,2,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
title("Modelo 1", line=-1.2, cex.main=1.5)

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
title("Modelo 2", line=-1.2, cex.main=1.5)

par(mar=c(3,3,2,1))
plot(ytest,ythat1,xlab=expression(italic(y)[test]),col="#B03060",pch=16,
     main=paste0("Error abs. medio: ",round(error1, 3)),cex.main=1.5,
     ylab=expression(hat(italic(y))[test]), xlim=c(-2.3,2.3),ylim=c(-1.5,1.5),type="p")
abline(0,1,col="gray10", lwd=2)

plot(ytest,ythat2,xlab=expression(italic(y)[test]), col="#B03060",pch=16,
     main=paste0("Error abs. medio: ",round(error2, 3)),cex.main=1.5,
     ylab=expression(hat(italic(y))[test]), xlim=c(-2.3,2.3),ylim=c(-1.5,1.5))
abline(0,1,col="gray10", lwd=2)

par(mar=c(0,2,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
title("Modelo 3", line=-1.2, cex.main=1.5)

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
title("Modelo 4", line=-1.2, cex.main=1.5)

par(mar=c(3,3,2,1))
plot(ytest,ythat3,xlab=expression(italic(y)[test]), col="#B03060",pch=16,
     main=paste0("Error abs. medio: ",round(error3, 3)),cex.main=1.5,
     ylab=expression(hat(italic(y))[test]), xlim=c(-2.3,2.3),ylim=c(-1.5,1.5))
abline(0,1,col="gray10", lwd=2)

plot(ytest,ythat4,xlab=expression(italic(y)[test]), col="#B03060",pch=16,
     main=paste0("Error abs. medio: ",round(error4, 3)),cex.main=1.5,
     ylab=expression(hat(italic(y))[test]), xlim=c(-2.3,2.3),ylim=c(-1.5,1.5))
abline(0,1,col="gray10", lwd=2)

## Pregunta 2 ####

### Parámetros de prueba de la muestra ####
ybar <- mean(ytrain)
NUCLEOS <- 11 ### Número de núcleos a usar con foreach

#### Modelo 1 ####
B <- 10000
i_n <- diag(1, n)
m1sig <- BETA1[,p+1]

##### Iniciar foreach ####
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)

ñ <- proc.time()
media1 <- foreach(s=1:B,.combine = "rbind")%dopar%{
  set.seed(100+s)
  muestra1 <- c(mvtnorm::rmvnorm(n = 1, mean = xtrain%*%BETA1[s,-(p+1)], sigma = i_n*m1sig[s]))
  media <- mean(muestra1)
  return(media)
}
stopCluster(cl)
proc.time()-ñ 

#### Modelo 2 ####
B <- 10000
i_n <- diag(1, n)
m2sig <- BETA2[,p+1]

##### Iniciar foreach ####
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)

ñ <- proc.time()
media2 <- foreach(s=1:B,.combine = "rbind")%dopar%{
  set.seed(100+s)
  muestra1 <- c(mvtnorm::rmvnorm(n = 1, mean = xtrain%*%BETA2[s,-(p+1)], sigma = i_n*m2sig[s]))
  media <- mean(muestra1)
  return(media)
}
stopCluster(cl)
proc.time()-ñ 

#### Modelo 3 ####
B <- 10000
i_n <- diag(1, n)
m3sig <- BETA3[,p+1]

##### Iniciar foreach ####
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)

ñ <- proc.time()
media3 <- foreach(s=1:B,.combine = "rbind")%dopar%{
  set.seed(100+s)
  muestra1 <- c(mvtnorm::rmvnorm(n = 1, mean = xtrain%*%BETA3[s,-c(p+1,p+2)], sigma = i_n*m3sig[s]))
  media <- mean(muestra1)
  return(media)
}
stopCluster(cl)
proc.time()-ñ 

#### Modelo 4 ####
B <- 5000
i_n <- diag(1, n)
DY <-abs(outer( (1:n),(1:n) ,"-"))
m4sig <- BETA4[,p+1]
m4rho <- BETA4[,p+2]

##### Iniciar foreach ####
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)

ñ <- proc.time()
media4 <- foreach(s=1:B,.combine = "rbind")%dopar%{
  set.seed(100+s)
  muestra1 <- c(mvtnorm::rmvnorm(n = 1, mean = xtrain%*%BETA4[s,-c(p+1,p+2)], sigma = m4rho[s]^DY*m4sig[s]))
  media <- mean(muestra1)
  return(media)
}
stopCluster(cl)
proc.time()-ñ 

### PPPs ####

cont1 <- 0
cont2 <- 0
cont3 <- 0
cont4 <- 0
B <- 10000

for(s in 1:B){
  cont1 <- ifelse(media1[s] < ybar, cont1 + 1, cont1)
  cont2 <- ifelse(media2[s] < ybar, cont2 + 1, cont2)
  cont3 <- ifelse(media3[s] < ybar, cont3 + 1, cont3)
  if (s <= 5000){
  cont4 <- ifelse(media4[s] < ybar, cont4 + 1, cont4)
  }
}

pppmedia <- c(cont1/B, cont2/B, cont3/B, cont4/(B/2));pppmedia

#### Gráfico predictiva posterior ####

par(mfrow = c(2, 2), oma = c(0,0, 2, 0)) 
# Grafico 1
valor_p <- cont1 / B
hist(media1,  col = "#9FFFE7",border = "#23BBC9", freq = F, main = paste("Distribución Predictiva Posterior- Modelo 1\n ppp media:", round(valor_p, 3)), xlab = "", ylab = "p(t|y)")
abline(v = mean(media1), col = "gray25", lwd = 3, lty=8)

# Grafico 2
valor_p <- cont2 / B
hist(media2, col = "#9FFFE7",border = "#23BBC9", freq = F, main = paste("Distribución Predictiva Posterior- Modelo 2\n ppp media:", round(valor_p, 3)), xlab = "", ylab = "p(t|y)")
abline(v = mean(media2), col = "gray25", lwd = 3,lty=8)

# Grafico 3
valor_p <- cont3 / B
hist(media3, col = "#9FFFE7",border = "#23BBC9", freq = F, main = paste("Distribución Predictiva Posterior- Modelo 3\n ppp media:", round(valor_p, 3)), xlab = "", ylab = "p(t|y)")
abline(v = mean(media3), col = "gray25", lwd = 3,lty=8)
legend_coord <- locator(1)
legend(x=legend_coord[1],y=legend_coord[2], legend = c("Media observada"), 
       col = "gray25", lty = 8, lwd = 2, bty = "n",xpd = TRUE, horiz = TRUE)
# Grafico 4
valor_p <- cont4 / (B/2)
hist(media4,  col = "#9FFFE7",border = "#23BBC9", freq = F, main = paste("Distribución Predictiva Posterior- Modelo 4\n ppp media:", round(valor_p, 3)),ylim = c(0,10), xlab = "", ylab = "p(t|y)")
abline(v = mean(media4), col = "gray25", lwd = 3,lty=8)
legend_coord <- locator(1)
legend(x=legend_coord[1],y=legend_coord[2], legend = c("Media observada"), 
       col = "gray25", lty = 8, lwd = 2, bty = "n",xpd = TRUE, horiz = TRUE)

## Pregunta 3 ####

### DIC #### 

#### Modelo 1 ####
lpyth_m1 <- mvtnorm::dmvnorm(x = ytrain, mean = xtrain%*%beta1hat, sigma = i_n*mean(BETA1[,p+1]), log = T)
pdic1 <- 2*(lpyth_m1-mean(LL1))
dic1 <- -2*lpyth_m1+2*pdic1

#### Modelo 2 ####
lpyth_m2 <- mvtnorm::dmvnorm(x = ytrain, mean = xtrain%*%beta2hat, sigma = i_n*mean(BETA2[,p+1]), log = T)
pdic2 <- 2*(lpyth_m2-mean(LL2))
dic2 <- -2*lpyth_m2+2*pdic2

#### Modelo 3 ####
lpyth_m3 <- mvtnorm::dmvnorm(x = ytrain, mean = xtrain%*%beta3hat, sigma = i_n*mean(BETA3[,p+1]), log = T)
pdic3 <- 2*(lpyth_m3-mean(LL3))
dic3 <- -2*lpyth_m3+2*pdic3

#### Modelo 4 ####
lpyth_m4 <- mvtnorm::dmvnorm(x = ytrain, mean = xtrain%*%beta4hat, sigma = mean(BETA4[,p+2])^DY*mean(BETA4[,p+1]), log = T)
pdic4 <- 2*(lpyth_m4-mean(LL4))
dic4 <- -2*lpyth_m4+2*pdic4


c(dic1, dic2, dic3, dic4)
