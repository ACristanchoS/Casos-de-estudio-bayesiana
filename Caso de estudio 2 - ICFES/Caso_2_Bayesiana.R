## Parcial 2 Bayesiana ##
setwd("C:/Users/EQUIPO/OneDrive/Documentos/Documentos Ander/UNAL/Estadística Bayesiana/Caso de estudio 2")
getwd()

library(readr)
library(dplyr)
library(readxl)
library(foreach)
library(doParallel)
library(xtable)
library(ggplot2)
library(sf)

missings <- function(x) return(sum(is.na(x)))

# Lectura de la base de datos 
icfes <- read.delim("SB11_20222.TXT", sep = ";")
icfes <- icfes %>% 
  filter(ESTU_NACIONALIDAD == "COLOMBIA",ESTU_PAIS_RESIDE == "COLOMBIA",ESTU_ESTADOINVESTIGACION == "PUBLICAR",COLE_COD_DEPTO_UBICACION != 88,!is.na(COLE_COD_MCPIO_UBICACION),!is.na(COLE_COD_DEPTO_UBICACION),!is.na(PUNT_GLOBAL)) %>% 
  select(PUNT_GLOBAL, COLE_DEPTO_UBICACION, COLE_COD_DEPTO_UBICACION, COLE_MCPIO_UBICACION, COLE_COD_MCPIO_UBICACION) %>%
  arrange(COLE_COD_DEPTO_UBICACION,COLE_COD_MCPIO_UBICACION)%>% as.data.frame()

y <- icfes$PUNT_GLOBAL

apply(icfes,2,missings)

pm <- as.data.frame(read_excel("pobreza monetaria.xls", sheet = 2,skip = 15,n_max = 25,))
colnames(pm)[1] <- "DPTO"
pm <- pm %>% 
  select(DPTO, "2018") %>%  as.data.frame()
COLE_COD_DEPTO_UBICACION <- c(5,8,11,13,15,17,18,19,20,27,23,25,41,44,47,50,52,54,63,66,68,70,73,76,NA)
pm <- cbind(pm,COLE_COD_DEPTO_UBICACION)

ed <- read.csv("estadísticas educación.csv")
ed <- ed %>% 
  filter(AÑO == "2022", CÓDIGO_DEPARTAMENTO != 88) %>% 
  select(CÓDIGO_MUNICIPIO, MUNICIPIO, CÓDIGO_DEPARTAMENTO, DEPARTAMENTO, COBERTURA_NETA_SECUNDARIA) %>% 
  arrange(CÓDIGO_MUNICIPIO, MUNICIPIO)%>% as.data.frame()

dptoshp<-st_read("MGN_DPTO_POLITICO.shp", quiet=TRUE)[-c(23),]
mpioshp <- st_read("MGN_MPIO_POLITICO.shp",quiet=TRUE)

## --------------------- Punto 1 ------------------------ ##
## Mapa de departamentos- Media global
mapdepto<-dptoshp %>% left_join(estadisticos, by=c("DPTO_CCDGO"="COLE_COD_DEPTO_UBICACION"))
box<-st_bbox(mapdepto)
ggplot()+
  geom_sf(data = mapdepto, aes(fill=yb), col="gray70", linetype="solid")+
  coord_sf(xlim=c(box$xmin, box$xmax), ylim=c(box$ymin,box$ymax),expand=F)+
  geom_sf_text(data=mapdepto,aes(label=round(yb, 2)),col="black",
               fontface="bold",size=2.5,fun.geometry=function(x) sf::st_centroid(x)) +
  scale_fill_gradient(low="#EBEBEB", high = "#FF024A")+
  labs(x="Longitud", y="Latitud", title="", fill="Media Global")+
  theme(panel.background = element_rect(fill = "white"))


## Mapa por departamento - Incidencia de pobreza
colnames(pm)<-c("Codigo", "Departamento", "2018", "Año")
mappobre <- dptoshp %>% left_join(pm, by = c("DPTO_CCDGO" = "COLE_COD_DEPTO_UBICACION"))
box <- st_bbox(mappobre)

ggplot() +
  geom_sf(data = mappobre, aes(fill = Año), col = "gray70", linetype = "solid") +
  coord_sf(xlim = c(box$xmin, box$xmax), ylim = c(box$ymin, box$ymax), expand = FALSE) +
  geom_sf_text(data = mappobre, aes(label = Año), col = "black",
               fontface = "bold", size = 2.5, fun.geometry = function(x) sf::st_centroid(x)) +
  scale_fill_gradient(low = "#EBEBEB", high = "#FF024A") +
  labs(x = "Longitud", y = "Latitud", title = "", fill = "Incidencia de \n pobreza (%)") +
  theme(panel.background = element_rect(fill = "white"))



## ------------------- PUNTO 2 --------------------- #
## Mapa de municipios - Puntaje global
mapmun <- mpioshp %>% left_join(estadisticosmcpio, by = c("MPIO_CCNCT" = "COLE_COD_MCPIO_UBICACION"))
box <- st_bbox(mapmun)
ggplot()+
  geom_sf(data = mapmun, aes(fill=yb), col="gray50", linetype="solid")+
  coord_sf(xlim=c(box$xmin, box$xmax), ylim=c(box$ymin,box$ymax),expand=F)+
  scale_fill_gradient(low="#FDEBE8", high = "#6F09FD")+
  labs(x="Longitud", y="Latitud", title="", fill="Media Global")+
  theme(panel.background = element_rect(fill = "white"))

## Mapa de muncipios cobertura neta en secundaria
ed$CÓDIGO_MUNICIPIO<-as.character(ed$CÓDIGO_MUNICIPIO)
mapcbsn <- mpioshp %>% left_join(cbns, by = c("MPIO_CCNCT" = "CÓDIGO_MUNICIPIO"))
box <- st_bbox(mapcbsn)
ggplot()+
  geom_sf(data = mapcbsn, aes(fill=COBERTURA_NETA_SECUNDARIA), col="gray50", linetype="solid")+
  coord_sf(xlim=c(box$xmin, box$xmax), ylim=c(box$ymin,box$ymax),expand=F)+
  scale_fill_gradient(low="#FDEBE8", high = "#6F09FD")+
  labs(x="Longitud", y="Latitud", title="", fill="Cobertura neta \n secundaria")+
  theme(panel.background = element_rect(fill = "white"))


# -------------------- MODELOS ---------------------- #

### Estadísticos suficientes ###
n <- length(y)
mean_y <- mean(icfes$PUNT_GLOBAL)
sum_y <- sum(icfes$PUNT_GLOBAL)

### Hiperparámetros ###
mu0 = 250
g20 = 50^2
nu0 = 1
s20 = 50^2

### Número de iteraciones
B <- 101000
S <- 10000 # Muestreo sistemático
cont <- 0

### Cada 10% anuncie
ncat <- floor(B/10)

### Matriz MCMC
THETA1 <- matrix(data = NA, nrow = S, ncol = 2)
LL1    <- matrix(data = NA, nrow = S, ncol = 1)

### Muestreador de Gibbs Modelo 1
set.seed(1010)
sig2 <- 1/rgamma(n = 1, shape = nu0/2, rate = nu0*s20/2)

### Para optimizar el código
op1 <- c(1/g20,mu0/g20,nu0*s20,(n-1)*var(y))
nun   <- nu0 + n

# cadena
ñ <- proc.time() 
set.seed(1010)
for(b in 1:B) {
  # actualizar theta
  t2n   <- 1/(op1[1] + n/sig2)      
  mun   <- t2n*(op1[2] + sum_y/sig2)
  theta <- rnorm(n = 1, mean = mun, sd = sqrt(t2n))
  # actualizar sigma^2
  s2n   <- (op1[3] + op1[4]+n*(mean_y-theta)^2)/nun
  sig2 <- 1/rgamma(n = 1, shape = nun/2, rate = nun*s2n/2)
  # almacenar
  if(b > 1000 && b%%10 == 0){
    cont = cont + 1
  THETA1[cont,] <- c(theta, sig2)
  LL1[cont] <- sum(dnorm(x = y, mean = theta, sd = sqrt(sig2), log = T))
  }
  # progreso
  if (b%%ncat == 0) 
    cat(100*round(b/B, 1), "% completado ... \n", sep = "")
}
colnames(THETA1) <- c("theta", "sig2")
colnames(LL1) <- c("ll1")
THETA1 <- as.data.frame(THETA1)
LL1    <- as.data.frame(LL1)
proc.time()-ñ 

neff_THETA1 <- coda::effectiveSize(THETA1)
EMC1 <- round(apply(X = THETA1, MARGIN = 2, FUN = sd)/sqrt(neff_THETA1), 3)
100*abs(EMC1/colMeans(THETA1))

# Muestreador de Gibbs Modelo 2
m <- 32
Y <- vector(mode = "list", length = m)
g <- rep(NA, n)
for (j in 1:m) {
  idx <- icfes$COLE_COD_DEPTO_UBICACION == unique(icfes$COLE_COD_DEPTO_UBICACION)[j]
  g[idx] <- j
  Y[[j]] <- y[idx]
}
# Estadísticos para cada departamento
estadisticos <- icfes %>% 
  group_by(COLE_COD_DEPTO_UBICACION) %>% 
  summarise(COLE_COD_DEPTO_UBICACION = unique(COLE_COD_DEPTO_UBICACION), 
            COLE_DEPTO_UBICACION = unique(COLE_DEPTO_UBICACION), 
            nj = n(), 
            yb = mean(PUNT_GLOBAL), 
            s2 = var(PUNT_GLOBAL))

nj <- estadisticos$nj
yb <- estadisticos$yb
s2 <- estadisticos$s2

# Hiperparámetros
mu0  <- 250
g20  <- 50^2
eta0 <- 1
t20  <- 50^2
nu0  <- 1
s20  <- 50^2
cont <- 0

## Valores iniciales ##
theta <- yb
sig2  <- mean(s2)
mu    <- mean(theta)
tau2  <- var(theta)
# almacenamiento
THETA2 <- matrix(data = NA, nrow = S, ncol = m+3)
LL2    <- matrix(data = NA, nrow = S, ncol = 1)

ñ <- proc.time()
set.seed(1010)
# cadena
for (b in 1:B) {
  # actualizar theta
  vtheta <- 1/(1/tau2 + nj/sig2)
  theta  <- rnorm(n = m, mean = vtheta*(mu/tau2 + nj*yb/sig2), sd = sqrt(vtheta))
  # actualizar sigma^2
  sig2 <- 1/rgamma(n = 1, shape = 0.5*(nu0 + n), rate = 0.5*(nu0*s20 + sum((nj-1)*s2 + nj*(yb - theta)^2)))
  # actualizar mu
  vmu <- 1/(1/g20 + m/tau2)
  mu  <- rnorm(n = 1, mean = vmu*(mu0/g20 + m*mean(theta)/tau2), sd = sqrt(vmu)) 
  # actualizar tau^2
  tau2 <- 1/rgamma(n = 1, shape = 0.5*(eta0 + m), rate = 0.5*(eta0*t20 + (m-1)*var(theta) + m*(mean(theta) - mu)^2))
  # almacenar valores
  
  
  if(b > 1000 && b%%10 == 0){
    cont = cont + 1
    THETA2[cont,] <- c(theta, sig2, mu, tau2)
    LL2[cont] <- sum(dnorm(x = y, mean = rep(theta, nj), sd = sqrt(sig2), log = T))
  }
  
  if (b%%ncat == 0) 
    cat(100*round(b/B, 1), "% completado ... \n", sep = "")
}
proc.time()-ñ 

colnames(THETA2) <- c(paste0("theta",1:m), "sig2", "mu", "tau2")
colnames(LL2) <- c("ll2")
THETA2 <- as.data.frame(THETA2)
LL2    <- as.data.frame(LL2)

neff_THETA2 <- coda::effectiveSize(THETA2)
EMC2 <- round(apply(X = THETA2, MARGIN = 2, FUN = sd)/sqrt(neff_THETA2), 3)
CV2 <- 100*abs(EMC2/colMeans(THETA2))
summary(CV2[1:m])
CV2[m+(1:3)]


# Muestreador de Gibbs Modelo 3

# Hiperparámetros
mu0  <- 250 
g20  <- 50^2
eta0 <- 1  
t20  <- 50^2
nu <- 1  
al0  <- 1
be0  <- 1/50^2 
cont <- 0

# valores iniciales
theta <- yb
sig2  <- s2  # sigma_j^2
mu    <- mean(theta)
tau2  <- var(theta)
ups2  <- 500^2  # sigma^2

THETA3 <- matrix(data = NA, nrow = S, ncol = 2*m+3)
LL3    <- matrix(data = NA, nrow = S, ncol = 1)

# cadena
ñ <- proc.time() 
set.seed(1010)
for (b in 1:B) {
  # actualizar theta
  vtheta <- 1/(1/tau2 + nj/sig2)
  theta  <- rnorm(n = m, mean = vtheta*(mu/tau2 + nj*yb/sig2), sd = sqrt(vtheta))
  # actualizar sigma_j^2
  sig2 <- 1/rgamma(n = m, shape = 0.5*(nu + nj), rate = 0.5*(nu*ups2 + (nj-1)*s2 + nj*(yb - theta)^2))
  # actualizar mu
  vmu <- 1/(1/g20 + m/tau2)
  mu  <- rnorm(n = 1, mean = vmu*(mu0/g20 + m*mean(theta)/tau2), sd = sqrt(vmu))
  # actualizar tau2
  tau2 <- 1/rgamma(n = 1, shape = 0.5*(eta0+m), rate = 0.5*(eta0*t20 + (m-1)*var(theta) + m*(mean(theta) - mu)^2))
  # actualizar sigma^2
  ups2 <- rgamma(n = 1, shape = 0.5*(al0 + m*nu), rate = 0.5*(be0 + nu*sum(1/sig2)))
  
  if(b > 1000 && b%%10 == 0){
    cont = cont + 1
    THETA3[cont,] <- c(theta, sig2, mu, tau2, ups2)
    LL3[cont] <- sum(dnorm(x = y, mean = rep(theta, nj), sd = sqrt(rep(sig2, nj)), log = T))
  }
  if (b%%ncat == 0) 
    cat(100*round(b/B, 1), "% completado ... \n", sep = "")
}
proc.time()-ñ 

colnames(THETA3) <- c(paste0("theta",1:m), paste0("sig2_",1:m), "mu", "tau2","ups2")
colnames(LL3) <- c("ll3")
THETA3 <- as.data.frame(THETA3)
LL3    <- as.data.frame(LL3)

neff_THETA3 <- coda::effectiveSize(THETA3)
EMC3 <- round(apply(X = THETA3, MARGIN = 2, FUN = sd)/sqrt(neff_THETA3), 3)
CV3 <- 100*abs(EMC3/colMeans(THETA3))
summary(CV3[1:m])
summary(CV3[m+(1:m)])
CV3[2*m+(1:3)]

# Muestreador de Gibbs Modelo 4
o <- 1112
YM <- vector(mode = "list", length = o)
h <- rep(NA, n)
for (j in 1:o) {
  idx <- icfes$COLE_COD_MCPIO_UBICACION == unique(icfes$COLE_COD_MCPIO_UBICACION)[j]
  h[idx] <- j
  YM[[j]] <- y[idx]
}
# Estadísticos para cada municipio
estadisticosmcpio <- icfes %>% 
  group_by(COLE_COD_MCPIO_UBICACION) %>% 
  summarise(COLE_COD_MCPIO_UBICACION = unique(COLE_COD_MCPIO_UBICACION),
           COLE_MCPIO_UBICACION = unique(COLE_MCPIO_UBICACION),
           njm = n(), 
           ybm = mean(PUNT_GLOBAL), 
           s2m = var(PUNT_GLOBAL),
           minm = min(PUNT_GLOBAL),
           maxm = max(PUNT_GLOBAL),
           ricm = quantile(PUNT_GLOBAL, 0.75)-quantile(PUNT_GLOBAL, 0.25),
           medianam = quantile(PUNT_GLOBAL, 0.5))

pnummcpio <- icfes %>%
  group_by(COLE_COD_DEPTO_UBICACION, COLE_MCPIO_UBICACION) %>% 
  summarise(ksy = n())

nmcpios <- pnummcpio %>%
  group_by(COLE_COD_DEPTO_UBICACION) %>% 
  summarise(nm = n())
  

minm <- estadisticosmcpio$minm
maxm <- estadisticosmcpio$maxm
ricm <- estadisticosmcpio$ricm
medianam <- estadisticosmcpio$medianam
nm <- nmcpios$nm
grupos <- rep(1:m,nm)
njm <- estadisticosmcpio$njm
ybm <- estadisticosmcpio$ybm
s2m <- estadisticosmcpio$s2m
yb <- tapply(ybm, grupos, mean)
s2  <- tapply(ybm, grupos, var)
s2[3] <- 0


## Hiperparámetros
xi0 <- 1
k20 <- 50^2
mu0  <- 250
g20  <- 50^2
eta0 <- 1
t20  <- 50^2
nu0  <- 1
s20  <- 50^2
cont <- 0

# valores iniciales
zeta <- ybm
kap2 <- mean(s2m)
theta <- yb
sig2  <- mean(s2)
mu    <- mean(theta)
tau2  <- var(theta)

THETA4 <- matrix(data = NA, nrow = S, ncol = o+m+4)
LL4    <- matrix(data = NA, nrow = S, ncol = 1)

ñ <- proc.time()
set.seed(1010)
# cadena
for (b in 1:B) {
  # actualizar zeta
  vzeta <- 1/(1/sig2+njm/kap2)
  zeta <- rnorm(n = o, mean = vzeta*(rep(theta,nm)/sig2+njm*ybm/kap2), sd = sqrt(vzeta))
  # actualizar kappa
  kap2 <- 1/rgamma(n = 1, shape = 0.5*(xi0+n), rate = 0.5*(xi0*k20+sum((njm-1)*s2m+njm*(ybm - zeta)^2)))
  # actualizar theta
  vtheta <- 1/(1/tau2 + nm/sig2)
  theta  <- rnorm(n = m, mean = vtheta*(mu/tau2 + nm*tapply(zeta, grupos, mean)/sig2), sd = sqrt(vtheta))
  # actualizar sigma
  sig2 <- 1/rgamma(n = 1, shape = 0.5*(nu0 + o), rate = 0.5*(nu0*s20 + sum((zeta-rep(theta,nm))^2)))
  # actualizar mu
  vmu <- 1/(1/g20 + m/tau2)
  mu  <- rnorm(n = 1, mean = vmu*(mu0/g20 + m*mean(theta)/tau2), sd = sqrt(vmu)) 
  # actualizar tau^2
  tau2 <- 1/rgamma(n = 1, shape = 0.5*(eta0 + m), rate = 0.5*(eta0*t20 + (m-1)*var(theta) + m*(mean(theta) - mu)^2))
  # almacenar valores
  
  if(b > 1000 && b%%10 == 0){
    cont = cont + 1
    THETA4[cont,] <- c(zeta, theta, kap2, sig2, mu, tau2)
    LL4[cont] <- sum(dnorm(x = y, mean = rep(zeta, njm), sd = sqrt(kap2), log = T)) 
  }
  
  if (b%%ncat == 0) 
    cat(100*round(b/B, 1), "% completado ... \n", sep = "")
}
proc.time()-ñ 

colnames(THETA4) <- c(paste0("zeta",1:o), paste0("theta",1:m), "kap2", "sig2", "mu","tau2")
colnames(LL4) <- c("ll4")
THETA4 <- as.data.frame(THETA4)
LL4    <- as.data.frame(LL4)

neff_THETA4 <- coda::effectiveSize(THETA4)
EMC4 <- round(apply(X = THETA4, MARGIN = 2, FUN = sd)/sqrt(neff_THETA4), 3)
CV4 <- 100*abs(EMC4/colMeans(THETA4))
summary(CV4[1:o])
summary(CV4[o+(1:m)])
CV4[o+m+(1:4)]

# Muestreador de Gibbs Modelo 5

# Hiperparámetros
xi0 <- 1
k20 <- 50^2
mu0  <- 250
g20  <- 50^2
eta0 <- 1
t20  <- 50^2
nu  <- 1
al0  <- 1
be0  <- 1/50^2 
cont <- 0

# valores iniciales
zeta <- ybm
kap2 <- mean(s2m)
theta <- yb
sig2  <- s2
mu    <- mean(theta)
tau2  <- var(theta)
ups2  <- 500^2

THETA5 <- matrix(data = NA, nrow = S, ncol = o+2*m+4)
LL5    <- matrix(data = NA, nrow = S, ncol = 1)

ñ <- proc.time()
set.seed(1010)
# cadena
for (b in 1:B) {
  # actualizar sigma^2_j
  sig2 <- 1/rgamma(n = m, shape = 0.5*(nu + nm), rate = 0.5*(nu*ups2 + tapply((zeta-rep(theta,nm))^2,grupos,sum)))
  # actualizar zeta
  vzeta <- 1/(1/rep(sig2,nm)+njm/kap2)
  zeta <- rnorm(n = o, mean = vzeta*(rep(theta,nm)/rep(sig2,nm)+njm*ybm/kap2), sd = sqrt(vzeta))
  # actualizar kappa
  kap2 <- 1/rgamma(n = 1, shape = 0.5*(xi0+n), rate = 0.5*(xi0*k20+sum((njm-1)*s2m+njm*(ybm - zeta)^2)))
  # actualizar theta
  vtheta <- 1/(1/tau2 + nm/sig2)
  theta  <- rnorm(n = m, mean = vtheta*(mu/tau2 + nm*tapply(zeta, grupos, mean)/sig2), sd = sqrt(vtheta))
  # actualizar mu
  vmu <- 1/(1/g20 + m/tau2)
  mu  <- rnorm(n = 1, mean = vmu*(mu0/g20 + m*mean(theta)/tau2), sd = sqrt(vmu)) 
  # actualizar tau^2
  tau2 <- 1/rgamma(n = 1, shape = 0.5*(eta0 + m), rate = 0.5*(eta0*t20 + (m-1)*var(theta) + m*(mean(theta) - mu)^2))
  # actualizar sigma^2
  ups2 <- rgamma(n = 1, shape = 0.5*(al0 + m*nu), rate = 0.5*(be0 + nu*sum(1/sig2)))
  # almacenar valores
  
  if(b > 1000 && b%%10 == 0){
    cont = cont + 1
    THETA5[cont,] <- c(zeta, theta, sig2, kap2, mu, tau2, ups2)
    LL5[cont] <- sum(dnorm(x = y, mean = rep(zeta, njm), sd = sqrt(kap2), log = T)) 
  }
  
  if (b%%ncat == 0) 
    cat(100*round(b/B, 1), "% completado ... \n", sep = "")
}
proc.time()-ñ 

colnames(THETA5) <- c(paste0("zeta",1:o), paste0("theta",1:m), paste0("sig2_",1:m), "kap2", "mu","tau2", "ups2")
colnames(LL5) <- c("ll5")
THETA5 <- as.data.frame(THETA5)
LL5    <- as.data.frame(LL5)

neff_THETA5 <- coda::effectiveSize(THETA5)
EMC5 <- round(apply(X = THETA5, MARGIN = 2, FUN = sd)/sqrt(neff_THETA5), 3)
CV5 <- 100*abs(EMC5/colMeans(THETA5))
summary(CV5[1:o])
summary(CV5[o+(1:m)])
summary(CV5[o+m+(1:m)])
CV5[o+2*m+(1:4)]

## ------------------ PUNTO 5 -------------------- #

## DIC para los 5 modelos
tbay1 <- colMeans(THETA1)
tbay2 <- colMeans(THETA2)
tbay3 <- colMeans(THETA3)
tbay4 <- colMeans(THETA4)
tbay5 <- colMeans(THETA5)

lpyth_m1 <- sum(dnorm(x = y, mean = tbay1[1], sd = sqrt(tbay1[2]), log = T))
pdic1 <- 2*(lpyth_m1-mean(LL1$ll1))
dic1 <- -2*lpyth_m1+2*pdic1

lpyth_m2 <- sum(dnorm(x = y, mean = rep(tbay2[1:32],nj), sd = sqrt(tbay2[33]), log = T))
pdic2 <- 2*(lpyth_m2-mean(LL2$ll2))
dic2 <- -2*lpyth_m2+2*pdic2

lpyth_m3 <- sum(dnorm(x = y, mean = rep(tbay3[1:32],nj), sd = rep(sqrt(tbay3[33:64]),nj), log = T))
pdic3 <- 2*(lpyth_m3-mean(LL3$ll3))
dic3 <- -2*lpyth_m3+2*pdic3

lpyth_m4 <- sum(dnorm(x = y, mean = rep(tbay4[1:1112],njm), sd = sqrt(tbay4[1145]), log = T))
pdic4 <- 2*(lpyth_m4-mean(LL4$ll4))
dic4 <- -2*lpyth_m4+2*pdic4

lpyth_m5 <- sum(dnorm(x = y, mean = rep(tbay5[1:1112],njm), sd = sqrt(tbay5[1177]), log = T))
pdic5 <- 2*(lpyth_m5-mean(LL5$ll5))
dic5 <- -2*lpyth_m5+2*pdic5

c(dic1, dic2, dic3, dic4, dic5)

# WAIC
WAIC <- c(0,0)
NUCLEOS <- 11 ### Número de núcleos a usar con foreach


### WAIC Modelo 1
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)

ñ <- proc.time()
WAICM1 <- foreach(i = 1:n,.combine = "+")%dopar%{
  # lppd
  tmp1 <- dnorm(x = y[i], mean = THETA1$theta, sd = sqrt(THETA1$sig2))
  WAIC[1] <- log(mean(tmp1))
  # pWAIC
  tmp2 <- dnorm(x = y[i], mean = THETA1$theta, sd =  sqrt(THETA1$sig2), log = T)
  WAIC[2] <- 2*(log(mean(tmp1)) - mean(tmp2))
  return(WAIC)
}
proc.time()-ñ 
stopCluster(cl) 

waic_m1 <- -2*WAICM1[1] + 2*WAICM1[2]

## WAIC Modelo 2
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)

ñ <- proc.time()
WAICM2 <- foreach(i = 1:n,.combine = "+")%dopar%{
  # lppd
  tmp1 <- dnorm(x = y[i], mean = THETA2[,g[i]], sd = sqrt(THETA2$sig2))
  WAIC[1] <- log(mean(tmp1))
  # pWAIC
  tmp2 <- dnorm(x = y[i], mean = THETA2[,g[i]], sd =  sqrt(THETA2$sig2), log = T)
  WAIC[2] <- 2*(log(mean(tmp1)) - mean(tmp2))
  return(WAIC)
}
proc.time()-ñ 
stopCluster(cl) 

waic_m2 <- -2*WAICM2[1] + 2*WAICM2[2]

## WAIC Modelo 3
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)

ñ <- proc.time()
WAICM3 <- foreach(i = 1:n,.combine = "+")%dopar%{
     # lppd
     tmp1 <- dnorm(x = y[i], mean = THETA3[,g[i]], sd = sqrt(THETA3[,m+g[i]]))
     WAIC[1] <- log(mean(tmp1))
     # pWAIC
     tmp2 <- dnorm(x = y[i], mean = THETA3[,g[i]], sd =  sqrt(THETA3[,m+g[i]]), log = T)
     WAIC[2] <- 2*(log(mean(tmp1)) - mean(tmp2))
     return(WAIC)
}
proc.time()-ñ 
stopCluster(cl) 

waic_m3 <- -2*WAICM3[1] + 2*WAICM3[2]

## WAIC Modelo 4
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)

ñ <- proc.time()
WAICM4 <- foreach(i = 1:n,.combine = "+")%dopar%{
    # lppd
    tmp1 <- dnorm(x = y[i], mean = THETA4[,h[i]], sd = sqrt(THETA4$kap2))
    WAIC[1] <- log(mean(tmp1))
    # pWAIC
    tmp2 <- dnorm(x = y[i], mean = THETA4[,h[i]], sd =  sqrt(THETA4$kap2), log = T)
    WAIC[2] <- 2*(log(mean(tmp1)) - mean(tmp2))
    return(WAIC)
}
proc.time()-ñ 
stopCluster(cl) 

waic_m4 <- -2*WAICM4[1] + 2*WAICM4[2]

## WAIC Modelo 5
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)

ñ <- proc.time()
WAICM5 <- foreach(i = 1:n,.combine = "+")%dopar%{
  # lppd
  tmp1 <- dnorm(x = y[i], mean = THETA5[,h[i]], sd = sqrt(THETA5$kap2))
  WAIC[1] <- log(mean(tmp1))
  # pWAIC
  tmp2 <- dnorm(x = y[i], mean = THETA5[,h[i]], sd =  sqrt(THETA5$kap2), log = T)
  WAIC[2] <- 2*(log(mean(tmp1)) - mean(tmp2))
  return(WAIC)
}
proc.time()-ñ 
stopCluster(cl) 

waic_m5 <- -2*WAICM5[1] + 2*WAICM5[2]

waics <- c(waic_m1, waic_m2, waic_m3, waic_m4, waic_m5)
waics

## ---------------- PUNTO 6 ----------------- ##

m1mu <- c(mean(THETA1$theta),quantile(THETA1$theta, c(0.025,0.975)))
m2mu <- c(mean(THETA2$mu),quantile(THETA2$mu, c(0.025,0.975)))
m3mu <- c(mean(THETA3$mu),quantile(THETA3$mu, c(0.025,0.975)))
m4mu <- c(mean(THETA4$mu),quantile(THETA4$mu, c(0.025,0.975)))
m5mu <- c(mean(THETA5$mu),quantile(THETA5$mu, c(0.025,0.975)))

xtable(rbind(m1mu, m2mu, m3mu, m4mu, m5mu))

## --------------- PUNTO 7 ------------------ ##
#ranking bayesiano
THETA <- THETA5[, 1113:1144]
ids2 <- estadisticos$COLE_DEPTO_UBICACION
that <- colMeans(THETA[, 1:m])
ic1 <- apply(X = THETA[, 1:m], MARGIN = 2, FUN = function(x) quantile(x, c(0.025, 0.975)))
ranking <- order(that)
ids2 <- ids2[ranking]
that <- that[ranking]
ic1 <- ic1[, ranking]
colo <- rep(2, m)
colo[which(ic1[1,] > 250)] <- 1
colo[which(ic1[2,] < 250)] <- 3
colo <- c("#00C603", "black", "#EF0000")[colo]

# Ajusta el tamaño del gráfico para mostrar todas las etiquetas
par(mar = c(3, 6, 2, 2))

plot(NA, NA, xlab = "Puntaje", ylab = "", main = "Ranking Bayesiano: Modelo 5", xlim = c(190, 300), ylim = c(1, m), cex.axis = 0.7,  yaxt = "n")
axis(side = 2, at = 1:m, labels = ids2, las = 1, cex.axis = 0.5)  # Ajusta el tamaño de los labels con cex.axis
# Las etiquetas se mostrarán horizontalmente

abline(v = 250, col = "gray", lwd = 3)
abline(h = 1:m, col = "lightgray", lwd = 1)

for (j in 1:m) {
  segments(x0 = ic1[1, j], y0 = j, x1 = ic1[2, j], y1 = j, col = colo[j])
  lines(x = that[j], y = j, type = "p", pch = 16, cex = 0.8, col = colo[j])
}


# Raking frecuentista
yb<-estadisticos$yb
s2<-estadisticos$s2
nj<-estadisticos$nj
# ranking frecuentista
ids2 <- estadisticos$COLE_DEPTO_UBICACION
that <- yb
ic1  <- NULL
for (j in 1:m){
  ic1  <- cbind(ic1, yb[j] + c(-1,1)*qnorm(p = 0.975)*sqrt(s2[j])/sqrt(nj[j]))}
ranking <- order(that) 
ids2 <- ids2[ ranking]
that <- that[ ranking]
ic1  <- ic1 [,ranking]
colo <- rep(2,m)
colo[which(ic1[1,]>250)] <- 1
colo[which(ic1[2,]<250)] <- 3
colo <- c("#00C603", "black", "#EF0000")[colo]
par(mar = c(3, 6, 2, 2))

# gráfico
plot(NA, NA, xlab = "Puntaje", ylab = "", main = "Ranking Frecuentista", xlim = c(190,300), ylim = c(1,m), cex.axis = 0.75, yaxt = "n")
axis(side = 2, at = 1:m, labels = ids2, las = 1, cex.axis = 0.5) 
abline(v = 250,  col = "gray", lwd = 3)
abline(h = 1:m, col = "lightgray", lwd = 1)
for (j in 1:m) {
  segments(x0 = ic1[1,j], y0 = j, x1 = ic1[2,j], y1 = j, col = colo[j])
  lines(x = that[j], y = j, type = "p", pch = 16, cex = 0.8, col = colo[j])
}
## ----#-------------- PUNTO 8 --------------------------- #

thetas_dep7 <- as.matrix(THETA5[,1113:1144])
H <- 5

xi <- matrix(data = NA, nrow = S, ncol = m)

set.seed(1010)
for (b in 1:S) {
  tmp <- stats::kmeans(x = thetas_dep7[b,], centers = H)
  xi[b,] <- tmp$cluster
}

deps <- unique(icfes$COLE_DEPTO_UBICACION)


A <- matrix(data = 0, nrow = m, ncol = m)
for (b in 1:S) {
  for (i in 1:(m-1)) {
    for (j in (i+1):m) {
      if (xi[b,i] == xi[b,j]) {
        A[i,j] <- A[i,j] + 1/S
      } 
    }
  }
}

A <- A + t(A)
diag(A) <- 1
colnames(A) <- deps
rownames(A) <- deps

thetapos <- colMeans(thetas_dep7)
indices <- rev(order(thetapos))
# se organizan las observaciones de acuerdo a la partición verdadera
A <- A[indices,indices]
names <- colnames(A)
# visualización de la matriz de incidencia
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
corrplot::corrplot(corr = A, is.corr = FALSE, addgrid.col = NA, method = "color", tl.pos = "lt", tl.cex = 0.7, tl.col = "black", 
                   title = "Matriz de incidencia por departamentos", mar=c(0,0,1,0))


# Mapa de cluster para departamentos

set.seed(1010)
kmeans <- kmeans(x = thetapos, centers = 5)
cluskmeans <- kmeans$cluster
estadisticos<-cbind(estadisticos,cluskmeans)
mapdepto<-dptoshp %>% left_join(estadisticos, by=c("DPTO_CCDGO"="COLE_COD_DEPTO_UBICACION"))
box<-st_bbox(mapdepto)
ggplot() +
  geom_sf(data = mapdepto, aes(fill = factor(cluskmeans)), col = "gray30", linetype = "solid") +
  coord_sf(xlim = c(box$xmin, box$xmax), ylim = c(box$ymin, box$ymax), expand = F) +
  scale_fill_manual(values = c("#FF8888", "#A6C5FF","#F52FFF","#F0FF48", "#90FF71"),
                    breaks = c("1", "2", "3", "4", "5"),
                    labels = c("Grupo 1", "Grupo 2", "Grupo 3", "Grupo 4", "Grupo 5")) +
  labs(x = "Longitud", y = "Latitud", title = "", fill = " ") +
  theme(panel.background = element_rect(fill = "white"))

## ------------------------ PUNTO 9 ------------------------- #

depsna <- c(81, 85, 86, 91, 94, 95, 97, 99)

depsx <- as.matrix(THETA5[1113:1136])
depsy <- as.matrix(THETA5[1137:1144])

dpreg1 <- matrix(data = NA, nrow = S, ncol = length(depsna))

yreg1 <- pm$`2018`

for(b in 1:S){
  Xm <- as.data.frame(cbind(yreg1,as.vector(depsx[b,])))
  colnames(Xm) <- c("2018","theta")
  reg <- lm(`2018` ~ 1 + theta, data = Xm)
  dpreg1[b,] <- reg$coefficient[1]+reg$coefficient[2]*depsy[b,]
}
dpreg1 <- as.data.frame(dpreg1)
colnames(dpreg1) <- depsna
tbayreg1 <- colMeans(dpreg1)
tab9 <- t(rbind(tbayreg1,apply(dpreg1,2, FUN = function(x)quantile(x,0.025)),apply(dpreg1,2, FUN = function(x)quantile(x,0.975))))
tab9 <- as.data.frame(tab9)
xtable(arrange(tab9,desc(tbayreg1)))

## -------------------- PUNTO 10 ------------------------ ##
zetas_10 <- as.matrix(THETA5[,1:1112])
H <- 8

xi <- matrix(data = NA, nrow = S, ncol = o)

set.seed(1010)
for (b in 1:S) {
  tmp <- stats::kmeans(x = zetas_10[b,], centers = H)
  xi[b,] <- tmp$cluster
}

mcpos <- estadisticosmcpio$COLE_MCPIO_UBICACION

ñ <- proc.time()

C <- matrix(data = 0, nrow = o, ncol = o)
for(b in 1:S) {
  for (i in 1:(o-1)) {
    for (j in (i+1):o) {
      if (xi[b,i] == xi[b,j]) {
        C[i,j] <- C[i,j] + 1/S
      } 
    }
  }
}
proc.time()-ñ 

C <- C + t(C)
diag(C) <- 1
colnames(C) <- mcpos
rownames(C) <- mcpos
zetapos <- colMeans(zetas_10)
indices2 <- order(zetapos)
# se organizan las observaciones de acuerdo a la partición verdadera
C <- C[indices2,indices2]
# visualización de la matriz de incidencia
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
corrplot::corrplot(corr = C, is.corr = FALSE, addgrid.col = NA, method = "color", tl.pos = "n", 
                   title = "Matriz de incidencia por municipios", mar=c(0,0,1,0))

# Mapa de cluster municipios
zetapos <- colMeans(zetas_10)
set.seed(1010)
kmeans <- kmeans(x = zetapos, centers = 8)
cluskmeans1 <- kmeans$cluster
estadisticosmcpio<-cbind(estadisticosmcpio,cluskmeans1)
mapmun <- mpioshp %>% left_join(estadisticosmcpio, by = c("MPIO_CCNCT" = "COLE_COD_MCPIO_UBICACION"))
mapmun$cluskmeans1<-as.factor(mapmun$cluskmeans1)
box <- st_bbox(mapmun)

colores <- c("#8DD3C7","#FFFFB3", "#FB8072", "#FDB462", "#FCCDE5", "#BC80BD", "#FFED6F","#20DE8B")
ggplot() +
  geom_sf(data = mapmun, aes(fill = cluskmeans1), col = "#717171", linetype = "solid") +
  coord_sf(xlim = c(box$xmin, box$xmax), ylim = c(box$ymin, box$ymax), expand = FALSE) +
  scale_fill_manual(values = colores, breaks = c("1", "2", "3", "4", "5", "6", "7", "8"),
                    labels = c("Grupo 1", "Grupo 2", "Grupo 3", "Grupo 4", "Grupo 5", "Grupo 6", "Grupo 7", "Grupo 8")) +
  labs(x = "Longitud", y = "Latitud", title = "", fill = " ") +
  theme(panel.background = element_rect(fill = "white"))

## ---------------- PUNTO 11 -------------------- #

mcps <- estadisticosmcpio %>% 
  left_join(ed, by = c("COLE_COD_MCPIO_UBICACION" = "CÓDIGO_MUNICIPIO")) %>% 
  as.data.frame()
which(is.na(mcps$COBERTURA_NETA_SECUNDARIA))
mcps[1098,]
mcps[582,]
mcpsna <- c(27086,94663)

mcpsx <- as.matrix(THETA5[1:1112])
mcpsx <- mcpsx[,-c(582,1098)]
mcpsy <- as.matrix(THETA5[,c(582,1098)])
mcpreg1 <- matrix(data = NA, nrow = S, ncol = length(mcpsna))
yreg2 <- mcps$COBERTURA_NETA_SECUNDARIA

for (i in 1:1112) {
if(is.na(yreg2[i]) == TRUE){
    yreg2 <- yreg2[-i]
}
}

for(b in 1:S){
  Xm <- as.data.frame(cbind(yreg2,as.vector(mcpsx[b,])))
  colnames(Xm) <- c("cob","zeta")
  reg <- lm(cob ~ 1 + zeta, data = Xm)
  mcpreg1[b,] <- reg$coefficient[1]+reg$coefficient[2]*mcpsy[b,]
}

mcpreg1 <- as.data.frame(mcpreg1)
colnames(mcpreg1) <- mcpsna
tbayreg2 <- colMeans(mcpreg1)
tab11 <- t(rbind(tbayreg2,apply(mcpreg1,2, FUN = function(x)quantile(x,0.025)),apply(mcpreg1,2, FUN = function(x)quantile(x,0.975))))
tab11 <- as.data.frame(tab11)
xtable(arrange(tab11,desc(tbayreg2)))


## ---------------- PUNTO 12 --------------------- ##

### Parámetros de prueba de las muestras
kappa_12 <- as.matrix(THETA5$kap2)
mcpogp <- rep(1:o,njm)

### Mínimo
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)
ñ <- proc.time()

minimo <- foreach(i=1:S,.combine = "rbind")%dopar%{
  set.seed(100+i)
  muestra1 <- rnorm(n = n, mean = rep(zetas_10[i,],njm), sd = sqrt(kappa_12[i]))
  minimo <- tapply(muestra1, mcpogp, min)
  return(minimo)
}
stopCluster(cl)
proc.time()-ñ 

cont <- 0
for(i in 1:S){
  cont <- ifelse(minimo[i,] < minm, cont + 1, cont)
}
pppmin <- cont/S

### Máximo
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)
ñ <- proc.time()

maximo <- foreach(i=1:S,.combine = "rbind")%dopar%{
  set.seed(100+i)
  muestra1 <- rnorm(n = n, mean = rep(zetas_10[i,],njm), sd = sqrt(kappa_12[i]))
  maximo <- tapply(muestra1, mcpogp, max)
  return(maximo)
}
stopCluster(cl)
proc.time()-ñ 

cont <- 0
for(i in 1:S){
  cont <- ifelse(maximo[i,] < maxm, cont + 1, cont)
}
pppmax <- cont/S

### Rango intercuartilico
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)
ñ <- proc.time()

ric <- foreach(i=1:S,.combine = "rbind")%dopar%{
  set.seed(100+i)
  muestra1 <- rnorm(n = n, mean = rep(zetas_10[i,],njm), sd = sqrt(kappa_12[i]))
  ric <- tapply(muestra1, mcpogp, FUN = function(x) {quantile(x, 0.75)-quantile(x,0.25)} )
  return(ric)
}
stopCluster(cl)
proc.time()-ñ 

cont <- 0
for(i in 1:S){
  cont <- ifelse(ric[i,] < ricm, cont + 1, cont)
}
pppric <- cont/S

### Media
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)
ñ <- proc.time()

media <- foreach(i=1:S,.combine = "rbind")%dopar%{
  set.seed(100+i)
  muestra1 <- rnorm(n = n, mean = rep(zetas_10[i,],njm), sd = sqrt(kappa_12[i]))
  media <- tapply(muestra1, mcpogp, mean)
  return(media)
}
stopCluster(cl)
proc.time()-ñ 

cont <- 0
for(i in 1:S){
  cont <- ifelse(media[i,] < ybm, cont + 1, cont)
}
pppmean <- cont/S

### Mediana
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)

ñ <- proc.time()
mediana <- foreach(i=1:S,.combine = "rbind")%dopar%{
  set.seed(100+i)
  muestra1 <- rnorm(n = n, mean = rep(zetas_10[i,],njm), sd = sqrt(kappa_12[i]))
  mediana <- tapply(muestra1, mcpogp, function(x) quantile(x, 0.5))
  return(mediana)
}
stopCluster(cl)
proc.time()-ñ 

cont <- 0
for(i in 1:S){
  cont <- ifelse(mediana[i,] < medianam, cont + 1, cont)
}
pppmediana <- cont/S

### Desviación estándar
cl <- makeCluster(NUCLEOS)
registerDoParallel(cl)

ñ <- proc.time()
sdev <- foreach(i=1:S,.combine = "rbind")%dopar%{
  set.seed(100+i)
  muestra1 <- rnorm(n = n, mean = rep(zetas_10[i,],njm), sd = sqrt(kappa_12[i]))
  sdev <- tapply(muestra1, mcpogp, sd)
  return(sdev)
}
stopCluster(cl)
proc.time()-ñ 

cont <- 0
for(i in 1:S){
  cont <- ifelse(sdev[i,] < sqrt(s2m), cont + 1, cont)
}
pppsd <- cont/S


boxplot(pppmin, pppmax, pppric, pppmean, pppmediana, pppsd, outline = FALSE,
        names = c("Mínimo", "Máximo", "RIC", "Media", "Mediana", "SD"),
        col = "white", border = c("#FF5733", "#33FF57", "#5733FF", "#FF3361", "#FFD733", "#33B4FF"), main = "Bondad de ajuste")
stripchart(list(pppmin, pppmax, pppric, pppmean, pppmediana, pppsd), 
           vertical = TRUE, method = "jitter", pch = 1, col = c("#FF5733", "#33FF57", "#5733FF", "#FF3361", "#FFD733", "#33B4FF"), cex = 0.7, add = TRUE)
