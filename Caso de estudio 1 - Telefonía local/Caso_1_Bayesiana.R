### Directorio de trabajo
setwd("C:/Users/EQUIPO/OneDrive/Documentos/Documentos Ander/UNAL/Estadística Bayesiana")
getwd()

### Librerías a utilizar
if (!require(invgamma)){install.packages("invgamma")
  library(invgamma)
}
if (!require(doParallel)){install.packages("doParallel")
  library(doParallel)
}
if (!require(foreach)){install.packages("foreach")
  library(foreach)
}
if (!require(dplyr)){install.packages("dplyr")
  library(dplyr)
}

### Datos
data <- read.csv(file = "Verizon.csv")

ILEC <- data %>% 
  filter(Group == "ILEC") %>% as.data.frame()

CLEC <- data %>% 
  filter(Group == "CLEC") %>% as.data.frame()

y1 <- ILEC$Time
y2 <- CLEC$Time

### Hiperparámetros, tamaños de muestra y estadísticos suficientes
n1 <- length(y1);n1
n2 <- length(y2);n2

a1 <- a2 <- 3
b1 <- b2 <- 17

s1 <- sum(y1);s1
s2 <- sum(y2);s2

ybar1 <- 8.41
ybar2 <- 16.51

## Análisis bayesiano ##

### Parámetros de la distribución posterior
a1p <- a1 + n1;a1p
a2p <- a2 + n2;a2p
b1p <- b1 + s1;b1p
b2p <- b2 + s2;b2p

### Visualización distribución previa y posterior
par(mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
theta <- seq(0,25, length = 10000)
plot(NA, NA, xlim = c(0,25),ylim = c(0,2), xlab = expression(lambda),
     ylab = expression(paste("p(",lambda," | ",y,")",sep = "")),
     main = "Distribución Posterior")
lines(theta, dinvgamma(theta, shape = a1p, rate = b1p), col = 2, lwd = 2)
lines(theta, dinvgamma(theta, shape = a2p, rate = b2p), col = 4, lwd = 2)
lines(theta, dinvgamma(theta, shape = a1, rate=b1),lwd = 2)
grid(nx=10)
legend("right", legend = c("ILEC", "CLEC", "Previa"), 
       col = c(2, 4, 1), lty = c(1,1,1), lwd = 2, bty = "n")


# Ajuste del modelo -------------------------------------------------------

### Montecarlo para aproximar eta
set.seed(100)
{lambda1mc <- rinvgamma(n = 10000, shape = a1p, rate = b1p)
lambda2mc <- rinvgamma(n = 10000, shape = a2p, rate = b2p)}

etamc <- lambda1mc-lambda2mc

round(mean(etamc),3)
cvetamc <- round(sd(etamc)/abs(mean(etamc)),3);cvetamc
ICetamc <- round(quantile(etamc,probs=c(0.025,0.975)),3);ICetamc

### Visualización distribución posterior de eta
par(mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
hist(x = etamc, freq = F, col = "gray90", border = "gray90", xlim = c(-20,2), ylim = c(0,0.15),
     xlab = expression(eta), ylab = expression(paste("p(",eta," | ",y,")",sep = "")), 
     main = bquote(~"Distribución posterior de  "~eta~"="~lambda[1]-lambda[2]))
lines(density(etamc), col = 4, lwd = 2)
abline(v = quantile(x = etamc, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
abline(v = mean(etamc), lty = 2, lwd = 2, col = 2)
grid(nx=10)
legend("right", legend = c("Posterior", "IC 95%", "Media"), 
       col = c(4, 3, 2), lty = c(1,2,2), lwd = 2, bty = "n")


# Análisis de sensitividad ------------------------------------------------

### Hiperparámetros de prueba
a11 <- a12 <- 3
b11 <- b12 <- 17

a21 <- a22 <- 2
b21 <- b22 <- 8.5

a31 <- a32 <- 3
b31 <- 16.8
b32 <- 33

a41 <- a42 <- 2
b41 <- 8.4
b42 <- 16.5

# Montecarlo para eta previa
set.seed(100)
{lambda1p1mc <- rinvgamma(n = 10000, shape = a11, rate = b11)
  lambda2p1mc <- rinvgamma(n = 10000, shape = a12, rate = b12)
  
  lambda1p2mc <- rinvgamma(n = 10000, shape = a21, rate = b21)
  lambda2p2mc <- rinvgamma(n = 10000, shape = a22, rate = b22)
  
  lambda1p3mc <- rinvgamma(n = 10000, shape = a31, rate = b31)
  lambda2p3mc <- rinvgamma(n = 10000, shape = a32, rate = b32)
  
  lambda1p4mc <- rinvgamma(n = 10000, shape = a41, rate = b41)
  lambda2p4mc <- rinvgamma(n = 10000, shape = a42, rate = b42)
}

# Muestras de eta para cada previa
etap1mc <- lambda1p1mc-lambda2p1mc
etap2mc <- lambda1p2mc-lambda2p2mc
etap3mc <- lambda1p3mc-lambda2p3mc
etap4mc <- lambda1p4mc-lambda2p4mc

# Medias de las muestras de eta para cada previa
round(mean(etap1mc),3)
round(mean(etap2mc),3)
round(mean(etap3mc),3)
round(mean(etap4mc),3)

# Coeficiente de variación de las muestras de eta para cada previa
cvetap1mc <- round(sd(etap1mc)/abs(mean(etap1mc)),3);cvetap1mc
cvetap2mc <- round(sd(etap2mc)/abs(mean(etap2mc)),3);cvetap2mc
cvetap3mc <- round(sd(etap3mc)/abs(mean(etap3mc)),3);cvetap3mc
cvetap4mc <- round(sd(etap4mc)/abs(mean(etap4mc)),3);cvetap4mc

# Intervalos de confianza de las muestras de eta para cada previa
ICetap1mc <- round(quantile(etap1mc,probs=c(0.025,0.975)),3);ICetap1mc
ICetap2mc <- round(quantile(etap2mc,probs=c(0.025,0.975)),3);ICetap2mc
ICetap3mc <- round(quantile(etap3mc,probs=c(0.025,0.975)),3);ICetap3mc
ICetap4mc <- round(quantile(etap4mc,probs=c(0.025,0.975)),3);ICetap4mc

### Distribución posterior para cada previa

### Hiperparámetros para la distribución posterior

#Previa con a_k=3 y b_k=17
a11p <- a11 + n1;a11p
a12p <- a12 + n2;a12p
b11p <- b11 + s1;b11p
b12p <- b12 + s2;b12p

#Previa con a_k=2 y b_k=8.5
a21p <- a21 + n1;a21p
a22p <- a22 + n2;a22p
b21p <- b21 + s1;b21p
b22p <- b22 + s2;b22p

#Previa con a_k=3, b_1=16.8 y b_2=33
a31p <- a31 + n1;a31p
a32p <- a32 + n2;a32p
b31p <- b31 + s1;b31p
b32p <- b32 + s2;b32p

#Previa con a_k=2, b_1=8.4 y b_2=16.5
a41p <- a41 + n1;a41p
a42p <- a42 + n2;a42p
b41p <- b41 + s1;b41p
b42p <- b42 + s2;b42p

# Montecarlo para aproximar la posterior de eta para cada previa
set.seed(100)
{lambda1p1pmc <- rinvgamma(n = 10000, shape = a11p, rate = b11p)
  lambda2p1pmc <- rinvgamma(n = 10000, shape = a12p, rate = b12p)
  
  lambda1p2pmc <- rinvgamma(n = 10000, shape = a21p, rate = b21p)
  lambda2p2pmc <- rinvgamma(n = 10000, shape = a22p, rate = b22p)
  
  lambda1p3pmc <- rinvgamma(n = 10000, shape = a31p, rate = b31p)
  lambda2p3pmc <- rinvgamma(n = 10000, shape = a32p, rate = b32p)
  
  lambda1p4pmc <- rinvgamma(n = 10000, shape = a41p, rate = b41p)
  lambda2p4pmc <- rinvgamma(n = 10000, shape = a42p, rate = b42p)
}

# Muestras de eta posterior para cada previa
etap1pmc <- lambda1p1pmc-lambda2p1pmc
etap2pmc <- lambda1p2pmc-lambda2p2pmc
etap3pmc <- lambda1p3pmc-lambda2p3pmc
etap4pmc <- lambda1p4pmc-lambda2p4pmc

# Medias de las muestras de eta posterior para cada previa
round(mean(etap1pmc),3)
round(mean(etap2pmc),3)
round(mean(etap3pmc),3)
round(mean(etap4pmc),3)

# Coeficiente de variación de las muestras de eta para cada previa
cvetap1pmc <- round(sd(etap1pmc)/abs(mean(etap1pmc)),3);cvetap1pmc
cvetap2pmc <- round(sd(etap2pmc)/abs(mean(etap2pmc)),3);cvetap2pmc
cvetap3pmc <- round(sd(etap3pmc)/abs(mean(etap3pmc)),3);cvetap3pmc
cvetap4pmc <- round(sd(etap4pmc)/abs(mean(etap4pmc)),3);cvetap4pmc

# Intervalos de confianza de las muestras de eta para cada previa
ICetap1pmc <- round(quantile(etap1pmc,probs=c(0.025,0.975)),3);ICetap1pmc
ICetap2pmc <- round(quantile(etap2pmc,probs=c(0.025,0.975)),3);ICetap2pmc
ICetap3pmc <- round(quantile(etap3pmc,probs=c(0.025,0.975)),3);ICetap3pmc
ICetap4pmc <- round(quantile(etap4pmc,probs=c(0.025,0.975)),3);ICetap4pmc

### Panel con las visualizaciones de las posteriores con cada previa
par(mfrow=c(2,2),mar = c(5,3,1.4,1.4), mgp = c(1.75,0.75,0))
{hist(x = etap1pmc, freq = F, col = "gray90", border = "gray90", xlim = c(-20,2), ylim = c(0,0.15),
     xlab = expression(eta), ylab = expression(paste("p(",eta," | ",y,")",sep = "")), 
     main = bquote(~"Previa con "~a[k]~" = 3 y"~b[k]~" = 17"))
lines(density(etap1pmc), col = 4, lwd = 2)
abline(v = quantile(x = etap1pmc, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
abline(v = mean(etap2pmc), lty = 2, lwd = 2, col = 2)
grid(nx=10)}

{hist(x = etap2pmc, freq = F, col = "gray90", border = "gray90", xlim = c(-20,2), ylim = c(0,0.15),
      xlab = expression(eta), ylab = expression(paste("p(",eta," | ",y,")",sep = "")), 
      main = bquote(~"Previa con "~a[k]~" = 2 y"~b[k]~" = 8.5"))
  lines(density(etap2pmc), col = 4, lwd = 2)
  abline(v = quantile(x = etap2pmc, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
  abline(v = mean(etap2pmc), lty = 2, lwd = 2, col = 2)
  grid(nx=10)}

{hist(x = etap3pmc, freq = F, col = "gray90", border = "gray90", xlim = c(-20,2), ylim = c(0,0.15),
      xlab = expression(eta), ylab = expression(paste("p(",eta," | ",y,")",sep = "")), 
      main = bquote(~"Previa con "~a[k]~" = 3,"~b[1]~" = 16.8 y "~b[2]~" = 33"))
  lines(density(etap3pmc), col = 4, lwd = 2)
  abline(v = quantile(x = etap3pmc, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
  abline(v = mean(etap3pmc), lty = 2, lwd = 2, col = 2)
  grid(nx=10)}

legend_coord <- locator(1)
legend(x=legend_coord[1],y=legend_coord[2], legend = c("Posterior", "IC 95%", "Media"), 
       col = c(4, 3, 2), lty = c(1,2,2), lwd = 2, bty = "n",xpd = TRUE, horiz = TRUE)

{hist(x = etap4pmc, freq = F, col = "gray90", border = "gray90", xlim = c(-20,2), ylim = c(0,0.15),
      xlab = expression(eta), ylab = expression(paste("p(",eta," | ",y,")",sep = "")), 
      main = bquote(~"Previa con "~a[k]~" = 2,"~b[1]~" = 8.4 y "~b[2]~" = 16.5"))
  lines(density(etap4pmc), col = 4, lwd = 2)
  abline(v = quantile(x = etap4pmc, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
  abline(v = mean(etap4pmc), lty = 2, lwd = 2, col = 2)
  grid(nx=10)}

legend_coord <- locator(1)
legend(x=legend_coord[1],y=legend_coord[2], legend = c("Posterior", "IC 95%", "Media"), 
       col = c(4, 3, 2), lty = c(1,2,2), lwd = 2, bty = "n",xpd = TRUE, horiz = TRUE)


# Bondad de ajuste --------------------------------------------------------


### Parámetros de prueba de las muestras
ybar1
ybar2

sd01 <- sd(y1);sd01
sd02 <- sd(y2);sd02

B = 10000 # Número de iteraciones

stats1 <- matrix(data = NA, nrow = B, ncol = 2)
stats2 <- matrix(data = NA, nrow = B, ncol = 2)

####### USO FOREACH

# Configura el número de núcleos a utilizar, en mi caso tengo 12, usaré 10
cl <- makeCluster(10)

# Activo el uso de los núcleos
registerDoParallel(cl)

# Foreach para la bondad de ajuste
stats1 <- foreach(i=1:B,.combine = "rbind")%dopar%{
  set.seed(100)
  muestra1 <- rexp(n = n1, rate = 1/lambda1mc[i])
  a <- mean(muestra1)
  b <- sd(muestra1)
  c <- cbind(a,b)
  return(c)
}

stats2 <- foreach(i=1:B,.combine = "rbind")%dopar%{
  set.seed(100)
  muestra2 <- rexp(n = n2, rate = 1/lambda2mc[i])
  a <- mean(muestra2)
  b <- sd(muestra2)
  c <- cbind(a,b)
  return(c)
}
# Detenemos el uso de varios núcleos
stopCluster(cl)

# Matrices de los estadísticos calculados en las 10000 muestras para cada grupo
stats1
stats2

par(mfrow=c(2,2),mar = c(5,3,1.4,1.4), mgp = c(1.75,0.75,0))
{hist(x = stats1[,1], freq = F, col = "gray90", border = "gray90", xlim = c(7.5,9), ylim = c(0,2),
     xlab = expression(t), ylab =  expression(paste("p(",t," | ",y,")",sep = "")) , 
     main = bquote(~"Bondad de ajuste para"~bar(y)~" grupo 1"))
lines(density(stats1[,1]), col = 4, lwd = 2)
abline(v = quantile(x = stats1[,1], probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
abline(v = ybar1, lty = 2, lwd = 2, col = 2)
grid(nx=10)}

{hist(x = stats1[,2], freq = F, col = "gray90", border = "gray90", xlim = c(7.5,15), ylim = c(0,2),
      xlab = expression(t), ylab =  expression(paste("p(",t," | ",y,")",sep = "")) , 
      main = bquote(~"Bondad de ajuste para"~s~" grupo 1"))
  lines(density(stats1[,2]), col = 4, lwd = 2)
  abline(v = quantile(x = stats1[,2], probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
  abline(v = sd01, lty = 2, lwd = 2, col = 2)
  grid(nx=10)}

{hist(x = stats2[,1], freq = F, col = "gray90", border = "gray90", xlim = c(0,20), ylim = c(0,0.2),
      xlab = expression(t), ylab =  expression(paste("p(",t," | ",y,")",sep = "")) , 
      main = bquote(~"Bondad de ajuste para"~bar(y)~" grupo 2"))
  lines(density(stats2[,1]), col = 4, lwd = 2)
  abline(v = quantile(x = stats2[,1], probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
  abline(v = ybar2, lty = 2, lwd = 2, col = 2)
  grid(nx=10)}



# Análisis frecuentista ################################# 


# Normalidad asintótica ---------------------------------------------------

### Parámetros y estadísticos 
mean1MLE <- mean(y1)
mean2MLE <- mean(y2)
var1MLE <- mean1MLE^2/n1
var2MLE <- mean2MLE^2/n2

meaneta <- mean1MLE-mean2MLE
sdeta <- sqrt(var1MLE+var2MLE)


round(meaneta,3)
cveta <- round(sdeta/abs(meaneta),3);cveta
ICeta <- round(c(meaneta-qnorm(0.975)*sdeta,meaneta+qnorm(0.975)*sdeta),3);ICeta


# Bootstrap frecuentista paramétrico --------------------------------------
B = 10000
lambda1bp <- NULL
lambda2bp <- NULL

### Muestras a partir del modelo exponencial
set.seed(100)
for(i in 1:B){
  bp1 <- rexp(n = n1, rate = 1/ybar1)
  bp2 <- rexp(n = n2, rate = 1/ybar2)
  lambda1bp[i] = mean(bp1)
  lambda2bp[i] = mean(bp2)
}

### Estimadores de eta
etabp<-lambda1bp-lambda2bp

round(mean(etabp),3)
cvetabp <- round(sd(etabp)/abs(mean(etabp)),3);cvetabp
ICetabp <- round(quantile(etabp,probs=c(0.025,0.975)),3);ICetabp

### Visualización de la distribución del Bootstrap paramétrico
par(mfrow=c(1,1),mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
hist(x = etabp, freq = F, col = "gray90", border = "gray90", xlim = c(-20,2), ylim = c(0,0.15),
     xlab = expression(hat(eta)), ylab =  expression(paste("p(",hat(eta)," | ",y,")",sep = "")), 
     main = bquote(~"Densidad de "~hat(eta)~"="~hat(lambda[1])-hat(lambda[2])~" utilizando Bootstrap paramétrico"))
lines(density(etabp), col = 4, lwd = 2)
abline(v = quantile(x = etabp, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
abline(v = mean(etabp), lty = 2, lwd = 2, col = 2)
grid(nx=10)
legend("right", legend = c("Densidad", "IC 95%", "Media"), 
       col = c(4, 3, 2), lty = c(1,2,2), lwd = 2, bty = "n")



# Bootstrap frecuentista no paramétrico -----------------------------------
B = 10000
lambda1bnp <- NULL
lambda2bnp <- NULL

### Muestras a partir de la muestra inicial
set.seed(100)
for(i in 1:B){
  bp1 <- sample(x = y1, size = n1, replace = TRUE) 
  bp2 <- sample(x = y2, size = n2, replace = TRUE) 
  lambda1bnp[i] = mean(bp1)
  lambda2bnp[i] = mean(bp2)
}

### Estimadores de eta
etabnp<-lambda1bnp-lambda2bnp

round(mean(etabnp),3)
cvetabnp <- round(sd(etabnp)/abs(mean(etabnp)),3);cvetabnp
ICetabnp <- round(quantile(etabnp,probs=c(0.025,0.975)),3);ICetabnp

### Visualización Bootstrap no paramétrico
par(mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
hist(x = etabnp, freq = F, col = "gray90", border = "gray90", xlim = c(-20,2), ylim = c(0,0.15),
     xlab = expression(hat(eta)), ylab =  expression(paste("p(",hat(eta)," | ",y,")",sep = "")) , 
     main = bquote(~"Densidad de "~hat(eta)~"="~hat(lambda[1])-hat(lambda[2])~" utilizando Bootstrap no paramétrico"))
lines(density(etabnp), col = 4, lwd = 2)
abline(v = quantile(x = etabnp, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
abline(v = mean(etabnp), lty = 2, lwd = 2, col = 2)
grid(nx=10)
legend("right", legend = c("Densidad", "IC 95%", "Media"), 
       col = c(4, 3, 2), lty = c(1,2,2), lwd = 2, bty = "n")
