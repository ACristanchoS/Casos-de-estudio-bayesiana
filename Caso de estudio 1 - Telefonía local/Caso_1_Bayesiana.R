## Parcial 1 Bayesiana ##
setwd("C:/Users/valer/Downloads")
getwd()

library(invgamma)
library(xtable)
library(ggplot2)
library(ggExtra)
library(foreach)
library(doParallel)

# Lectura de la base de datos 
Y<-read.csv("Verizon.csv")
y1<-Y[Y$Group == "ILEC",1]
y2<-Y[Y$Group == "CLEC",1]

# Tamaños de muestra
n1<-length(y1) #1664 
n2<-length(y2) #23

# Estadisticos suficientes
s1<-sum(y1) #13996.92
s2<-sum(y2) #379.71

# Análisis Bayesiano ################################# 

# ------ Ajuste de los modelos Gamma- Inversa-Exponencial -----#
## Previa Gamma inversa (a,b)
a<-3
b<-17
## Parametros de la distribucion posterior (Actualizacion de parametros)
ap1<- a + n1
bp1<- b + s1
ap2<- a + n2
bp2<- b + s2

theta <- seq(0, 25, length = 5000)  
df <- data.frame(theta = theta,
                 ILEC = dinvgamma(theta, shape = a, rate = b),
                 CLEC = dinvgamma(theta, shape = ap1, rate = bp1),
                 PREVIA = dinvgamma(theta, shape = ap2, rate = bp2))

# Crear el gráfico en ggplot2
ggplot(df, aes(x = theta)) +
  geom_line(aes(y = ILEC, color = "ILEC"), size = 1) +
  geom_line(aes(y = CLEC, color = "CLEC"), size = 1) +
  geom_line(aes(y = PREVIA, color = "PREVIA"), size = 1) +
  labs(x = expression(lambda),
       y = expression(paste("p(", lambda, " | ", y, ")", sep = "")),
       title = "Distribución Posterior") +
  theme_minimal() +
  scale_color_manual(values = c( 1,"purple", "red")) +
  theme(legend.position = "topright", line_types=c("solid", "dashed","solid")) +
  guides(color = guide_legend(title = NULL)) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

## Distribucion posterior Monte Carlo

B<-10000
set.seed(123)
lambda1<-rinvgamma(n = B, shape = ap1, rate = bp1)
lambda2<-rinvgamma(n = B, shape = ap2, rate = bp2)

etha<-lambda1-lambda2 #Distribucion posterior eta por Monte Carlo
tab <- as.matrix(c(mean(etha), quantile(x = etha, probs = c(.025,.975)), sd(etha)/abs(mean(etha))))
rownames(tab) <- c("Estimación","Q2.5%","Q97.5%","Coef. Variacion")
xtable(tab, digits=3)

## Grafico Distribucion posterior de eta

data <- data.frame(etha)

ggplot(data, aes(x = etha, y = ..density..)) +
  geom_histogram(aes(y = ..density..), colour = "gray94", fill = "gray95") +
  geom_density(lwd = 1.3, color = "turquoise3") +
  labs(x = expression(eta), y = expression(paste("P( ", eta ,"| y)")),
       title = expression(paste("Distribución posterior de ", eta))) +
  geom_vline(
    aes(xintercept = quantile(x = etha, probs = 0.025)),
    linetype = "dashed",
    size = 1,
    color ="violetred"
  ) +
  geom_vline(
    aes(xintercept = quantile(x = etha, probs = 0.975)),
    linetype = "dashed",
    size = 1,
    color = "violetred"
  )+
  geom_vline(
    aes(xintercept = mean(etha)),
    linetype = "solid",
    size = 1,
    color = "orangered")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#------------------------- Punto 2 ------------------------ #


## Analisis de sensitividad

a<-c(3,2,3,2)
b1<-c(17,8.5,16.8,8.4)
b2<-c(17,8.5,33,16.5)

## Distribuciones de eta posteriores
set.seed(123)
etap1pmc<-rinvgamma(10000, shape = 3+n1, rate = 17+s1)-rinvgamma(10000, shape = 3+n2, rate = 17+s2)
etap2pmc<-rinvgamma(10000, shape = 2+n1, rate = 8.5+s1)-rinvgamma(10000, shape = 2+n2, rate = 8.5+s2)
etap3pmc<-rinvgamma(10000, shape = 3+n1, rate = 16.8+s1)-rinvgamma(10000, shape = 3+n2, rate = 33+s2)
etap4pmc<-rinvgamma(10000, shape = 2+n1, rate = 8.4+s1)-rinvgamma(10000, shape = 2+n2, rate = 16.5+s2)

## Estadísticos de interés
mean_apriori1<-b1/(a-1)
mean_apriori2<-b2/(a-1)
cv_apriori<-(1/sqrt(a-2))*100
mean_eta <- c(mean(etap1pmc),mean(etap2pmc),mean(etap3pmc),mean(etap4pmc))
cv_eta <- 100*c(sqrt(var(etap1pmc))/abs(mean(etap1pmc)),sqrt(var(etap2pmc))/abs(mean(etap2pmc)),sqrt(var(etap3pmc))/abs(mean(etap3pmc)),sqrt(var(etap4pmc))/abs(mean(etap4pmc)))
IC_eta <- rbind(quantile(etap1pmc, probs = c(0.025,0.975)),quantile(etap2pmc, probs = c(0.025,0.975)),quantile(etap3pmc, probs = c(0.025,0.975)),quantile(etap4pmc, probs = c(0.025,0.975)))
out <- cbind(mean_apriori1, cv_apriori, mean_apriori2, cv_apriori, mean_eta,cv_eta,IC_eta)
colnames(out)<-c("Media previa 1","CV previa 1","Media previa 2","CV previa 2","Media Posterior", "CV Posterior", "Q2.5%", "Q97.5%")
xtable(out, digits = 3)

## Visualización del análisis de sensitividad
par(mfrow=c(2,2),mar = c(5,3,1.4,1.4), mgp = c(1.75,0.75,0))
{hist(x = etap1pmc, freq = F, col = "gray90", border = "gray90", xlim = c(-20,2), ylim = c(0,0.15),
      xlab = expression(eta), ylab = expression(paste("p(",eta," | ",y,")",sep = "")), 
      main = bquote(~"Previa con "~a[k]~" = 3 y"~b[k]~" = 17"))
  lines(density(etap1pmc), col = 4, lwd = 2)
  abline(v = quantile(x = etap1pmc, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
  abline(v = mean(etap2pmc), lty = 1, lwd = 2, col = 2)
  grid(nx=10)}

{hist(x = etap2pmc, freq = F, col = "gray90", border = "gray90", xlim = c(-20,2), ylim = c(0,0.15),
      xlab = expression(eta), ylab = expression(paste("p(",eta," | ",y,")",sep = "")), 
      main = bquote(~"Previa con "~a[k]~" = 2 y"~b[k]~" = 8.5"))
  lines(density(etap2pmc), col = 4, lwd = 2)
  abline(v = quantile(x = etap2pmc, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
  abline(v = mean(etap2pmc), lty = 1, lwd = 2, col = 2)
  grid(nx=10)}

{hist(x = etap3pmc, freq = F, col = "gray90", border = "gray90", xlim = c(-20,2), ylim = c(0,0.15),
      xlab = expression(eta), ylab = expression(paste("p(",eta," | ",y,")",sep = "")), 
      main = bquote(~"Previa con "~a[k]~" = 3,"~b[1]~" = 16.8 y "~b[2]~" = 33"))
  lines(density(etap3pmc), col = 4, lwd = 2)
  abline(v = quantile(x = etap3pmc, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
  abline(v = mean(etap3pmc), lty = 1, lwd = 2, col = 2)
  grid(nx=10)}

legend_coord <- locator(1)
legend(x=legend_coord[1],y=legend_coord[2], legend = c("Posterior", "IC 95%", "Media"), 
       col = c(4, 3, 2), lty = c(1,2,1), lwd = 2, bty = "n",xpd = TRUE, horiz = TRUE)

{hist(x = etap4pmc, freq = F, col = "gray90", border = "gray90", xlim = c(-20,2), ylim = c(0,0.15),
      xlab = expression(eta), ylab = expression(paste("p(",eta," | ",y,")",sep = "")), 
      main = bquote(~"Previa con "~a[k]~" = 2,"~b[1]~" = 8.4 y "~b[2]~" = 16.5"))
  lines(density(etap4pmc), col = 4, lwd = 2)
  abline(v = quantile(x = etap4pmc, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
  abline(v = mean(etap4pmc), lty = 1, lwd = 2, col = 2)
  grid(nx=10)}

legend_coord <- locator(1)
legend(x=legend_coord[1],y=legend_coord[2], legend = c("Posterior", "IC 95%", "Media"), 
       col = c(4, 3, 2), lty = c(1,2,1), lwd = 2, bty = "n",xpd = TRUE, horiz = TRUE)

# --------------------- Punto 3 ---------------------- #
## Bondad de ajuste
ybarra1<-s1/n1
ybarra2<-s2/n2
sd11<-sd(y1)
sd22<-sd(y2)

media1<-NULL
media2<-NULL
sd1<-NULL
sd2<-NULL
set.seed(123)
for (i in 1:10000) {
  y1_i<-rexp(n1, rate=1/lambda1[i])
  media1[i]<-mean(y1_i)
  sd1[i]<-sd(y1_i)
  y2_i<-rexp(n2, rate=1/lambda2[i])
  media2[i]<-mean(y2_i)
  sd2[i]<-sd(y2_i)
}

par(mfrow=c(2,2),mar = c(5,3,1.4,1.4), mgp = c(1.75,0.75,0))
{hist(x = media1, freq = F, col = "gray90", border = "gray90", xlim = c(7,12), ylim = c(0, ceiling(max(density(media1)$y))), xlab = "t", ylab =  expression(paste("P( t | y)")), main = "Bondad de ajuste con media")
  grid(nx=10)
  lines(density(media1), col = 4, lwd = 2)
  abline(v = ybarra1, col = 14, lwd = 2, lty = 1)
  abline(v = quantile(x = media1, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)}

{hist(sd1, freq = F, xlim = c(7,15), ylim = c(0,1.5), col = "gray90", border = "gray90",  xlab = "t", ylab =  expression(paste("P( t | y)")), main = "Bondad de ajuste con sd" )
  grid(nx=10)
  lines(density(sd1), col = 4, lwd = 2)
  abline(v = sd11, col = 14, lwd = 2, lty = 1)
  abline(v = quantile(x = sd1, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)}

{hist(x = media2, freq = F, col = "gray90", border = "gray90", xlim = c(0,48), ylim = c(0, 0.15), xlab = "t", ylab =  expression(paste("P( t | y)")), main = "Bondad de ajuste con media")
  grid(nx=10)
  lines(density(media2), col = 4, lwd = 2)
  abline(v = ybarra2, col = 14, lwd = 2, lty = 1)
  abline(v = quantile(x = media2, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)}

legend_coord <- locator(1)
legend(x=legend_coord[1],y=legend_coord[2], legend = c("Posterior", "IC 95%", "Media"), 
       col = c(4, 3, 14), lty = c(1,2,2), lwd = 2, bty = "n",xpd = TRUE, horiz = TRUE)

{hist(sd2, freq = F, xlim = c(1,40), ylim = c(0,0.12), col = "gray90", border = "gray90",  xlab = "t", ylab =  expression(paste("P( t | y)")),main = "Bondad de ajuste con sd" )
  grid(nx=10)
  lines(density(sd2), col = 4, lwd = 2)
  abline(v = sd22, col = 14, lwd = 2, lty = 1)
  abline(v = quantile(x = sd2, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
}

legend_coord <- locator(1)
legend(x=legend_coord[1],y=legend_coord[2], legend = c("Posterior", "IC 95%", "t obs"), 
       col = c(4, 3, 14), lty = c(1,2,1), lwd = 2, bty = "n",xpd = TRUE, horiz = TRUE)

## ppp

media_1<-c(mean(media1>ybarra1), ybarra1, quantile(media1, probs = c(0.025, 0.975)))
media_2<-c(mean(media2>ybarra2), ybarra2, quantile(media2, probs = c(0.025, 0.975)))
sd_1<-c(mean(sd1>sd11), sd11, quantile(sd1, probs = c(0.025, 0.975)))
sd_2<-c(mean(sd2>sd22), sd22, quantile(sd2, probs = c(0.025, 0.975)))

ppp<-rbind(media_1,media_2, sd_1, sd_2)

colnames(ppp)<-c("PPP", "Estimación", "Lim. inf.", "Lim. Sup")
xtable(ppp, digits = 3)

## Dispersograma
disp1<-data.frame(media1,sd1)
disp2<-data.frame(media2,sd2)

## Dispersograma ILEC
disper1 <- ggplot(disp1, aes(x = media1, y = sd1)) +
  geom_point() +
  labs(
    x = "Media",
    y = "Desviación estándar",
    title = "Distribuciones predictivas de los estadisticos de prueba ILEC"
  )+theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(x = ybarra1, y = sd11, color = "magenta2", size = 4)+
  geom_vline(xintercept = ybarra1, color = "salmon", linetype = "dashed") +
  geom_hline(yintercept = sd11, color = "salmon", linetype = "dashed")
  
disper_1 <- ggMarginal(disper1, type = "histogram", fill="salmon")
print(disper_1)

## Dispersograma CLEC
disper2 <- ggplot(disp2, aes(x = media2, y = sd2)) +
  geom_point() +
  labs(
    x = "Media",
    y = "Desviación estándar",
    title = "Distribuciones predictivas de los estadisticos de prueba CLEC"
  )+theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(x = ybarra2, y = sd22, color = "palevioletred1", size = 4)+
  geom_vline(xintercept = ybarra2, color = "violetred", linetype = "dashed") +
  geom_hline(yintercept = sd22, color = "violetred", linetype = "dashed")

disper_2 <- ggMarginal(disper2, type = "histogram", fill="violetred")
print(disper_2)

#------------------------- Análisis frecuentista -------------------------# 


# Normalidad asintótica ---------------------------------------------------

### Parámetros y estadísticos 
ybarra1
ybarra2
var1MLE <- ybarra1^2/n1
var2MLE <- ybarra2^2/n2

mean_etaMLE <- ybarra1-ybarra2
sd_etaMLE <- sqrt(var1MLE+var2MLE)


round(mean_etaMLE,3)
cveta <- round(sd_etaMLE/abs(mean_etaMLE),3);cveta
ICeta <- round(c(mean_etaMLE-qnorm(0.975)*sd_etaMLE,mean_etaMLE+qnorm(0.975)*sd_etaMLE),3);ICeta
nasin<-cbind(round(mean_etaMLE,3),cveta,t(ICeta))

### Aproximación de la distribución de etaMLE
theta <- seq(-25,10, length = 100000)  
par(mfrow=c(1,1),mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
hist(x = qnorm(theta, mean = mean_etaMLE, sd = sd_etaMLE), freq = F, col = "gray90", border = "gray90", xlim = c(-20,5), ylim = c(0,0.15),
     xlab = expression(hat(eta)), ylab =  expression(paste("p(",hat(eta)," | ",y,")",sep = "")), 
     main = bquote(~"Distribución aproximada de "~hat(eta)~"="~hat(lambda[1])-hat(lambda[2])))
grid(nx = 30, ny = 20, col = "gray95")
lines(theta, dnorm(theta, mean = mean_etaMLE, sd = sd_etaMLE), col=4, lwd=2)
abline(v = c(qnorm(0.025, mean = mean_etaMLE, sd = sd_etaMLE), qnorm(0.975, mean = mean_etaMLE, sd = sd_etaMLE)), lty = 2, lwd = 2, col = 3)
abline(v = mean_etaMLE, lwd = 2, col = 2)
legend("right", legend = c("Densidad", "IC 95%", "Media"), 
       col = c(4, 3, 2), lty = c(1,2,1), lwd = 2, bty = "n")

# Bootstrap frecuentista paramétrico --------------------------------------
B = 10000
lambda1bp <- NULL
lambda2bp <- NULL

### Muestras a partir del modelo exponencial
set.seed(100)
for(i in 1:B){
  bp1 <- rexp(n = n1, rate = 1/ybarra1)
  bp2 <- rexp(n = n2, rate = 1/ybarra2)
  lambda1bp[i] = mean(bp1)
  lambda2bp[i] = mean(bp2)
}

### Estimadores de eta
etabp<-lambda1bp-lambda2bp

round(mean(etabp),3)
cvetabp <- round(sd(etabp)/abs(mean(etabp)),3);cvetabp
ICetabp <- round(quantile(etabp,probs=c(0.025,0.975)),3);ICetabp
bp <-cbind(round(mean(etabp),3),cvetabp,t(ICetabp))

### Visualización de la distribución del Bootstrap paramétrico
par(mfrow=c(1,1),mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
hist(x = etabp, freq = F, col = "gray90", border = "gray90", xlim = c(-20,2), ylim = c(0,0.15),
     xlab = expression(hat(eta)), ylab =  expression(paste("p(",hat(eta)," | ",y,")",sep = "")), 
     main = bquote(~"Distribución de "~hat(eta)~"="~hat(lambda[1])-hat(lambda[2])~" utilizando Bootstrap paramétrico"))
lines(density(etabp), col = 4, lwd = 2)
abline(v = quantile(x = etabp, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
abline(v = mean(etabp), lwd = 2, col = 2)
grid(nx=10)
legend("right", legend = c("Densidad", "IC 95%", "Media"), 
       col = c(4, 3, 2), lty = c(1,2,1), lwd = 2, bty = "n")



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

bnp <-cbind(round(mean(etabnp),3),cvetabp,t(ICetabnp))

### Visualización Bootstrap no paramétrico
par(mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
hist(x = etabnp, freq = F, col = "gray90", border = "gray90", xlim = c(-20,2), ylim = c(0,0.15),
     xlab = expression(hat(eta)), ylab =  expression(paste("p(",hat(eta)," | ",y,")",sep = "")) , 
     main = bquote(~"Distribución de "~hat(eta)~"="~hat(lambda[1])-hat(lambda[2])~" utilizando Bootstrap no paramétrico"))
lines(density(etabnp), col = 4, lwd = 2)
abline(v = quantile(x = etabnp, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = 3)
abline(v = mean(etabnp), lwd = 2, col = 2)
grid(nx=10)
legend("right", legend = c("Densidad", "IC 95%", "Media"), 
       col = c(4, 3, 2), lty = c(1,2,1), lwd = 2, bty = "n")

xtable(rbind(nasin,bp,bnp),digits = 3)


# ----------------------------Simulación ----------------------------- # 
ns <- c(10,20,50,100)
K <- 100000
Bay <- 10000
B <- 1000

ybarra1 <- mean(y1)
ybarra2 <- mean(y2)
etaobs <- ybarra1-ybarra2

a<-3
b<-17

as <- a + ns

propbay <- NULL
propmle <- NULL
propbp <- NULL
propbnp <- NULL

cont <- c(0,0,0,0)

# Configura el número de núcleos a utilizar, en mi caso tengo 12, usaré 10

NUCLEOS <- 10

cl <- makeCluster(NUCLEOS)

# Activo el uso de los núcleos
registerDoParallel(cl)

## Modelo bayesiano
ñ <- proc.time() 
set.seed(123)
for (j in 1:length(ns)) {
  contbay <- foreach (i = 1:K,.combine = "+")%dopar%{
    set.seed(123+i)
    muestra1 <- rexp(n = ns[j], rate = 1/ybarra1)
    muestra2 <- rexp(n = ns[j], rate = 1/ybarra2)
    
    ## Modelo bayesiano
    bs1<- b + sum(muestra1)
    bs2<- b + sum(muestra2)
    
    etasb<-1/rgamma(n = Bay, shape = as[j], rate = bs1)-1/rgamma(n = Bay, shape = as[j], rate = bs2)
    icbay<-quantile(x = etasb, probs = c(.025,.975))
    
    if(etaobs>=icbay[1] && etaobs<=icbay[2]){
      1
    }
    else{0}
  }
  propbay[j] <- contbay*100/K
}

registerDoParallel(cl)

## Normalidad Asintótica
set.seed(123)
for (j in 1:length(ns)) {
  contmle <- foreach (i = 1:K,.combine = "+")%dopar%{
    set.seed(123+i)
    muestra1 <- rexp(n = ns[j], rate = 1/ybarra1)
    muestra2 <- rexp(n = ns[j], rate = 1/ybarra2)
    
    mean1 <- mean(muestra1)
    mean2 <- mean(muestra2)
    etamle <- mean1-mean2
    sdetamle <- sqrt(mean1^2/ns[j]+mean2^2/ns[j])
    
    icmle <- c(etamle-1.959964*sdetamle,etamle+1.959964*sdetamle)
    
    if(etaobs>=icmle[1] && etaobs<=icmle[2]){
      1
    }
    else{0}
  }
  propmle[j] <- contmle*100/K
} 

registerDoParallel(cl)

## Bootstrap paramétrico
set.seed(123)
for (j in 1:length(ns)) {
  
  etasimbp <- NULL
  
  contbp <- foreach (i = 1:K,.combine = "+")%dopar%{
    set.seed(123+i)
    muestra1 <- rexp(n = ns[j], rate = 1/ybarra1)
    muestra2 <- rexp(n = ns[j], rate = 1/ybarra2)
    set.seed(123+i)
    for(l in 1:B){
      etasimbp[l] = mean(rexp(n = ns[j], rate = 1/mean(muestra1)))-mean(rexp(n = ns[j], rate = 1/mean(muestra2)))
    }
    icbp<-quantile(x = etasimbp, probs = c(.025,.975))
    
    if(etaobs>=icbp[1] && etaobs<=icbp[2]){
      1
    }
    else{0}
  }
  propbp[j] <- contbp*100/K
}

registerDoParallel(cl)

## Bootstrap no paramétrico
set.seed(123)
for (j in 1:length(ns)) {
  
  etasimbnp <- NULL
  
  contbnp <- foreach (i = 1:K,.combine = "+")%dopar%{
    set.seed(123+i)
    muestra1 <- rexp(n = ns[j], rate = 1/ybarra1)
    muestra2 <- rexp(n = ns[j], rate = 1/ybarra2)
    
    for(q in 1:B){
      etasimbnp[q] <- mean(sample(x = muestra1, size = ns[j], replace = TRUE))-mean(sample(x = muestra2, size = ns[j], replace = TRUE))
    }
    icbnp<-quantile(x = etasimbnp, probs = c(.025,.975))
    
    if(etaobs>=icbnp[1] && etaobs<=icbnp[2]){
      1
    }
    else{0}
  }
  propbnp[j] <- contbnp*100/K
}
props <- cbind(propbay,propmle,propbp,propbnp)

proc.time()-ñ   
# Detenemos el uso de varios núcleos
stopCluster(cl)

xtable(props,digits = 3)
