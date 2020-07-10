############# SEGUNDA ENTREGA - ANÁLISIS DE SERIES TEMPORALES #############
# Alumno: Fernando Martínez 

library(tseries)
library(forecast)
library(ggplot2)
library(gridExtra)
library(car)
library(nortest)
library(AdequacyModel)
library(lmtest)
library(quantmod)
library(dygraphs)
library(lessR)
library(lessR)

# Limpio la memoria
rm( list=ls() )
gc()

# Cargo la base con las series temporales
setwd("C:\\Users\\user\\Google Drive\\Austral\\Series Temporales\\Clase 3\\TP")
base=readRDS("series.rds")

# Grafico las 10 series
for (i in 1:10) {
  plot(base[[i]], ylab=i)
}

# Elijo la serie 1
ts<-base[[1]]
ts.plot(ts, main = "Serie 1", ylab = "")


# Planteo el test Dickey-Fuller para confirmar que es estacionaria
adf.test(ts) # p-value = 0.01. Es estacionaria

# Otra forma de comprobar la estacionariedad es detectar si todos los coeficientes de autocorrelación son iguales a cero

# Cargo la siguiente función de incorrelación que realiza un test de Ljung-Box o Box-Pierce para distintos lags
Incorrelation <- function(ts, type = c("Ljung-Box","Box-Pierce"), fitdf = 0){
  p_ljung_box = NULL
  s_ljung_box = NULL
  for(i in 0:(length(ts)/4)){
    p_ljung_box[i] = Box.test(ts,lag = i,type = type,fitdf = fitdf)$p.value
    s_ljung_box[i] = Box.test(ts,lag = i,type = type,fitdf = fitdf)$statistic
  }
  table = data.frame(j = 1:(length(ts)/4),
                     P_Value = p_ljung_box,
                     Statistic = s_ljung_box)
  return(table)
}

# Planteo el test de Ljung-Box. Si rechazo H0 significa que hay coeficientes de autocorrelación distintos a cero
Incorrelation(ts,"Ljung-Box")

inco_wn = Incorrelation(ts,"Ljung-Box")

# Grafico los p-value para disntitos lags
autoplot(ts(inco_wn$P_Value)) +
  ggtitle("Test de Ljung-Box", subtitle = "P-Value") +
  ylab("")
# Ningún p-value supera 0.05 por lo que no rechazo H0 y puedo considerar que es una serie estacionaria 

# Grafico FAS, FAC y FACP 
acf(ts,type = "covariance",plot = T)
ggAcf(ts) + ggtitle("FAC Serie 1")
ggAcf(ts,type = "partial") + ggtitle("FACP Serie 1")

##### MODELOS #####

## Planteo un primer modelo ARMA. AR: orden 1 - MA: orden 0

modelo1 <- arma(ts,order = c(1,0), include.intercept = T)
coef(modelo1)
summary(modelo1) 
# AIC = 562,3
# El coeficiente del AR1 es significativo (p-value menor a 0.5), mientras que el intercepto no lo es

# Ploteo la serie junto con el modelo
ts.plot(ts)
ARMA_fit1 <- ts - residuals(modelo1)
points(ARMA_fit1, type = "l", col = 2, lty = 2)

### Planteo un segundo modelo ARMA. AR: orden 1 - MA: orden 1

modelo2 <- arma(ts,order = c(1,1), include.intercept = T)
coef(modelo2)
summary(modelo2) 
# AIC = 563.12
# El coeficiente del AR1 es significativo (p-value menor a 0.5), mientras que el coeficiente del MA1 y el intercepto no lo son

# Ploteo la serie junto con el modelo
ts.plot(ts)
ARMA_fit2 <- ts - residuals(modelo2)
points(ARMA_fit2, type = "l", col = 2, lty = 2)

# Cargo la funcion AIC_Matrix
AIC_Matrix <- function(ts,p.order = 1, q.order = 1){
  require(forecast)
  aic_matrix = matrix(data = NA, nrow = p.order, ncol = q.order)
  for(i in 0:p.order){
    for(j in 0:q.order){
      aic = arima(ts,order = c(i,0,j))$aic
      aic_matrix[i,j] = aic
    }
  }
  rownames(aic_matrix) = c(1:p.order)
  colnames(aic_matrix) = c(1:q.order)
  return(aic_matrix)
}

AIC_Matrix(ts,p.order = 5,q.order = 5)

# Según esta función, la que mejor criterio de información (menor AIC) tiene es un ARMA de orden (2,1).

## Planteo un tercer modelo ARMA. AR: orden 2 - MA: orden 1

modelo3 <- arma(ts,order = c(2,1), include.intercept = T)
coef(modelo3)
summary(modelo3)
# AIC = 554.55
# salvo el intercepto, todos los coeficientes son significativos

# Analizo los residuos
residuos3 <- resid(modelo3)
Histogram(residuos3,density = T)

# Cargo la siguiente función que realiza test de normalidad
Normality_Test <- function(ts,type = c("JB", "AD", "SW")){
  require(tseries)
  require(nortest)
  if(type == "JB"){
    p_val = jarque.bera.test(ts)$p.value
    stat  = jarque.bera.test(ts)$statistic
  } else if(type == "AD"){
    p_val = ad.test(ts)$p.value
    stat  = ad.test(ts)$statistic
  } else {
    p_val = shapiro.test(ts)$p.value
    stat  = shapiro.test(ts)$statistic
  }
  
  table = data.frame(P_Value = p_val,
                     Statistic = stat)
  return(table)
}

# Verifico la normalidad de los residuos
Normality_Test(na.omit(residuos3),type = "JB") # p-value: 0.7281481
Normality_Test(na.omit(residuos3),type = "AD") # p-value: 0.588623
Normality_Test(na.omit(residuos3),type = "SW") # p-value: 0.6360564

# Con ninguno de los tres test rechazo el supuesto de normalidad. Pruebo la aleatoriedad de los residuos

Incorrelation(residuos3,"Ljung-Box")
inco_wn = Incorrelation(residuos3,"Ljung-Box")
# Grafico los p-value para disntitos lags
autoplot(ts(inco_wn$P_Value)) +
  ggtitle("Test de Ljung-Box", subtitle = "P-Value") +
  ylab("")

# todos los p-value son mayores a 0.05. No rechazo el supuesto de incorrelación (independencia)

# Grafico la FAC y FACP para los residuos buscando que queden dentro del rango
p1=ggAcf(residuos3) + ggtitle("FAC Residuos")
p2=ggAcf(residuos3,type = "partial") + ggtitle("FACP Residuos")
grid.arrange(p1, p2, ncol=1)

# Ploteo la serie junto con el modelo
ts.plot(ts)
ARMA_fit3 <- ts - residuals(modelo3)
points(ARMA_fit3, type = "l", col = 2, lty = 2)




