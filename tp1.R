clea#Limpia la memoria
rm( list=ls() )
gc()
options(scipen = 999) ### turn off scientific notation

#---------------------------------------------- 0 ----------------------------------------------
#Tomado de scipt de clase 4 - Sebastian Calcagno - Fernando Martinez
# FUNCION INCORRELACION

# Cargo la siguiente funci?n de incorrelaci?n que realiza un test de Ljung-Box o Box-Pierce para distintos lags

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

#---------------------------------------------- 0 ----------------------------------------------
# FUNCION TEST DE NORMALIDAD
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

#---------------------------------------------- 0 ----------------------------------------------

#Librerias
library(ggplot2)
library(fpp2)
library(tseries)
library(gridExtra) 
library(tseries)
library(forecast)
library(astsa)
library(PerformanceAnalytics)
library(quantmod)
library(Quandl)
library(rlang)
library(dygraphs)
library(urca)
library(dplyr)

## Carga de datos. 
df1 <- read.csv('data_rev_dir.csv', sep=',')
df2 <- read.csv('AST/data_rev_pos.csv', sep=',')

## Exploración incial de los datos
str(df1)
str(df2)

#Renombrando las variables
names(df1)[1] <- "Fecha"
names(df1)[2] <- "Revenue"
names(df2)[1] <- "Fecha"
names(df2)[2] <- "Revenue"

head(df1)
head(df2)

#Modificando la variable de fecha
dates1 <- df1$Fecha
dates1 <- as.Date(dates1, format = "%Y-%m-%d")
df1$Fecha <- dates1
class(df1$Fecha)
head(df1)

dates2 <- df2$Fecha
dates2 <- as.Date(dates2, format = "%Y-%m-%d")
df2$Fecha <- dates2
class(df2$Fecha)
head(df2)

#Se generan las series de tiempo agrupando por mes
#Formula tomada de: https://rpubs.com/DiegoUsai/565288

df1_y <- df1 %>% 
    group_by(month = lubridate::floor_date(Fecha, "month")) %>%
    summarize(summary_variable = sum(Revenue))
head(df1)

df1_y %>% summary()

ts1 <- ts(df1_y$summary_variable, start = c(2018, 10), frequency = 12)
autoplot(ts1, ts.colour = 'dark blue')+ ggtitle("Revenue Directo") + ylab("")

start(ts1)
end(ts1)

df2_y <- df2 %>% 
    group_by(month = lubridate::floor_date(Fecha, "month")) %>%
    summarize(summary_variable = sum(Revenue))
head(df2_y)

df2_y %>% summary()

ts2 <- ts(df2_y$summary_variable, start = c(2018, 10), frequency = 12)
autoplot(ts2, ts.colour = 'dark blue')+ ggtitle("Revenue POS") + ylab("")

start(ts2)
end(ts2)

#Descomposicion de las series

ts1 %>%  decompose(type = "additive") %>% autoplot()
ts2 %>%  decompose(type = "additive") %>% autoplot()

# Grafico la ACF y PACF

par(mfrow=c(1,2))
g1<-acf(ts1, main="Funcion de Autocorrelacion para Revenue Directo")
g2<-pacf(ts1, main="Funcion de Autocorrelacion Parcial para Revenue Directo")

par(mfrow=c(1,2))
g3<-acf(ts2, main="Funcion de Autocorrelacion para Revenue POS")
g4<-pacf(ts2, main="Funcion de Autocorrelacion Parcial para Revenue POS")

# VERIFICO LA ESTACIONARIEDAD

# 1) LJUNG - BOX

# Planteo el test de Ljung-Box. Si rechazo H0 significa que 
# hay coeficientes de autocorrelaci?n distintos a cero
Incorrelation(ts1,"Ljung-Box")
inco_wn1 = Incorrelation(ts1,"Ljung-Box")
autoplot(ts(inco_wn1$P_Value)) + ggtitle("Test de Ljung-Box para Revenue Directo", subtitle = "P-Value") + ylab("") # Grafico los p-value para disntitos lags

Incorrelation(ts2,"Ljung-Box")
inco_wn2 = Incorrelation(ts2,"Ljung-Box")
autoplot(ts(inco_wn2$P_Value)) + ggtitle("Test de Ljung-Box para Revenue POS", subtitle = "P-Value") + ylab("") # Grafico los p-value para disntitos lags

# 2) DICKEY-FULLER

#H0: RAIZ UNITARIA (NO ESTACIONARIO) 
#H1: NO HAY RAIZ UNITARIA (ESTACIONARIO)

#Directo:
adf.test(ts1)
# p-value = 0.5408 -  No Rechazo H0, no estacionario  REVISAR ESTO REVISAR

#POS
adf.test(ts2)
# p-value = 0.07845 - No Rechazo H0, no estacionario EVISAR ESTO REVISAR


# Si los estad?st?cos son menores a los valores cr?ticos, se rechaza H0. La serie es estacionaria
#Directo
summary(ur.df(ts1, type = "none", selectlags = c("AIC")))
summary(ur.df(ts1, type = "drift", selectlags = c("AIC")))
summary(ur.df(ts1, type = "trend", selectlags = c("AIC")))

#POS
summary(ur.df(ts2, type = "none", selectlags = c("AIC")))
summary(ur.df(ts2, type = "drift", selectlags = c("AIC")))
summary(ur.df(ts2, type = "trend", selectlags = c("AIC")))

# 3) KPSS

#H0: ESTACIONARIO
#H1: NO ESTACIONARIO

#Directo
kpss.test(ts1)
# p-value = 0.1

#POS
kpss.test(ts2)
# p-value = 0.1



##############
## FORECAST ##
##############

# ALTERNATIVA 1: Dividiendo en train y test como un problema supervisado normal
train <- window(ts1,start = 2018, end = c(2020,12))
test <- window(ts1, start = 2021)
# Ploteamos el conjunto de train y el conjunto de test
par(mfrow=c(1,1))
ts.plot(train,test, gpars = list(col = c("black", "red")))

# Entrenamos los siguientes modelos. M?s informaci?n sobre los modelos: https://otexts.com/fpp2/simple-methods.html
modelo1 <- meanf(train,h = 7) # Average Method
modelo2 <- naive(train,h = 7) # Naive
modelo3 <- rwf(train,h = 7) # Drift
modelo4 <- snaive(train,h = 7) # Sasonal Naive

autoplot(modelo1)
autoplot(modelo2)
autoplot(modelo3)
autoplot(modelo4)

# Defino el modelo a analizar (del 1 al 4)
modelo <- modelo4 # del 1 al 4

# Predicciones y ploteo
summary(modelo)
autoplot(modelo)

# Ploteamos el conjunto de test con las prediciones
ts.plot(test,modelo$mean , gpars = list(col = c("black", "red")))

# Defino los residuos. Teor?a sobre el an?lisis de los residuos: https://otexts.com/fpp2/residuals.html
residuos <- resid(modelo)

# Verifico la normalidad de los residuos
# H0: NORMALIDAD
Normality_Test(na.omit(residuos),type = "JB") 
Normality_Test(na.omit(residuos),type = "AD")
Normality_Test(na.omit(residuos),type = "SW") 

# Planteo el test de Ljung-Box. Si rechazo H0 significa que hay coeficientes de autocorrelaci?n distintos a cero: hay autocorrelaci?n
Incorrelation(residuos,"Ljung-Box")
inco_wn = Incorrelation(residuos,"Ljung-Box")

# Grafico los p-value para disntitos lags
autoplot(ts(inco_wn$P_Value)) +
  ggtitle("Test de Ljung-Box", subtitle = "P-Value") +
  ylab("")

# Chequeo los residuos
checkresiduals(modelo)

# Estos modelos simples no toman en consideraci?n la estructura autorregresiva de la serie, 
# por lo que el ajuste puede ser bastante bueno y los residuos pueden estar correlacionados. 
# El ajuste se observa mediante un esquema de validaci?n con tus predicciones. En estos modelos, 
# se prioriza el ajuste por encima del an?lisis de los residuos. Esto es diametralmente distinto en 
# los modelos ARIMA donde si en los residuos me queda cierta componente autorregresiva, debo volver 
# a especificar el modelo porque no es correcto.

####### Planteo un modelo ARIMA #######
modelo5 <- auto.arima(train, seasonal = FALSE) # ARIMA
summary(modelo5) # ARIMA(1,0,3) 


# Predicciones y ploteo
pred1 <- forecast(modelo5, 7)
autoplot(pred1)
pred1$mean

# Ploteamos el conjunto de test con las prediciones
ts.plot(test,pred1$mean , gpars = list(col = c("black", "red")))

# Defino los residuos
residuos <- resid(modelo5)

# Verifico la normalidad de los residuos
# H0: NORMALIDAD
Normality_Test(na.omit(residuos),type = "JB") 
Normality_Test(na.omit(residuos),type = "AD")
Normality_Test(na.omit(residuos),type = "SW") 
# Con lo tres test no rechazo el supuesto de normalidad

# Planteo el test de Ljung-Box. Si rechazo H0 significa que hay coeficientes de autocorrelaci?n distintos a cero: hay autocorrelaci?n
Incorrelation(residuos,"Ljung-Box")
inco_wn = Incorrelation(residuos,"Ljung-Box")

# Grafico los p-value para disntitos lags
autoplot(ts(inco_wn$P_Value)) +
  ggtitle("Test de Ljung-Box", subtitle = "P-Value") +
  ylab("")

# Chequeo los residuos
checkresiduals(modelo5)

# Se cumple la normalidad, pero hay correlaci?n en los residuos. Esto indica que queda informaci?n 
# en los residuos que puede usarse para calcular los pron?sticos, es decir que se puede mejorar el modelo.
# Por ejemplo podr?a incorporarse el lag 4 al modelo.

#### Esto No ####
####### Planteo un modelo ARIMA para la serie diferenciada #######
traind <- window(tsd,start = 1992,end = c(2007,4))
testd <- window(tsd, start = 2008)
traind

modelo6 <- auto.arima(traind, seasonal = FALSE) # ARIMA
summary(modelo6) # ARIMA(3,0,1)

# Predicciones y ploteo
pred2 <- forecast(modelo6, 10)
autoplot(pred2)
pred2$mean

# Realizamos la transformaci?n inversa
predaux <- c(pred2$x, pred2$mean) # Concateno la base de train con las predicciones
predaux1 <- cumsum(c(ausbeer1[1], predaux)) # Llevamos la base a los valores originales
tspred <- ts(predaux1 , start= 1992, freq = 4)   # Definimos la serie de tiempo
tspred
tspred_test <- window(tspred, start = 2008, freq = 4) # Ventana temporal de la serie de test
tspred_test 

# Ploteamos el conjunto de test con las prediciones
ts.plot(test,tspred_test, gpars = list(col = c("black", "red")))

# Defino los residuos
residuos <- resid(modelo6)

# Verifico la normalidad de los residuos
# H0: NORMALIDAD
Normality_Test(na.omit(residuos),type = "JB") 
Normality_Test(na.omit(residuos),type = "AD")
Normality_Test(na.omit(residuos),type = "SW") 
# Con lo tres test no rechazo el supuesto de normalidad

# Planteo el test de Ljung-Box. Si rechazo H0 significa que hay coeficientes de autocorrelaci?n distintos a cero: hay autocorrelaci?n
Incorrelation(residuos,"Ljung-Box")
inco_wn = Incorrelation(residuos,"Ljung-Box")

# Grafico los p-value para disntitos lags
autoplot(ts(inco_wn$P_Value)) +
  ggtitle("Test de Ljung-Box", subtitle = "P-Value") +
  ylab("")

# Chequeo los residuos
checkresiduals(modelo6)

# Sigo teniendo correlaci?n al igual que el caso anterior.

####### Planteo un modelo SARIMA ####### https://www.youtube.com/watch?v=WjeGUs6mzXg 
# ver tambi?n: https://otexts.com/fpp2/seasonal-arima.html. Especialmente la parte de comparaci?n de modelos.

modelo7 <- auto.arima(train, seasonal = TRUE) # SARIMA
summary(modelo7) # SARIMA(0,0,0)(0,1,0)[12]

# Predicciones y ploteo
pred3 <- forecast(modelo7, 7)
autoplot(pred3)
pred3$mean

# Ploteamos el conjunto de test con las prediciones
ts.plot(test,pred3$mean , gpars = list(col = c("black", "red")))

# Defino los residuos
residuos <- resid(modelo7)

# Verifico la normalidad de los residuos
# H0: NORMALIDAD
Normality_Test(na.omit(residuos),type = "JB") 
Normality_Test(na.omit(residuos),type = "AD")
Normality_Test(na.omit(residuos),type = "SW") 
# Rechazo el supuesto de normalidad. No son normales.

# Planteo el test de Ljung-Box. Si rechazo H0 significa que hay coeficientes de autocorrelaci?n distintos a cero: hay autocorrelaci?n
Incorrelation(residuos,"Ljung-Box")
inco_wn = Incorrelation(residuos,"Ljung-Box")

# Grafico los p-value para disntitos lags
autoplot(ts(inco_wn$P_Value)) +
  ggtitle("Test de Ljung-Box", subtitle = "P-Value") +
  ylab("")
# todos los p-value son mayores a 0.05. No hay autocorrelaci?n.

# Chequeo los residuos
checkresiduals(modelo7)

# Los residuos est?n incorrelacionados, pero no son normales. Que los resiudos no sean normales afecta la 
# precisi?n del intervalo de predicci?n. Pero no necesariamente hay que descartar el modelo. Puede llegar a
# corregirse mediante una trasnformaci?n de Box-Cox

#######  Verifico la precisi?n de los modelos ####### 
# Video sobre evaluaci?n de modelos: https://www.youtube.com/watch?v=kgBDQ3baESw&list=PLvcbYUQ5t0UHOLnBzl46_Q6QKtFgfMGc3&index=11
# Video sobre las medidas de precisi?n: https://www.youtube.com/watch?v=KzHJXdFJSIQ
accuracy(modelo1, test)
accuracy(modelo2, test)
accuracy(modelo3, test)
accuracy(modelo4, test)
accuracy(pred1, test)
accuracy(tspred_test , test)
accuracy(pred3, test)

# Info adicional sobre las medidas de precisi?n: https://stats.stackexchange.com/questions/194453/interpreting-accuracy-results-for-an-arima-model-fit

# ALTERNATIVA 2: Incorporando Cross-Validation para series de tiempo

# TEOR?A:
# https://campus.datacamp.com/courses/forecasting-in-r/benchmark-methods-and-forecast-accuracy?ex=9
# https://medium.com/@soumyachess1496/cross-validation-in-time-series-566ae4981ce4

# EJEMPLOS PR?CTICOS:
# https://rpubs.com/arquez9512/628219
# https://www.r-bloggers.com/2019/07/time-series-forecast-cross-validation-by-ellis2013nz/

?tsCV

# Error calculado mediante Time Series Cross-Validation para snaive
e1 <- tsCV(ts1, snaive, h=7) 
e1
# Calculo el RSME
sqrt(mean(e1^2, na.rm=TRUE)) # 10106342
# Calculo el MSE. Se calcula como el promedio de los errores
mse1 <- colMeans(e1^2, na.rm = T)
mse1
# Plotear los valores MSE contra el horizonte de predicci?n
data.frame(h = 1:7, MSE = mse1) %>%
  ggplot(aes(x = h, y = MSE)) + geom_point()

# Error calculado mediante Time Series Cross-Validation para ARIMA. Fuente: https://pkg.robjhyndman.com/forecast/reference/tsCV.html
far2 <- function(x, h){forecast(Arima(x, order=c(0,0,0)), h=h)}
e2 <- tsCV(ts1, far2, h=10)
# Calculo el RSME
sqrt(mean(e2^2, na.rm=TRUE)) # 58.1004
# Calculo el MSE
mse2 <- colMeans(e1^2, na.rm = T)
mse2
# Plotear los valores MSE contra el horizonte de predicci?n
data.frame(h = 1:7, MSE = mse2) %>%
  ggplot(aes(x = h, y = MSE)) + geom_point()

