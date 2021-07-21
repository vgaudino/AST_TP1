#Limpia la memoria
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
df1 <- read.csv('../AST/data_rev_dir.csv', sep=',')
df2 <- read.csv('../AST/data_rev_pos.csv', sep=',')

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

