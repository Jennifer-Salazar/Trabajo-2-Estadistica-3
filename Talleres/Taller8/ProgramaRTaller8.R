library(forecast)
library(TSA)
library(lmtest)
#Creando función usuario para obtener test Box-Pierce y Ljung-Box
BP.LB.test=function(serie,maxlag,type){
aux=floor(maxlag/6) 
X.squared=c(rep(NA,aux))
df=c(rep(NA,aux))
p.value=c(rep(NA,aux))
for(i in 1:aux){
test=Box.test(serie,lag=(6*i),type=type)
X.squared[i]=test[[1]]
df[i]=test[[2]]
p.value[i]=test[[3]]
}
lag=6*c(1:aux)
teste=as.data.frame(cbind(X.squared,df,p.value))
rownames(teste)=lag
teste
}

#Creando función usuario crit.inf.resid() para calcular C_n^*(p)
crit.inf.resid=function(residuales,n.par,AIC="TRUE"){
if(AIC=="TRUE"){
#Calcula AIC
CI=log(mean(residuales^2))+2*n.par/length(residuales)
}
if(AIC=="FALSE"){
#Calcula BIC
CI=log(mean(residuales^2))+n.par*log(length(residuales))/length(residuales)
}
CI
}

#datos 1
#Leer AR2SIMUL.txt
datos1=scan(file.choose())
datos1=ts(datos1,freq=1)

layout(rbind(c(1,1,2,2),c(0,3,3,0)))
plot(datos1)
abline(h=mean(datos1))
acf(datos1,ci.type="ma",lag.max=25,main="ACF muestral datos1")
pacf(datos1,lag.max=25,main="PACF muestral datos1")

#Identificación modelos sobre datos 1
eacf(datos1,ar.max=16,ma.max=16) #EACF
auto.arima(datos1,ic="aic")
auto.arima(datos1,ic="bic")

##Modelos sobre datos 1
#Modelo AR(2) con media cero
mod1=Arima(datos1,order=c(2,0,0),include.mean=F,method="ML")
df1=length(datos1)-2 #Número parámetros es p=2
coeftest(mod1,df=df1)

yhat1=mod1$fitted

#ARMA(2,2) con media cero
mod2=Arima(datos1,order=c(2,0,2),include.mean=F,method="ML")
df2=length(datos1)-4 #Número parámetros es p=4
coeftest(mod2,df=df2)

yhat2=mod2$fitted

#Medidas de ajuste
aic1=exp(crit.inf.resid(residuals(mod1),n.par=2))
aic2=exp(crit.inf.resid(residuals(mod2),n.par=4))
bic1=exp(crit.inf.resid(residuals(mod1),n.par=2,AIC="FALSE"))
bic2=exp(crit.inf.resid(residuals(mod2),n.par=4,AIC="FALSE"))

criterios=data.frame(AIC=c(aic1,aic2),BIC=c(bic1,bic2),row.names=c("AR(2)","ARMA(2,2)"))
criterios

#Gráficos de los ajustes
win.graph(width=10,height=5)
layout(matrix(c(1,1,2,2),ncol=4))
plot(datos1)
lines(yhat1,col=2)
legend("topleft",legend=c("Datos 1","Ajuste AR(2)"),col=1:2,lwd=2)

plot(datos1)
lines(yhat2,col=2)
legend("topleft",legend=c("Datos 1","Ajuste ARMA(2,2)"),col=1:2,lwd=2)

#Gráficos residuos
win.graph() 
layout(rbind(c(1,1,2,2),c(3,3,4,4)))
plot(residuals(mod1),ylim=c(-2.5,3.1))
abline(h=c(-2*sqrt(mod1$sigma2),0,2*sqrt(mod1$sigma2)),col=2)

plot(residuals(mod2),ylim=c(-2.5,3.1))
abline(h=c(-2*sqrt(mod2$sigma2),0,2*sqrt(mod2$sigma2)),col=2)

plot(as.numeric(yhat1),residuals(mod1),ylim=c(-2.5,3.1))
abline(h=c(-2*sqrt(mod1$sigma2),0,2*sqrt(mod1$sigma2)),col=2)

plot(as.numeric(yhat2),residuals(mod2),ylim=c(-2.5,3.1))
abline(h=c(-2*sqrt(mod2$sigma2),0,2*sqrt(mod2$sigma2)),col=2)


#ACF y PACF residuales
win.graph() 
layout(rbind(c(1,1,2,2),c(3,3,4,4)))
acf(as.numeric(residuals(mod1)),ci.type="ma",lag.max=25,main="ACF AR(2) sobre datos1",ci.col=2)
acf(as.numeric(residuals(mod2)),ci.type="ma",lag.max=25,main="ACF ARMA(2,2) sobre datos1",ci.col=2)
pacf(as.numeric(residuals(mod1)),lag.max=25,main="PACF modelo AR(2) sobre datos1",ci.col=2)
pacf(as.numeric(residuals(mod2)),lag.max=25,main="PACF modelo ARMA(2,2) sobre datos1",ci.col=2)

#Tests Ljun-Box
BP.LB.test(residuals(mod1),maxlag=25,type="Ljung")
BP.LB.test(residuals(mod2),maxlag=25,type="Ljung")

#tests de normalidad
shapiro.test(residuals(mod1))
shapiro.test(residuals(mod2))

#Gráficos de probabilidad normal
win.graph(width=10,height=5)
layout(matrix(c(1,1,2,2),ncol=4))
qqnorm(residuals(mod1));qqline(residuals(mod1),col=2)
legend("topleft",legend="Modelo AR(2) sobre datos1")
qqnorm(residuals(mod2));qqline(residuals(mod2),col=2)
legend("topleft",legend="Modelo ARMA(2,2) sobre datos1")


#datos 2
#Leer MA2SIMULb.txt
datos2=scan(file.choose())
datos2=ts(datos2,freq=1)

win.graph()
layout(rbind(c(1,1,2,2),c(0,3,3,0)))
plot(datos2)
abline(h=mean(datos2))
acf(datos2,ci.type="ma",lag.max=25,main="ACF muestral datos2")
pacf(datos2,lag.max=25,main="PACF muestral datos2")

#Identificación con algunas funciones R de modelos sobre datos 2
eacf(datos2,ar.max=16,ma.max=16) #EACF
auto.arima(datos2,ic="aic")
auto.arima(datos2,ic="bic")

##Modelos sobre datos 2
#Modelo MA(2) con media cero
mod1b=Arima(datos2,order=c(0,0,2),include.mean=F,method="ML")
df1b=length(datos2)-2 #Número parámetros es p=2
coeftest(mod1b,df=df1b)

yhat1b=mod1b$fitted

#ARMA(3,2) con media cero
mod2b=Arima(datos2,order=c(3,0,2),include.mean=F,method="ML")
df2b=length(datos2)-5 #Número parámetros es p=5
coeftest(mod2b,df=df2b)

yhat2b=mod2b$fitted

#Medidas de ajuste
aic1b=exp(crit.inf.resid(residuals(mod1b),n.par=2))
aic2b=exp(crit.inf.resid(residuals(mod2b),n.par=5))
bic1b=exp(crit.inf.resid(residuals(mod1b),n.par=2,AIC="FALSE"))
bic2b=exp(crit.inf.resid(residuals(mod2b),n.par=5,AIC="FALSE"))

criteriosb=data.frame(AIC=c(aic1b,aic2b),BIC=c(bic1b,bic2b),row.names=c("MA(2)","ARMA(3,2)"))
criteriosb

#Gráficos de los ajustes
win.graph(width=10,height=5)
layout(matrix(c(1,1,2,2),ncol=4))
plot(datos2)
lines(yhat1b,col=2)
legend("topleft",legend=c("datos2","Ajuste MA(2)"),col=1:2,lwd=2)

plot(datos2)
lines(yhat2b,col=2)
legend("topleft",legend=c("datos2","Ajuste ARMA(3,2)"),col=1:2,lwd=2)

#Gráficos residuos
win.graph() 
layout(rbind(c(1,1,2,2),c(3,3,4,4)))
plot(residuals(mod1b),ylim=c(-2.5,3.1))
abline(h=c(-2*sqrt(mod1b$sigma2),0,2*sqrt(mod1b$sigma2)),col=2)

plot(residuals(mod2b),ylim=c(-2.5,3.1))
abline(h=c(-2*sqrt(mod2b$sigma2),0,2*sqrt(mod2b$sigma2)),col=2)

plot(as.numeric(yhat1b),residuals(mod1b),ylim=c(-2.5,3.1))
abline(h=c(-2*sqrt(mod1b$sigma2),0,2*sqrt(mod1b$sigma2)),col=2)

plot(as.numeric(yhat2b),residuals(mod2b),ylim=c(-2.5,3.1))
abline(h=c(-2*sqrt(mod2b$sigma2),0,2*sqrt(mod2b$sigma2)),col=2)

#ACF y PACF residuales
win.graph() 
layout(rbind(c(1,1,2,2),c(3,3,4,4)))
acf(as.numeric(residuals(mod1b)),ci.type="ma",lag.max=25,main="ACF MA(2) sobre datos2",ci.col=2)
acf(as.numeric(residuals(mod2b)),ci.type="ma",lag.max=25,main="ACF ARMA(3,2) sobre datos2",ci.col=2)
pacf(as.numeric(residuals(mod1b)),lag.max=25,main="PACF modelo MA(2) sobre datos2",ci.col=2)
pacf(as.numeric(residuals(mod2b)),lag.max=25,main="PACF modelo ARMA(3,2) sobre datos2",ci.col=2)

#Tests Ljun-Box
BP.LB.test(residuals(mod1b),maxlag=25,type="Ljung")
BP.LB.test(residuals(mod2b),maxlag=25,type="Ljung")

#tests de normalidad
shapiro.test(residuals(mod1b))
shapiro.test(residuals(mod2b))

#Gráficos de probabilidad normal
win.graph(width=10,height=5)
layout(matrix(c(1,1,2,2),ncol=4))
qqnorm(residuals(mod1b));qqline(residuals(mod1b),col=2)
legend("topleft",legend="Modelo MA(2) sobre datos2")
qqnorm(residuals(mod2b));qqline(residuals(mod2b),col=2)
legend("topleft",legend="Modelo ARMA(3,2) sobre datos2")


##datos 3
#Leer ARMA1.1SIMUL.txt
datos3=scan(file.choose())
datos3=ts(datos3,freq=1)

win.graph()
layout(rbind(c(1,1,2,2),c(0,3,3,0)))
plot(datos3)
abline(h=mean(datos3))
acf(datos3,ci.type="ma",lag.max=25,main="ACF muestral datos3")
pacf(datos3,lag.max=25,main="PACF muestral datos3")

#Identificación con algunas funciones R de modelos sobre datos 3
eacf(datos3,ar.max=16,ma.max=16)
auto.arima(datos3,ic="aic")
auto.arima(datos3,ic="bic")

##modelos sobre datos 3
#Modelo ARMA(1,1) con media cero
mod1c=Arima(datos3,order=c(1,0,1),include.mean=F,method="ML")
df1c=length(datos3)-2 #Número parámetros es p=2
coeftest(mod1c,df=df1c)

yhat1c=mod1c$fitted

#ARMA(2,2) con media cero
mod2c=Arima(datos3,order=c(2,0,2),include.mean=F,method="ML")
df2c=length(datos3)-4 #Número parámetros es p=4
coeftest(mod2c,df=df2c)

yhat2c=mod2c$fitted

#Medidas de ajuste
aic1c=exp(crit.inf.resid(residuals(mod1c),n.par=2))
aic2c=exp(crit.inf.resid(residuals(mod2c),n.par=4))
bic1c=exp(crit.inf.resid(residuals(mod1c),n.par=2,AIC="FALSE"))
bic2c=exp(crit.inf.resid(residuals(mod2c),n.par=4,AIC="FALSE"))

criteriosc=data.frame(AIC=c(aic1c,aic2c),BIC=c(bic1c,bic2c),row.names=c("ARMA(1,1)","ARMA(2,2)"))
criteriosc

#Gráficos de los ajustes
win.graph(width=10,height=5)
layout(matrix(c(1,1,2,2),ncol=4))
plot(datos3)
lines(yhat1c,col=2)
legend("topleft",legend=c("datos3","Ajuste ARMA(1,1)"),col=1:2,lwd=2)

plot(datos3)
lines(yhat2c,col=2)
legend("topleft",legend=c("datos3","Ajuste ARMA(2,2)"),col=1:2,lwd=2)

#Gráficos residuos
win.graph() 
layout(rbind(c(1,1,2,2),c(3,3,4,4)))
plot(residuals(mod1c),ylim=c(-2.5,3.1))
abline(h=c(-2*sqrt(mod1c$sigma2),0,2*sqrt(mod1c$sigma2)),col=2)

plot(residuals(mod2c),ylim=c(-2.5,3.1))
abline(h=c(-2*sqrt(mod2c$sigma2),0,2*sqrt(mod2c$sigma2)),col=2)

plot(as.numeric(yhat1c),residuals(mod1c),ylim=c(-2.5,3.1))
abline(h=c(-2*sqrt(mod1c$sigma2),0,2*sqrt(mod1c$sigma2)),col=2)

plot(as.numeric(yhat2c),residuals(mod2c),ylim=c(-2.5,3.1))
abline(h=c(-2*sqrt(mod2c$sigma2),0,2*sqrt(mod2c$sigma2)),col=2)

#ACF y PACF residuales
win.graph() 
layout(rbind(c(1,1,2,2),c(3,3,4,4)))
acf(as.numeric(residuals(mod1c)),ci.type="ma",lag.max=25,main="ACF ARMA(1,1) sobre datos3",ci.col=2)
acf(as.numeric(residuals(mod2c)),ci.type="ma",lag.max=25,main="ACF ARMA(2,2) sobre datos3",ci.col=2)
pacf(as.numeric(residuals(mod1c)),lag.max=25,main="PACF modelo ARMA(1,1) sobre datos3",ci.col=2)
pacf(as.numeric(residuals(mod2c)),lag.max=25,main="PACF modelo ARMA(2,2) sobre datos3",ci.col=2)

#Tests Ljung-Box
BP.LB.test(residuals(mod1c),maxlag=25,type="Ljung")
BP.LB.test(residuals(mod2c),maxlag=25,type="Ljung")

#tests de normalidad
shapiro.test(residuals(mod1c))
shapiro.test(residuals(mod2c))

#Gráficos de probabilidad normal
win.graph(width=10,height=5)
layout(matrix(c(1,1,2,2),ncol=4))
qqnorm(residuals(mod1c));qqline(residuals(mod1c),col=2)
legend("topleft",legend="Modelo ARMA(1,1) sobre datos3")
qqnorm(residuals(mod2c));qqline(residuals(mod2c),col=2)
legend("topleft",legend="Modelo ARMA(2,2)sobre datos3")


