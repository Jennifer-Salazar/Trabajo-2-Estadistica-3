library(forecast)
library(car)
library(TSA)
library(FitAR)
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

#NOTA: Esta programación corre en la versión R para windows 3.3.3 o superior
#Creando una función usuario para obtener test Durbin-Watson 
#para autocorrelación de orden 1, positiva y negativa
pruebaDW1=function(modelo){
dwneg=durbinWatsonTest(modelo,max.lag=1,method="normal",alternative="negative")
dwpos=durbinWatsonTest(modelo,max.lag=1,method="normal",alternative="positive")
res=data.frame(1,dwneg$r,dwneg$dw,dwpos$p,dwneg$p)
names(res)=c("lag","rho estimado","Estadístico D-W","VP rho>0","VP rho<0")
res
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

#Función para calcular la amplitud de los I.P
amplitud=function(LIP,LSP){
a=LSP-LIP
am=mean(a)
am
}

#Función para calcular la cobertura de los I.P 
cobertura=function(real,LIP,LSP){
I=ifelse(real>=LIP & real<=LSP,1,0)
p=mean(I)
p
}


#LECTURA DE LOS DATOS
yt=scan()
589    561    640    656    727    697    640    599    568    577    553    582
600    566    653    673    742    716    660    617    583    587    565    598
628    618    688    705    770    736    678    639    604    611    594    634
658    622    709    722    782    756    702    653    615    621    602    635
677    635    736    755    811    798    735    697    661    667    645    688
713    667    762    784    837    817    767    722    681    687    660    698
717    696    775    796    858    826    783    740    701    706    677    711
734    690    785    805    871    845    801    764    725    723    690    734
750    707    807    824    886    859    819    783    740    747    711    751
804    756    860    878    942    913    869    834    790    800    763    800
826    799    890    900    961    935    894    855    809    810    766    805
821    773    883    898    957    924    881    837    784    791    760    802
828    778    889    902    969    947    908    867    815    812    773    813
834    782    892    903    966    937    896    858    817    827    797    843

yt=ts(yt,frequency=12, start=c(1962,1)) #serie con todos los datos

#Grafico la serie de tiempo
plot(yt,main="Produccion de leche mensual 1962.1-1975.12")

#Grafico de Descomposición en estacionalidad, tendencia y error
plot(decompose(yt))

#Boxplots según meses año calendario y periodograma
boxplot(yt~cycle(yt),names=month.abb)
periodogram(diff(yt),lwd=3);abline(v=c(1:6)/12,col=2,lty=2)

#FAC sobre la serie de producción de leche
win.graph(width=4.8,height=4.8,pointsize=8) 
acf(as.numeric(yt),ci.type="ma",lag.max=48,main="ACF de la serie de producción mendual de leche")

#Definiendo variables para ajuste sin los últimos 12 meses
n=length(yt)-12
t=1:n; t2=t^2;t3=t^3
yt2=ts(yt[t],frequency=12, start=c(1962,1)) #serie recortada
mes=seasonaldummy(yt2)
I1=mes[,1]
I2=mes[,2]
I3=mes[,3]
I4=mes[,4]
I5=mes[,5]
I6=mes[,6]
I7=mes[,7]
I8=mes[,8]
I9=mes[,9]
I10=mes[,10]
I11=mes[,11]

X3=cbind(t,t2,t3,mes) #matriz de regresión que será usada en modelos de regresión
                      #con errores ARMA, en los períodos de ajuste

#Definiendo variables para pronósticos
tnuevo=(n+1):length(yt)
t2nuevo=tnuevo^2
t3nuevo=tnuevo^3
mesnuevo=seasonaldummy(yt2,h=12)
#Separando una a una las 11 indicadoras para los tiempos de pronóstico
I1n=mesnuevo[,1]
I2n=mesnuevo[,2]
I3n=mesnuevo[,3]
I4n=mesnuevo[,4]
I5n=mesnuevo[,5]
I6n=mesnuevo[,6]
I7n=mesnuevo[,7]
I8n=mesnuevo[,8]
I9n=mesnuevo[,9]
I10n=mesnuevo[,10]
I11n=mesnuevo[,11]

X3nuevo=cbind(t=tnuevo,t2=t2nuevo,t3=t3nuevo,mes=mesnuevo) #matriz de regresión con valores de predictores
                                                           #para pronósticos en modelos con errores ARMA

ytnuevo=ts(yt[tnuevo],freq=12,start=c(1975,1)) #Serie de los últimos 12 valores observados
                                               #para realizar la validación cruzada


#Ajuste modelo 3 con tendencia cúbica  y errores R.B
modelo3=lm(yt2~t+t2+t3+I1+I2+I3+I4+I5+I6+I7+I8+I9+I10+I11); summary(modelo3)

yhat3=ts(fitted(modelo3),freq=12,start=c(1962,1))

plot(yt)
lines(yhat3,col=2)
legend("topleft",legend=c("Real","Ajuste modelo3"),lwd=2,col=c(1,2))

#Medidas de ajuste para modelo 3 con exp(C^*_n(p))
aic3=exp(crit.inf.resid(residuales=residuals(modelo3),n.par=15));aic3
bic3=exp(crit.inf.resid(residuales=residuals(modelo3),n.par=15,AIC="FALSE"));bic3

#Pronósticos modelo 3
pred3=predict(modelo3,newdata=data.frame(t=tnuevo,t2=t2nuevo,t3=t3nuevo,I1=I1n,I2=I2n,I3=I3n,I4=I4n,I5=I5n,I6=I6n,I7=I7n,I8=I8n,I9=I9n,I10=I10n,I11=I11n),interval="prediction")
pred3=ts(pred3,freq=12,start=c(1975,1))
pred3

#Amplitud media de los I.P
ampmod3=amplitud(LIP=pred3[,2],LSP=pred3[,3])
ampmod3

#cobertura de los I.P
cobmod3=cobertura(real=ytnuevo,LIP=pred3[,2],LSP=pred3[,3])
cobmod3

#precisión pronósticos puntuales
accuracy(pred3[,1],ytnuevo)

#Gráfico de residuales vs tiempo modelo3
win.graph(width=4.8,height=4.8,pointsize=8) 
plot(t,residuals(modelo3),type="o"); abline(h=0,col=2)
abline(h=c(-2*summary(modelo3)$sigma,2*summary(modelo3)$sigma),col=2)

#Gráfico de residuales vs predichos modelo3
win.graph(width=4.8,height=4.8,pointsize=8)
plot(fitted(modelo3), residuals(modelo3), type="p"); abline(h=0,col=2)
abline(h=c(-2*summary(modelo3)$sigma,2*summary(modelo3)$sigma),col=2)


#ACF o FAC residuales modelo3
win.graph(width=4.8,height=4.8,pointsize=8) 
acf(residuals(modelo3),ci.type="ma",lag.max=48,main="ACF modelo3",ci.col=2)

#PACF sobre residuales modelo3
win.graph(width=4.8,height=4.8,pointsize=8) 
pacf(residuals(modelo3),lag.max=48,main="PACF modelo3",ci.col=2)


#Usando la función usuario BP.LB.test()para obtener test B-P y L-B de ruido blanco
#Para los residuales modelo 3
BP.LB.test(residuals(modelo3),maxlag=48,type="Box")
BP.LB.test(residuals(modelo3),maxlag=48,type="Ljung")

#Usando la función usuario pruebaDW1 sobre residuales modelo 3
pruebaDW1(modelo3)


#######################MODELOS CON TENDENCIA CÚBICA ESTACIONALIDAD Y ERRORES ARMA###################### 

#Identificación con criterio AIC
auto.arima(residuals(modelo3),ic="aic")

#lo siguiente puede arrojar un modelo SARMA
et3=ts(residuals(modelo3),frequency=12, start=c(1962,1)) #serie de tiempo residuos modelo 3
auto.arima(et3,ic="aic")

SelectModel(et3, lag.max = 36, Criterion="AIC",ARModel = "AR")
SelectModel(et3, lag.max = 36, Criterion="BIC",ARModel = "AR")


#Selección con criterio BIC, TABLERO DE 12x12
plot(armasubsets(et3,nar=12,nma=12,y.name='AR',ar.method='ml'))

#Selección con EACF
eacf(residuals(modelo3),ar.max=24,ma.max=24)

####MODELOS CON ERRORES ARMA.################################### 
#Modelo con errores AR(14)
modelo3b=Arima(yt2,order=c(14,0,0),xreg=X3,method="ML") 
dfb=n-29 #n-Total de parámetros 
coeftest(modelo3b,df=dfb)

yhat3b=modelo3b$fitted #para obtener valores ajustados
                                                                   
plot(yt)
lines(yhat3b,col=2) 
legend("topleft",legend=c("Real","ajustada modelo3b"),col=c(1,2),lty=1)

aic3b=exp(crit.inf.resid(residuales=residuals(modelo3b),n.par=29));aic3b
bic3b=exp(crit.inf.resid(residuales=residuals(modelo3b),n.par=29,AIC="FALSE"));bic3b

#Gráfico de residuales de ajuste vs tiempo
win.graph(width=4.8,height=4.8,pointsize=8) 
plot(t,residuals(modelo3b),type="o"); abline(h=0,col=2)
abline(h=c(-2*sqrt(modelo3b$sigma2),2*sqrt(modelo3b$sigma2)),col=2)

#Gráfico de residuales de ajuste vs ajustados
win.graph(width=4.8,height=4.8,pointsize=8)
plot(yhat3b, residuals(modelo3b), type="p"); abline(h=0,col=2)
abline(h=c(-2*sqrt(modelo3b$sigma2),2*sqrt(modelo3b$sigma2)),col=2)

#FAC sobre residuales de ajuste modelo3b
win.graph(width=4.8,height=4.8,pointsize=8) 
acf(as.numeric(residuals(modelo3b)),ci.type="ma",lag.max=48,main="ACF modelo3b",ci.col=2)

#PACF sobre residuales de ajuste modelo3b
win.graph(width=4.8,height=4.8,pointsize=8) 
pacf(as.numeric(residuals(modelo3b)),lag.max=48,main="PACF modelo3b",ci.col=2)

#Box-pierce & Ljung-Box
BP.LB.test(residuals(modelo3b),maxlag=48,type="Box")
BP.LB.test(residuals(modelo3b),maxlag=48,type="Ljung")

shapiro.test(residuals(modelo3b))
win.graph(width=4.8,height=4.8,pointsize=8)
qqnorm(residuals(modelo3b),main="Gráfico de normalidad residuos modelo3b")
qqline(residuals(modelo3b),col=2)

#Pronósticos e I.P del 95%
pred3b=ts(as.data.frame(forecast(modelo3b,xreg=X3nuevo,level=95)),freq=12,start=c(1975,1)) #Matriz de regresión en pronósticos es X3nuevo

#Precisión pronósticos puntuales
accuracy(pred3b[,1],ytnuevo)

#Precisión pronósticos por I.P
ampmod3b=amplitud(LIP=pred3b[,2],LSP=pred3b[,3])
ampmod3b

cobmod3b=cobertura(real=ytnuevo,LIP=pred3b[,2],LSP=pred3b[,3])
cobmod3b


#Modelo con errores ARMA(2,10) pero sólo con phi1, phi2 y theta10
modelo3c=Arima(yt2,order=c(2,0,10),fixed=c(NA,NA,rep(0,9),NA,rep(NA,15)),xreg=X3,method="ML")
dfc=n-18 #n-Total parámetros
coeftest(modelo3c,df=dfc)

#Calculando valores ajustados
yhat3c=modelo3c$fitted
                                                                  
plot(yt)
lines(yhat3c,col=2)
legend("topleft",legend=c("Real","ajustada modelo3c"),col=c(1,2),lty=1)

aic3c=exp(crit.inf.resid(residuales=residuals(modelo3c),n.par=18));aic3c
bic3c=exp(crit.inf.resid(residuales=residuals(modelo3c),n.par=18,AIC="FALSE"));bic3c

#Gráfico de residuales de ajuste vs tiempo
win.graph(width=4.8,height=4.8,pointsize=8) 
plot(t,residuals(modelo3c),type="o"); abline(h=0,col=2)
abline(h=c(-2*sqrt(modelo3c$sigma2),2*sqrt(modelo3c$sigma2)),col=2)

#Gráfico de residuales de ajuste vs ajustados
win.graph(width=4.8,height=4.8,pointsize=8)
plot(yhat3c, residuals(modelo3c), type="p"); abline(h=0,col=2)
abline(h=c(-2*sqrt(modelo3c$sigma2),2*sqrt(modelo3c$sigma2)),col=2)

#FAC sobre residuales de ajuste modelo3c
win.graph(width=4.8,height=4.8,pointsize=8) 
acf(as.numeric(residuals(modelo3c)),ci.type="ma",lag.max=48,main="ACF modelo3c",ci.col=2)

#PACF sobre residuales de ajuste modelo3c
win.graph(width=4.8,height=4.8,pointsize=8) 
pacf(as.numeric(residuals(modelo3c)),lag.max=48,main="PACF modelo3c",ci.col=2)

#Box-pierce & Ljung-Box
BP.LB.test(residuals(modelo3c),maxlag=48,type="Box")
BP.LB.test(residuals(modelo3c),maxlag=48,type="Ljung")

shapiro.test(residuals(modelo3c))
win.graph(width=4.8,height=4.8,pointsize=8)
qqnorm(residuals(modelo3c),main="Gráfico de normalidad residuos modelo3c")
qqline(residuals(modelo3c),col=2)

#Pronósticos e I.P del 95%
pred3c=ts(as.data.frame(forecast(modelo3c,xreg=X3nuevo,level=95)),freq=12,start=c(1975,1)) #matriz de regresión en pronósticos es X3nuevo
pred3c

#Precisión pronósticos puntuales
accuracy(pred3c[,1],ytnuevo)

#Precisión pronósticos por I.P
ampmod3c=amplitud(LIP=pred3c[,2],LSP=pred3c[,3])
ampmod3c

cobmod3c=cobertura(real=ytnuevo,LIP=pred3c[,2],LSP=pred3c[,3])
cobmod3c

#Modelo con errores ARMA(3,2)
modelo3d=Arima(yt2,order=c(3,0,2),xreg=X3,method="ML")
dfd=n-20 #n-Total parámetros
coeftest(modelo3d,df=dfd)

#Cálculo valores ajustados
yhat3d=modelo3d$fitted

plot(yt)
lines(yhat3d,col=2)
legend("topleft",legend=c("Real","ajustada modelo3d"),col=c(1,2),lty=1)

aic3d=exp(crit.inf.resid(residuales=residuals(modelo3d),n.par=20));aic3d
bic3d=exp(crit.inf.resid(residuales=residuals(modelo3d),n.par=20,AIC="FALSE"));bic3d


#Gráfico de residuales de ajuste vs tiempo
win.graph(width=4.8,height=4.8,pointsize=8) 
plot(t,residuals(modelo3d),type="o"); abline(h=0,col=2)
abline(h=c(-2*sqrt(modelo3d$sigma2),2*sqrt(modelo3d$sigma2)),col=2)

#Gráfico de residuales de ajuste vs ajustados
win.graph(width=4.8,height=4.8,pointsize=8)
plot(yhat3d, residuals(modelo3d), type="p"); abline(h=0,col=2)
abline(h=c(-2*sqrt(modelo3d$sigma2),2*sqrt(modelo3d$sigma2)),col=2)

#FAC sobre residuales de ajuste modelo3d
win.graph(width=4.8,height=4.8,pointsize=8) 
acf(as.numeric(residuals(modelo3d)),ci.type="ma",lag.max=48,main="ACF modelo3d",ci.col=2)

#PACF sobre residuales de ajuste modelo3d
win.graph(width=4.8,height=4.8,pointsize=8) 
pacf(as.numeric(residuals(modelo3d)),lag.max=48,main="PACF modelo3d",ci.col=2)

#Box-pierce & Ljung-Box
BP.LB.test(residuals(modelo3d),maxlag=48,type="Box")
BP.LB.test(residuals(modelo3d),maxlag=48,type="Ljung")

shapiro.test(residuals(modelo3d))
win.graph(width=4.8,height=4.8,pointsize=8)
qqnorm(residuals(modelo3d),main="Gráfico de normalidad residuos modelo3d")
qqline(residuals(modelo3d),col=2)

#Pronósticos e I.P del 95%
pred3d=ts(as.data.frame(forecast(modelo3d,xreg=X3nuevo,level=95)),freq=12,start=c(1975,1)) #matriz de regresión en pronósticos es X3nuevo
pred3d

#Precisión pronósticos puntuales
accuracy(pred3d[,1],ytnuevo)

#Precisión pronósticos por I.P
ampmod3d=amplitud(LIP=pred3d[,2],LSP=pred3d[,3])
ampmod3d

cobmod3d=cobertura(real=ytnuevo,LIP=pred3d[,2],LSP=pred3d[,3])
cobmod3d

#Modelo con errores ARMA(4,0)xARMA(2,0)[12]
modelo3e=Arima(yt2,order=c(4,0,0),seasonal=list(order=c(2,0,0)),xreg=X3,method="ML")
dfe=n-21 #n-Total parámetros
coeftest(modelo3e, df=dfe)

#Cálculo valores ajustados
yhat3e=modelo3e$fitted

plot(yt)
lines(yhat3e,col=2)
legend("topleft",legend=c("Real","ajustada modelo3e"),col=c(1,2),lty=1)

aic3e=exp(crit.inf.resid(residuales=residuals(modelo3e),n.par=21));aic3e
bic3e=exp(crit.inf.resid(residuales=residuals(modelo3e),n.par=21,AIC="FALSE"));bic3e

#Gráfico de residuales de ajuste vs tiempo
win.graph(width=4.8,height=4.8,pointsize=8) 
plot(t,residuals(modelo3e),type="o"); abline(h=0,col=2)
abline(h=c(-2*sqrt(modelo3e$sigma2),2*sqrt(modelo3e$sigma2)),col=2)

#Gráfico de residuales de ajuste vs ajustados
win.graph(width=4.8,height=4.8,pointsize=8)
plot(yhat3e, residuals(modelo3e), type="p"); abline(h=0,col=2)
abline(h=c(-2*sqrt(modelo3e$sigma2),2*sqrt(modelo3e$sigma2)),col=2)

#FAC sobre residuales de ajuste modelo3e
win.graph(width=4.8,height=4.8,pointsize=8) 
acf(as.numeric(residuals(modelo3e)),ci.type="ma",lag.max=48,main="ACF modelo3e",ci.col=2)

#PACF sobre residuales de ajuste modelo3e
win.graph(width=4.8,height=4.8,pointsize=8) 
pacf(as.numeric(residuals(modelo3e)),lag.max=48,main="PACF modelo3e",ci.col=2)

#Box-pierce & Ljung-Box
BP.LB.test(residuals(modelo3e),maxlag=48,type="Box")
BP.LB.test(residuals(modelo3e),maxlag=48,type="Ljung")

shapiro.test(residuals(modelo3e))
win.graph(width=4.8,height=4.8,pointsize=8)
qqnorm(residuals(modelo3e),main="Gráfico de normalidad residuos modelo3e")
qqline(residuals(modelo3e),col=2)

#Pronósticos e I.P del 95%
pred3e=ts(as.data.frame(forecast(modelo3e,xreg=X3nuevo,level=95)),freq=12,start=c(1975,1)) #matriz de regresión en pronósticos es X3nuevo
pred3e

#Precisión pronósticos puntuales
accuracy(pred3e[,1],ytnuevo)

#Precisión pronósticos por I.P
ampmod3e=amplitud(LIP=pred3e[,2],LSP=pred3e[,3])
ampmod3e

cobmod3e=cobertura(real=ytnuevo,LIP=pred3e[,2],LSP=pred3e[,3])
cobmod3e

#Gráfico comparativo de pronósticos ex-post
plot(ytnuevo,type="b",pch=19,col=1,lwd=2,xaxt="n",ylim=c(750,1000))
axis(1,at=time(ytnuevo),labels=c("I-75","II-75","III-75","IV-75","V-75","VI-75","VII-75","VIII-75","IX-75","X-75","XI-75","XII-75"),cex.axis=0.7)
lines(pred3[,1],col=2,type="b",pch=2,lwd=2)
lines(pred3b[,1],col=4,type="b",pch=4,lwd=2)
lines(pred3c[,1],col=5,type="b",pch=5,lwd=2)
lines(pred3d[,1],col=6,type="b",pch=6,lwd=2)
lines(pred3e[,1],col=7,type="b",pch=7,lwd=2)
legend("topright",legend=c("Real","pron. modelo3","pron. modelo3b","pron. modelo3c","pron. modelo3d","pron. modelo3e"),col=c(1,2,4,5,6,7),pch=c(19,2,4,5,6,7),lwd=2)

#Tabla con todas las medidas de ajuste
tablamedidas=data.frame(p=c(15,29,18,20,21),AIC=c(aic3,aic3b,aic3c,aic3d,aic3e),BIC=c(bic3,bic3b,bic3c,bic3d,bic3e),row.names=c("modelo3","modelo3b","modelo3c","modelo3d","modelo3e"))
tablamedidas

#Tabla con todas las medidas de precisión de pronósticos
tablaprec=data.frame(AmplitudI.P=c(ampmod3,ampmod3b,ampmod3c,ampmod3d,ampmod3e),Cobertura=c(cobmod3,cobmod3b,cobmod3c,cobmod3d,cobmod3d)*100,
RMSE=c(accuracy(pred3[,1],ytnuevo)[2],accuracy(pred3b[,1],ytnuevo)[2],accuracy(pred3c[,1],ytnuevo)[2],accuracy(pred3d[,1],ytnuevo)[2],accuracy(pred3e[,1],ytnuevo)[2]),
MAE=c(accuracy(pred3[,1],ytnuevo)[3],accuracy(pred3b[,1],ytnuevo)[3],accuracy(pred3c[,1],ytnuevo)[3],accuracy(pred3d[,1],ytnuevo)[3],accuracy(pred3e[,1],ytnuevo)[3]),
MAPE=c(accuracy(pred3[,1],ytnuevo)[5],accuracy(pred3b[,1],ytnuevo)[5],accuracy(pred3c[,1],ytnuevo)[5],accuracy(pred3d[,1],ytnuevo)[5],accuracy(pred3e[,1],ytnuevo)[5]),
row.names=c("modelo3","modelo3b","modelo3c","modelo3d","modelo3e"))
tablaprec