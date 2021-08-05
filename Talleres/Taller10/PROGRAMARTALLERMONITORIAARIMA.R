library(forecast)
library(timsac)
library(TSA)
library(lmtest)

#Creando función usuario para obtener test Box-Pierce y Ljung-Box
BP.LB.test=function(serie,maxlag,type="Box"){
aux=floor(maxlag/6);
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

ARIMA1.1.0=ts(scan(file.choose(),skip=1),freq=1)#leer ARIMA110.SIMUL.txt
layout(rbind(c(1,1,2,2),c(3,3,4,4)))
plot(ARIMA1.1.0)
acf(ARIMA1.1.0,ci.type="ma")
plot(diff(ARIMA1.1.0))
abline(h=mean(diff(ARIMA1.1.0)))
abline(h=c(-2*sd(diff(ARIMA1.1.0)),2*sd(diff(ARIMA1.1.0))),lty=2)
acf(diff(ARIMA1.1.0),ci.type="ma")

####################################################
#Identificación y ajuste con los primeros 191 datos
####################################################
n=length(ARIMA1.1.0)-10
t=1:n
yt=ts(ARIMA1.1.0[t],freq=1) #serie con sólo 191 datos

layout(rbind(c(1,1,2,2),c(3,3,4,4)))
plot(diff(yt))
abline(h=mean(diff(ARIMA1.1.0)))
abline(h=c(-2*sd(diff(ARIMA1.1.0)),2*sd(diff(ARIMA1.1.0))),lty=2)
plot(time(diff(yt)),as.numeric(diff(yt)),xlab="Time")
abline(h=mean(diff(ARIMA1.1.0)))
abline(h=c(-2*sd(diff(ARIMA1.1.0)),2*sd(diff(ARIMA1.1.0))),lty=2)
acf(diff(yt),ci.type="ma")
pacf(diff(yt))
eacf(diff(yt))

win.graph()
plot(armasubsets(diff(yt),nar=12,nma=12,y.name='test',ar.method='ml'))

autoarmafit(diff(yt)) #identificando ARMA's sobre serie diferencia
#Identificación con auto.arima usando criterio de información AICC (criterio que usa por defecto)
auto.arima(diff(yt)) #auto.arima sobre serie diferencia identifica ARIMA(p,0,q)
auto.arima(yt) #auto.arima sobre serie sin diferenciar identifica ARIMA(p,d,q)

mean(diff(yt)) #verificando media

modelo1=Arima(yt,order=c(1,1,0),method="ML")
coeftest(modelo1) #valores P según N(0,1)

layout(rbind(c(1,2),c(3,4)))
plot.ts(residuals(modelo1))
abline(h=c(-2*sqrt(modelo1$sigma2),0,2*sqrt(modelo1$sigma2)),col=2)
plot(modelo1$fitted,residuals(modelo1))
abline(h=c(-2*sqrt(modelo1$sigma2),0,2*sqrt(modelo1$sigma2)),col=2)
acf(residuals(modelo1),ci.type="ma")
pacf(residuals(modelo1))

BP.LB.test(residuals(modelo1),maxlag=36,type="Ljung")
shapiro.test(residuals(modelo1))


modelo2=Arima(yt,order=c(1,1,0),include.drift=TRUE,method="ML") 
coeftest(modelo2)

layout(rbind(c(1,2),c(3,4)))
plot.ts(residuals(modelo2))
abline(h=c(-2*sqrt(modelo2$sigma2),0,2*sqrt(modelo2$sigma2)),col=2)
plot(modelo2$fitted,residuals(modelo2))
abline(h=c(-2*sqrt(modelo2$sigma2),0,2*sqrt(modelo2$sigma2)),col=2)
acf(residuals(modelo2),ci.type="ma")
pacf(residuals(modelo2))

BP.LB.test(residuals(modelo2),maxlag=36,type="Ljung")
shapiro.test(residuals(modelo2))

modelo3=Arima(yt,order=c(2,1,1),method="ML")
coeftest(modelo3)

layout(rbind(c(1,2),c(3,4)))
plot.ts(residuals(modelo3))
abline(h=c(-2*sqrt(modelo3$sigma2),0,2*sqrt(modelo3$sigma2)),col=2)
plot(modelo3$fitted,residuals(modelo3))
abline(h=c(-2*sqrt(modelo3$sigma2),0,2*sqrt(modelo3$sigma2)),col=2)
acf(residuals(modelo3),ci.type="ma")
pacf(residuals(modelo3))

BP.LB.test(residuals(modelo3),maxlag=36,type="Ljung")
shapiro.test(residuals(modelo3))

modelo3b=Arima(yt,order=c(2,1,1),include.drift=TRUE,method="ML")
coeftest(modelo3b)

layout(rbind(c(1,2),c(3,4)))
plot.ts(residuals(modelo3b))
abline(h=c(-2*sqrt(modelo3b$sigma2),0,2*sqrt(modelo3b$sigma2)),col=2)
plot(modelo3b$fitted,residuals(modelo3b))
abline(h=c(-2*sqrt(modelo3b$sigma2),0,2*sqrt(modelo3b$sigma2)),col=2)
acf(residuals(modelo3b),ci.type="ma")
pacf(residuals(modelo3b))

BP.LB.test(residuals(modelo3b),maxlag=36,type="Ljung")
shapiro.test(residuals(modelo3b))

modelo4=Arima(yt,order=c(1,1,1),method="ML")
coeftest(modelo4)

layout(rbind(c(1,2),c(3,4)))
plot.ts(residuals(modelo4))
abline(h=c(-2*sqrt(modelo4$sigma2),0,2*sqrt(modelo4$sigma2)),col=2)
plot(modelo4$fitted,residuals(modelo4))
abline(h=c(-2*sqrt(modelo4$sigma2),0,2*sqrt(modelo4$sigma2)),col=2)
acf(residuals(modelo4),ci.type="ma")
pacf(residuals(modelo4))

BP.LB.test(residuals(modelo4),maxlag=36,type="Ljung")
shapiro.test(residuals(modelo4))

modelo4b=Arima(yt,order=c(1,1,1),include.drift=TRUE,method="ML")
coeftest(modelo4b)

layout(rbind(c(1,2),c(3,4)))
plot.ts(residuals(modelo4b))
abline(h=c(-2*sqrt(modelo4b$sigma2),0,2*sqrt(modelo4b$sigma2)),col=2)
plot(modelo4b$fitted,residuals(modelo4b))
abline(h=c(-2*sqrt(modelo4b$sigma2),0,2*sqrt(modelo4b$sigma2)),col=2)
acf(residuals(modelo4b),ci.type="ma")
pacf(residuals(modelo4b))

BP.LB.test(residuals(modelo4b),maxlag=36,type="Ljung")
shapiro.test(residuals(modelo4b))

modelo5=Arima(yt,order=c(2,1,2),method="ML")
coeftest(modelo5)

layout(rbind(c(1,2),c(3,4)))
plot.ts(residuals(modelo5))
abline(h=c(-2*sqrt(modelo5$sigma2),0,2*sqrt(modelo5$sigma2)),col=2)
plot(modelo5$fitted,residuals(modelo5))
abline(h=c(-2*sqrt(modelo5$sigma2),0,2*sqrt(modelo5$sigma2)),col=2)
acf(residuals(modelo5),ci.type="ma")
pacf(residuals(modelo5))

BP.LB.test(residuals(modelo5),maxlag=36,type="Ljung")
shapiro.test(residuals(modelo5))

modelo5b=Arima(yt,order=c(2,1,2),include.drift=T,method="ML")
coeftest(modelo5b)

layout(rbind(c(1,2),c(3,4)))
plot.ts(residuals(modelo5b))
abline(h=c(-2*sqrt(modelo5b$sigma2),0,2*sqrt(modelo5b$sigma2)),col=2)
plot(modelo5b$fitted,residuals(modelo5b))
abline(h=c(-2*sqrt(modelo5b$sigma2),0,2*sqrt(modelo5b$sigma2)),col=2)
acf(residuals(modelo5b),ci.type="ma")
pacf(residuals(modelo5b))

BP.LB.test(residuals(modelo5b),maxlag=36,type="Ljung")
shapiro.test(residuals(modelo5b))

win.graph(width=16,height=8)
layout(rbind(c(1,2,3,4),c(5,6,7,8)))
plot(ARIMA1.1.0)
lines(modelo1$fitted,col=2,lwd=2)
legend("bottomleft",legend=c("Original","Ajuste modelo 1"),col=1:2,lty=1)
plot(ARIMA1.1.0)
lines(modelo2$fitted,col=2,lwd=2)
legend("bottomleft",legend=c("Original","Ajuste modelo 2"),col=1:2,lty=1)
plot(ARIMA1.1.0)
lines(modelo3$fitted,col=2,lwd=2)
legend("bottomleft",legend=c("Original","Ajuste modelo 3"),col=1:2,lty=1)
plot(ARIMA1.1.0)
lines(modelo3b$fitted,col=2,lwd=2)
legend("bottomleft",legend=c("Original","Ajuste modelo 3b"),col=1:2,lty=1)
plot(ARIMA1.1.0)
lines(modelo4$fitted,col=2,lwd=2)
legend("bottomleft",legend=c("Original","Ajuste modelo 4"),col=1:2,lty=1)
plot(ARIMA1.1.0)
lines(modelo4b$fitted,col=2,lwd=2)
legend("bottomleft",legend=c("Original","Ajuste modelo 4b"),col=1:2,lty=1)
plot(ARIMA1.1.0)
lines(modelo5$fitted,col=2,lwd=2)
legend("bottomleft",legend=c("Original","Ajuste modelo 5"),col=1:2,lty=1)
plot(ARIMA1.1.0)
lines(modelo5b$fitted,col=2,lwd=2)
legend("bottomleft",legend=c("Original","Ajuste modelo 5b"),col=1:2,lty=1)

#Cálculo de medidas de bondad de ajuste
aic1=exp(crit.inf.resid(residuals(modelo1),n.par=1))
aic2=exp(crit.inf.resid(residuals(modelo2),n.par=2))
aic3=exp(crit.inf.resid(residuals(modelo3),n.par=3))
aic3b=exp(crit.inf.resid(residuals(modelo3b),n.par=4))
aic4=exp(crit.inf.resid(residuals(modelo4),n.par=2))
aic4b=exp(crit.inf.resid(residuals(modelo4b),n.par=3))
aic5=exp(crit.inf.resid(residuals(modelo5),n.par=4))
aic5b=exp(crit.inf.resid(residuals(modelo5b),n.par=5))

bic1=exp(crit.inf.resid(residuals(modelo1),n.par=1,AIC="FALSE"))
bic2=exp(crit.inf.resid(residuals(modelo2),n.par=2,AIC="FALSE"))
bic3=exp(crit.inf.resid(residuals(modelo3),n.par=3,AIC="FALSE"))
bic3b=exp(crit.inf.resid(residuals(modelo3b),n.par=4,AIC="FALSE"))
bic4=exp(crit.inf.resid(residuals(modelo4),n.par=2,AIC="FALSE"))
bic4b=exp(crit.inf.resid(residuals(modelo4b),n.par=3,AIC="FALSE"))
bic5=exp(crit.inf.resid(residuals(modelo5),n.par=4,AIC="FALSE"))
bic5b=exp(crit.inf.resid(residuals(modelo5b),n.par=5,AIC="FALSE"))

criterios=data.frame(p=c(1,2,3,4,2,3,4,5),AIC=c(aic1,aic2,aic3,aic3b,aic4,aic4b,aic5,aic5b),BIC=c(bic1,bic2,bic3,bic3b,bic4,bic4b,bic5,bic5b),row.names=c("modelo1","modelo2","modelo3","modelo3b","modelo4","modelo4b","modelo5","modelo5b"))
criterios

######################################################
#Predicciones de las últimas 10 observaciones y medidas de precisión pronósticos
######################################################
tnuevo=(n+1):length(ARIMA1.1.0) #valor índice de tiempo en los pronósticos
ynuevo=ts(ARIMA1.1.0[tnuevo],freq=1,start=192) #valor serie en las últimos 10 tiempo

pred1=forecast(modelo1,h=10,level=95)
pred1=cbind(Predic=pred1$mean,LIP95=pred1$lower,LSP95=pred1$upper)
pred2=forecast(modelo2,h=10,level=95)
pred2=cbind(Predic=pred2$mean,LIP95=pred2$lower,LSP95=pred2$upper)
pred3=forecast(modelo3,h=10,level=95)
pred3=cbind(Predic=pred3$mean,LIP95=pred3$lower,LSP95=pred3$upper)
pred3b=forecast(modelo3b,h=10,level=95)
pred3b=cbind(Predic=pred3b$mean,LIP95=pred3b$lower,LSP95=pred3b$upper)
pred4=forecast(modelo4,h=10,level=95)
pred4=cbind(Predic=pred4$mean,LIP95=pred4$lower,LSP95=pred4$upper)
pred4b=forecast(modelo4b,h=10,level=95)
pred4b=cbind(Predic=pred4b$mean,LIP95=pred4b$lower,LSP95=pred4b$upper)
pred5=forecast(modelo5,h=10,level=95)
pred5=cbind(Predic=pred5$mean,LIP95=pred5$lower,LSP95=pred5$upper)
pred5b=forecast(modelo5b,h=10,level=95)
pred5b=cbind(Predic=pred5b$mean,LIP95=pred5b$lower,LSP95=pred5b$upper)

plot(ynuevo,col="brown",lwd=2,type="b",pch=19,ylim=c(-73,-65))
lines(pred1[,1],col=2,lwd=2,type="b",pch=1)
lines(pred2[,1],col=3,lwd=2,type="b",pch=2)
lines(pred3[,1],col=4,lwd=2,type="b",pch=3)
lines(pred3b[,1],col=5,lwd=2,type="b",pch=4)
lines(pred4[,1],col=6,lwd=2,type="b",pch=5)
lines(pred4b[,1],col=7,lwd=2,type="b",pch=6)
lines(pred5[,1],col=8,lwd=2,type="b",pch=7)
lines(pred5b[,1],col=9,lwd=2,type="b",pch=8)
legend("topright",legend=c("Original","pred1","Pred2","Pred3","pred3b","Pred4","Pred4b","Pred5","Pred5b"),col=c("brown",2:9),pch=c(19,1:8),lty=1,lwd=2)


amp1=amplitud(LIP=pred1[,2],LSP=pred1[,3])
amp2=amplitud(LIP=pred2[,2],LSP=pred2[,3])
amp3=amplitud(LIP=pred3[,2],LSP=pred3[,3])
amp3b=amplitud(LIP=pred3b[,2],LSP=pred3b[,3])
amp4=amplitud(LIP=pred4[,2],LSP=pred4[,3])
amp4b=amplitud(LIP=pred4b[,2],LSP=pred4b[,3])
amp5=amplitud(LIP=pred5[,2],LSP=pred5[,3])
amp5b=amplitud(LIP=pred5b[,2],LSP=pred5b[,3])
cob1=cobertura(real=ynuevo,LIP=pred1[,2],LSP=pred1[,3])
cob2=cobertura(real=ynuevo,LIP=pred2[,2],LSP=pred2[,3])
cob3=cobertura(real=ynuevo,LIP=pred3[,2],LSP=pred3[,3])
cob3b=cobertura(real=ynuevo,LIP=pred3b[,2],LSP=pred3b[,3])
cob4=cobertura(real=ynuevo,LIP=pred4[,2],LSP=pred4[,3])
cob4b=cobertura(real=ynuevo,LIP=pred4b[,2],LSP=pred4b[,3])
cob5=cobertura(real=ynuevo,LIP=pred5[,2],LSP=pred5[,3])
cob5b=cobertura(real=ynuevo,LIP=pred5b[,2],LSP=pred5b[,3])

#Tabla con todas las medidas de precisión de pronósticos
tablaprec=data.frame(AmplitudI.P=c(amp1,amp2,amp3,amp3b,amp4,amp4b,amp5,amp5b),Cobertura=c(cob1,cob2,cob3,cob3b,cob4,cob4b,cob5,cob5b)*100,
RMSE=c(accuracy(pred1[,1],ynuevo)[2],accuracy(pred2[,1],ynuevo)[2],accuracy(pred3[,1],ynuevo)[2],accuracy(pred3b[,1],ynuevo)[2],accuracy(pred4[,1],ynuevo)[2],accuracy(pred4b[,1],ynuevo)[2],accuracy(pred5[,1],ynuevo)[2],accuracy(pred5b[,1],ynuevo)[2]),
MAE=c(accuracy(pred1[,1],ynuevo)[3],accuracy(pred2[,1],ynuevo)[3],accuracy(pred3[,1],ynuevo)[3],accuracy(pred3b[,1],ynuevo)[3],accuracy(pred4[,1],ynuevo)[3],accuracy(pred4b[,1],ynuevo)[3],accuracy(pred5[,1],ynuevo)[3],accuracy(pred5b[,1],ynuevo)[3]),
MAPE=c(accuracy(pred1[,1],ynuevo)[5],accuracy(pred2[,1],ynuevo)[5],accuracy(pred3[,1],ynuevo)[5],accuracy(pred3b[,1],ynuevo)[5],accuracy(pred4[,1],ynuevo)[5],accuracy(pred4b[,1],ynuevo)[5],accuracy(pred5[,1],ynuevo)[5],accuracy(pred5b[,1],ynuevo)[5]),
row.names=c("modelo1","modelo2","modelo3","modelo3b","modelo4","modelo4b","modelo5","modelo5b"))
round(tablaprec,digits=2)





