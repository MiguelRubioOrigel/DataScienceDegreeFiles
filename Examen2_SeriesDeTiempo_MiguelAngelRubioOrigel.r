library(forecast)
library(tsoutliers)
library(tseries)
library(quantmod)
library(fUnitRoots)
library(fGarch)

ds <- read.csv("/uploads/iteso-2examen-2021.csv")
A <- ts(ds$V4)
head(A)
ts.plot(A)
####################################################################PRUEBAS EN NIVELES
adf.test(A)
#	Augmented Dickey-Fuller Test

#data:  A
#Dickey-Fuller = -5.5866, Lag order = 5, p-value = 0.01
#alternative hypothesis: stationary

#Warning message:
#In adf.test(A) : p-value smaller than printed p-value
#Existe suficiente evidencia estadistica para rechazar H0 (por lo tanto, la evidencia indica la serie SI es estacionaria)
kpss.test(A)

#	KPSS Test for Level Stationarity

#data:  A
#KPSS Level = 0.074001, Truncation lag parameter = 4, p-value =0.1

#Warning message:
#In kpss.test(A) : p-value greater than printed p-value
#No hay suficiente evidencia estadistica para rechazar H0 (por lo tanto, la evidencia indica la serie SI es estacionaria)
####################################################################IDENTIFICACION GARCH(p,q)
mean(A)
#1.907147
B <- (A-mean(A))^2
adf.test(B)
#	Augmented Dickey-Fuller Test

#data:  B
#Dickey-Fuller = -3.8565, Lag order = 5, p-value = 0.01744
#alternative hypothesis: stationary
#Existe suficiente evidencia estadistica al 95% de confianza para rechazar H0 (por lo tanto, la evidencia indica la serie SI es estacionaria)
kpss.test(B)
#	KPSS Test for Level Stationarity

#data:  B
#KPSS Level = 0.52296, Truncation lag parameter = 4, p-value =0.0365
#Existe suficiente evidencia estadistica al 95% de confianza para rechazar H0 (por lo tanto, la evidencia indica la serie NO es estacionaria)
#Por lo tanto, hay una contradiccion entre Dicky-Fuller y KPSS

pacf(B)
#indica rezagos 1, 2 y 5 son significativos, por lo tanto ARCH = 5
acf(B)
#confirma rezagos 1, 2 y 5 son significativos, por lo tanto ARCH = 5

#Entonces p + q = 5
archTest <- function(rtn,m=10){
  # Perform Lagrange Multiplier Test for ARCH effect of a time series
  # rtn: time series
  # m: selected AR order
  # Source:, 
  y=(rtn-mean(rtn))^2
  T=length(rtn)
  atsq=y[(m+1):T]
  x=matrix(0,(T-m),m)
  for (i in 1:m){
    x[,i]=y[(m+1-i):(T-i)]
  }
  md=lm(atsq~x)
  summary(md)
}

archTest(A,20)
#Call:
#lm(formula = atsq ~ x)

#Residuals:
#    Min      1Q  Median      3Q     Max 
#-18.261  -4.389  -2.374   0.945  72.172 

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  3.170790   1.532993   2.068  0.04022 * 
#x1           0.148477   0.079081   1.878  0.06228 . 
#x2           0.103157   0.080033   1.289  0.19930   
#x3           0.045683   0.080468   0.568  0.57103   
#x4          -0.036111   0.081147  -0.445  0.65692   
#x5           0.238050   0.083143   2.863  0.00476 **
#x6           0.034330   0.084930   0.404  0.68660   
#x7           0.015108   0.085320   0.177  0.85968   
#x8          -0.007067   0.086447  -0.082  0.93495   
#x9           0.087659   0.099939   0.877  0.38174   
#x10         -0.077519   0.100752  -0.769  0.44280   
#x11         -0.026865   0.100773  -0.267  0.79013   
#x12          0.016559   0.100748   0.164  0.86965   
#x13         -0.047605   0.100834  -0.472  0.63749   
#x14         -0.095884   0.102345  -0.937  0.35025   
#x15          0.081648   0.103000   0.793  0.42914   
#x16         -0.011377   0.100336  -0.113  0.90987   
#x17          0.011886   0.100441   0.118  0.90595   
#x18         -0.026811   0.100574  -0.267  0.79013   
#x19         -0.013003   0.099249  -0.131  0.89593   
#x20          0.082036   0.098909   0.829  0.40812   

#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 11.33 on 159 degrees of freedom
#Multiple R-squared:  0.1291,	Adjusted R-squared:  0.01959 
#F-statistic: 1.179 on 20 and 159 DF,  p-value: 0.279

#archTest confirma p+q = 5
m1=garchFit(~1+garch(5,0),data=A,trace=F)
summary(m1)

#Title:
# GARCH Modelling 

#Call:
# garchFit(formula = ~1 + garch(5, 0), data = A, trace = F) 

#Mean and Variance Equation:
# data ~ 1 + garch(5, 0)
#<environment: 0x5558a065f958>
# [data = A]

#Conditional Distribution:
# norm 

#Coefficient(s):
#      mu     omega    alpha1    alpha2    alpha3    alpha4    alpha5  
#2.044175  0.669424  0.165208  0.549883  0.112158  0.130788  0.092646  

#Std. Errors:
# based on Hessian 

#Error Analysis:
#        Estimate  Std. Error  t value Pr(>|t|)    
#mu       2.04418     0.10585   19.312  < 2e-16 ***
#omega    0.66942     0.30405    2.202 0.027688 *  
#alpha1   0.16521     0.11498    1.437 0.150759    
#alpha2   0.54988     0.15526    3.542 0.000398 ***
#alpha3   0.11216     0.09141    1.227 0.219813    
#alpha4   0.13079     0.07208    1.814 0.069601 .  
#alpha5   0.09265     0.08694    1.066 0.286574    

#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Log Likelihood:
# -434.0522    normalized:  -2.170261 

#Description:
# Sun May 16 00:58:07 2021 by user:  


#Standardised Residuals Tests:
#                                Statistic p-Value  
# Jarque-Bera Test   R    Chi^2  3.026388  0.2202055
# Shapiro-Wilk Test  R    W      0.9903461 0.2006356
# Ljung-Box Test     R    Q(10)  6.907568  0.7341413
# Ljung-Box Test     R    Q(15)  7.231055  0.950918 
# Ljung-Box Test     R    Q(20)  17.77342  0.6023315
# Ljung-Box Test     R^2  Q(10)  4.096374  0.9428947
# Ljung-Box Test     R^2  Q(15)  14.26269  0.5057069
# Ljung-Box Test     R^2  Q(20)  16.49212  0.685663 
# LM Arch Test       R    TR^2   5.722815  0.929405 

#Information Criterion Statistics:
#     AIC      BIC      SIC     HQIC 
#4.410522 4.525963 4.408180 4.457239 

#Parece que es significativo hasta el alpha 2, aunque me preocupa que parece haber residuales significativos en el rango Q(15) y Q(20) en los test de Jung Box

m2=garchFit(~1+garch(2,3),data=A,trace=F)
summary(m2)

#Title:
# GARCH Modelling 

#Call:
# garchFit(formula = ~1 + garch(2, 3), data = A, trace = F) 

#Mean and Variance Equation:
# data ~ 1 + garch(2, 3)
#<environment: 0x55589cdd0b60>
# [data = A]

#Conditional Distribution:
# norm 

#Coefficient(s):
#        mu       omega      alpha1      alpha2       beta1       beta2  
#2.04093190  0.38735427  0.18568502  0.52339581  0.19322447  0.15558738  
#     beta3  
#0.00000001  

#Std. Errors:
# based on Hessian 

#Error Analysis:
#        Estimate  Std. Error  t value Pr(>|t|)    
#mu     2.041e+00   1.063e-01   19.204  < 2e-16 ***
#omega  3.874e-01   2.312e-01    1.676  0.09380 .  
#alpha1 1.857e-01   1.286e-01    1.444  0.14872    
#alpha2 5.234e-01   1.844e-01    2.838  0.00453 ** 
#beta1  1.932e-01   2.687e-01    0.719  0.47207    
#beta2  1.556e-01   1.634e-01    0.952  0.34088    
#beta3  1.000e-08   2.643e-01    0.000  1.00000    

#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Log Likelihood:
# -434.2531    normalized:  -2.171265 

#Description:
# Sun May 16 01:10:20 2021 by user:  


#Standardised Residuals Tests:
#                                Statistic p-Value  
# Jarque-Bera Test   R    Chi^2  3.168283  0.2051238
# Shapiro-Wilk Test  R    W      0.9898671 0.1708594
# Ljung-Box Test     R    Q(10)  7.076063  0.7182452
# Ljung-Box Test     R    Q(15)  7.360131  0.9468753
# Ljung-Box Test     R    Q(20)  17.70589  0.6067757
# Ljung-Box Test     R^2  Q(10)  3.942378  0.949909 
# Ljung-Box Test     R^2  Q(15)  14.82052  0.4644208
# Ljung-Box Test     R^2  Q(20)  16.60476  0.6784686
# LM Arch Test       R    TR^2   5.981802  0.9169966

#Information Criterion Statistics:
#     AIC      BIC      SIC     HQIC 
#4.412531 4.527972 4.410189 4.459248

#GARCH(2,3) obtiene resultados muy similares a GARCH(5,0)

m3=garchFit(~1+garch(2,0),data=A,trace=F)
summary(m3)

#Title:
# GARCH Modelling 

#Call:
# garchFit(formula = ~1 + garch(2, 0), data = A, trace = F) 

#Mean and Variance Equation:
# data ~ 1 + garch(2, 0)
#<environment: 0x5558a7df4810>
# [data = A]

#Conditional Distribution:
# norm 

#Coefficient(s):
#     mu    omega   alpha1   alpha2  
#2.00010  1.85754  0.35851  0.46644  

#Std. Errors:
# based on Hessian 

#Error Analysis:
#        Estimate  Std. Error  t value Pr(>|t|)    
#mu        2.0001      0.1205   16.593  < 2e-16 ***
#omega     1.8575      0.4590    4.047 5.19e-05 ***
#alpha1    0.3585      0.1390    2.579  0.00992 ** 
#alpha2    0.4664      0.1597    2.920  0.00350 ** 

#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Log Likelihood:
# -443.7394    normalized:  -2.218697 

#Description:
# Sun May 16 01:14:47 2021 by user:  


#Standardised Residuals Tests:
#                                Statistic p-Value   
# Jarque-Bera Test   R    Chi^2  3.638044  0.1621843 
# Shapiro-Wilk Test  R    W      0.9916986 0.3113353 
# Ljung-Box Test     R    Q(10)  5.433401  0.8604126 
# Ljung-Box Test     R    Q(15)  6.087467  0.9782408 
# Ljung-Box Test     R    Q(20)  16.46734  0.6872412 
# Ljung-Box Test     R^2  Q(10)  11.55396  0.3160166 
# Ljung-Box Test     R^2  Q(15)  23.0616   0.08284079
# Ljung-Box Test     R^2  Q(20)  24.27568  0.2305571 
# LM Arch Test       R    TR^2   15.42878  0.2188208 

#Information Criterion Statistics:
#     AIC      BIC      SIC     HQIC 
#4.477394 4.543360 4.476614 4.504089 

#GARCH(2,0) no es buen modelo, nos baja mucho el p-value del Q(10) por lo que hay rezagos significativos en los primeros 10 que se estan ignorando, el AIC tambien empeora

m4=garchFit(~1+garch(2,2),data=A,trace=F)
summary(m4)

#Title:
# GARCH Modelling 

#Call:
# garchFit(formula = ~1 + garch(2, 2), data = A, trace = F) 

#Mean and Variance Equation:
# data ~ 1 + garch(2, 2)
#<environment: 0x5558a0efa9b8>
# [data = A]

#Conditional Distribution:
# norm 

#Coefficient(s):
#     mu    omega   alpha1   alpha2    beta1    beta2  
#2.04481  0.38438  0.18746  0.51893  0.19591  0.15493  

#Std. Errors:
# based on Hessian 

#Error Analysis:
#        Estimate  Std. Error  t value Pr(>|t|)    
#mu        2.0448      0.1056   19.364  < 2e-16 ***
#omega     0.3844      0.2297    1.673  0.09424 .  
#alpha1    0.1875      0.1205    1.556  0.11982    
#alpha2    0.5189      0.1740    2.983  0.00285 ** 
#beta1     0.1959      0.1827    1.072  0.28368    
#beta2     0.1549      0.1281    1.209  0.22664    

#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Log Likelihood:
# -434.1218    normalized:  -2.170609 

#Description:
# Sun May 16 01:24:06 2021 by user:  

#Standardised Residuals Tests:
#                                Statistic p-Value  
# Jarque-Bera Test   R    Chi^2  3.166056  0.2053524
# Shapiro-Wilk Test  R    W      0.989827  0.1685615
# Ljung-Box Test     R    Q(10)  7.069628  0.7188557
# Ljung-Box Test     R    Q(15)  7.348606  0.9472446
# Ljung-Box Test     R    Q(20)  17.68606  0.6080805
# Ljung-Box Test     R^2  Q(10)  3.973977  0.9485133
# Ljung-Box Test     R^2  Q(15)  14.76138  0.4687385
# Ljung-Box Test     R^2  Q(20)  16.52679  0.6834527
# LM Arch Test       R    TR^2   6.018227  0.9151605

#Information Criterion Statistics:
#     AIC      BIC      SIC     HQIC 
#4.401218 4.500168 4.399487 4.441262

#GARCH(2,2) baja el p-Value del Q(10) respecto del GARCH(2,3) por lo que hay rezagos significativos que se estan ignorando en ese rango

m5=garchFit(~1+garch(2,1),data=A,trace=F)
summary(m5)

#Title:
# GARCH Modelling 

#Call:
# garchFit(formula = ~1 + garch(2, 1), data = A, trace = F) 

#Mean and Variance Equation:
# data ~ 1 + garch(2, 1)
#<environment: 0x5595addd69b8>
# [data = A]

#Conditional Distribution:
# norm 

#Coefficient(s):
#     mu    omega   alpha1   alpha2    beta1  
#2.05202  0.33661  0.19259  0.42619  0.43374  

#Std. Errors:
# based on Hessian 

#Error Analysis:
#        Estimate  Std. Error  t value Pr(>|t|)    
#mu       2.05202     0.10716   19.149  < 2e-16 ***
#omega    0.33661     0.19277    1.746   0.0808 .  
#alpha1   0.19259     0.12018    1.602   0.1091    
#alpha2   0.42619     0.17618    2.419   0.0156 *  
#beta1    0.43374     0.09983    4.345 1.39e-05 ***

#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Log Likelihood:
# -434.5537    normalized:  -2.172768 

#Description:
# Sun May 16 04:36:59 2021 by user:  


#Standardised Residuals Tests:
#                                Statistic p-Value  
# Jarque-Bera Test   R    Chi^2  2.794299  0.2473009
# Shapiro-Wilk Test  R    W      0.9901973 0.1909092
# Ljung-Box Test     R    Q(10)  7.292691  0.6975543
# Ljung-Box Test     R    Q(15)  7.741908  0.9337125
# Ljung-Box Test     R    Q(20)  17.75424  0.6035939
# Ljung-Box Test     R^2  Q(10)  4.818476  0.9029693
# Ljung-Box Test     R^2  Q(15)  15.00119  0.4513312
# Ljung-Box Test     R^2  Q(20)  17.24698  0.6368826
# LM Arch Test       R    TR^2   7.017931  0.8564264

#Information Criterion Statistics:
#     AIC      BIC      SIC     HQIC 
#4.395537 4.477994 4.394327 4.428906

#Hasta ahorita el mejor modelo que encuentro es GARCH(2,3), ya que tiene los mejores resultados en L-Jung Box para el Q(10)
#Sin embargo, me preocupa que el GARCH(2,1) si muestre significativo el coaeficiente beta1

res1 <- residuals(m1,standardize=T)
res2 <- residuals(m2,standardize=T)
res3 <- residuals(m3,standardize=T)
res4 <- residuals(m4,standardize=T)
res5 <- residuals(m5,standardize=T)
tdx <- c(1:200)
plot(tdx,res1,xlab="m1",ylab="stand-resi",type="l")
plot(tdx,res2,xlab="m2",ylab="stand-resi",type="l")
plot(tdx,res3,xlab="m3",ylab="stand-resi",type="l")
plot(tdx,res4,xlab="m4",ylab="stand-resi",type="l")
plot(tdx,res5,xlab="m5",ylab="stand-resi",type="l")

acf(res1^2,lag=20)
pacf(res1^2,lag=20)

acf(res2^2,lag=20)
pacf(res2^2,lag=20)

acf(res3^2,lag=20)
pacf(res3^2,lag=20)
#el modelo 3 muestra el rezago 5 como significativo, lo cual ya era obvio desde que obtuvimos p+q=5

acf(res4^2,lag=20)
pacf(res4^2,lag=20)

acf(res5^2,lag=20)
pacf(res5^2,lag=20)
#Todos los modelos muestran el rezago 15 como significativo

m6=garchFit(~1+garch(2,3),data=A,trace=F,cond.dist="std")
m7=garchFit(~1+garch(2,3),data=A,trace=F,cond.dist="sstd")

plot(m2) #opcion 13 muestra apego a distribucion normal
plot(m6) #opcion 13 no muestra apego a distribucion T-Student
plot(m7) #opcion 13 no muestra apego a distribucion T-Student con sesgo

#CONCLUSIONES:  El mejor modelo me parece que es un GARCH(2,3) ya que el L-Jung box test para R^2 en los primeros 10 rezagos se acerca mas al 95%, 
#aunque me queda claro que el residual 15 parece significativo en todos los modelos, incluido este. 
#Sin embargo, el que GARCH(2,1) si muestre el beta1 significativo me preocupa, no se si eso es mas importante que los L-Jung Box test…
#Los residuales parecen apegarse mas a una distribución normal, excepto por la segunda cola que están fuera de la diagonal, 
#aun asi supera el apego a la diagonal en las T-Student.