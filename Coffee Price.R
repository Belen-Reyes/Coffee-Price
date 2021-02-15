# COFFEE PRICE

# Packages ----

library(tseries)
library(ggplot2)
library(rugarch)
library(nortest)
library(dplyr)
library(FinTS)
library(fGarch)
library(stats)
library(forecast)
library(PerformanceAnalytics)
library(moments)
library(gridExtra)
library(dlm)
library(timeSeries)
library(rjags)
library(coda)
library(lmtest)
library(astsa)

# Data ----
data <- read.csv('C:/Users/Belen/Documents/Belén/Séptimo semestre/Bayesiana/Proyecto Final/Coffee.csv', sep=';')
head(data)
cafe <- data[,2]
cafe <- ts(cafe, frequency = 12, start = c(1990,1), end = c(2020,11))
# NA's
table(is.na(cafe))
```
# Descriptive Analysis ----
summary(cafe)

kurtosis_cafe <- kurtosis(cafe)
skewness_cafe <- skewness(cafe)
# Histogram
hist(cafe, breaks =20, main = "Price histogram", probability= T, col = "#528792", 
     border ="white", ylab = "Probability", xlab= "Price (USD)")
lines(density(cafe), col = "#CACACA", lwd = 2)
legend('topright', c('Kurtosis: 3.5', 'Skewness: 0.73'))
# Time series
plot(cafe, main = "Monthly coffee prices 1990-2020", col= "#6CB4AB", 
     ylab = "Price (USD)", xlab= "Year", lwd = 2)
# BoxPlot
boxplot(cafe~cycle(cafe), names = month.abb, col ="#B2D6A9", 
        main="Coffee prices according to the month", xlab = "Month", ylab = "Price")
# Time Series Decomposition ----
cafe_decomp<-decompose(cafe, type= "multiplicative")
trend <- cafe_decomp$trend
seasonal <- cafe_decomp$seasonal
season <- ts(seasonal[1:60], frequency = 12, start = c(1990,1), end = c(1994,12))
# Trend and seasonal
par(mfrow=c(1,2))
plot(trend, main = "Trend", col= "#EBB00B", ylab = " ", xlab= "Year", lwd =2)
plot(season, main = "Seasonal", col= "#EF7E16", ylab = " ", xlab= "Year", lwd =2)
par(mfrow=c(1,1))
ts.plot(cbind(cafe, trend, seasonal*trend), col=c("#BFBFBF", "#D1A4C2","#AF4551"),
        lty=c(1,1,2), lwd =c(1,2,2), 
        main ="Coffee prices (trend and seasonal)", xlab = "Year", 
        ylab = "Price")
legend("topleft", c("Original Series","Trend", "Seasonal"), 
       col=c("#BFBFBF", "#C29CBE","#CD4864"), lty=c(1,1,2), lwd =c(1,2,2), 
       cex=.7, bty = "n")
# Model Fit ----
# Stationarity 
adf.test(cafe)
pp.test(cafe)
kpss.test(cafe)
t.test(cafe)

# Number of differences
nsdiffs(cafe) 
ndiffs(cafe)

cafe_dif <- diff(cafe)
cafe_dif2 <- diff(cafe_dif, lag = 12)

# Stationarity tests
adf.test(cafe_dif2)
pp.test(cafe_dif2)
kpss.test(cafe_dif2)

x_dif2 <- 1:length(cafe_dif2)
y_dif2 <- as.numeric(cafe_dif2)
bptest(y_dif2~x_dif2)
t.test(cafe_dif2)

# Correlograms
par(mfrow=c(2,2), cex = 0.7)

acf(cafe, lag.max = 100, main = 'Coffee ACF', ylab = '', xlab = '')
pacf(cafe, lag.max = 100,  main = 'Coffee PACF', ylab = '', xlab = '')

acf(cafe_dif2, lag.max = 100,  main = 'Coffee ACF (with differences)', ylab = '')
pacf(cafe_dif2, lag.max = 100, main = 'Coffee PACF (with differences)', ylab = '')

# Model without transfomations ----

auto.arima(cafe)
auto.arima(cafe_dif2)

# BoxCox Transformation ----
lambda <- BoxCox.lambda(cafe)

cafe_boxcox <- cafe^lambda
nsdiffs(cafe_boxcox) 
ndiffs(cafe_boxcox)
cafe_boxcox_dif <- diff(cafe_boxcox)
cafe_boxcox_dif2 <- diff(cafe_boxcox_dif, lag = 12)

# Correlograms 
par(mfrow=c(1,2))
acf(cafe_boxcox_dif2, lag.max = 100,  main = 'Coffee ACF BoxCox', ylab = '')
pacf(cafe_boxcox_dif2, lag.max = 100, main = 'Coffee PACF BoxCox', ylab = '')

auto.arima(cafe_boxcox_dif2)
# Model Comparison ----

sarima_1 <- arima(cafe,order=c(1,1,2), seasonal = list(order=c(0,1,0),period = 12))
sarima_2 <- arima(cafe,order=c(2,1,0), seasonal = list(order=c(2,1,0),period = 12))
sarima_3 <- arima(cafe,order=c(0,1,2), seasonal = list(order=c(1,1,0),period = 12))
sarima_4 <- arima(cafe,order=c(0,1,2), seasonal = list(order=c(0,1,0),period = 12))
sarima_5 <- arima(cafe,order=c(0,1,2), seasonal = list(order=c(0,1,1),period = 12))
BC_sarima1 <- arima(cafe_boxcox,order=c(2,1,1), seasonal = list(order=c(0,1,0),period = 12))
BC_sarima2 <- arima(cafe_boxcox,order=c(0,1,1), seasonal = list(order=c(2,1,0),period = 12))

  # Sarima_1 
  summary(sarima_1)
  confint(sarima_1)
  BIC(sarima_1)
  AIC(sarima_1)

    # Assumption validation
    
    # Normality
    ad.test(sarima_1$residuals)
    shapiro.test(sarima_1$residuals)
    par(mfrow=c(1,1))
    qqnorm(sarima_1$residuals)
    qqline(sarima_1$residuals)
    # Homoscedasticity
    x_1 <- 1:length(sarima_1$residuals)
    y_1 <- as.numeric(sarima_1$residuals)
    bptest(y_1~x_1)
    plot(sarima_1$residuals)
    # Zero mean
    t.test(sarima_1$residuals)
    # No correlation
    checkresiduals(sarima_1$residuals)
    tsdiag(sarima_1)

  # Sarima_2 
  summary(sarima_2)
  confint(sarima_2)
  BIC(sarima_2)
  AIC(sarima_2)

    # Assumption validation
    
    # Normality
    ad.test(sarima_2$residuals)
    shapiro.test(sarima_2$residuals)
    qqnorm(sarima_2$residuals)
    qqline(sarima_2$residuals)
    # Homocedasticity
    x_2 <- 1:length(sarima_2$residuals)
    y_2 <- as.numeric(sarima_2$residuals)
    bptest(y_2~x_2)
    plot(sarima_2$residuals)
    # Zero mean
    t.test(sarima_2$residuals)
    # No correlation
    checkresiduals(sarima_2$residuals)
    tsdiag(sarima_2)

  # Sarima_3 
  summary(sarima_3)
  confint(sarima_3)
  BIC(sarima_3)
  AIC(sarima_3)
  
    # Assumption validation
    
    # Normality
    ad.test(sarima_3$residuals)
    shapiro.test(sarima_3$residuals)
    qqnorm(sarima_3$residuals)
    qqline(sarima_3$residuals)
    # Homocedasticity
    x_3 <- 1:length(sarima_3$residuals)
    y_3 <- as.numeric(sarima_3$residuals)
    bptest(y_3~x_3)
    plot(sarima_3$residuals)
    # Zero mean
    t.test(sarima_3$residuals)
    #No correlation
    checkresiduals(sarima_3$residuals)
    tsdiag(sarima_3)

  # Sarima_4 
  summary(sarima_4)
  confint(sarima_4)
  BIC(sarima_4)
  AIC(sarima_4)
  
    # Assumption validation
    
    # Normality
    ad.test(sarima_4$residuals)
    shapiro.test(sarima_4$residuals)
    qqnorm(sarima_4$residuals)
    qqline(sarima_4$residuals)
    # Homocedasticity
    x_4 <- 1:length(sarima_4$residuals)
    y_4 <- as.numeric(sarima_4$residuals)
    bptest(y_4~x_4)
    plot(sarima_4$residuals)
    # Zero mean
    t.test(sarima_4$residuals)
    #No Correlation
    checkresiduals(sarima_4$residuals)
    tsdiag(sarima_4)

  # Sarima_5 
  summary(sarima_5)
  confint(sarima_5)
  BIC(sarima_5)
  AIC(sarima_5)
  
    # Assumption validation
    
    # Normality
    ad.test(sarima_5$residuals)
    shapiro.test(sarima_5$residuals)
    qqnorm(sarima_5$residuals)
    qqline(sarima_5$residuals)
    # Homocedasticity
    x_5 <- 1:length(sarima_5$residuals)
    y_5 <- as.numeric(sarima_5$residuals)
    bptest(y_5~x_5)
    plot(sarima_5$residuals)
    # Zero mean
    t.test(sarima_5$residuals)
    #No correlation
    checkresiduals(sarima_5$residuals)
    tsdiag(sarima_5)

  # BC_sarima1 
  summary(BC_sarima1)
  confint(BC_sarima1)
  BIC(BC_sarima1)
  AIC(BC_sarima1)
  
    # Assumption validation
    
    # Normality
    ad.test(BC_sarima1$residuals)
    shapiro.test(BC_sarima1$residuals)
    qqnorm(BC_sarima1$residuals)
    qqline(BC_sarima1$residuals)
    # Homocedasticity
    BCx_1 <- 1:length(BC_sarima1$residuals)
    BCy_1 <- as.numeric(BC_sarima1$residuals)
    bptest(BCy_1~BCx_1)
    plot(BC_sarima1$residuals)
    # Zero mean
    t.test(BC_sarima1$residuals)
    #No correlation
    checkresiduals(BC_sarima1$residuals)
    tsdiag(BC_sarima1)

  # BC_sarima2 
  summary(BC_sarima2)
  confint(BC_sarima2)
  BIC(BC_sarima2)
  AIC(BC_sarima2)
  
    # Assumption validation
    
    # Normality
    ad.test(BC_sarima2$residuals)
    shapiro.test(BC_sarima2$residuals)
    qqnorm(BC_sarima2$residuals)
    qqline(BC_sarima2$residuals)
    # Homocedasticity
    BCx_2 <- 1:length(BC_sarima2$residuals)
    BCy_2 <- as.numeric(BC_sarima2$residuals)
    bptest(BCy_2~BCx_2)
    plot(BC_sarima2$residuals)
    # Zero mean
    t.test(BC_sarima2$residuals)
    #No correlation
    checkresiduals(BC_sarima2$residuals)
    tsdiag(BC_sarima2)
    

# Model SARIMA (0,1,2)(0,1,1)[12] and forecast ----
# Fit
ajuste_1 <- cafe-residuals(sarima_5)
    
ts.plot(cafe, main = 'SARIMA(0,1,2)(0,1,1)[12]', xlab = 'Year', ylab = 'Price', 
        col = '#08829E', lwd =1)
points(ajuste_1, type = 'l', col = '#A1CAD3', lwd =1)
# Prediction
plot(forecast(sarima_5,h=12), include = 60, col = c('#4D8C99','#357887'), 
     shadecols= c('#E7EAEA','#BEDADF'), ylab='Price', xlab = 'Year',
     main = "Forecast SARIMA(0,1,2)(0,1,1)[12]", lwd =2)
# Bayesian Model ----
# JAGS
n <- length(cafe)
data <- list(
  y = as.integer(cafe),
  n = length(cafe)
)
# Initial values
inits <- function(){ list(
  
  theta_1 = rnorm(1,0,0.1),
  theta_2 = rnorm(1,0,0.1),
  Theta_1 = rnorm(1,0,0.1),
  tau = runif(1,0,0.001),
  tau_z = runif(1,0,0.001))
  
}

params <- c('theta_1', 'theta_2', 'Theta_1', 'sigma2','sigma2z','y.pred')

modelo_SARIMA <-  "model{

  for(i in 1:n){
    z[i] ~ dnorm(0,1/sigma2z)
  }
  
  for(i in 15:n){
    y[i]~ dnorm(f[i],1/sigma2)
    f[i] <- alpha + y[i-1]+y[i-12]-y[i-13]+z[i]-theta_1*z[i-1]-theta_2*z[i-2]-
            Theta_1*(z[i-12]-theta_1*z[i-13]-theta_2*z[i-14])
  }
  
  # Predictions 
  for(i in 15:(12+15)){
               y.pred[i] ~dnorm(f.pred[i],1/sigma2)
               f.pred[i] <- alpha + y.pred[i-1]+y.pred[i-12]-y.pred[i-13] 
               + z.pred[i]- theta_1*z.pred[i-1]-theta_2*z.pred[i-2] - 
               Theta_1*(z.pred[i-12]-theta_1*z.pred[i-13]-theta_2*z.pred[i-14])
  }
  
  for(i in 15:(12+15)){
              z.pred[i] ~ dnorm(0,1/sigma2z)
  }
  
  for(i in 1:14){
              z.pred[i] <- z[n-14+i] 
              y.pred[i] <- y[n-14+i]
  }
  
  
  for(i in 1:14){
    y[i]~ dnorm(f[i],1/sigma2)
    f[i] <- alpha -theta_1*z[i]-theta_2*z[i]-
            Theta_1*(z[i]-theta_1*z[i]-theta_2*z[i])
  }

# Informative prior
  alpha ~ dnorm(0,0.01)
  theta_1 ~ dnorm(0,1)
  theta_2 ~ dnorm(0,1)
  Theta_1 ~ dnorm(0,1)
  tau ~ dunif(0.00001,0.001)
  tau_z ~ dunif(0.00001,0.001)
  sigma2 <- 1/tau
  sigma2z <- 1/tau_z
  
}"

set.seed(100)
fit <- jags.model(textConnection(modelo_SARIMA), data, inits = inits, n.chains=3)

update(fit,5000)
sample.sarima <- coda.samples(fit,params,n.iter = 5000, thin = 1)
#plot(sample.sarima)

# Values 
summary(sample.sarima)

#gelman.diag(sample.sarima, confidence = 0.95, transform = F)
#gelman.plot(sample.sarima)
# JAGS Forecast ----
 
aux<-summary(sample.sarima)
pred = aux$statistics[19:30,1]
sd.error= aux$statistics[19:30,3]
cuantil = aux$quantiles[19:30,c(1,5)]

# Confidence intervals 2 SE (Naive SD) : 
lim_sup = pred[1:12] + 2*sd.error[1:12]
lim_inf = pred[1:12] - 2*sd.error[1:12]

x <- seq.Date(as.Date('2016-01-01'),as.Date('2020-11-01'), by = 'month')
x_1 <- seq.Date(as.Date('2020-12-01'),as.Date('2021-11-01'), by = 'month')
upper <- cuantil[,2]
lower <- cuantil[,1]
e1<-ggplot()+
  geom_line(aes(x = x, y = cafe[313:371], color = "Data"), size = 0.8)+
  geom_ribbon(aes(x = x_1, ymax =upper, ymin = lower, fill = '.'), alpha = 0.4)+
  geom_line(aes(x = x_1, y = pred, color = "Forecast"), size = 0.8)+
  scale_fill_manual(values =c('#c6dbef','#c6dbef','#c6dbef','#c6dbef',
                              '#c6dbef','#c6dbef','#c6dbef','#c6dbef'), name = 'Bands')+
  scale_color_manual(values = c('#51B1D4', '#00779B'), name = 'Series')+
  labs(y = 'Price', x = 'Year')+
  ggtitle('Coffee price forecast') +  
  theme_test() +
  theme(plot.title = element_text(size= 12, hjust = 0.5),
        axis.title.x = element_text(size = 10, color = 'grey20'),
        axis.title.y = element_text(size = 10, color = 'grey20'))  +
  theme(legend.position = "none")

e2<-ggplot()+
  geom_line(aes(x = x, y = cafe[313:371], color = "Data"), size = 0.8)+
  geom_ribbon(aes(x = x_1, ymax =lim_sup, ymin = lim_inf, fill = '.'), alpha = 0.6)+
  geom_line(aes(x = x_1, y = pred, color = "Forecast"), size = 0.8)+
  scale_fill_manual(values =c('#c6dbef','#c6dbef','#c6dbef','#c6dbef',
                              '#c6dbef','#c6dbef','#c6dbef','#c6dbef'), name = 'Bands')+
  scale_color_manual(values = c('#51B1D4', '#00779B'), name = 'Series')+
  labs(y = 'Price', x = 'Year')+
  ggtitle('Coffee price prediction') +  
  theme_test() +
  theme(plot.title = element_text(size= 11, hjust = 0.5),
        axis.title.x = element_text(size = 10, color = 'grey20'),
        axis.title.y = element_text(size = 10, color = 'grey20'))  +
  theme(legend.position = "none")

grid.arrange(e1, e2, ncol=2)

# Linear Dynamic Models ----
# JAGS 

Nnew = 12 ### Number of predictions
n <- length(cafe)
data <- list(
  Y = as.integer(cafe),
  N = length(cafe),
  Nnew = Nnew
)

inits <- function(){ list(
  
  theta_1 = rnorm(1,0,0.1),
  theta_2 = rnorm(1,0,0.1),
  Theta_1 = rnorm(1,0,0.1),
  tau = runif(1,0,0.001),
  tau_z = runif(1,0,0.001))
}

params <- c("sd.q", "sd.r", "mu","phi1","phi2","Ynew")

modelo_DINAMICO <-  "model {  
   # priors on parameters
   mu ~ dnorm(0, 0.01); 
   tau.pro ~ dunif(0.00001,0.001); 
   sd.q <- 1/sqrt(tau.pro);
   tau.obs ~ dunif(0.00001,0.001);
   sd.r <- 1/sqrt(tau.obs); 
   phi1 ~ dnorm(0,1);
   phi2 ~ dnorm(0,1);
### Modelo   
   X[1] <- mu;
   X[2] <- mu;
   predY[1] <- X[1];
   predY[2] <- X[2];
   Y[1] ~ dnorm(X[1], tau.obs);
   Y[2] ~ dnorm(X[2], tau.obs);

   for(i in 3:N) {
      predX[i] <- phi1*X[i-1]+phi2*X[i-2]; 
      X[i] ~ dnorm(predX[i],tau.pro); # Process variation
      predY[i] <- X[i];
      Y[i] ~ dnorm(X[i], tau.obs); # Observation variation
   }
### Prediccion
   predXnew[1] <- phi1*X[N]+phi2*X[N-1]; 
   predXnew[2] <- phi1*Xnew[2-1]+phi2*X[N]; 
   for(i in 3:Nnew) {
      predXnew[i] <- phi1*Xnew[i-1]+phi2*Xnew[i-2]; 
   }
   for(i in 1:Nnew) {
      Xnew[i] ~ dnorm(predXnew[i],tau.pro); # Process variation
      predYnew[i] <- Xnew[i];
      Ynew[i] ~ dnorm(Xnew[i], tau.obs); # Observation variation
   }
} "

set.seed(100)

fit <- jags.model(textConnection(modelo_DINAMICO), data, inits = inits, n.chains=3)

update(fit,10000)
sample.dinamico <- coda.samples(fit,params,n.iter = 10000, thin = 1)
#plot(sample.dinamico)

summary(sample.dinamico)
#gelman.diag(sample.dinamico, confidence = 0.95, transform = F)
#gelman.plot(sample.dinamico)


# LDM JAGS Forecast ----
N = length(cafe)
cafe[(N-Nnew):N]
(cafepred = summary(sample.dinamico)$statistics[1:Nnew,1])
(cuantil = summary(sample.dinamico)$quantiles[1:Nnew,c(1,5)])

x <- seq.Date(as.Date('2016-01-01'),as.Date('2020-11-01'), by = 'month')
x_1 <- seq.Date(as.Date('2020-12-01'),as.Date('2021-11-01'), by = 'month')
upper <- cuantil[,2]
lower <- cuantil[,1]
ggplot()+
  geom_line(aes(x = x, y = cafe[313:371], color = "Data"), size = 0.8)+
  geom_ribbon(aes(x = x_1, ymax = upper, ymin = lower, fill = '95%'), alpha = 0.5)+
  geom_line(aes(x = x_1, y = cafepred, color = "Forecast"), size = 0.8)+
  scale_fill_manual(values =c('#c6dbef','#9ecae1'), name = 'Bands')+
  scale_color_manual(values = c('#51B1D4', '#3896B8'), name = 'Series')+
  labs(y = 'Price', x = 'Year')+
  ggtitle('Coffee price forecast LDM (JAGS)') +  theme_test() +
  theme(plot.title = element_text(size= 12, hjust = 0.5))+
  theme(legend.position = "none")


# Fit with ldm package ----
library(dlm)
# Parameterization of the two variances as a function of the maximum likelihood

build <- function(parm){
  dlmModPoly(order = 1, dV = exp(parm[1]), dW = exp(parm[2]))
}
ajuste <- dlmMLE(cafe, rep(0,2), lower=c(1e-6,0), build, hessian = TRUE)
ajuste$convergence
unlist(build(ajuste$par))[c('V','W')]
avar <- solve(ajuste$hessian)
sqrt(diag(avar))

dlmcafe <- build(ajuste$par)
cafeFit <- dlmFilter(cafe, dlmcafe)
# Fit 
plot(cafe, col = '#00779B', lwd = 1.5, main = 'Ajuste del modelo', xlab = 'Año', 
     ylab = 'Precio') 
lines(dropFirst(dlmBSample(cafeFit)), col = "#51B1D4", type = "o")

# Forecast with ldm package ----
set.seed(120)
n <- 12
num <- 100
y1 <- dlmForecast(cafeFit, nAhead=n,sampleNew=num)
summary(y1)

# Example
x <- c()
y<-c()
for(j in 1:12){
  for(i in 1:100){
    x[i]<-y1$newObs[[i]][j,1]
  }
  
  y[j]<-median(x)
}
x <- seq.Date(as.Date('2016-01-01'),as.Date('2020-11-01'), by = 'month')
x_1 <- seq.Date(as.Date('2020-12-01'),as.Date('2021-11-01'), by = 'month')

p1 = ggplot()+
  geom_line(aes(x =x, y = cafe[313:371], color = "Data"), size = 0.6)+
  geom_line(aes(x = x_1, y = y1$newObs[[3]], color = "Estimated"), size = 0.6)+
  scale_color_manual(values = c('#51B1D4', '#00779B', '#51B1D4'), name = 'Series')+
  labs(y = 'Price', x = 'Year')+
  ggtitle(' ') +  theme_test() +
  theme(plot.title = element_text(size= 12, hjust = 0.5),
        axis.title.x = element_text(size = 10, color = 'grey20'),
        axis.title.y = element_text(size = 10, color = 'grey20'))+
  theme(legend.position = "none")

p2 = ggplot()+
  geom_line(aes(x =x, y = cafe[313:371], color = "Data"), size = 0.6)+
  geom_line(aes(x = x_1, y = y1$newObs[[4]], color = "Estimated"), size = 0.6)+
  scale_color_manual(values = c('#51B1D4', '#00779B', '#51B1D4'), name = 'Series')+
  labs(y = 'Price', x = 'Year')+
  ggtitle('Forecast example') +  theme_test() +
  theme(plot.title = element_text(size= 12, hjust = 0.5),
        axis.title.x = element_text(size = 10, color = 'grey20'),
        axis.title.y = element_text(size = 10, color = 'grey20'))+
  theme(legend.position = "none")

p3 = ggplot()+
  geom_line(aes(x =x, y = cafe[313:371], color = "Data"), size = 0.6)+
  geom_line(aes(x = x_1, y = y1$newObs[[5]], color = "Estimated"), size = 0.6)+
  scale_color_manual(values = c('#51B1D4', '#00779B', '#51B1D4'), name = 'Series')+
  labs(y = 'Price', x = 'Year')+
  ggtitle(' ') +  theme_test() +
  theme(plot.title = element_text(size= 12, hjust = 0.5),
        axis.title.x = element_text(size = 10, color = 'grey20'),
        axis.title.y = element_text(size = 10, color = 'grey20'))+
  theme(legend.position = "none")

grid.arrange(p1, p2, p3, ncol=3)

# FORECAST
x <- c()
y<-c()
for(j in 1:12){
  for(i in 1:100){
    x[i]<-y1$newObs[[i]][j,1]
  }
  
  y[j]<-median(x)
}
x <- seq.Date(as.Date('2016-01-01'),as.Date('2020-11-01'), by = 'month')
x_1 <- seq.Date(as.Date('2020-12-01'),as.Date('2021-11-01'), by = 'month')
upper <- cuantil[,2]
lower <- cuantil[,1]
ggplot()+
  geom_line(aes(x = x, y = cafe[313:371], color = "Data"), size = 0.8)+
  geom_ribbon(aes(x = x_1, ymax = y1$newObs[[1]], ymin = y1$newObs[[2]], fill = '.'), alpha = 0.4)+
  geom_ribbon(aes(x = x_1, ymax = y1$newObs[[3]], ymin = y1$newObs[[4]], fill = '.'), alpha = 0.4)+
  geom_ribbon(aes(x = x_1, ymax = y1$newObs[[5]], ymin = y1$newObs[[6]], fill = '.'), alpha = 0.4)+
  geom_ribbon(aes(x = x_1, ymax = y1$newObs[[7]], ymin = y1$newObs[[8]], fill = '.'), alpha = 0.4)+
  geom_ribbon(aes(x = x_1, ymax = y1$newObs[[9]], ymin = y1$newObs[[10]], fill = '.'), alpha = 0.4)+
  geom_ribbon(aes(x = x_1, ymax = y1$newObs[[11]], ymin = y1$newObs[[12]], fill  = '.'), alpha = 0.4)+
  geom_ribbon(aes(x = x_1, ymax = y1$newObs[[13]], ymin = y1$newObs[[14]], fill  = '.'), alpha = 0.4)+
  geom_line(aes(x = x_1, y = y, color = "Forecast"), size = 0.8)+
  scale_fill_manual(values =c('#c6dbef','#c6dbef','#c6dbef','#c6dbef','#c6dbef','#c6dbef','#c6dbef','#c6dbef'), name = 'Bandas')+
  scale_color_manual(values = c('#51B1D4', '#00779B'), name = 'Series')+
  labs(y = 'Price', x = 'Year')+
  ggtitle('Coffee price forecats ldm (meadian)') +  theme_test() +
  theme(plot.title = element_text(size= 12, hjust = 0.5),
        axis.title.x = element_text(size = 10, color = 'grey20'),
        axis.title.y = element_text(size = 10, color = 'grey20'))+
  theme(legend.position = "none")


# Training and validation sets ----
aux<- cafe[1:359]
aux2<- cafe[360:371]

build <- function(parm){
  dlmModPoly(order = 1, dV = exp(parm[1]), dW = exp(parm[2]))
}
ajuste <- dlmMLE(aux, rep(0,2), lower=c(1e-6,0), build, hessian = TRUE)
ajuste$convergence
unlist(build(ajuste$par))[c('V','W')]
avar <- solve(ajuste$hessian)
sqrt(diag(avar))

dlmaux <- build(ajuste$par)
auxFit <- dlmFilter(aux, dlmaux)
set.seed(120)
n <- 12
num <- 100
y1 <- dlmForecast(auxFit, nAhead=n,sampleNew=num)
summary(y1)

# Plot
x<-c()
y<-c()
for(j in 1:n){
  for(i in 1:num){
    x[i]<-y1$newObs[[i]][j,1]
  }
  
  y[j]<-median(x)}

x_3 <- seq.Date(as.Date('2016-01-01'),as.Date('2019-11-01'), by = 'month')
x_2 <-seq.Date(as.Date('2019-12-01'),as.Date('2020-11-01'), by = 'month')

ggplot()+
  geom_line(aes(x = x_3, y = cafe[313:359], color = "Series"), size = 0.8)+
  geom_ribbon(aes(x = x_2, ymax = y1$newObs[[1]], ymin = y1$newObs[[2]], 
                  fill = '.'), alpha = 0.4)+
  geom_ribbon(aes(x = x_2, ymax = y1$newObs[[3]], ymin = y1$newObs[[4]], 
                  fill = '.'), alpha = 0.4)+
  geom_ribbon(aes(x = x_2, ymax = y1$newObs[[5]], ymin = y1$newObs[[6]], 
                  fill = '.'), alpha = 0.4)+
  geom_ribbon(aes(x = x_2, ymax = y1$newObs[[7]], ymin = y1$newObs[[8]], 
                  fill = '.'), alpha = 0.4)+
  geom_ribbon(aes(x = x_2, ymax = y1$newObs[[9]], ymin = y1$newObs[[10]], 
                  fill = '.'), alpha = 0.4)+
  geom_ribbon(aes(x = x_2, ymax = y1$newObs[[11]], ymin = y1$newObs[[12]], 
                  fill  = '.'), alpha = 0.4)+
  geom_ribbon(aes(x = x_2, ymax = y1$newObs[[13]], ymin = y1$newObs[[14]], 
                  fill  = '.'), alpha = 0.4)+
  geom_line(aes(x = x_2, y = aux2, color = "Estimated"), size = 0.8)+
  geom_line(aes(x = x_2, y = y, color = "Series"), size = 0.8)+
  scale_color_manual(values = c('#51B1D4', '#00779B', '#51B1D4'), name = 'Series')+
  scale_fill_manual(values =c('#c6dbef','#c6dbef','#c6dbef','#c6dbef','#c6dbef',
                              '#c6dbef','#c6dbef','#c6dbef'), name = 'Bands')+
  labs(y = 'Price', x = 'Year')+
  ggtitle('Comparison of the model with the original series') +  theme_test() +
  theme(plot.title = element_text(size= 12, hjust = 0.5),
        axis.title.x = element_text(size = 10, color = 'grey20'),
        axis.title.y = element_text(size = 10, color = 'grey20'))+
  theme(legend.position = "none")
