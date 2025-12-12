

# install.packages("xts")
# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("lubridate")
# install.packages("dynlm")
# install.packages("dynsim")
# install.packages("forecast")
# install.packages("rlang")
# install.packages("devtools")
# install.packages("pkgload")
# install.packages("vctrs")
# install.packages("tseries")
# install.packages("astsa")
# install.packages("plotly")
# install.packages("hrbrthemes")
# install.packages("tibbletime")
# install.packages("anomalize")
# install.packages("timetk")
# install.packages("fpp2")
# install.packages(adress, repos = NULL, type = "source")

library(readr)
library(readxl)
library(timetk)
library(tibbletime)
library(anomalize)
library(xts)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(dynlm)
library(forecast)
library(dynsim)
library(devtools)
library(tsdl)
library(tseries)
library(astsa)
library(plotly)
library(fpp2)
library(hrbrthemes)
library(gridExtra)
library(tseries)
library(uroot)
library(car)
library(bestNormalize)
library(stats)
library(fpp)
library(zoom)


### EDA ####

UNRATE <- read_excel("UNRATE.xls")

head(UNRATE)
tail(UNRATE)
UNRATE <- as.data.frame(UNRATE)
class(UNRATE)



unrate2 <- ts(UNRATE[,2], start=c(1948,1), frequency = 12)


head(unrate2,120)


sum(is.na(unrate2))

summary(unrate2)
sd(unrate2)

autoplot(unrate2, main = "Unemployement Rate", col ="red", ylab = "Rate %")

str(UNRATE)

UNRATE$observation_date <- paste(UNRATE$observation_date, "01", sep = "-")
UNRATE$observation_date <- as.Date(UNRATE$observation_date, format = "%Y-%m-%d")

p <- UNRATE %>%
  ggplot(aes(x = observation_date, y = UNRATE)) +
  geom_area(fill = "cyan4", alpha=0.5) +
  geom_line(color="cyan4") +
  ylab("Unemployement Rate")+ labs(title= "Time Series Plot of Monthly Unemployement Rate") +
  theme_ipsum()
p <- ggplotly(p)
p


cycle(unrate2)

boxplot(unrate2~cycle(unrate2), xlab = "Month", ylab = "Unemployement Rate")


unrate3 <- UNRATE %>% dplyr::mutate(year = lubridate::year(observation_date), month = lubridate::month(observation_date))
unrate3
str(unrate3)

unrate3$month <- as.factor(unrate3$month)
unrate3$year <- as.factor(unrate3$year)

bp <- ggplot(unrate3, aes(x= month, y= UNRATE, fill=month))+
  geom_boxplot()+
  labs(title = "Boxplot Across Months", x= "Month", y= "Unemployement Rate")
bp


bp2 <- ggplot(unrate3, aes(x = year, y=UNRATE, fill=year ))+
  geom_boxplot()+
  labs(title = "Boxplot Across Years", x="Year", y="Unemployement Rate")
bp2


a <- ggplot(unrate3, aes(as.numeric(as.character(month)), UNRATE))+
  geom_line(aes(color=year))+
  xlab("Month")

a


par(mfrow = c(2,2))


unrate4 <- unrate3[1:288,]
unrate5 <- unrate3[289:576,]
unrate6 <- unrate3[577:898,]




par(mfrow = c(2,2))


b <- ggplot(unrate4, aes(as.numeric(as.character(month)), UNRATE))+
  geom_line(aes(color=year))+
  xlab("Month")+
  ylab("UNRATE")

c <- ggplot(unrate5, aes(as.numeric(as.character(month)), UNRATE))+
  geom_line(aes(color=year))+
  xlab("Month")+
  ylab("UNRATE")

d <- ggplot(unrate6, aes(as.numeric(as.character(month)), UNRATE))+
  geom_line(aes(color=year))+
  xlab("Month")+
  ylab("UNRATE")

grid.arrange(a,b,c,d)


### Anomaly detection ####

UNRATE$observation_date <- as.Date(UNRATE$observation_date, format ="%Y-%m-%d")
UNRATE$rate <- UNRATE$UNRATE
UNRATE$month <- UNRATE$observation_date



unrate <- UNRATE %>% select(month, rate)

unrate_tibble <- as_tibble(unrate)
class(unrate_tibble)
head(unrate_tibble)

f <- unrate_tibble %>%
  time_decompose(rate, method = "stl", frequency = "auto", trend = "auto") %>%
  anomalize(remainder, method = "gesd", alpha = 0.05, max_anoms = 0.2) %>%
  plot_anomaly_decomposition()



g <- unrate_tibble %>%
  time_decompose(rate) %>%
  anomalize(remainder) %>%
  time_recompose() %>%
  plot_anomalies(time_recomposed = T, ncol = 3, alpha_dots = 0.5)

grid.arrange(f,g)




### Anomaly Extraction ####


anomalous <- unrate_tibble %>%
  time_decompose(rate) %>%
  anomalize(remainder) %>%
  time_recompose() %>%
  filter(anomaly == "Yes")

print(anomalous, n = 100)


unrate$rate <- tsclean(unrate$rate)

unrate_tibble_clean <- as_tibble(unrate)
class(unrate_tibble_clean)


h <- unrate_tibble_clean %>%
  time_decompose(rate, method = "stl", frequency = "auto", trend = "auto") %>%
  anomalize(remainder, method = "gesd", alpha = 0.05, max_anoms = 0.2) %>%
  plot_anomaly_decomposition()

h

j <- unrate_tibble_clean %>%
  time_decompose(rate) %>%
  anomalize(remainder) %>%
  time_recompose() %>%
  plot_anomalies(time_recomposed = T, ncol = 3, alpha_dots = 0.5)

j


### Data splitting ####
autoplot(unrate2)


qqPlot(unrate2, ylab = "Rate")

lambda <- 0.15
unrate_transformed <- BoxCox(unrate2, lambda)


qqPlot(unrate_transformed)

qqPlot(unrate_transformed$x.t, ylab = "Transformed Rate")

trainindex <- 1:(length(unrate_transformed)-100)
train <- unrate_transformed[trainindex]
test <- unrate_transformed[-trainindex]

train
train_ts <- ts(train, start=c(1948,1), frequency = 12)


plot(train, type = "l")


test
test_ts <- ts(test, start=c(2014,6), frequency = 12)

autoplot(test_ts)


fit1 <- ses(train_ts, alpha=0.2, initial="simple", h=100)
fit2 <- ses(train_ts, alpha=0.6, initial="simple", h=100)
fit3 <- ses(train_ts, h=100)

autoplot(train_ts) +
  autolayer(fit1,PI=F,series="alpha=0.2")+
  autolayer(fit2, PI=F,series="alpha=0.6") +
  autolayer(fit3, PI=F,series="alpha=0.99")+
  autolayer(test_ts,series="actual",color="black")+
  ggtitle("Forecast from SES'S")+
  theme_minimal()



summary(fit1)
summary(fit2)
summary(fit3)


plot(fit1)


# According to error measurements Fit3 is better fit.






### Nonstationay Test ####


p1 <- ggAcf(train)
p2 <- ggPacf(train)

grid.arrange(p1,p2, nrow = 1)

acf2(as.vector(train), main = "ACF and PACF graphs of Train Dataset")

train_ts <- ts(train, start=c(1948,1), frequency = 12)


kpss.test(train_ts, null = c("Level"))
kpss.test(train_ts, null = c("Trend"))

PP.test(train_ts)

HEGY.test(train_ts, itsd = c(1,1,0), regvar = 0, 
          selectlags = list(mode = "signf", Pmax = NULL))

# Regular unit root exists
train_ts_d <- diff(train_ts)
autoplot(train_ts_d, col="red", main = "Differenced Unemployment Rate Graph", ylab = "Rate")
kpss.test(train_ts_d, null = c("Level"))
kpss.test(train_ts_d, null = c("Trend"))
ADF.test(train_ts_d, itsd = c(1,1,0))
HEGY.test(train_ts_d, itsd = c(1,1,0), regvar = 0, 
          selectlags = list(mode ="signf", Pmax = NULL))

### SARIMA Modelling ####




acf2(as.vector(train_ts_d), max.lag = 60, main = "P/ACF Plots of Transformed Dataset")

sarima_model <- Arima(train_ts, order = c(3,1,1), 
                      seasonal = list(order = c(2,0,2), period = 12))


sarima_model

sarima_model_auto <- auto.arima(train_ts)


sarima_model_auto


sarima.for(train_ts, n.ahead = 60, 0,1,0,1,0,1,12)
plot(unrate2)


fore_sarima <- forecast(sarima_model, h =12)
df_sarima <- as.data.frame(fore_sarima)



sarima_model_auto %>% forecast() %>% autoplot( ylab= "Rate")



fr1 <- ets(train_ts, model = "MMM")
fore1 <- forecast(fr1, h = 12)

fr2 <- ets(train_ts, model = "ZZZ")
fore2 <- forecast(fr2, h =12)

fr3 <- tbats(train_ts)
fore3 <- forecast(fr3, h = 12)

fr4 <- nnetar(train_ts)
fore4 <- forecast(fr4, h = 12, PI = T)

autoplot(fore1)+
  autolayer(fitted(fore1))+
  autolayer(test_ts)

autoplot(fore2)+
  autolayer(fitted(fore2))+
  autolayer(test_ts)




autoplot(fore3)+
  autolayer(fitted(fore3))+ 
  autolayer(test_ts)

autoplot(fore4)+
  autolayer(test_ts)


accuracy(fore1, test_ts)
accuracy(fore2, test_ts)
accuracy(fore3, test_ts)
accuracy(fore4, test_ts)


Fit <- nnetar(unrate2)

forecast <- forecast(Fit, h = 60, PI =T)

autoplot(forecast, ylab = "Rate")


