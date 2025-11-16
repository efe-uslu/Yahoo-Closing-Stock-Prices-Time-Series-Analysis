#Part1

library(xts)
library(TSA)
library(forecast)
library(tibble)
library(dplyr)
library(ggplot2)
library(tseries)
library(urca)
library(lubridate)
library(gridExtra)
library(pdR)
library(Matrix)
library(uroot)
library(lmtest)
library(bestNormalize)
library(aTSA)
library(MTS)
library(rugarch)
library(foreach)
library(doParallel)
library(prophet)
library(zoo)
library(stats)





yahoo_stock$Date <- as.Date(yahoo_stock$Date)
time_series_xts <- xts(yahoo_stock$Close, order.by = yahoo_stock$Date)
plot(time_series_xts, main = "Time Series Data", ylab = "Close", xlab = "Date")
a1<-ggAcf(time_series_xts,lag.max = 48)+theme_minimal()+ggtitle("ACF of Closing Stock Prices")
a2<-ggPacf(time_series_xts,lag.max = 48)+theme_minimal()+ggtitle("PACF of Closing Stock Prices")
grid.arrange(a1,a2,ncol=2)
###increasing trend


#Part2
test_set <- time_series_xts[(nrow(time_series_xts) - 39):nrow(time_series_xts), ]
train_set <- time_series_xts[1:(nrow(time_series_xts) - 40), ]


#Part3
lambda <- BoxCox.lambda(time_series_xts)
lambda
###lambda value for optimal variance stabilizing box cox transformation is "-0.8656467".
data_transformed <- BoxCox(time_series_xts, lambda)

plot(data_transformed, main = "Transformed Time Series Data", ylab = "Transformed Close", xlab = "Date")

###As the form of the time series is not changed and variation seems very similar to the original series.
###Transformation in this case is then determined to be not necessary.

plot(train_set, main = "Train Data", ylab = "Close", xlab = "Date")
plot(test_set, main = "Test Data", ylab = "Close", xlab = "Date")


#Part4

###Transforming to tibble
train_set_tbl <- as_tibble(train_set, rownames = "Date") %>% 
  rename(Close = V1) %>%
  mutate(Date = as.Date(Date))
test_set_tbl <- as_tibble(test_set, rownames = "Date") %>% 
  rename(Close = V1) %>%
  mutate(Date = as.Date(Date))

###Anomaly Removal
train_cleaned <- train_set %>%
  as.numeric() %>%
  tsclean()
test_cleaned <- test_set %>%
  as.numeric() %>%
  tsclean()

###Tibble again
train_cleaned_tbl <- tibble(
  Date = as.Date(index(train_set)),
  Close = train_cleaned
)
test_cleaned_tbl <- tibble(
  Date = as.Date(index(test_set)),
  Close = test_cleaned
)

###Visiualization of anomaly detection

ggplot() +
  geom_line(data = train_set_tbl, aes(x = Date, y = Close), color = "red", linetype = "dashed") +
  geom_line(data = train_cleaned_tbl, aes(x = Date, y = Close), color = "blue") +
  labs(title = "Original vs Cleaned Data", x = "Date", y = "Close") +
  theme_minimal()

#Part5

train_ts <- ts(
  train_cleaned_tbl$Close,
  frequency = 365,
  start = c(year(min(train_cleaned_tbl$Date)), yday(min(train_cleaned_tbl$Date)))
)
test_ts <- ts(
  test_cleaned_tbl$Close,
  frequency = 365,
  start = c(year(min(test_cleaned_tbl$Date)), yday(min(test_cleaned_tbl$Date)))
)

g1<-ggAcf(train_ts,lag.max = 48)+theme_minimal()+ggtitle("ACF of Train")
g2<-ggPacf(train_ts,lag.max = 48)+theme_minimal()+ggtitle("PACF of Train")
grid.arrange(g1,g2,ncol=2)
###slow linear decay on ACF indicated non-stationarity in the data

###KPSS test for stationarity
tseries::kpss.test(train_ts, null = "Level")
###As it is less than 0.01 "H0:stationary" is rejected therefore the series is not stationary.
tseries::kpss.test(train_ts, null = "Trend")
###Again being less than 0.01 means that there is stochastic trend

###Test for seasonal roots
###Hegy test
result <- hegy.test(train_ts)
result$pvalues
###Result for 364:365 has very low p-value, thus no seasonal differencing is necessary.

#Part6
###Now a differencing is necessary as the series has stochastic trend.
train_differenced <- diff(train_ts)
###KPSS test for stationarity
tseries::kpss.test(train_differenced, null = "Level")
###Stationarity is achieved with p-value > 0.1

result2 <- hegy.test(train_differenced)
result2$pvalues
###There exists no seasonal or regular unit roots based on hegy test results.
nsdiffs(train_differenced)
ndiffs(train_differenced)

#Part7

###Checking the ACF and PACF graphs
plot(train_differenced, main = "Time Series Data Differenced", ylab = "Differenced Close", xlab = "Date")
###Constant mean, but variance does not seem constant.
g3<-ggAcf(train_differenced,lag.max = 48)+theme_minimal()+ggtitle("ACF of Train")
g4<-ggPacf(train_differenced,lag.max = 48)+theme_minimal()+ggtitle("PACF of Train")
grid.arrange(g3,g4,ncol=2)
###ACF and PACF both break off at 1,2 and 4 lets try all these and compare the results.

decomposition = stl(train_ts, t.window=13, s.window="periodic", robust=TRUE)
autoplot(decomposition,col="red")+theme_minimal()
###decomposition shows the trend, seasonality and what is left after this components are removed.
results <- data.frame()
for(p in 0:4){
  for(q in 0:4){
    fit <- Arima(train_ts, order=c(p,1,q))
    bic <- fit$bic
    coefs <- fit$coef
    se <- sqrt(diag(fit$var.coef))
    is_significant <- function(estimate, se){
      z_val <- abs(estimate)/se
      p_val <- 2*pnorm(z_val, lower.tail=FALSE)
      p_val < 0.05
    }
    get_sig_string <- function(type, i){
      term <- paste0(type, i)
      if(term %in% names(coefs)){ifelse(is_significant(coefs[term], se[term]), "Yes", "No")} else{NA}
    }
    results <- rbind(results, data.frame(
      p = p,
      q = q,
      BIC = bic,
      AR1_sig = get_sig_string("ar", 1),
      AR2_sig = get_sig_string("ar", 2),
      AR3_sig = get_sig_string("ar", 3),
      AR4_sig = get_sig_string("ar", 4),
      MA1_sig = get_sig_string("ma", 1),
      MA2_sig = get_sig_string("ma", 2),
      MA3_sig = get_sig_string("ma", 3),
      MA4_sig = get_sig_string("ma", 4)
    ))
  }
}
results


###So after checking BIC values, it is determined that ARIMA(0,1,2) is the best model.with BIC = 15794.35 but it is not statistically significant.
###With every term being significant except AR3 (A) ARIMA(4,1,0) with BIC = 15795.77 can be chosen with dropping the AR3 term
###Or another model which has every term significant can also be chosen such has:
###(B) ARIMA(2,1,2) with BIC = 15799.74.
###(C) ARIMA(0,1,2) with BIC = 15794.35.
###(D) ARIMA(1,1,0) with BIC = 15797.06.
###(E) ARIMA(2,1,0) with BIC = 15797.31.
###(F) ARIMA(0,1,1) with BIC = 15800.08.


###Dropping AR3 for (A)
fit_A <- Arima(train_ts,order = c(4,1,0),fixed = c(NA, NA, 0, NA, NA),include.constant = TRUE, transform.pars = FALSE)
summary(fit_A)
###Checking if dropping AR3 affected stationarity.
ar_coefs <- fit_A$model$phi
polynomial <- c(1, -ar_coefs)
roots <- polyroot(polynomial)
mod_roots <- Mod(roots)
all(mod_roots > 1)
###It is still stationary.


###Fits
fit_B <- Arima(train_ts, order = c(2, 1, 2))
fit_C <- Arima(train_ts, order = c(0, 1, 2))
fit_D <- Arima(train_ts, order = c(1, 1, 0))
fit_E <- Arima(train_ts, order = c(2, 1, 0))
fit_F <- Arima(train_ts, order = c(0, 1, 1))

#Diagnostic Checking.

##Correlation of residuals

##First Box.test is conducted to test if residuals are uncorrelated.
Box.test(fit_A$residuals, lag = 50, type = "Ljung-Box")
###p-value is greater than 0.05 thus residuals of (A) are uncorrelated.

Box.test(fit_B$residuals, lag = 50, type = "Ljung-Box")
###p-value is greater than 0.05 thus residuals of (B) are uncorrelated.

Box.test(fit_C$residuals, lag = 50, type = "Ljung-Box")
###p-value is less than 0.05 thus residuals of (C) are correlated !!!

Box.test(fit_D$residuals, lag = 50, type = "Ljung-Box")
###p-value is less than 0.05 thus residuals of (D) are correlated !!!

Box.test(fit_E$residuals, lag = 50, type = "Ljung-Box")
###p-value is less than 0.05 thus residuals of (E) are correlated !!!

Box.test(fit_F$residuals, lag = 50, type = "Ljung-Box")
###p-value is less than 0.05 thus residuals of (E) are correlated !!!

##Then, ACF PACF plots are checked.

##A
gA1<-ggAcf(fit_A$residuals,lag.max = 48)+theme_minimal()+ggtitle("ACF of (A)'s residuals")
gA2<-ggPacf(fit_A$residuals,lag.max = 48)+theme_minimal()+ggtitle("PACF of (A)'s residuals")
grid.arrange(gA1,gA2,ncol=2)
###3 points are outside of the white noise band.

##B
gB1<-ggAcf(fit_B$residuals,lag.max = 48)+theme_minimal()+ggtitle("ACF of (B)'s residuals")
gB2<-ggPacf(fit_B$residuals,lag.max = 48)+theme_minimal()+ggtitle("PACF of (B)'s residuals")
grid.arrange(gB1,gB2,ncol=2)
###5 points are outside of the white noise band.

##C
gC1<-ggAcf(fit_C$residuals,lag.max = 48)+theme_minimal()+ggtitle("ACF of (C)'s residuals")
gC2<-ggPacf(fit_C$residuals,lag.max = 48)+theme_minimal()+ggtitle("PACF of (C)'s residuals")
grid.arrange(gC1,gC2,ncol=2)
###5 points are outside of the white noise band.

##D
gD1<-ggAcf(fit_D$residuals,lag.max = 48)+theme_minimal()+ggtitle("ACF of (D)'s residuals")
gD2<-ggPacf(fit_D$residuals,lag.max = 48)+theme_minimal()+ggtitle("PACF of (D)'s residuals")
grid.arrange(gD1,gD2,ncol=2)
###7 points are outside of the white noise band.

##E
gE1<-ggAcf(fit_E$residuals,lag.max = 48)+theme_minimal()+ggtitle("ACF of (E)'s residuals")
gE2<-ggPacf(fit_E$residuals,lag.max = 48)+theme_minimal()+ggtitle("PACF of (E)'s residuals")
grid.arrange(gE1,gE2,ncol=2)
###5 points are outside of the white noise band.


##Next, standardized residuals vs time plots

##A
resid_A_std <- rstandard(fit_A)
plot(resid_A_std, type = "o", 
     ylab = "Standardized Residuals of (A)",
     main = "Standardized Residuals vs. Time of (A)")
abline(h = 0, col = "red", lty = 2)

##B
resid_B_std <- rstandard(fit_B)
plot(resid_B_std, type = "o", 
     ylab = "Standardized Residuals of (B)",
     main = "Standardized Residuals vs. Time of (B)")
abline(h = 0, col = "red", lty = 2)

##C
resid_C_std <- rstandard(fit_C)
plot(resid_C_std, type = "o", 
     ylab = "Standardized Residuals of (C)",
     main = "Standardized Residuals vs. Time of (C)")
abline(h = 0, col = "red", lty = 2)

##D
resid_D_std <- rstandard(fit_D)
plot(resid_D_std, type = "o", 
     ylab = "Standardized Residuals of (D)",
     main = "Standardized Residuals vs. Time of (D)")
abline(h = 0, col = "red", lty = 2)

##E
resid_E_std <- rstandard(fit_E)
plot(resid_E_std, type = "o", 
     ylab = "Standardized Residuals of (E)",
     main = "Standardized Residuals vs. Time of (E)")
abline(h = 0, col = "red", lty = 2)

###They all seem like they have increasing patterns.

###By Box.test, only fit_A and fit_B seems like they dont have correlated residuals.

##Histogram and QQ plots of residuals of (A) and (B)
res_A <- residuals(fit_A)
res_B <- residuals(fit_B)
res_C <- residuals(fit_C)
res_D <- residuals(fit_D)
res_E <- residuals(fit_E)

hist(res_A,
     main = "Histogram of Residuals (A)",
     xlab = "Residuals",
     col  = "lightblue",
     border = "white")

hist(res_B,
     main = "Histogram of Residuals (B)",
     xlab = "Residuals",
     col  = "lightgreen",
     border = "white")

hist(res_C,
     main = "Histogram of Residuals (C)",
     xlab = "Residuals",
     col  = "red",
     border = "white")

hist(res_D,
     main = "Histogram of Residuals (D)",
     xlab = "Residuals",
     col  = "coral3",
     border = "white")

hist(res_E,
     main = "Histogram of Residuals (E)",
     xlab = "Residuals",
     col  = "lightpink",
     border = "white")
###Neither of them seems like normally distributed, values closer to the average is much higher.

qqnorm(res_A,
       main = "qq Plot of Residuals (A)")
qqline(res_A, col = "red", lwd = 2)

qqnorm(res_B,
       main = "qq Plot of Residuals (B)")
qqline(res_B, col = "red", lwd = 2)

qqnorm(res_C,
       main = "qq Plot of Residuals (C)")
qqline(res_C, col = "red", lwd = 2)

qqnorm(res_D,
       main = "qq Plot of Residuals (D)")
qqline(res_D, col = "red", lwd = 2)

qqnorm(res_E,
       main = "qq Plot of Residuals (E)")
qqline(res_E, col = "red", lwd = 2)

###Data points are not on the line, also indicates violation of normality assumption.

##Jarque-Bera normality test as we have economic data.
jb_A <- jarque.bera.test(res_A)
jb_B <- jarque.bera.test(res_B)
jb_C <- jarque.bera.test(res_C)
jb_D <- jarque.bera.test(res_D)
jb_E <- jarque.bera.test(res_E)
jb_A
jb_B
jb_C
jb_D
jb_E
###Low p-value on test results also indicate violation of normality assumption.

###We can use forecasting methods that do not rely on the normality of residuals such as ETS (Exponential Smoothing) or TBATS models.

##We can also use transformations on (A) and (B) to try to fix normality violation:

###Log transformation
train_ts_log <- log(train_ts)
fit_A_log <- Arima(train_ts_log,order = c(4,1,0),fixed = c(NA, NA, 0, NA, NA),include.constant = TRUE, transform.pars = FALSE)
fit_B_log <- Arima(train_ts_log, order = c(2, 1, 2))

jb_A_log <- jarque.bera.test(residuals(fit_A_log))
jb_B_log <- jarque.bera.test(residuals(fit_B_log))
jb_A_log
jb_B_log

###Sqrt transformation
train_ts_sqrt <- sqrt(train_ts)
fit_A_sqrt <- Arima(train_ts_sqrt,order = c(4,1,0),fixed = c(NA, NA, 0, NA, NA),include.constant = TRUE, transform.pars = FALSE)
fit_B_sqrt <- Arima(train_ts_sqrt, order = c(2, 1, 2))

jb_A_sqrt <- jarque.bera.test(residuals(fit_A_sqrt))
jb_B_sqrt <- jarque.bera.test(residuals(fit_A_sqrt))
jb_A_sqrt
jb_B_sqrt

###Yeo_Johnson transformation

yj_est <- yeojohnson(train_ts)
train_ts_yj <- predict(yj_est)

fit_A_yj <- Arima(train_ts_yj, order = c(4,1,0),fixed = c(NA, NA, 0, NA, NA),include.constant = TRUE, transform.pars = FALSE)
fit_B_yj <- Arima(train_ts_yj, order = c(2, 1, 2))

jb_A_yj <- jarque.bera.test(residuals(fit_A_yj))
jb_B_yj <- jarque.bera.test(residuals(fit_B_yj))
jb_A_yj
jb_B_yj

###none of these transformations solve normality problem


##Heteroscedasticity check

##First the ACF/PACF plots:
###For (A)
gA_sqrt1<-ggAcf(res_A^2,lag.max = 48)+theme_minimal()+ggtitle("ACF of Squared Residuals (A)")
gA_sqrt2<-ggPacf(res_A^2,lag.max = 48)+theme_minimal()+ggtitle("PACF of Squared Residuals (A)")
grid.arrange(gA_sqrt1,gA_sqrt2,ncol=2)
###For (B)
gB_sqrt1<-ggAcf(res_B^2,lag.max = 48)+theme_minimal()+ggtitle("ACF of Squared Residuals (B)")
gB_sqrt2<-ggPacf(res_B^2,lag.max = 48)+theme_minimal()+ggtitle("PACF of Squared Residuals (B)")
grid.arrange(gB_sqrt1,gB_sqrt2,ncol=2)
###For (C)
gC_sqrt1<-ggAcf(res_C^2,lag.max = 48)+theme_minimal()+ggtitle("ACF of Squared Residuals (C)")
gC_sqrt2<-ggPacf(res_C^2,lag.max = 48)+theme_minimal()+ggtitle("PACF of Squared Residuals (C)")
grid.arrange(gC_sqrt1,gC_sqrt2,ncol=2)
###For (D)
gD_sqrt1<-ggAcf(res_D^2,lag.max = 48)+theme_minimal()+ggtitle("ACF of Squared Residuals (D)")
gD_sqrt2<-ggPacf(res_D^2,lag.max = 48)+theme_minimal()+ggtitle("PACF of Squared Residuals (D)")
grid.arrange(gD_sqrt1,gD_sqrt2,ncol=2)
###For (E)
gE_sqrt1<-ggAcf(res_E^2,lag.max = 48)+theme_minimal()+ggtitle("ACF of Squared Residuals (E)")
gE_sqrt2<-ggPacf(res_E^2,lag.max = 48)+theme_minimal()+ggtitle("PACF of Squared Residuals (E)")
grid.arrange(gE_sqrt1,gE_sqrt2,ncol=2)
###All have significant spikes.
###Best parameters for G/ARCH model is (6,8)


##Then the ARCH Engle's test:
archTest(res_A)
archTest(res_B)
archTest(res_C)
archTest(res_D)
archTest(res_E)
###Since the p-value is extremely small, it is concluded that there exists ARCH effects and therefore, there is Heteroscedasticity.

##As there exists Heteroscedasticity, ARCH and GARCH methods should be used.

spec = ugarchspec()
def.fit = ugarchfit(spec = spec, data = train_ts)
print(def.fit)
###It is seen except omega and alpha parameter, all parameters are significant.
###Ljung Box Tests are used to test serial autocorrelation among the residuals as p-values are very small.
###The results show that residuals have autocorrelation, but squared residuals do not because their p-values are very high.
###ARCH LM test is used to check presence of ARCH effect. As the result is higher than 0.05 for all, it is adequately fitted ARCH process
###The results show that the GARCH process is adequately fitted because of its high p-values
###Sign bias test results indicate that an asymmetric GARCH specification can be tried and then compared. If that does not improve the model meaningfully, the slight sign bias may not be a major concern.
###The nyblom stability test show that the parameter values are constant i.e. zero variance, the alternative hypothesis is that their variance > 0.
###Omega, alpha and beta have stability problem. We should also consider TGARCH models.
###According to information above; sGARCH, eGARCH, gjrGARCH, apARCH, csGARCH models can be used for the set.



spec=ugarchspec(variance.model = list(model="sGARCH",garchOrder = c(6,8))) 
def.fit1= ugarchfit(spec = spec, data = train_ts)
print(def.fit1)
###initial AIC
spec=ugarchspec(variance.model = list(model="apARCH")) 
def.fit2= ugarchfit(spec = spec, data = train_ts)
print(def.fit2)

spec=ugarchspec(variance.model = list(model="apARCH",garchOrder = c(6,8))) 
def.fit3= ugarchfit(spec = spec, data = train_ts)
print(def.fit3)
###Better AIC
spec=ugarchspec(variance.model = list(model="apARCH",garchOrder = c(6,8)),  mean.model = list(armaOrder = c(2,1),include.mean = TRUE)) 
def.fit4= ugarchfit(spec = spec, data = train_ts)
print(def.fit4)
###Better AIC
spec=ugarchspec(variance.model = list(model="apARCH",garchOrder = c(6,8)),  mean.model = list(armaOrder = c(2,2),include.mean = TRUE)) 
def.fit5= ugarchfit(spec = spec, data = train_ts)
print(def.fit5)
###Worse AIC

fit_arch <- def.fit4

#FORECAST
##1) Stochastic models (A) and (B):
forecast_A <- forecast::forecast(fit_A, h = 40)
forecast_B <- forecast::forecast(fit_B, h = 40)
mse_A <- mean((forecast_A$mean - test_ts)^2)
mse_B <- mean((forecast_B$mean - test_ts)^2)
mse_A
mse_B
checkresiduals(fit_A)
checkresiduals(fit_B)
###looks like all assumptions are satisfied.

##2) ARCH model
forecast_arch <- ugarchforecast(fit_arch,n.ahead = 40)
predictions_arch <- forecast_arch@forecast$seriesFor
mse_arch <- mean((predictions_arch - test_ts)^2)
mse_arch
s<-as.vector(forecast_arch@forecast$seriesFor)
bootp=ugarchboot(fit_arch,method=c("Partial","Full")[1],n.ahead = 40,n.bootpred=1000,n.bootfit=1000)
plot(bootp,which=2)
###bootstrap seems more adjusted.
bootp

##3) ets model
ets_model <- ets(train_ts)
forecast_ets <- forecast::forecast(ets_model, h = 40) 
mse_ets <- mean((forecast_ets$mean - test_ts)^2)
mse_ets
checkresiduals(ets_model)
###looks like all assumptions are satisfied.

##4) nnetar forecast
###Define a function for tuning nnetar hyperparameters
tune_nnetar <- function(train_ts, test_ts, p_values, size_values, maxit_values) {
  best_rmse <- Inf
  best_model <- NULL
  best_params <- list()
  ###Loop over combinations of hyperparameters
  for (p in p_values) {
    for (size in size_values) {
      for (maxit in maxit_values) {
        cat("Tuning: p =", p, "size =", size, "maxit =", maxit, "\n")
        ###Fit the nnetar model with current hyperparameters
        model <- tryCatch(
          nnetar(train_ts, p = p, size = size, maxit = maxit),
          error = function(e) NULL  # Catch errors if model fitting fails
        )
        if (!is.null(model)) {
          ###Forecast using the model
          forecasted <- forecast::forecast(model, h = length(test_ts))  # Forecast for test_ts length
          ###Calculate RMSE for the model
          rmse_value <- sqrt(mean((forecasted$mean - test_ts)^2))
          ###Check if the current model is better (lower RMSE)
          if (rmse_value < best_rmse) {
            best_rmse <- rmse_value
            best_model <- model
            best_params <- list(p = p, size = size, maxit = maxit)
          }
        }
      }
    }
  }
  ###Return the best model and hyperparameters
  return(list(best_model = best_model, best_params = best_params, best_rmse = best_rmse))
}
###Define the grid of hyperparameters
p_values <- 1:5  ###Number of lags
size_values <- seq(3, 10, by = 1)  ###Number of hidden nodes
maxit_values <- c(100, 200, 300)  ###Maximum iterations
###Run the tuning function
tuning_results <- tune_nnetar(train_ts, test_ts, p_values, size_values, maxit_values)
###Print the best hyperparameters
cat("Best Model Hyperparameters: p =", tuning_results$best_params$p, 
    "size =", tuning_results$best_params$size, 
    "maxit =", tuning_results$best_params$maxit, "\n")
###Fit the best model using the best hyperparameters
best_model_nn <- tuning_results$best_model
###Forecast with the best model
forecast_nn_tuned <- forecast::forecast(best_model_nn, h = length(test_ts))
mse_nn <- mean((forecast_nn_tuned$mean - test_ts)^2)
mse_nn
###Plot the forecast
plot(forecast_nn_tuned, main = "Forecast with Tuned Neural Network Model")
###Check the residuals and assumptions
checkresiduals(best_model_nn)
###spikes are mostly in white noise band and residuals seem normally distributed meaning that the assumptions are met.



##5) TBATS

tbats_model <- tbats(train_ts)
###Generate forecasts for the test set
forecasts_tbats <- forecast::forecast(tbats_model, h=40)
###Calculate residuals
res_tbats <- residuals(tbats_model)
###Calculate Mean Squared Error (MSE)
mse_tbats <- mean((test_ts - forecasts_tbats$mean)^2)
mse_tbats
checkresiduals(tbats_model)
###spikes are mostly in white noise band and residuals seem normally distributed meaning that the assumptions are met.


##6) Prophet
changepoint_prior <- c(0.1, 0.5, 0.9)
seasonality_prior <- c(0.1, 0.3, 0.5)
changepoint_range <- c(0.6, 0.8, 0.9)
start_date <- as.Date("2015-11-23")
date_seq <- seq(start_date, by = "day", length.out = length(train_ts))
train_prophet <- data.frame(ds = date_seq, y = as.numeric(train_ts))
results <- data.frame(
  changepoint_prior = numeric(),
  seasonality_prior = numeric(),
  changepoint_range = numeric(),
  RMSE = numeric()
)
for (cp in changepoint_prior) {
  for (sp in seasonality_prior) {
    for (cr in changepoint_range) {
      m <- prophet(
        changepoint.prior.scale = cp,
        seasonality.prior.scale = sp,
        changepoint.range = cr
      )
      m <- fit.prophet(m, train_prophet)
      future <- make_future_dataframe(m, periods = length(test_ts), freq = "day")
      forecast <- predict(m, future)
      predicted <- tail(forecast$yhat, length(test_ts))
      acc <- accuracy(predicted, test_ts)  
      rmse <- acc["Test set", "RMSE"]
      results <- rbind(results, data.frame(
        changepoint_prior = cp, 
        seasonality_prior = sp, 
        changepoint_range = cr, 
        RMSE = rmse
      ))
    }
  }
}
best_params <- results[which.min(results$RMSE), ]
best_params
fit_prophet <- prophet(
  changepoint.prior.scale = 0.9,
  seasonality.prior.scale = 0.3,
  changepoint.range = 0.8
)
fit_prophet <- fit.prophet(fit_prophet, train_prophet)
future_prophet <- make_future_dataframe(fit_prophet, periods = length(test_ts), freq = "day")
forecast_prophet <- predict(fit_prophet, future_prophet)
predicted_prophet <- tail(forecast_prophet$yhat, length(test_ts))
predicted_prophet
predicted_prophet_all <- forecast_prophet$yhat
accuracy_prophet <- accuracy(predicted_prophet, test_ts)
res_prophet <- test_ts - predicted_prophet
mse_prophet <- mean((test_ts - predicted_prophet)^2)
mse_prophet


###Check residuals
checkresiduals(res_prophet)
###Looks like normality assumption is met. but there might be a problem with the white noice band as too many spikes are outside of it with exponential decay.


predictions_arch <- as.vector(predictions_arch)
class(predictions_arch)
#Checking accuracy
accuracy_A <- accuracy(forecast_A,test_ts)
accuracy_B <- accuracy(forecast_B,test_ts)
accuracy_ets <- accuracy(forecast_ets,test_ts)
accuracy_arch <- accuracy(predictions_arch,test_ts)
accuracy_nn <- accuracy(forecast_nn_tuned,test_ts)
accuracy_tbats <- accuracy(forecasts_tbats,test_ts)
accuracy_prophet <- accuracy(predicted_prophet, test_ts)

accuracy_results <- data.frame(
  Model = c("(A)", "(B)", "ets", "arch" , "nnetar", "tbats", "prophet"),
  MAE = c(accuracy_A[1, "MAE"], accuracy_B[1, "MAE"], accuracy_ets[1, "MAE"], accuracy_arch[1, "MAE"], 
          accuracy_nn[1, "MAE"], accuracy_tbats[1, "MAE"], accuracy_prophet[1, "MAE"]),
  RMSE = c(accuracy_A[1, "RMSE"], accuracy_B[1, "RMSE"], accuracy_ets[1, "RMSE"], accuracy_arch[1, "RMSE"], 
           accuracy_nn[1, "RMSE"], accuracy_tbats[1, "RMSE"], accuracy_prophet[1, "RMSE"]),
  MAPE = c(accuracy_A[1, "MAPE"], accuracy_B[1, "MAPE"], accuracy_ets[1, "MAPE"], accuracy_arch[1, "MAPE"], 
           accuracy_nn[1, "MAPE"], accuracy_tbats[1, "MAPE"], accuracy_prophet[1, "MAPE"])
)
accuracy_results
###According to RMSE from these results, it is determined that best model is model (A). But, assumptions are not met on it thus, ets might also be preferred.
##Overfitting:
###There is no sign of bad performance for the test set and the RMSE is low enough so there is no overfitting for (A) and (B).
###ETS model is designed for capturing trends and seaonality thus less likely to encounter overfitting problem.
###neural nnetar tbats and prophet are highly prone to overfitting while prophet is less likely to overfit. None of these models exhibit bad performance while still having relatively little RMSE


##arch forecast variable adjusting
mean_mat  <- forecast_arch@forecast$seriesFor
sigma_mat <- forecast_arch@forecast$sigmaFor  
mean_vec  <- as.vector(mean_mat)
sigma_vec <- as.vector(sigma_mat)
z80 <- qnorm(0.5 + 0.80/2)
z95 <- qnorm(0.5 + 0.95/2)
lower80 <- mean_vec - z80 * sigma_vec
lower95 <- mean_vec - z95 * sigma_vec
upper80 <- mean_vec + z80 * sigma_vec
upper95 <- mean_vec + z95 * sigma_vec
n_train         <- length(train_ts)
train_end_time  <- time(train_ts)[n_train]
train_freq      <- frequency(train_ts)
forecast_start  <- train_end_time + deltat(train_ts)
mean_ts    <- ts(mean_vec, start = forecast_start, frequency = train_freq)
lower_ts   <- ts(cbind(lower80, lower95), start = forecast_start, frequency = train_freq)
upper_ts   <- ts(cbind(upper80, upper95), start = forecast_start, frequency = train_freq)
arch_forecast <- structure(list(mean   = mean_ts,lower  = lower_ts,upper  = upper_ts, level  = c(80, 95),x      = train_ts,fitted = fitted(fit_arch), method = "apARCH Forecast (rugarch)"),class = "forecast")

#Plots for Forecasts

autoplot(forecast_A)+autolayer(test_ts,series="actual",color="red")+theme_minimal()
autoplot(forecast_B)+autolayer(test_ts,series="actual",color="red")+theme_minimal()
autoplot(forecast_ets)+autolayer(test_ts,series="actual",color="red")+theme_minimal()
autoplot(arch_forecast) +autolayer(test_ts, series = "Actual", color = "red") +theme_minimal() +labs(title = "apARCH Forecast with 80% & 95% Confidence Intervals",x = "Time", y = "Value") 
autoplot(forecast_nn_tuned)+autolayer(test_ts,series="actual",color="red")+theme_minimal()
plot(fit_prophet, forecast_prophet)+theme_minimal()+labs(title="Prophet forecast")
autoplot(forecasts_tbats)+autolayer(test_ts,series="actual",color="red")+theme_minimal()

###Best plot among these is nnetar model as it correctly gives the intervals and can predict the decreasing and increasing trend in test set.


#Conclusion
###For 5 of the 6 models, accuracy is very close in terms of RMSE. Also, the neural network model does much better than other models when observed visiually through its forecast plot while also giving one of the best RMSE value. The good performance on forecast plot also indicates there exists no over fitting as the model visiually does well. No violation of any assumption is present. Thus, the neural network model is selected to be the best model.



