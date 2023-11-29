# rm(list = ls()) 
#detach("package:tseries", unload = TRUE)
#detach("package:tidyverse", unload = TRUE)
#detach("package:itsmr", unload = TRUE)

library("tseries")
library("tidyverse")
library("itsmr")

# Load the data
df <- readRDS("maternal_mortality_rate_2021.rds")
df <- df[1:135,] # Dados até 15 de maio
# Transform the data to time series object (day frequency)
df_ts <- ts(df$rate, start = c(2021, 1), frequency = 365)

h <- 12 # forecast window

# Split data in train and test
df_train <- df_ts[1:(length(df_ts) - h)]
df_train <- ts(df_train)

plot.ts(df_train)
#--------------------------- Fit ARMA ------------------------------------------
#Possible combinations including: 0, 1 parameter, 2 parameters, ... 6 parameters
#1 + 6 + 15 + 20 + 15 + 6 + 1 = 64 models 
#64*64 = 4096 models considering, ar 0 up to ar 6 and ma 0 up to ma 6

# All coef combinations
list_combinations <- list(
  "c(NA, NA, NA, NA, NA, NA)" = c(NA, NA, NA, NA, NA, NA), 
  "c(0, NA, NA, NA, NA, NA)" = c(0, NA, NA, NA, NA, NA),
  "c(NA, 0, NA, NA, NA, NA)" = c(NA, 0, NA, NA, NA, NA),
  "c(NA, NA, 0, NA, NA, NA)" = c(NA, NA, 0, NA, NA, NA), 
  "c(NA, NA, NA, 0, NA, NA)" = c(NA, NA, NA, 0, NA, NA), 
  "c(NA, NA, NA, NA, 0, NA)" = c(NA, NA, NA, NA, 0, NA), 
  "c(NA, NA, NA, NA, NA, 0)" = c(NA, NA, NA, NA, NA, 0),
  "c(0, 0, NA, NA, NA, NA)" = c(0, 0, NA, NA, NA, NA), 
  "c(0, NA, 0, NA, NA, NA)" = c(0, NA, 0, NA, NA, NA),
  "c(0, NA, NA, 0, NA, NA)" = c(0, NA, NA, 0, NA, NA),
  "c(0, NA, NA, NA, 0, NA)" = c(0, NA, NA, NA, 0, NA), 
  "c(0, NA, NA, NA, NA, 0)" = c(0, NA, NA, NA, NA, 0), 
  "c(NA, 0, 0, NA, NA, NA)" = c(NA, 0, 0, NA, NA, NA), 
  "c(NA, 0, NA, 0, NA, NA)" = c(NA, 0, NA, 0, NA, NA), 
  "c(NA, 0, NA, NA, 0, NA)" = c(NA, 0, NA, NA, 0, NA), 
  "c(NA, 0, NA, NA, NA, 0)" = c(NA, 0, NA, NA, NA, 0), 
  "c(NA, NA, 0, 0, NA, NA)" = c(NA, NA, 0, 0, NA, NA), 
  "c(NA, NA, 0, NA, 0, NA)" = c(NA, NA, 0, NA, 0, NA), 
  "c(NA, NA, 0, NA, NA, 0)" = c(NA, NA, 0, NA, NA, 0),
  "c(NA, NA, NA, 0, 0, NA)" = c(NA, NA, NA, 0, 0, NA),
  "c(NA, NA, NA, 0, NA, 0)" = c(NA, NA, NA, 0, NA, 0),
  "c(NA, NA, NA, NA, 0, 0)" = c(NA, NA, NA, NA, 0, 0),
  "c(0, 0, 0, NA, NA, NA)" = c(0, 0, 0, NA, NA, NA),
  "c(0, 0, NA, 0, NA, NA)" = c(0, 0, NA, 0, NA, NA),
  "c(0, 0, NA, NA, 0, NA)" = c(0, 0, NA, NA, 0, NA),
  "c(0, 0, NA, NA, NA, 0)" = c(0, 0, NA, NA, NA, 0),
  "c(0, NA, 0, 0, NA, NA)" = c(0, NA, 0, 0, NA, NA),
  "c(0, NA, 0, NA, 0, NA)" = c(0, NA, 0, NA, 0, NA),
  "c(0, NA, 0, NA, NA, 0)" = c(0, NA, 0, NA, NA, 0),
  "c(0, NA, NA, 0, 0, NA)" = c(0, NA, NA, 0, 0, NA),
  "c(0, NA, NA, 0, NA, 0)" = c(0, NA, NA, 0, NA, 0),
  "c(0, NA, NA, NA, 0, 0)" = c(0, NA, NA, NA, 0, 0),
  "c(NA, 0, 0, 0, NA, NA)" = c(NA, 0, 0, 0, NA, NA),
  "c(NA, 0, 0, NA, 0, NA)" = c(NA, 0, 0, NA, 0, NA),
  "c(NA, 0, 0, NA, NA, 0)" = c(NA, 0, 0, NA, NA, 0),
  "c(NA, 0, NA, 0, 0, NA)" = c(NA, 0, NA, 0, 0, NA),
  "c(NA, 0, NA, 0, NA, 0)" = c(NA, 0, NA, 0, NA, 0),
  "c(NA, 0, NA, NA, 0, 0)" = c(NA, 0, NA, NA, 0, 0),
  "c(NA, NA, 0, 0, 0, NA)" = c(NA, NA, 0, 0, 0, NA),
  "c(NA, NA, 0, 0, NA, 0)" = c(NA, NA, 0, 0, NA, 0),
  "c(NA, NA, 0, NA, 0, 0)" = c(NA, NA, 0, NA, 0, 0),
  "c(NA, NA, NA, 0, 0, 0)" = c(NA, NA, NA, 0, 0, 0),
  "c(0, 0, 0, 0, NA, NA)" = c(0, 0, 0, 0, NA, NA),
  "c(0, 0, 0, NA, 0, NA)" = c(0, 0, 0, NA, 0, NA),
  "c(0, 0, 0, NA, NA, 0)" = c(0, 0, 0, NA, NA, 0),
  "c(0, 0, NA, 0, 0, NA)" = c(0, 0, NA, 0, 0, NA),
  "c(0, 0, NA, 0, NA, 0)" = c(0, 0, NA, 0, NA, 0),
  "c(0, 0, NA, NA, 0, 0)" = c(0, 0, NA, NA, 0, 0),
  "c(0, NA, 0, 0, 0, NA)" = c(0, NA, 0, 0, 0, NA),
  "c(0, NA, 0, 0, NA, 0)" = c(0, NA, 0, 0, NA, 0),
  "c(0, NA, 0, NA, 0, 0)" = c(0, NA, 0, NA, 0, 0),
  "c(0, NA, NA, 0, 0, 0)" = c(0, NA, NA, 0, 0, 0),
  "c(NA, 0, 0, 0, 0, NA)" = c(NA, 0, 0, 0, 0, NA),
  "c(NA, 0, 0, 0, NA, 0)" = c(NA, 0, 0, 0, NA, 0),
  "c(NA, 0, 0, NA, 0, 0)" = c(NA, 0, 0, NA, 0, 0),
  "c(NA, 0, NA, 0, 0, 0)" = c(NA, 0, NA, 0, 0, 0),
  "c(NA, NA, 0, 0, 0, 0)" = c(NA, NA, 0, 0, 0, 0),
  "c(0, 0, 0, 0, 0, NA)" = c(0, 0, 0, 0, 0, NA),
  "c(0, 0, 0, 0, NA, 0)" = c(0, 0, 0, 0, NA, 0),
  "c(0, 0, 0, NA, 0, 0)" = c(0, 0, 0, NA, 0, 0),
  "c(0, 0, NA, 0, 0, 0)" = c(0, 0, NA, 0, 0, 0),
  "c(0, NA, 0, 0, 0, 0)" = c(0, NA, 0, 0, 0, 0),
  "c(NA, 0, 0, 0, 0, 0)" = c(NA, 0, 0, 0, 0, 0),
  "c(0, 0, 0, 0, 0, 0)" = c(0, 0, 0, 0, 0, 0)
  )

# Start variables
valid_models <- list()
valid_models_aic <- list()
valid_models_bic <- list()
  
for (k in 1:length(list_combinations)) {
  ar <- list_combinations[[k]]
  for (j in 1:length(list_combinations)) {
    tryCatch({ # don't stop on errors
      fit_arima <- forecast::Arima(
        df_train, 
        order = c(6, 1, 6), 
        fixed = append(ar, list_combinations[[j]])
      )
      coef_df <- lmtest::coeftest(fit_arima)
      coef_df <- as.data.frame(coef_df[,])
      coef_pvalues <- coef_df$`Pr(>|z|)` # p-values
      
      coef_names <- rownames(coef_df)
      coef_ar <- coef_df[grepl("ar", coef_names, fixed = TRUE), ]
      coef_ar_values <- coef_ar$Estimate # ar coef
      coef_ma <- coef_df[grepl("ma", coef_names, fixed = TRUE), ]
      coef_ma_values <- coef_ma$Estimate # ma coef
      
      if (
        all(coef_pvalues < 0.05) & # 5% significance
        length(coef_ar_values) != 0 & # At leat 1 ar coef
        all(coef_ar_values > -0.9 & coef_ar_values < 0.9) & # Causal model
        all(coef_ma_values > -1.0 & coef_ma_values < 1.0) # Invertible model
        ){
        valid_models <- append(
        valid_models, 
        list(append(ar, list_combinations[[j]]))
        )
        valid_models_aic <- append(
          valid_models_aic, 
          list(fit_arima$aic)
        )
        valid_models_bic <- append(
          valid_models_bic, 
          list(fit_arima$bic)
        )
      }
    }, error = function(e) {
      
    }, warning = function(w) {
      
    })
  }
}

valid_models # 62 of 4096 
valid_models_aic # min: AIC -1516.579
valid_models_bic # min: BIC -1504.638

# BEST MODEL BY AIC
which.min(valid_models_aic) # position of min in the vector
valid_models_aic[[which.min(valid_models_aic)]] # AIC 
best_arma_combination_aic <- valid_models[[which.min(valid_models_aic)]] # best

fit_arima <- forecast::Arima(
  df_train, 
  order = c(6, 1, 6), 
  fixed = best_arma_combination_aic # NA  0  0  0  0 NA  0 NA  0  0  0 NA
)
lmtest::coeftest(fit_arima);cat("AIC", fit_arima$aic);cat("|BIC", fit_arima$bic)

cpgram(fit_arima$resid, main = "Cumulative Periodogram of Residuals")
car::qqPlot(fit_arima$res)
tseries::jarque.bera.test(residuals(fit_arima))
plotc(df_train,fit_arima$fitted)
Box.test(residuals(fit_arima), type = "Box-Pierce")
Box.test(residuals(fit_arima), type = "Ljung-Box")
forecast::checkresiduals(fit_arima, test = "LB")

# BEST MODEL BY BIC
which.min(valid_models_bic) # position of min in the vector
valid_models_bic[[which.min(valid_models_bic)]] # BIC
best_arma_combination_bic <- valid_models[[which.min(valid_models_bic)]] # best

fit_arima <- forecast::Arima(
  df_train, 
  order = c(6, 1, 6), 
  fixed = best_arma_combination_bic # 0  0  0  0 NA  0 NA  0  0  0 NA  0
)
lmtest::coeftest(fit_arima);cat("AIC", fit_arima$aic);cat("|BIC", fit_arima$bic)

cpgram(fit_arima$resid, main = "Cumulative Periodogram of Residuals")
car::qqPlot(fit_arima$res)
tseries::jarque.bera.test(residuals(fit_arima))
plotc(df_train,fit_arima$fitted)
Box.test(residuals(fit_arima), type = "Box-Pierce")
Box.test(residuals(fit_arima), type = "Ljung-Box")
forecast::checkresiduals(fit_arima, test = "LB")
