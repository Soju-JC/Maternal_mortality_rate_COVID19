# rm(list = ls()) 

# Options for karma():
# resid = 1 : standardized residual (or Person’s)
# resid = 2 : deviance residuals 2
# resid = 3 : quantile residuals

# diag = 0 : do not plot any graph (useful for simulations)
# diag = 1 : plot graphs
# diag = 2 : save graphs on ps files

# link = "logit"
# link = "probit"
# link = "cloglog"

# h is the prediction window

# Example
# karma(
#   y, # data
#   ar = p, # AR orders (accepts individual coef)
#   ma = q, # MA orders (accepts individual coef)
#   h = 7, # prediction window = 7
#   diag = 0, # do not plot
#   resid = 1, # standardized residual
#   link = "logit", # g(*)
#   prec_start = 5 # First guess for precision
# )
#-------------------------------------------------------------------------------
#detach("package:tseries", unload = TRUE)
#detach("package:tidyverse", unload = TRUE)
#detach("package:itsmr", unload = TRUE)

library("tseries")
library("tidyverse")
library("itsmr")
library("tsoutliers")

# Load the data
df <- readRDS("maternal_mortality_rate_2021.rds")
df <- df[1:135,] # up to day 15/May
df_ts <- ts(df$rate, start = c(2021, 1), frequency = 365)

h <- 12 # forecast window

# Split data in train and test
df_train <- df_ts[1:(length(df_ts) - h)]
df_train <- ts(df_train)

df_test <- df_ts[(length(df_train) + 1):(length(df_ts))]
df_test <- ts(df_test)

plot.ts(df_train)
#-------------------------------------------------------------------------------
# Required functions
# https://github.com/fabiobayer/KARMA
source('supporting_scripts/kum-mu-phi.r')
# Modified version including precision guess value
source('supporting_scripts/karma.modified.fit.r') 
# Modified version including precision guess value
source('supporting_scripts/karma.modified.r') 

#d <- 0 # No difference transformation
d <- 1 # One difference transformation

if(d == 0){
  y <- df_train
} else if (d > 0) {
  y <- diff(df_train, differences = d)
  # Transform the data into double-bounded (0-1)
  a = min(y)
  b = max(y)
  y = (y - a)/(b - a)
} else {warning("d NOT SUPPORTED!")}

n_data = length(y)

# Correction for values on the extremes(when d > 0 with a = min() and b = max())
#(According to Cribari-Neto, Zeileis, 2009, "Beta Regression in R")
extremes = FALSE
if (min(y) == 0 || max(y) == 1) {
  extremes = TRUE
  # Smithson and Verkuilen (2006)
  y = (y*(n_data-1)+0.5)/n_data
}

#-------------------------- Fit KARMA ------------------------------------------
#Possible combinations including: 0, 1 parameter, 2 parameters, ... 6 parameters
#1 + 6 + 15 + 20 + 15 + 6 + 1 = 64 models 
#64*64 = 4096 models considering, ar 0 up to ar 6 and ma 0 up to ma 6

source('supporting_scripts/model_orders.r')

list_combinations <- list_combinations_BK6 # all combinations up to order 6

# Start variables
valid_models_ar_coef <- list()
valid_models_ma_coef <- list()
valid_models_aic <- list()
valid_models_bic <- list()
  
for (k in 1:length(list_combinations)) {
  ar <- list_combinations[[k]]
  for (j in 1:length(list_combinations)) {
    tryCatch({ # don't stop on errors
      fit_karma <- karma(
        y,
        ar = ar,
        ma = list_combinations[[j]],
        h = h,
        diag = 0,
        #resid = 1,
        link = "logit",
        prec_start = 5
      )
      coef_df <- as.data.frame(fit_karma$model)
      coef_pvalues <- coef_df$`Pr(>|z|)` # p-values
      
      coef_names <- rownames(coef_df)
      coef_ar <- coef_df[grepl("phi", coef_names, fixed = TRUE), ]
      coef_ar_values <- coef_ar$Estimate # ar coef
      coef_ma <- coef_df[grepl("theta", coef_names, fixed = TRUE), ]
      coef_ma_values <- coef_ma$Estimate # ma coef
      
      # Verify the valid models
      if (
        all(coef_pvalues < 0.05) & # 5% significance
        length(coef_ar_values) != 0 & # At leat 1 ar coef
        all(coef_ar_values > -0.9 & coef_ar_values < 0.9) & # Causal model
        all(coef_ma_values > -1.0 & coef_ma_values < 1.0) 
        ){
        valid_models_ar_coef <- append(
          valid_models_ar_coef, 
          list(ar)
        )
        valid_models_ma_coef <- append(
          valid_models_ma_coef, 
          list(list_combinations[[j]])
        )
        valid_models_aic <- append(
          valid_models_aic, 
          list(fit_karma$aic)
        )
        valid_models_bic <- append(
          valid_models_bic, 
          list(fit_karma$bic)
        )
      }
    }, error = function(e) {
      
    }, warning = function(w) {
      
    })
  }
}

valid_models_ar_coef # - of 4096 
valid_models_ma_coef # - of 4096
valid_models_aic # min: AIC - | model selected: AIC -
valid_models_bic # min: BIC - | ~

# BEST MODEL BY AIC
which.min(valid_models_aic) # position of min in the vector
valid_models_aic[[which.min(valid_models_aic)]] # AIC 
best_arma_combination_ar_aic <- 
  valid_models_ar_coef[[which.min(valid_models_aic)]] # best
best_arma_combination_ma_aic <- 
  valid_models_ma_coef[[which.min(valid_models_aic)]] # best

fit_karma_best_aic <- karma(
  y,
  ar = best_arma_combination_ar_aic, # 
  ma = best_arma_combination_ma_aic, # 
  h = h,
  diag = 1,
  resid = 3,
  link = "logit",
  prec_start = 5
)

cpgram(fit_karma_best_aic$resid3, main = "Cumulative Periodogram of Residuals")
car::qqPlot(fit_karma_best_aic$resid3)
tseries::jarque.bera.test(fit_karma_best_aic$resid3)
Box.test(fit_karma_best_aic$resid3, type = "Box-Pierce")
Box.test(fit_karma_best_aic$resid3, type = "Ljung-Box")
forecast::checkresiduals(fit_karma_best_aic$resid3, test = "LB")
acf(fit_karma_best_aic$resid3)

# BEST MODEL BY BIC
which.min(valid_models_bic) # position of min in the vector
valid_models_bic[[which.min(valid_models_bic)]] # BIC
best_arma_combination_ar_bic <- 
  valid_models_ar_coef[[which.min(valid_models_bic)]] # best
best_arma_combination_ma_bic <- 
  valid_models_ma_coef[[which.min(valid_models_bic)]] # best

fit_karma_best_bic <- karma(
  y,
  ar = best_arma_combination_ar_bic, # 
  ma = best_arma_combination_ma_bic, # 
  h = h,
  diag = 1,
  resid = 3,
  link = "logit",
  prec_start = 5
)

cpgram(fit_karma_best_bic$resid3, main = "Cumulative Periodogram of Residuals")
car::qqPlot(fit_karma_best_bic$resid3)
tseries::jarque.bera.test(fit_karma_best_bic$resid3)
Box.test(fit_karma_best_bic$resid3, type = "Box-Pierce")
Box.test(fit_karma_best_bic$resid3, type = "Ljung-Box")
forecast::checkresiduals(fit_karma_best_bic$resid3, test = "LB")
acf(fit_karma_best_bic$resid3)
