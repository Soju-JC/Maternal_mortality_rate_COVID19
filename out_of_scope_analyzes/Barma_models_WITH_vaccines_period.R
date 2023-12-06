# rm(list = ls()) 

# Options for barma():
# resid = 1 : standardized residual
# resid = 2 : standardized residual 2
# resid = 3 : standardized weighted residual
# resid = 4 : deviance residual

# diag = 0 : do not plot any graph (useful for simulations)
# diag = 1 : plot graphs
# diag = 2 : save graphs on ps files

# link = "logit"
# link = "probit"
# link = "cloglog"

# h is the prediction window

#Example
#barma(y, ar = p, ma = q, h = h, diag = 0, resid = 1, link = "logit")

#detach("package:tseries", unload = TRUE)
#detach("package:tidyverse", unload = TRUE)
#detach("package:itsmr", unload = TRUE)

library("tseries")
library("tidyverse")
library("itsmr")
library("tsoutliers")

# Load the data
df <- readRDS("maternal_mortality_rate_2021.rds")

# Transform the data to time series object (day frequency)
df_ts <- ts(df$rate, start = c(2021, 1), frequency = 365)

h <- 7 # forecast window

# Split data in train and test
df_train <- df_ts[1:(length(df_ts) - h)]
df_train <- ts(df_train)

#-------------------------------------------------------------------------------
########################## Análise de intervenção ##############################
#-------------------------------------------------------------------------------
# outliers_excess_ts <- tsoutliers::tso(df_train) #Detectando outiliers.
# plot(outliers_excess_ts) #Gráfico de detecção.
# outliers_excess_ts$outliers #Informações gerais dos pontos encontrados.
# outliers_idx <- outliers_excess_ts$outliers$ind #Posição dos outliers na série.
# 
# #length of our time series
# n <- length(df_train) # 128 if up to 135
# 
# (col_int <- outliers(c("TC", "AO"), outliers_idx)) #Colunas de interesse.
# 
# #Esqueleto da série identificando a posição do outlier.
# esq_posi <- outliers.effects(col_int, n) 
# # efeito_outlier_series <- as.matrix(
# #   rowSums(esq_posi[, c("LS81")])
# # )
# efeito_outlier_series <- esq_posi
# 
# #Coeficientes dos efeitos dos outliers.
# omega_hat <- unlist(outliers_excess_ts$outliers["coefhat"])
# 
# #Calculando vetor que representa o efeito do outlier.
# out_effect <- efeito_outlier_series
# out_effect <- as.data.frame(out_effect)
# colnames(out_effect) <- c("TC149", "AO171")
# out_effect[out_effect$AO171 == 1, ] <- as.numeric(omega_hat[2])
# out_effect$TC149 <- out_effect$TC149 * as.numeric(omega_hat[1])
# ao_effect_ts <- ts(out_effect$AO171)
# ts_effect_ts <- ts(out_effect$TC149)
# 
# #Substraindo o efeito da intervenção.
# df_train_clean <- df_train - ao_effect_ts # extract aditive effect
# 
# #ts_effect_ts[171] <- 0 # This observation has been treated as aditive
# df_train_clean <- df_train_clean - ts_effect_ts # extract transient effect
# 
# plot(cbind("Original" = df_train,
#            "Without outliers" = df_train_clean,
#            "Additive effect" = ao_effect_ts,
#            "Transient change effect" = ts_effect_ts),
#      main = "Time series and outlier effects")
# 
# 
# plot.ts(df_train_clean)
#-------------------------------------------------------------------------------
# Required functions
source("supporting_scripts/barma.r")
source("supporting_scripts/barma.fit.r")


#d <- 0 # No difference transformation
d <- 1 # One difference transformation

# Apply difference transformation to make data stationary (when applied)
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
  y = (y*(n_data-1)+0.5)/n_data
}

#-------------------------- Fit BARMA ------------------------------------------
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
      fit_barma <- barma(
        y,
        ar = ar,
        ma = list_combinations[[j]],
        h = h,
        diag = 0,
        #resid = 1,
        link = "logit"
      )
      coef_df <- as.data.frame(fit_barma$model)
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
        all(coef_ma_values > -1.0 & coef_ma_values < 1.0) # Invertible model
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
          list(fit_barma$aic)
        )
        valid_models_bic <- append(
          valid_models_bic, 
          list(fit_barma$bic)
        )
      }
    }, error = function(e) {
      
    }, warning = function(w) {
      
    })
  }
}

valid_models_ar_coef # 217 of 4096 (14: 61)
valid_models_ma_coef # 217 of 4096
valid_models_aic # min: AIC -156.8501
valid_models_bic # min: BIC -131.2524

# BEST MODEL BY AIC
which.min(valid_models_aic) # position of min in the vector
valid_models_aic[[which.min(valid_models_aic)]] # AIC 
best_arma_combination_ar_aic <- 
  valid_models_ar_coef[[which.min(valid_models_aic)]] # best
best_arma_combination_ma_aic <- 
  valid_models_ma_coef[[which.min(valid_models_aic)]] # best

fit_barma_best_aic <- barma(
  y,
  ar = best_arma_combination_ar_aic, # 1 2 5
  ma = best_arma_combination_ma_aic, # 1 3 5 6
  h = h,
  diag = 1,
  resid = 1,
  link = "logit"
)

cpgram(fit_barma_best_aic$resid1, main = "Cumulative Periodogram of Residuals")
car::qqPlot(fit_barma_best_aic$resid1)
tseries::jarque.bera.test(fit_barma_best_aic$resid1)
Box.test(fit_barma_best_aic$resid1, type = "Box-Pierce")
Box.test(fit_barma_best_aic$resid1, type = "Ljung-Box")
forecast::checkresiduals(fit_barma_best_aic$resid1, test = "LB")
acf(fit_barma_best_aic$resid1)

# BEST MODEL BY BIC
which.min(valid_models_bic) # position of min in the vector
valid_models_bic[[which.min(valid_models_bic)]] # BIC
best_arma_combination_ar_bic <- 
  valid_models_ar_coef[[which.min(valid_models_bic)]] # best
best_arma_combination_ma_bic <- 
  valid_models_ma_coef[[which.min(valid_models_bic)]] # best

fit_barma_best_bic <- barma(
  y,
  ar = best_arma_combination_ar_bic, # 1 2 5
  ma = best_arma_combination_ma_bic, # 1 3 5 6
  h = h,
  diag = 1,
  resid = 1,
  link = "logit"
)

cpgram(fit_barma_best_bic$resid1, main = "Cumulative Periodogram of Residuals")
car::qqPlot(fit_barma_best_bic$resid1)
tseries::jarque.bera.test(fit_barma_best_bic$resid1)
Box.test(fit_barma_best_bic$resid1, type = "Box-Pierce")
Box.test(fit_barma_best_bic$resid1, type = "Ljung-Box")
forecast::checkresiduals(fit_barma_best_bic$resid1, test = "LB")
acf(fit_barma_best_bic$resid1)
