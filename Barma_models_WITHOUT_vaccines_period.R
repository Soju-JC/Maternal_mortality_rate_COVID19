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
#df <- df[1:120,] # Dados fim de abril
df <- df[1:135,] # Dados até 15 de maio
df_ts <- ts(df$rate, start = c(2021, 1), frequency = 365)

h <- 12 # forecast window

# Split data in train and test
df_train <- df_ts[1:(length(df_ts) - h)]
df_train <- ts(df_train)

df_test <- df_ts[(length(df_train) + 1):(length(df_ts))]
df_test <- ts(df_test)

plot.ts(df_train)

#-------------------------------------------------------------------------------
########################## Análise de intervenção ##########################
#-------------------------------------------------------------------------------
outliers_excess_ts <- tsoutliers::tso(df_train)
plot(outliers_excess_ts) #Gráfico de detecção.
outliers_excess_ts$outliers #Informações gerais dos pontos encontrados.
outliers_idx <- outliers_excess_ts$outliers$ind #Posição dos outliers na série.
# plot(forecast::tsclean(df_train))
# plot(df_train)
#length of our time series
n <- length(df_train)

(col_int <- outliers(c("AO", "LS", "AO", "AO"), outliers_idx[c(1, 3, 4, 6)])) #Colunas de interesse.

#Esqueleto da série identificando a posição do outlier.
esq_posi <- outliers.effects(col_int, n)
# efeito_outlier_series <- as.matrix(
#   rowSums(esq_posi[, c("LS81")])
# )
efeito_outlier_series <- esq_posi

#Coeficientes dos efeitos dos outliers.
omega_hat <- unlist(outliers_excess_ts$outliers["coefhat"])

#Calculando vetor que representa o efeito do outlier.
out_effect <- efeito_outlier_series
out_effect <- as.data.frame(out_effect)
colnames(out_effect) <- c("AO65", "LS81", "AO101", "AO115")
out_effect[out_effect$AO65 == 1, "AO65"] <- as.numeric(omega_hat[1])
out_effect[out_effect$LS81 == 1, "LS81"] <- as.numeric(omega_hat[3])
out_effect[out_effect$AO101 == 1, "AO101"] <- as.numeric(omega_hat[4])
out_effect[out_effect$AO115 == 1, "AO115"] <- as.numeric(omega_hat[6])
#out_effect[out_effect$LS81 == 1, "LS81"] <- as.numeric(omega_hat[2])
#out_effect$TC149 <- out_effect$TC149 * as.numeric(omega_hat[1])
ao1_effect_ts <- ts(out_effect$AO65)
ls1_effect_ts <- ts(out_effect$LS81)
ao2_effect_ts <- ts(out_effect$AO101)
ao3_effect_ts <- ts(out_effect$AO115)
#ls_effect_ts <- ts(out_effect$LS81)
#
# #Substraindo o efeito da intervenção.
df_train_clean <- df_train - ao1_effect_ts # extract aditive effect
df_train_clean <- df_train_clean - ls1_effect_ts # extract level shift effect
df_train_clean <- df_train_clean - ao2_effect_ts # extract aditive effect
df_train_clean <- df_train_clean - ao3_effect_ts # extract aditive effect
#df_train_clean <- df_train_clean - ao2_effect_ts # extract aditive effect
#ts_effect_ts[171] <- 0 # This observation has been treated as aditive
#df_train_clean <- df_train_clean - ls_effect_ts # extract transient effect
#
plot(cbind("Original" = df_train,
           "Without outliers" = df_train_clean,
           "Additive effect 1" = ao1_effect_ts,
           "Level Shift effect 1" = ls1_effect_ts,
           #"Additive effect 2" = ao1_effect_ts,
           #"Additive effect 3" = ao1_effect_ts,
           #"Additive effect" = ao2_effect_ts,
           #"Transient change effect" = ls_effect_ts),
           main = "Time series and outlier effects"))
# # #
# # #
# # # plot.ts(df_train_clean)
df_train <- df_train_clean
#-------------------------------------------------------------------------------

# Required functions
source("barma.r")
source("barma.fit.r")


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

# # All coef combinations up to order 6
list_combinations <- list(
  "c(NA)" = c(NA),
  "c(1)" = c(1),
  "c(2)" = c(2),
  "c(3)" = c(3),
  "c(4)" = c(4),
  "c(5)" = c(5),
  "c(6)" = c(6),
  "c(1, 2)" = c(1, 2),
  "c(1, 3)" = c(1, 3),
  "c(1, 4)" = c(1, 4),
  "c(1, 5)" = c(1, 5),
  "c(1, 6)" = c(1, 6),
  "c(2, 3)" = c(2, 3),
  "c(2, 4)" = c(2, 4),
  "c(2, 5)" = c(2, 5),
  "c(2, 6)" = c(2, 6),
  "c(3, 4)" = c(3, 4),
  "c(3, 5)" = c(3, 5),
  "c(3, 6)" = c(3, 6),
  "c(4, 5)" = c(4, 5),
  "c(4, 6)" = c(4, 6),
  "c(5, 6)" = c(5, 6),
  "c(1, 2, 3)" = c(1, 2, 3),
  "c(1, 2, 4)" = c(1, 2, 4),
  "c(1, 2, 5)" = c(1, 2, 5),
  "c(1, 2, 6)" = c(1, 2, 6),
  "c(1, 3, 4)" = c(1, 3, 4),
  "c(1, 3, 5)" = c(1, 3, 5),
  "c(1, 3, 6)" = c(1, 3, 6),
  "c(1, 4, 5)" = c(1, 4, 5),
  "c(1, 4, 6)" = c(1, 4, 6),
  "c(1, 5, 6)" = c(1, 5, 6),
  "c(2, 3, 4)" = c(2, 3, 4),
  "c(2, 3, 5)" = c(2, 3, 5),
  "c(2, 3, 6)" = c(2, 3, 6),
  "c(2, 4, 5)" = c(2, 4, 5),
  "c(2, 4, 6)" = c(2, 4, 6),
  "c(2, 5, 6)" = c(2, 5, 6),
  "c(3, 4, 5)" = c(3, 4, 5),
  "c(3, 4, 6)" = c(3, 4, 6),
  "c(3, 5, 6)" = c(3, 5, 6),
  "c(4, 5, 6)" = c(4, 5, 6),
  "c(1, 2, 3, 4)" = c(1, 2, 3, 4),
  "c(1, 2, 3, 5)" = c(1, 2, 3, 5),
  "c(1, 2, 3, 6)" = c(1, 2, 3, 6),
  "c(1, 2, 4, 5)" = c(1, 2, 4, 5),
  "c(1, 2, 4, 6)" = c(1, 2, 4, 6),
  "c(1, 2, 5, 6)" = c(1, 2, 5, 6),
  "c(1, 3, 4, 5)" = c(1, 3, 4, 5),
  "c(1, 3, 4, 6)" = c(1, 3, 4, 6),
  "c(1, 3, 5, 6)" = c(1, 3, 5, 6),
  "c(1, 4, 5, 6)" = c(1, 4, 5, 6),
  "c(2, 3, 4, 5)" = c(2, 3, 4, 5),
  "c(2, 3, 4, 6)" = c(2, 3, 4, 6),
  "c(2, 3, 5, 6)" = c(2, 3, 5, 6),
  "c(2, 4, 5, 6)" = c(2, 4, 5, 6),
  "c(3, 4, 5, 6)" = c(3, 4, 5, 6),
  "c(1, 2, 3, 4, 5)" = c(1, 2, 3, 4, 5),
  "c(1, 2, 3, 4, 6)" = c(1, 2, 3, 4, 6),
  "c(1, 2, 3, 5, 6)" = c(1, 2, 3, 5, 6),
  "c(1, 2, 4, 5, 6)" = c(1, 2, 4, 5, 6),
  "c(1, 3, 4, 5, 6)" = c(1, 3, 4, 5, 6),
  "c(2, 3, 4, 5, 6)" = c(2, 3, 4, 5, 6),
  "c(1, 2, 3, 4, 5, 6)" = c(1, 2, 3, 4, 5, 6)
  )

# All coef combinations up to order 5
# list_combinations <- list(
#   "c(NA)" = c(NA),
#   "c(1)" = c(1),
#   "c(2)" = c(2),
#   "c(3)" = c(3),
#   "c(4)" = c(4),
#   "c(5)" = c(5),
#   "c(1, 2)" = c(1, 2),
#   "c(1, 3)" = c(1, 3),
#   "c(1, 4)" = c(1, 4),
#   "c(1, 5)" = c(1, 5),
#   "c(2, 3)" = c(2, 3),
#   "c(2, 4)" = c(2, 4),
#   "c(2, 5)" = c(2, 5),
#   "c(3, 4)" = c(3, 4),
#   "c(3, 5)" = c(3, 5),
#   "c(4, 5)" = c(4, 5),
#   "c(1, 2, 3)" = c(1, 2, 3),
#   "c(1, 2, 4)" = c(1, 2, 4),
#   "c(1, 2, 5)" = c(1, 2, 5),
#   "c(1, 3, 4)" = c(1, 3, 4),
#   "c(1, 3, 5)" = c(1, 3, 5),
#   "c(1, 4, 5)" = c(1, 4, 5),
#   "c(2, 3, 4)" = c(2, 3, 4),
#   "c(2, 3, 5)" = c(2, 3, 5),
#   "c(2, 4, 5)" = c(2, 4, 5),
#   "c(3, 4, 5)" = c(3, 4, 5),
#   "c(1, 2, 3, 4)" = c(1, 2, 3, 4),
#   "c(1, 2, 3, 5)" = c(1, 2, 3, 5),
#   "c(1, 2, 4, 5)" = c(1, 2, 4, 5),
#   "c(1, 3, 4, 5)" = c(1, 3, 4, 5),
#   "c(2, 3, 4, 5)" = c(2, 3, 4, 5),
#   "c(1, 2, 3, 4, 5)" = c(1, 2, 3, 4, 5)
# )

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
        resid = 1,
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

valid_models_ar_coef # 78 of 4096 
valid_models_ma_coef # 78 of 4096
valid_models_aic # min: AIC -146.8761 | model selected: AIC -142.48082
valid_models_bic # min: BIC -121.6399 | ~

# BEST MODEL BY AIC
which.min(valid_models_aic) # position of min in the vector
valid_models_aic[[which.min(valid_models_aic)]] # AIC 
best_arma_combination_ar_aic <- 
  valid_models_ar_coef[[which.min(valid_models_aic)]] # best
best_arma_combination_ma_aic <- 
  valid_models_ma_coef[[which.min(valid_models_aic)]] # best

fit_barma_best_aic <- barma(
  y,
  ar = best_arma_combination_ar_aic, # 1, 2, 5, 6
  ma = best_arma_combination_ma_aic, # 3, 4, 5
  h = h,
  diag = 1,
  resid = 3,
  link = "logit"
)

cpgram(fit_barma_best_aic$resid3, main = "Cumulative Periodogram of Residuals")
car::qqPlot(fit_barma_best_aic$resid3)
tseries::jarque.bera.test(fit_barma_best_aic$resid3)
Box.test(fit_barma_best_aic$resid3, type = "Box-Pierce")
Box.test(fit_barma_best_aic$resid3, type = "Ljung-Box")
forecast::checkresiduals(fit_barma_best_aic$resid3, test = "LB")
acf(fit_barma_best_aic$resid3)

# BEST MODEL BY BIC
which.min(valid_models_bic) # position of min in the vector
valid_models_bic[[which.min(valid_models_bic)]] # BIC
best_arma_combination_ar_bic <- 
  valid_models_ar_coef[[which.min(valid_models_bic)]] # best
best_arma_combination_ma_bic <- 
  valid_models_ma_coef[[which.min(valid_models_bic)]] # best

fit_barma_best_bic <- barma(
  y,
  ar = best_arma_combination_ar_bic, # 1, 2, 5, 6
  ma = best_arma_combination_ma_bic, # 3, 4, 5
  h = h,
  diag = 1,
  resid = 3,
  link = "logit"
)

cpgram(fit_barma_best_bic$resid3, main = "Cumulative Periodogram of Residuals")
car::qqPlot(fit_barma_best_bic$resid3)
tseries::jarque.bera.test(fit_barma_best_bic$resid3)
Box.test(fit_barma_best_bic$resid3, type = "Box-Pierce")
Box.test(fit_barma_best_bic$resid3, type = "Ljung-Box")
forecast::checkresiduals(fit_barma_best_bic$resid3, test = "LB")
acf(fit_barma_best_bic$resid3)
