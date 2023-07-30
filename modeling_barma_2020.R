# rm(list = ls()) 
library("e1071") # To calculate asymetry and kustosis.
library("tseries")
library("tidyverse")

# Load the data
df <- readRDS("maternal_mortality_rate_2020.rds")

# Transform the data to time series object (day frequency)
df_ts <- ts(df$rate, start = c(2020, 1), frequency = 366)

h <- 7 # forecast window

# Split data in train and test
df_train <- df_ts[1:(length(df_ts) - h)]
df_train <- ts(df_train)

df_test <- df_ts[(length(df_train) + 1):(length(df_ts))]
df_test <- ts(df_test)

# Data description.
summary(df_train)  # Resume
var(df_train)      # Variance.
skewness(df_train) # Skewness.
kurtosis(df_train) # Kurtosis.

# Apply difference transformation to make data stationary
df_train_diff <- diff(df_train)

# Check out autocorrelation and partial autocorrelation functions
acf(df_train) # Before difference transformation
acf(df_train_diff) # After difference transformation
pacf(df_train) # Before difference transformation
pacf(df_train_diff) # After difference transformation

# Perform ADF and KPSS tests (indicates stationary data for diff data)
adf.test(df_train) # Not stationary
kpss.test(df_train) # Not stationary
adf.test(df_train_diff) # Stationary data
kpss.test(df_train_diff) # Stationary data

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

best_aic <- Inf
best_fit <- NULL
max_order = 3 # check up to order 3

#-------------------------------------------------------------------------------
# Select model considering every AR and MA combination (better, but take longer)
#-------------------------------------------------------------------------------
# Generate all combinations of values from 1 to max_order without repetition
ar_combinations <- lapply(1:max_order, function(x) combn(1:max_order, x, simplify = FALSE))
ar_combinations <- unlist(ar_combinations, recursive = FALSE)

ma_combinations <- lapply(1:max_order, function(x) combn(1:max_order, x, simplify = FALSE))
ma_combinations <- unlist(ma_combinations, recursive = FALSE)

# Select the model with the orders combination that gives the lower AIC
for (p in ar_combinations) {
  for (q in ma_combinations) {
    fit <- suppressWarnings(
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
      barma(y, ar = p, ma = q, h = h, diag = 0, resid = 1, link = "logit")
    )
    
    if (fit$conv != 0){
      # Skip this iteration if optimization algorithm do not converge
      next
    }
    else if (any(is.na(fit$forecast))){
      # Skip this iteration if exists missing values among the forecasts
      next
    }
    else if (!(all(abs(fit$phi) < 0.9) & all(abs(fit$theta) < 0.9))){
      # Skip this iteration if any estimates out of the interval (-0.9, 0.9)
      next
    }
    else if (!(all(fit$pvalues < 0.05))){
      # Skip this iteration if exists non significant parameters at 5%
      next
    }
    else if (fit$aic < best_aic){
      best_aic <- fit$aic
      best_fit <- fit
      barma_model_order <- list(ar = p, ma = q)
    }
  }
}

# #-------------------------------------------------------------------------------
# # Select model considering only one phi and one theta (worst, but faster)
# #-------------------------------------------------------------------------------
# # Generate all combinations of AR and MA orders up to max_order
# orders <- expand.grid(ar = 1:max_order, ma = 1:max_order)
# 
# # Select the model with the orders combination that gives the lower AIC
# for (i in 1:nrow(orders)) {
#   p <- orders$ar[i]
#   q <- orders$ma[i]
# 
#   fit <- suppressWarnings(
#     # Options for barma():
#     # resid = 1 : standardized residual
#     # resid = 2 : standardized residual 2
#     # resid = 3 : standardized weighted residual
#     # resid = 4 : deviance residual
#     # diag = 0 : do not plot any graph (useful for simulations)
#     # diag = 1 : plot graphs
#     # diag = 2 : save graphs on ps files
#     # link = "logit"
#     # link = "probit"
#     # link = "cloglog"
#     # h is the prediction window
#     barma(y, ar = p, ma = q, h = h, diag = 0, resid = 1, link = "logit")
#   )
# 
#   if (fit$conv != 0){
#     # Skip this iteration if optimization algorithm do not converge
#     next
#   }
# 
#   else if (any(is.na(fit$forecast))){
#     # Skip this iteration if exists missing values among the forecasts
#     next
#   }
#   
#   else if (!(all(abs(fit$phi) < 0.9) && all(abs(fit$theta) < 0.9))){
#     # Skip this iteration if parameters are not in the desired range (-1, 1)
#     next
#   }
#   
#   else if (fit$aic < best_aic){
#     best_aic <- fit$aic
#     best_fit <- fit
#     barma_model_order <- list(ar = p, ma = q)
#   }
# }


# best_fit <- barma(y, ar = c(3), ma = c(1), h = h,
#                   diag = 0, resid = 1, link = "logit")
# barma_model_order <- list(ar = c(3), ma = c(1))

barma_forecast <- best_fit$forecast
barma_fitted <- best_fit$fitted
res <- best_fit$resid1 

#################### Undo transformations and differences ####################
if (extremes){
  barma_forecast <- ((barma_forecast*n_data)-0.5)/(n_data-1)
  barma_fitted <- ((barma_fitted*n_data)-0.5)/(n_data-1)
}

if (d > 0){
  barma_forecast <- barma_forecast*(b-a)+a
  barma_forecast <- diffinv(barma_forecast, 
                            differences = d,
                            xi = df_train[length(df_train)])
  barma_forecast <- barma_forecast[-1]
  barma_fitted <- barma_fitted*(b-a)+a
}

# Fitted has max(ar,ma) NAs at the beginning
barma_nas = max(max(barma_model_order$ar), max(barma_model_order$ma))
if (any(is.na(barma_model_order$ar))) barma_nas = max(barma_model_order$ma)
if (any(is.na(barma_model_order$ma))) barma_nas = max(barma_model_order$ar)

# Ignore the NA's as the beginning
barma_fitted <- barma_fitted[(barma_nas+1):length(barma_fitted)]
if (d > 0) {
  barma_fitted <- df_train[
    (barma_nas+1):(barma_nas+length(barma_fitted))
  ]+barma_fitted
}

date_start_forecast = end(df_ts)
date_start_forecast[2] = date_start_forecast[2] - (h-1)

ts_barma_forecast <- ts(barma_forecast, 
                        start = date_start_forecast,
                        frequency = 366)
ts_barma_fitted <- ts(barma_fitted,
                      start = c(2020, barma_nas+2), 
                      frequency = 366)

################################################################################
######################## Diagnostics and visualizations ########################
################################################################################

# Create a data frame with the observed and predicted values
# Maternal Mortality Ratio (deaths per 100,000 live births)
df_plot <- data.frame(
  date = time(df_train[1:length(barma_fitted)]),
  observed = as.numeric(df_train[1:length(barma_fitted)])*100000,
  predicted = as.numeric(barma_fitted)*100000
)

# Add a column with the date in Date format
df_plot$date_formatted <- as.Date(paste0("2020-", df_plot$date),
                                  format = "%Y-%j")

#-------------- Save high resolution png plot ---------------
png("In_sample_plot_2020.png", units = "in", width = 8, height = 4, res = 300)

# In-sample plot
ggplot(df_plot, aes(x = date_formatted)) +
  geom_line(aes(y = observed), 
            color = "black",
            size = 0.8) +
  geom_line(aes(y = predicted), 
            color = "red", 
            linetype = "twodash",
            size = 1.0) +
  labs(x = "\nData",
       y = "Razão de mortalidade materna 2020\n (mortes por 100.000 nascidos vivos)\n",
       title = "Previsão intra-amostra") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank()) +
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b %d") +
  scale_y_continuous(labels = scales::comma) +
  annotate("text", 
           x = max(df_plot$date_formatted) - 50, 
           y = max(df_plot$observed) + 25, 
           label = "Observado",
           hjust = 0,
           vjust = 1, 
           color = "black") +
  annotate("text", 
           x = max(df_plot$date_formatted) - 50,
           y = max(df_plot$observed),
           label = "Previsão", 
           hjust = 0, 
           vjust = 1,
           color = "red")
dev.off()

#-------------------------- Fit ARIMA ----------------------------

# best_fit_arima <- auto.arima(df_train)
best_fit_arima <- stats::arima(df_train,
                               order = c(2, 1, 0)) 
lmtest::coeftest(best_fit_arima)

# Predict using ARIMA model
arima_forecast <- forecast(best_fit_arima, h = 7)$mean
arima_forecast <- as.numeric(arima_forecast)

# Jarque–Bera test
# H0: joint hypothesis of skewness=0 and excess kurtosis=0
tseries::jarque.bera.test(residuals(best_fit_arima))

############################ Portmanteau tests #################################
#----------------------------- Kwan and Sim -----------------------------------
#------------------  Q4 test (Goodness-of-fit article) ------------------------

# Lag m vector
vm <- c(seq(from = 3, to = 30, by = 1))
f <- 2 #degrees of freedom
#res_aux <-  res # βARMA residuals
#res <-  res_aux # βARMA residuals
#res <- best_fit_arima$residuals # ARIMA residuals

# Initializing response vector
QKW4 <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
Q4 <- matrix(rep(NA,30), nrow = 1, ncol = 30) 

# P-values
source("Kwan_Chest.R")

for(m in vm)
{
  # Kwan and Sim test(4) 
  QKW4[,m] <- Kwan.sim.chest(res, 
                             lag = m, 
                             fitdf = f, 
                             type = "correlation",
                             test = 4)$p.value 
  # New portmanteau test statistic (Goodness-of-fit article)
  Q4[,m] <- Kwan.sim.chest(res, 
                           lag = m, 
                           fitdf = f,
                           type = "partial", 
                           test = 4)$p.value         
}

# Graphics for portmanteau tests (Are showing that the model is suitable)
nlag <- 30

# Kwan-Sim (QKW4)
plot(1L:nlag, 
     QKW4, 
     xlab = "lag",
     ylab = "p-values",
     las = 1,
     ylim = c(0,1),
     lwd = 1.5)
abline(h = 0.05, lty = 2, col = "blue", lwd = 1.5)

# New portmanteau test statistic (Q4) 
plot(1L:nlag,
     Q4,
     xlab = "lag",
     ylab = "p-values",
     las = 1,
     ylim = c(0,1),
     lwd = 1.5)
abline(h = 0.05, lty = 2, col = "blue", lwd = 1.5)

#-------------------------- Model comparison ----------------------------

# Initializing the result matrix
table_metrics <- matrix(rep(NA, 42), nrow = 6, ncol = 7)                           

# Fill table with MAE, RMSE, and sMAPE metrics
# ARIMA(1, 1, 1) and βARMA(3, 1) com phi1 e phi2 fixos
for(i in 1:7)
{
  # MAE
  table_metrics[1,i] <- round(
    mean(abs((df_test[1:i] - barma_forecast[1:i])))*100,
    4
  )
  table_metrics[2,i] <- round(
    mean(abs((df_test[1:i] - arima_forecast[1:i])))*100,
    4
  )
  
  # RMSE
  table_metrics[3,i] <- round(
    sqrt(mean((df_test[1:i] - barma_forecast[1:i])^2))*100,
    4
  )
  table_metrics[4,i] <- round(
    sqrt(mean((df_test[1:i] - arima_forecast[1:i])^2))*100,
    4
  )
  
  # sMAPE
  table_metrics[5,i] <- round(
    mean(2 * abs(df_test[1:i] - barma_forecast[1:i])/
           (abs(df_test[1:i]) + abs(barma_forecast[1:i]))),
    4
  )
  table_metrics[6,i] <- round(
    mean(2 * abs(df_test[1:i] - arima_forecast[1:i])/
           (abs(df_test[1:i]) + abs(arima_forecast[1:i]))),
    4
  )
}

# Set the row names to indicate the metric and method
rownames(table_metrics) <- c("MAE(BARMA)x100", 
                             "MAE(ARIMA)x100", 
                             "RMSE(BARMA)",
                             "RMSE(ARIMA)",
                             "sMAPE(BARMA)", 
                             "sMAPE(ARIMA)")
colnames(table_metrics) <- c("h=1", 
                             "h=2", 
                             "h=3", 
                             "h=4", 
                             "h=5", 
                             "h=6", 
                             "h=7")

table_metrics
