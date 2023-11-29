# rm(list = ls()) 
detach("package:e1071", unload = TRUE)
detach("package:tseries", unload = TRUE)
detach("package:tidyverse", unload = TRUE)
detach("package:itsmr", unload = TRUE)
detach("package:tsoutliers", unload = TRUE)

library("e1071") # To calculate asymetry and kustosis.
library("tseries")
library("tidyverse")
library("itsmr")
library("tsoutliers")

# Load the data
df <- readRDS("maternal_mortality_rate_2021.rds")
df <- df[1:135,] # Dados fim de abril
#df <- df[1:135,] # Dados até 15 de maio
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
# outliers_excess_ts <- tsoutliers::tso(df_train)
# plot(outliers_excess_ts) #Gráfico de detecção.
# outliers_excess_ts$outliers #Informações gerais dos pontos encontrados.
# outliers_idx <- outliers_excess_ts$outliers$ind #Posição dos outliers na série.
# # plot(forecast::tsclean(df_train))
# # plot(df_train)
# #length of our time series
# n <- length(df_train)
# 
# (col_int <- outliers(c("AO", "LS", "AO", "AO"), outliers_idx[c(1, 3, 4, 6)])) #Colunas de interesse.
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
# colnames(out_effect) <- c("AO65", "LS81", "AO101", "AO115")
# out_effect[out_effect$AO65 == 1, "AO65"] <- as.numeric(omega_hat[1])
# out_effect[out_effect$LS81 == 1, "LS81"] <- as.numeric(omega_hat[3])
# out_effect[out_effect$AO101 == 1, "AO101"] <- as.numeric(omega_hat[4])
# out_effect[out_effect$AO115 == 1, "AO115"] <- as.numeric(omega_hat[6])
# #out_effect[out_effect$LS81 == 1, "LS81"] <- as.numeric(omega_hat[2])
# #out_effect$TC149 <- out_effect$TC149 * as.numeric(omega_hat[1])
# ao1_effect_ts <- ts(out_effect$AO65)
# ls1_effect_ts <- ts(out_effect$LS81)
# ao2_effect_ts <- ts(out_effect$AO101)
# ao3_effect_ts <- ts(out_effect$AO115)
# #ls_effect_ts <- ts(out_effect$LS81)
# #
# # #Substraindo o efeito da intervenção.
# df_train_clean <- df_train - ao1_effect_ts # extract aditive effect
# df_train_clean <- df_train_clean - ls1_effect_ts # extract level shift effect
# df_train_clean <- df_train_clean - ao2_effect_ts # extract aditive effect
# df_train_clean <- df_train_clean - ao3_effect_ts # extract aditive effect
# #df_train_clean <- df_train_clean - ao2_effect_ts # extract aditive effect
# #ts_effect_ts[171] <- 0 # This observation has been treated as aditive
# #df_train_clean <- df_train_clean - ls_effect_ts # extract transient effect
# #
# plot(cbind("Original" = df_train,
#            "Without outliers" = df_train_clean,
#            "Additive effect 1" = ao1_effect_ts,
#            "Level Shift effect 1" = ls1_effect_ts,
#            #"Additive effect 2" = ao1_effect_ts,
#            #"Additive effect 3" = ao1_effect_ts,
#            #"Additive effect" = ao2_effect_ts,
#            #"Transient change effect" = ls_effect_ts),
#      main = "Time series and outlier effects"))
# # # #
# # # #
# # # # plot.ts(df_train_clean)
# df_train <- df_train_clean
################################################################################

# Data description.
summary(df_train)  # Resume
var(df_train)      # Variance.
skewness(df_train) # Skewness.
kurtosis(df_train) # Kurtosis.

# Apply difference transformation to make data stationary
df_train_diff <- diff(df_train)
plotc(df_train, df_train_diff) # diff = 1 vs diff = 0

# Longmemory tests (d = 0.3527578 seems to highlight a period of 7)
# d.GPH <- fracdiff::fdGPH(df_train,bandw.exp = 0.7)
# for(i in 1:10)
#   print(cbind(beta=i*0.1,GPH=fracdiff::fdGPH(df_train,bandw.exp = i*0.1)$d))
# df_train_diff <- fracdiff::diffseries(df_train,fracdiff::fdGPH(df_train,bandw.exp = 0.7)$d)


# Check out autocorrelation and partial autocorrelation functions
# BEFORE difference transformation
par(mfrow=c(2,2))
plot(df_train,
     xlab = "Time",
     ylab = "Ratio of deaths", 
     pch = 18)
TSA::periodogram(df_train,
                 xlab = "Frequency",
                 ylab = "Periodogram")

acf(df_train,
    xlab = "Lag (days)",
    ylab = "ACF",
    main = "",
    lag.max = 70, 
    pch = 18)
pacf(df_train,
     xlab = "Lag (days)",
     ylab = "PACF",
     main = "",
     lag.max =70, 
     pch = 18)

# Check out autocorrelation and partial autocorrelation functions
# AFTER difference transformation
par(mfrow=c(2,2))
plot.ts(df_train_diff,
     xlab = "Time",
     ylab = "Ratio of deaths", 
     pch = 18)
TSA::periodogram(df_train_diff,
                 xlab = "Frequency",
                 ylab = "Periodogram")

acf(df_train_diff,
    xlab = "Lag (days)",
    ylab = "ACF",
    main = "",
    lag.max = 25, 
    pch = 18,
    ci.type = "ma")
acf(df_train_diff,
     xlab = "Lag (days)",
     ylab = "PACF",
     main = "",
     lag.max = 25, 
     pch = 18,
     type = "partial")

# sigma^2 is lower with diff = 1
forecast::Arima(df_train, order = c(0, 0, 0)) # non-stationary
forecast::Arima(df_train, order = c(0, 1, 0)) # diff optmized
forecast::Arima(df_train, order = c(0, 2, 0)) # overdiferenciated


# Perform ADF and KPSS tests (indicates stationary data for diff data)
adf.test(df_train, alternative = "stationary") # Not stationary
kpss.test(df_train) # Not stationary
adf.test(df_train_diff, alternative = "stationary") # Stationary data
kpss.test(df_train_diff) # Stationary data

#-------------------------------------------------------------------------------
# Required functions
source("barma.r")
source("barma.fit.r")


d <- 1 # No difference transformation
#d <- 1 # One difference transformation

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

#--------------------------------- Fit BARMA -----------------------------------
# BEST MODEL BARMA AIC 
fit_barma_best <- barma( 
  y,
  ar = c(1, 2, 5, 6),
  ma = c(3, 4, 5),
  h = h,
  diag = 1,
  resid = 1,
  link = "logit" 
)
fit_barma_best_order <- list(
  ar = c(1, 2, 5, 6),
  ma = c(3, 4, 5)
)

# # MODEL 2 BARMA AIC (BIC selected the same)
# fit_barma_best <- barma(
#   y,
#   ar = c(1, 2, 5),
#   ma = c(1, 3, 5, 6),
#   h = h,
#   diag = 1,
#   resid = 3,
#   link = "logit"
# )
# fit_barma_best_order <- list(
#   ar = c(1, 2, 5),
#   ma = c(1, 3, 5, 6)
# )

# # BARMAX
# fit_barma_best <- barma( 
#   y,
#   ar = c(4, 5),
#   ma = c(1, 3, 4, 5, 6),
#   h = h,
#   diag = 1,
#   resid = 1,
#   link = "logit",
#   X = as.matrix(ao1_effect_ts),
#   X_hat = as.matrix(rep(0, length(df_test))) 
# )
# fit_barma_best_order <- list(
#   ar = c(4, 5),
#   ma = c(1, 3, 4, 5, 6)
# )

# test <- if_else(ao1_effect_ts == 0, 0, 1)
# fit_barma_best <- barma( 
#   y,
#   ar = c(1, 2, 5, 6),
#   ma = c(3, 4, 5),
#   h = h,
#   diag = 1,
#   resid = 1,
#   link = "logit",
#   X = as.matrix(test),
#   X_hat = as.matrix(rep(0, length(df_test))) 
# )
# fit_barma_best_order <- list(
#   ar = c(1, 2, 5, 6),
#   ma = c(3, 4, 5)
# )

best_fit <- fit_barma_best
barma_model_order <- fit_barma_best_order
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

# Outlier
#barma_fitted <- barma_fitted + ao1_effect_ts[1:length(barma_fitted)]

date_start_forecast = end(df_ts)
date_start_forecast[2] = date_start_forecast[2] - (h-1)

ts_barma_forecast <- ts(barma_forecast, 
                        start = date_start_forecast,
                        frequency = 365)
ts_barma_fitted <- ts(barma_fitted,
                      start = c(2021, barma_nas+2), 
                      frequency = 365)

#-------------------------- Fit ARIMA ----------------------------
# BEST MODEL BY AIC: -1516.579
fit_arima <- forecast::Arima( 
  df_train,
  order = c(6, 1, 6),
  fixed = c(NA, 0, 0, 0, 0, NA, # ar
            0, NA, 0, 0, 0, NA) # ma
)

best_fit_arima <- fit_arima

# # BEST MODEL BY BIC: -1504.638 
# fit_arima <- forecast::Arima(
#   df_train,
#   order = c(6, 1, 6),
#   fixed = c(0, 0, 0, 0, NA, 0, # ar
#             NA, 0, 0, 0, NA, 0) # ma
# )
# 
# best_fit_arima <- fit_arima
# 
# # MODEL 2 BY AIC
# fit_arima <- forecast::Arima(
#   df_train,
#   order = c(6, 1, 6),
#   fixed = c(NA, NA, NA, 0, 0, NA, # ar
#             0, NA, 0, NA, 0, NA) # ma
# )
# 
# best_fit_arima <- fit_arima
# 
# # MODEL 2 BY BIC
# fit_arima <- forecast::Arima(
#   df_train,
#   order = c(6, 1, 6),
#   fixed = c(0, 0, 0, 0, NA, 0, # ar
#             NA, 0, 0, 0, NA, 0) # ma
# )
# 
# best_fit_arima <- fit_arima
 

# #BEST MODEL BY auto.arima
# #Seems to perform better, but in reality, after some investigation, this actually
# #produces fake good predictions. In the sense that it beats other models
# #sometimes, but just because the prediction line is flat all along, not getting
# #much variability.
# fit_arima_auto <- forecast::auto.arima(
#   df_train,
#   max.p = 6,
#   max.q = 6,
#   max.order = 6,
#   max.d = 1,
#   nmodels = 4096,
#   stepwise = FALSE,
#   allowdrift = TRUE,
#   lambda=FALSE
# )
# best_fit_arima <- fit_arima_auto

lmtest::coeftest(best_fit_arima) 
forecast::checkresiduals(best_fit_arima, test = "LB")
Box.test(residuals(best_fit_arima), type = "Box-Pierce")
Box.test(residuals(best_fit_arima), type = "Ljung-Box")

# Jarque–Bera test
# H0: joint hypothesis of skewness=0 and excess kurtosis=0
tseries::jarque.bera.test(residuals(best_fit_arima))

plotc(df_train,best_fit_arima$fitted)
forecast::accuracy(best_fit_arima$fitted, df_train)
tsdiag(best_fit_arima, gof.lag = 25)
test(best_fit_arima$residuals)
car::qqPlot(best_fit_arima$res)

# Forecasts
arima_forecast <- forecast::forecast(best_fit_arima, h = h)$mean
arima_forecast <- as.numeric(arima_forecast)

################################################################################
######################## Diagnostics and visualizations ########################
################################################################################
# Fitted
# Create a data frame with the observed and predicted values
# Maternal Mortality Ratio (deaths per 100,000 live births)
#BARMA
df_plot <- data.frame(
  date = time(df_train[1:length(barma_fitted)]),
  observed = as.numeric(df_train[1:length(barma_fitted)])*100000,
  predicted = as.numeric(barma_fitted)*100000
)

# df_plot2 <- data.frame(
#   date = time(df_train[1:length(barma_fitted)]),
#   observed = as.numeric(df_train[1:length(barma_fitted)]),
#   predicted = as.numeric(barma_fitted)
# )

# Add a column with the date in Date format
df_plot$date_formatted <- as.Date(paste0("2021-", df_plot$date),
                                  format = "%Y-%j")

# ARIMA
df_plot_arma <- data.frame(
  date = time(df_train[1:length(best_fit_arima$fitted)]),
  observed = as.numeric(df_train[1:length(best_fit_arima$fitted)])*100000,
  predicted = as.numeric(best_fit_arima$fitted)*100000
)

df_plot$predicted2 <- df_plot_arma[1:nrow(df_plot), "predicted"]

#-------------- Save high resolution png plot ---------------
#png("new_In_sample_plot_2021.png", units = "in", width = 8, height = 4, res = 300)

# In-sample plot
insample_plot1 <- 
  ggplot(df_plot, aes(x = date_formatted)) +
  geom_line(aes(y = observed), 
            color = "#646464",
            size = 0.8) +
  geom_point(aes(y = observed), 
             color = "#646464") +
  geom_line(aes(y = predicted), 
            color = "#E40808", 
            linetype = "solid",
            size = 0.9) +
  geom_point(aes(y = predicted), 
             color = "#E40808") +
  labs(x = "\nData (dia/mês)",
       y = "Razão de mortalidade materna em 2021\n (mortes por 100.000 nascidos vivos)\n",
       title = "Previsão intra-amostra") +
  theme_bw() +  
  theme(plot.title = element_text(color = "black", 
                                  face = "bold"),  
        axis.title = element_text(color = "black", 
                                  face = "bold"),  
        axis.text = element_text(color = "black", 
                                 face = "bold"),  
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%d/%m",
               limits = c(min(df_plot$date_formatted), max(df_plot$date_formatted))) +
  scale_y_continuous(labels = scales::comma) +
  coord_cartesian(xlim = c(min(df_plot$date_formatted), max(df_plot$date_formatted)),
                  ylim = c(min(df_plot$observed), max(df_plot$observed))) +
  annotate("text", 
           x = max(df_plot$date_formatted) - 112, 
           y = max(df_plot$observed) + 0, 
           label = "Observado",
           hjust = 0,
           vjust = 1, 
           color = "#646464") +
  annotate("text", 
           x = max(df_plot$date_formatted) - 112,
           y = max(df_plot$observed) -25,
           label = "Previsão βARMA", 
           hjust = 0, 
           vjust = 1,
           color = "#E40808") 

insample_plot2 <- 
  ggplot(df_plot, aes(x = date_formatted)) +
  geom_line(aes(y = observed), 
            color = "#646464",
            size = 0.8) +
  geom_point(aes(y = observed), 
             color = "#646464") +
  geom_line(aes(y = predicted2), 
            color = "blue", 
            linetype = "solid",
            size = 0.9) +
  geom_point(aes(y = predicted2), 
             color = "blue") +
  labs(x = "\nData (dia/mês)",
       y = "Razão de mortalidade materna em 2021\n (mortes por 100.000 nascidos vivos)\n",
       title = "Previsão intra-amostra") +
  theme_bw() +  
  theme(plot.title = element_text(color = "black", 
                                  face = "bold"),  
        axis.title = element_text(color = "black", 
                                  face = "bold"),  
        axis.text = element_text(color = "black", 
                                 face = "bold"),  
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%d/%m",
               limits = c(min(df_plot$date_formatted), max(df_plot$date_formatted))) +
  scale_y_continuous(labels = scales::comma) +
  coord_cartesian(xlim = c(min(df_plot$date_formatted), max(df_plot$date_formatted)),
                  ylim = c(min(df_plot$observed), max(df_plot$observed))) +
  annotate("text", 
           x = max(df_plot$date_formatted) - 112, 
           y = max(df_plot$observed) + 0, 
           label = "Observado",
           hjust = 0,
           vjust = 1, 
           color = "#646464") +
  annotate("text", 
           x = max(df_plot$date_formatted) - 112,
           y = max(df_plot$observed) -25,
           label = "Previsão ARMA", 
           hjust = 0, 
           vjust = 1,
           color = "darkblue")
#dev.off()


# BARMA Residuals
df_resid <- data.frame(Index = 1:length(res), Residuals = res)
# ARMA residuals
df_resid2 <- data.frame(Index = 1:length(fit_arima$residuals), 
                        Residuals = as.numeric(fit_arima$residuals))
tseries::jarque.bera.test(df_resid$Residuals)

resid_pondera1 <- 
ggplot(df_resid, aes(x = Index, y = Residuals)) +
  geom_point(shape = 4) +  
  geom_hline(yintercept = -3, 
             linetype = "dashed", 
             color = "darkblue", 
             size = 1) +  
  geom_hline(yintercept = 3, 
             linetype = "dashed", 
             color = "darkblue", 
             size = 1) +  
  geom_hline(yintercept = -2, 
             linetype = "dotted", 
             color = "lightblue", 
             size = 1) +  
  geom_hline(yintercept = 2, 
             linetype = "dotted", 
             color = "lightblue", 
             size = 1) +  
  coord_cartesian(ylim = c(-5, 5)) +  
  labs(x = "Índice", 
       y = "Resíduo", 
       title = "Resíduo Padronizado") +
  theme_bw() +  
  theme(plot.title = element_text(color = "black", 
                                  face = "bold"),  
        axis.title = element_text(color = "black", 
                                  face = "bold"),  
        axis.text = element_text(color = "black", 
                                 face = "bold"),  
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1)) 


p_qqplot <- 
ggplot(df_resid, aes(sample = Residuals)) +
  ggplot2::stat_qq(shape = 4, 
                   color = "black") +  
  qqplotr::stat_qq_line(color = "darkblue", 
                        size = 1) +  
  qqplotr::stat_qq_band(bandType = "pointwise", 
                        mapping = aes(fill = "Bootstrap"), 
                        conf = 0.95, 
                        alpha = 0.3) +  
  scale_fill_manual(values = c("Bootstrap" = "#0066cc")) +
  guides(fill = FALSE) +
  labs(x = "Quantis teóricos (normal)", 
       y = "Quantis amostrais", 
       title = "Gráfico Q-Q") +
  theme_bw() +  
  theme(plot.title = element_text(color = "black", 
                                  face = "bold"),
        axis.title = element_text(color = "black", 
                                  face = "bold"), 
        axis.text = element_text(color = "black", 
                                 face = "bold"), 
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))

p_densidade <- 
ggplot(df_resid, aes(x=Residuals)) +
  geom_density(aes(color = "Densidade dos resíduos"), 
               alpha=.2, fill="#0066cc",
               size = 1.1,
               linetype = "solid") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(df_resid$Residuals), 
                            sd = sd(df_resid$Residuals)), 
                aes(color = "Aproximação da distribuição normal"),
                size = 1.1,
                linetype = "dashed") +
  labs(title='Distribuição dos resíduos', 
       x='Resíduos', 
       y='Densidade') +
  theme_bw() +  
  theme(plot.title = element_text(color = "black", 
                                  face = "bold"),
        axis.title = element_text(color = "black", 
                                  face = "bold"), 
        axis.text = element_text(color = "black", 
                                 face = "bold"), 
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1),
        legend.position = c(0.20, 0.85),
        legend.background = element_rect(fill = "transparent"),  
        legend.title = element_blank()) +  
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 0.5)) +
  scale_color_manual(
    values = c("Densidade dos resíduos" = "#0066cc", 
               "Aproximação da distribuição normal" = "darkblue"))

p_densidade_arma <- 
  ggplot(df_resid2, aes(x=Residuals)) +
  geom_density(aes(color = "Densidade dos resíduos (modelo ARMA)"), 
               alpha=.2, fill="#0066cc",
               size = 1.1,
               linetype = "solid") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(df_resid2$Residuals), 
                            sd = sd(df_resid2$Residuals)), 
                aes(color = "Aproximação da distribuição normal"),
                size = 1.1,
                linetype = "dashed") +
  labs(title='Distribuição dos resíduos', 
       x='Resíduos', 
       y='Densidade') +
  theme_bw() +  
  theme(plot.title = element_text(color = "black", 
                                  face = "bold"),
        axis.title = element_text(color = "black", 
                                  face = "bold"), 
        axis.text = element_text(color = "black", 
                                 face = "bold"), 
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1),
        legend.position = c(0.20, 0.85),
        legend.background = element_rect(fill = "transparent"),  
        legend.title = element_blank()) +  
  #coord_cartesian(xlim = c(-3, 3), ylim = c(0, 0.5)) +
  scale_color_manual(
    values = c("Densidade dos resíduos" = "#0066cc", 
               "Aproximação da distribuição normal" = "darkblue"))


acf_funct_barma <- 
forecast::ggAcf(res) +
  labs(title = "Função de autocorrelação dos resíduos") +
  theme_bw() +  
  theme(
    plot.title = element_text(color = "black", face = "bold"),  
    axis.title = element_text(color = "black", face = "bold"),  
    axis.text = element_text(color = "black", face = "bold"),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) 

acf_funct_arma <- 
  forecast::ggAcf(df_resid2$Residuals) +
  labs(title = "Função de autocorrelação dos resíduos (modelo ARMA)") +
  theme_bw() +  
  theme(
    plot.title = element_text(color = "black", face = "bold"),  
    axis.title = element_text(color = "black", face = "bold"),  
    axis.text = element_text(color = "black", face = "bold"),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) 

ggplot(df_resid, aes(x = Index)) +
  geom_line(aes(y = Residuals), 
            color = "#646464",
            size = 0.8) +
  geom_point(aes(y = Residuals), 
             color = "#646464")
############################ Portmanteau tests #################################
#----------------------------- Kwan and Sim -----------------------------------
#------------------  Q4 test (Goodness-of-fit article) ------------------------

# res <- best_fit_arima$residuals # ARIMA residuals
# res_aux <-  res # βARMA residuals
# res <-  res_aux # βARMA residuals
# res <- best_fit_arima$residuals # ARIMA residuals

# Initializing response vector
QKW4 <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
Q4 <- matrix(rep(NA,30), nrow = 1, ncol = 30) 

# P-values
source("Kwan_Chest.R")

# Lag m vector
m <- 8
vm <- c(seq(from = m, to = 30, by = 1))

#res <- res33
f <- 7
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
# plot(1L:nlag, 
#      QKW4, 
#      xlab = "lag",
#      ylab = "p-values",
#      las = 1,
#      ylim = c(0,1),
#      lwd = 1.5)
# abline(h = 0.05, lty = 2, col = "blue", lwd = 1.5)
# 
# # New portmanteau test statistic (Q4) 
# plot(1L:nlag,
#      Q4,
#      xlab = "lag",
#      ylab = "p-values",
#      las = 1,
#      ylim = c(0,1),
#      lwd = 1.5)
# abline(h = 0.05, lty = 2, col = "blue", lwd = 1.5)

df_q4 <- data.frame(Lag = 1:nlag, PValue = as.numeric(Q4))
df_qQKW4 <- data.frame(Lag = 1:nlag, PValue = as.numeric(QKW4))

pormt_test <- 
ggplot(df_q4, aes(x = Lag, y = PValue)) +
  geom_point(shape = 18, color = "darkblue", size = 3) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
  coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
  labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Q4") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),  
        axis.title = element_text(color = "black", face = "bold"),  
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))

pormt_test_arima <- 
ggplot(df_qQKW4, aes(x = Lag, y = PValue)) +
  geom_point(shape = 18, color = "darkblue", size = 3) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
  coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
  labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Q4") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),  
        axis.title = element_text(color = "black", face = "bold"),  
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))

# p <- ggfortify::ggcpgram(df_resid$Residuals)
# p +  # Set y-axis limits
#   labs(x = "Lag", y = "Valores p", title = "Periodograma acumulado dos resíduos") +
#   theme_bw() +  # Use a theme with a boxed layout
#   theme(plot.title = element_text(color = "black", face = "bold"),
#         axis.title = element_text(color = "black", face = "bold"),
#         axis.text = element_text(color = "black"),
#         panel.border = element_rect(color = "black", fill = NA, size = 1))

#-------------------------------------------------------------------------------
# Add ARMA and BARMA predictions to data frame
df_plot_test <- data.frame(
  date = time(df_test),
  observed = as.numeric(df_test)*100000,
  predicted = as.numeric(barma_forecast)*100000
)
df_plot_test$arma_predicted <- as.numeric(arima_forecast)*100000

# Add a column with the date in Date format
# Make sure the year is correct
date_formatted <- df[(nrow(df) - h + 1):nrow(df), "date"]
df_plot_test$date_formatted <- date_formatted$date

# Define the start and end dates for May
start_date <- as.Date(min(df_plot_test$date_formatted))
end_date <- as.Date(max(df_plot_test$date_formatted))

# Out-of-sample plot
predi_both <- 
ggplot(df_plot_test, aes(x = date_formatted)) +
  geom_line(aes(y = observed, color = "Observado", linetype = "Observado"), 
            size = 0.8) +
  geom_line(aes(y = predicted, color = "βARMA", linetype = "βARMA"), 
            size = 0.8) +
  geom_line(aes(y = arma_predicted, color = "ARMA", linetype = "ARMA"),
            size = 0.8) +  # Add this line
  labs(x = "\nData (dia/mês)",
       y = "Razão de mortalidade materna em 2021\n (mortes por 100.000 nascidos vivos)\n",
       title = "Previsão fora da amostra",
       color = "") +
  theme_bw() +
  theme(plot.title = element_text(color = "black", 
                                  face = "bold"),  
        axis.title = element_text(color = "black", 
                                  face = "bold"),  
        axis.text = element_text(color = "black", 
                                 face = "bold"),  
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_x_date(date_breaks = "1 day",
               date_labels = "%d/%m") +
  # scale_x_date(date_breaks = "1 day",
  #              date_labels = "%b %d") +
  scale_y_continuous(labels = scales::comma) +
  coord_cartesian(ylim = c(0, 400)) +
  scale_color_manual(values = c("Observado" = "black", 
                                "βARMA" = "red",
                                "ARMA" = "blue"),
                     breaks = c("Observado", "ARMA", "βARMA")) +
  scale_linetype_manual(values=c("Observado" = "solid", 
                                 "βARMA" = "dashed", 
                                 "ARMA" = "dotted"),
                        breaks = c("Observado", "ARMA", "βARMA"), 
                        name = "")

library("patchwork")
png("pred_in_plot_2021.png", units = "in", width = 14, height = 5, res = 300)
# out and in sample plots
(insample_plot1 +
    insample_plot2 +
    #predi_both +
    plot_layout(ncol = 2))
dev.off()

png("pred_out_plot_2021.png", units = "in", width = 8, height = 5, res = 300)
# out and in sample plots
predi_both  
dev.off()

png("barma_diag_2021.png", units = "in", width = 13, height = 7, res = 300)
(resid_pondera1 +
    acf_funct_barma +
    p_qqplot +
    p_densidade +
    plot_layout(ncol = 2))
dev.off()

# portmanteau Q4
pormt_test
pormt_test_arima


resid_pondera1

# acf
acf_funct_barma

# normality
p_densidade
p_qqplot
tseries::jarque.bera.test(df_resid$Residuals)

cpgram(residuals(best_fit_arima), main = "Periodograma acumulado dos resíduos")
Box.test(residuals(best_fit_arima), type = "Box-Pierce")
Box.test(residuals(best_fit_arima), type = "Ljung-Box")

# Jarque–Bera test
# H0: joint hypothesis of skewness=0 and excess kurtosis=0
tseries::jarque.bera.test(residuals(best_fit_arima))
#-------------------------- Model comparison ----------------------------
# Initializing the result matrix
table_metrics <- matrix(rep(NA, 6*h), nrow = 6, ncol = h)                           

# Fill table with MAE, RMSE, and sMAPE metrics
for(i in 1:h)
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
                             "h=7",
                             "h=8",
                             "h=9",
                             "h=10",
                             "h=11",
                             "h=12"#,
                             # "h=13",
                             # "h=14"
                             )

table_metrics
