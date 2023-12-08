# rm(list = ls()) 
detach("package:e1071", unload = TRUE)
detach("package:tseries", unload = TRUE)
detach("package:tidyverse", unload = TRUE)
detach("package:itsmr", unload = TRUE)
detach("package:tsoutliers", unload = TRUE)

library("e1071") # asymmetry and kustosis.
library("tseries") # time series
library("tidyverse") # data manipulation
library("itsmr") # time series
library("tsoutliers") # outliers

################################################################################
#----------------------------- Load the data -----------------------------------
################################################################################
df <- readRDS("maternal_mortality_rate_2021.rds")
df <- df[1:135,] # up to May 15
df_ts <- ts(df$rate, start = c(2021, 1), frequency = 365)

h <- 7 # forecast window

# Split data in train and test
df_train <- df_ts[1:(length(df_ts) - h)]
df_train <- ts(df_train)

df_test <- df_ts[(length(df_train) + 1):(length(df_ts))]
df_test <- ts(df_test)

plot.ts(df_train)

################################################################################
########################### Intervention Analysis ##############################
################################################################################

# outliers_excess_ts <- tsoutliers::tso(df_train_clean)
# plot(outliers_excess_ts) #Detection graph.
# outliers_excess_ts$outliers #General information about the points found.
# outliers_idx <- outliers_excess_ts$outliers$ind #Position of outliers in the series.
# # plot(forecast::tsclean(df_train))
# # plot(df_train)
# #length of our time series
# n <- length(df_train)
# 
# (col_int <- outliers(c("AO", "LS", "AO", "AO"), outliers_idx[c(1, 3, 4, 6)])) #Columns of interest.
# 
# #Series skeleton identifying the position of the outlier.
# esq_posi <- outliers.effects(col_int, n)
# # efeito_outlier_series <- as.matrix(
# #   rowSums(esq_posi[, c("LS81")])
# # )
# efeito_outlier_series <- esq_posi
# 
# #Outlier effect coefficients.
# omega_hat <- unlist(outliers_excess_ts$outliers["coefhat"])
# 
# #Calculating vector that represents the effect of the outlier.
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
# # #Subtracting the effect of the intervention.
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
#            main = "Time series and outlier effects"))
# # 
# # 
# # # plot.ts(df_train_clean)
# df_train <- df_train_clean

################################################################################
#------------------------------- Check up --------------------------------------
################################################################################

# Data description.
summary(df_train)  # Resume
var(df_train)      # Variance.
skewness(df_train) # Skewness.
kurtosis(df_train) # Kurtosis.

# Apply difference transformation to make data stationary
df_train_diff <- diff(df_train)
plotc(df_train, df_train_diff) # diff = 1 vs diff = 0

# #Longmemory tests (d = 0.3527578 seems to highlight a period of 7)
# d.GPH <- fracdiff::fdGPH(df_train,bandw.exp = 0.7)
# for(i in 1:10)
#   print(cbind(beta=i*0.1,GPH=fracdiff::fdGPH(df_train,bandw.exp = i*0.1)$d))
# df_train_diff <- fracdiff::diffseries(
#   df_train,fracdiff::fdGPH(df_train,bandw.exp = 0.7)$d
#   )


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
source("supporting_scripts/barma.r") # barma
source("supporting_scripts/barma.fit.r") # barma
source('supporting_scripts/kum-mu-phi.r') # karma
source('supporting_scripts/karma.modified.fit.r') # karma
source('supporting_scripts/karma.modified.r') # karma
source("supporting_scripts/Kwan_Chest.R") # portmanteau
source("supporting_scripts/dufour.r") # portmanteau (BARMA specific)

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
  c = sd(y) 
  # Smithson and Verkuilen (2006)
  #y = (y - a)/(b - a)
  y = (y - a + c) / (b + c - a + c)
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

summary(y)
plot.ts(y, ylim = c(0, 1))
################################################################################
#--------------------------------- Fit BARMA -----------------------------------
################################################################################
#---------------------------------
# min AIC and BIC (Portmontau failed)
fit_barma_best <- barma( 
  y,
  ar = c(1, 2, 4, 5),
  ma = c(3, 4, 5),
  h = h,
  diag = 0,
  resid = 1,
  link = "logit" 
)
fit_barma_best_order <- list(
  ar = c(1, 2, 4, 5),
  ma = c(3, 4, 5)
)

# second_min AIC and BIC (resid 3) ---------- SELECTED βARMA MODEL
fit_barma_best <- barma( 
  y,
  ar = c(1, 2, 5, 6),
  ma = c(3, 5),
  h = h,
  diag = 0,
  resid = 1,
  link = "logit" 
)
fit_barma_best_order <- list(
  ar = c(1, 2, 5, 6),
  ma = c(3, 5)
)

barma_forecast <- fit_barma_best$forecast
barma_fitted <- fit_barma_best$fitted
res_barma <- fit_barma_best$resid3 

#----------------------- Undo transformations and differences ------------------
if (extremes){
  barma_forecast <- ((barma_forecast*n_data)-0.5)/(n_data-1)
  barma_fitted <- ((barma_fitted*n_data)-0.5)/(n_data-1)
}

if (d > 0){
  barma_forecast <- barma_forecast*(b+c-a+c)+(a-c)#barma_forecast*(b-a)+a
  # cumsum(c(df_train[length(df_train)], barma_forecast))
  barma_forecast <- diffinv(barma_forecast, 
                            differences = d,
                            xi = df_train[length(df_train)])
  barma_forecast <- barma_forecast[-1]
  barma_fitted <- barma_fitted*(b+c-a+c)+(a-c)#barma_fitted*(b-a)+a
}

# Fitted has max(ar,ma) NAs at the beginning
barma_nas = max(max(fit_barma_best_order$ar), max(fit_barma_best_order$ma))
if (any(is.na(fit_barma_best_order$ar))) barma_nas = max(fit_barma_best_order$ma)
if (any(is.na(fit_barma_best_order$ma))) barma_nas = max(fit_barma_best_order$ar)

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
                        frequency = 365)
ts_barma_fitted <- ts(barma_fitted,
                      start = c(2021, barma_nas+2), 
                      frequency = 365)

#---------------------------- Diagnostics and Graphs ---------------------------

# BARMA Residuals
df_resid_barma <- data.frame(Index = 1:length(res_barma), Residuals = res_barma)

# H0: joint hypothesis of skewness=0 and excess kurtosis=0
tseries::jarque.bera.test(df_resid_barma$Residuals) # Jarque–Bera
# H0: data came from a normally distributed population
shapiro.test(df_resid_barma$Residuals) # Shapiro-Wilk

# Residuals x Index
resid_index_barma <- 
  ggplot(df_resid_barma, aes(x = Index, y = Residuals)) +
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
       title = "Resíduo ponderado padronizado (βARMA)") +
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

# Q-Q plot
barma_qqplot <- 
  ggplot(df_resid_barma, aes(sample = Residuals)) +
  ggplot2::stat_qq(shape = 4, 
                   color = "black") +  
  qqplotr::stat_qq_line(color = "darkblue", 
                        size = 1) +  
  qqplotr::stat_qq_band(bandType = "pointwise", 
                        mapping = aes(fill = "Bootstrap"), 
                        conf = 0.95, 
                        alpha = 0.3) +  
  scale_fill_manual(values = c("Bootstrap" = "#0066cc")) +
  guides(fill = "none") +
  labs(x = "Quantis teóricos (normal)", 
       y = "Quantis amostrais", 
       title = "Gráfico Q-Q dos resíduos (βARMA)") +
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

# Density
barma_densidade <- 
  ggplot(df_resid_barma, aes(x=Residuals)) +
  geom_density(aes(color = "Densidade dos resíduos"), 
               alpha=.2, fill="#0066cc",
               size = 1.1,
               linetype = "solid") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(df_resid_barma$Residuals), 
                            sd = sd(df_resid_barma$Residuals)), 
                aes(color = "Aproximação da distribuição normal"),
                size = 1.1,
                linetype = "dashed") +
  labs(title='Distribuição dos resíduos (βARMA)', 
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

# ACF
acf_func_barma <- 
  forecast::ggAcf(df_resid_barma$Residuals, lag.max = 15) +
  labs(title = "Função de autocorrelação dos resíduos (βARMA)") +
  theme_bw() +  
  theme(
    plot.title = element_text(color = "black", face = "bold"),  
    axis.title = element_text(color = "black", face = "bold"),  
    axis.text = element_text(color = "black", face = "bold"),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) 

#PACF
pacf_func_barma <- 
  forecast::ggPacf(df_resid_barma$Residuals, lag.max = 15) +
  labs(title = "Função de autocorrelação parcial dos resíduos (βARMA)") +
  theme_bw() +  
  theme(
    plot.title = element_text(color = "black", face = "bold"),  
    axis.title = element_text(color = "black", face = "bold"),  
    axis.text = element_text(color = "black", face = "bold"),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) 

# Cumulative periodogram
acum_periodogram_barma <- 
  ggfortify::ggcpgram(df_resid_barma$Residuals, taper = 0.4,) +
  labs(x = "Lag", y = "Valores p", title = "Periodograma acumulado dos resíduos (βARMA)") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),
        axis.title = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

cpgram(df_resid_barma$Residuals, main = "Cumulative Periodogram of Residuals (βARMA)")

#-------------------------------- Portmanteau ----------------------------------  
# Ruey Tsay (2005) Analysis of Financial Time Series, 2nd ed. (Wiley, ch. 2)
# (My suggestion) calculate to m = p' + q' + log(T)... | (m > p + q)

# Initializing response vector
KwanS_P <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
Q4_P <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
LB_P <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
DUF_P <- matrix(rep(NA,30), nrow = 1, ncol = 30) 

KwanS_stat <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
Q4_stat <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
LB_stat <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
DUF_stat <- matrix(rep(NA,30), nrow = 1, ncol = 30) 

# df
barma_fitdf <-  4 + 3 # p'+q'

# Lag m vector
# min(10, T/5) robjhyndman
m <- barma_fitdf + round(log(length(df_resid_barma$Residuals))) # p'+q'+log(T)
vm <- c(seq(from = m, to = 30, by = 1))

for(m in vm)
{
  # Kwan and Sim test(4) 
  KwanS_P[,m] <- Kwan.sim.chest(df_resid_barma$Residuals, 
                             lag = m, 
                             fitdf = barma_fitdf, 
                             type = "correlation",
                             test = 4)$p.value 
  # New portmanteau test statistic (Goodness-of-fit article)
  Q4_P[,m] <- Kwan.sim.chest(df_resid_barma$Residuals, 
                           lag = m, 
                           fitdf = barma_fitdf,
                           type = "partial", 
                           test = 4)$p.value 
  # Ljung-Box
  LB_P[,m] <- Box.test(df_resid_barma$Residuals, 
                     lag=m, 
                     fitdf=barma_fitdf, 
                     type="Ljung-Box")$p.value 
  
  # Robust Dufour and Roy nonparametric test (Residuals ranks)
  DUF_P[,m] <- Dufour.test(df_resid_barma$Residuals, 
                         lag=m, 
                         fitdf=barma_fitdf)$p.value 
  
  # Kwan and Sim test(4) 
  KwanS_stat[,m] <- Kwan.sim.chest(df_resid_barma$Residuals, 
                              lag = m, 
                              fitdf = barma_fitdf, 
                              type = "correlation",
                              test = 4)$statistic
  # New portmanteau test statistic (Goodness-of-fit article)
  Q4_stat[,m] <- Kwan.sim.chest(df_resid_barma$Residuals, 
                           lag = m, 
                           fitdf = barma_fitdf,
                           type = "partial", 
                           test = 4)$statistic 
  # Ljung-Box
  LB_stat[,m] <- Box.test(df_resid_barma$Residuals, 
                     lag=m, 
                     fitdf=barma_fitdf, 
                     type="Ljung-Box")$statistic 
  
  # Robust Dufour and Roy nonparametric test (Residuals ranks)
  DUF_stat[,m] <- Dufour.test(df_resid_barma$Residuals, 
                         lag=m, 
                         fitdf=barma_fitdf)$statistic
}

nlag <- 30

df_KwanS_P <- data.frame(Lag = 1:nlag, PValue = as.numeric(KwanS_P))
df_q4_P <- data.frame(Lag = 1:nlag, PValue = as.numeric(Q4_P))
df_LB_P <- data.frame(Lag = 1:nlag, PValue = as.numeric(LB_P))
df_DUF_P <- data.frame(Lag = 1:nlag, PValue = as.numeric(DUF_P))

df_KwanS_stat <- data.frame(Lag = 1:nlag, stat = as.numeric(KwanS_stat))
df_q4_stat <- data.frame(Lag = 1:nlag, stat = as.numeric(Q4_stat))
df_LB_stat <- data.frame(Lag = 1:nlag, stat = as.numeric(LB_stat))
df_DUF_stat <- data.frame(Lag = 1:nlag, stat = as.numeric(DUF_stat))

# Kwan and Sim test
pormt_KwanS_barma <- 
  ggplot(df_KwanS_P[!is.na(df_KwanS_P$PValue),], aes(x = Lag, y = PValue)) +
  geom_point(shape = 18, color = "darkblue", size = 3) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
  coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
  labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Kwan and Sim 4 (βARMA)") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),  
        axis.title = element_text(color = "black", face = "bold"),  
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))

# Q4 (BARMA specific)
pormt_q4_barma <- 
  ggplot(df_q4_P[!is.na(df_q4_P$PValue),], aes(x = Lag, y = PValue)) +
  geom_point(shape = 18, color = "darkblue", size = 3) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
  coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
  labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Q4 (βARMA)") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),  
        axis.title = element_text(color = "black", face = "bold"),  
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))

# Ljung-Box
pormt_LB_barma <- 
  ggplot(df_LB_P[!is.na(df_LB_P$PValue),], aes(x = Lag, y = PValue)) +
  geom_point(shape = 18, color = "darkblue", size = 3) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
  coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
  labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Ljung-Box (βARMA)") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),  
        axis.title = element_text(color = "black", face = "bold"),  
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))

# Robust Dufour and Roy 
#Dufour, J. M., & Roy, R. (1986). Generalized portmanteau statistics and tests 
#of randomness. Communications in Statistics - Theory and Methods
pormt_DUF_barma <- 
  ggplot(df_DUF_P[!is.na(df_DUF_P$PValue),], aes(x = Lag, y = PValue)) +
  geom_point(shape = 18, color = "darkblue", size = 3) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
  coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limit
  labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Dufour e Roy com ranks (βARMA)") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),  
        axis.title = element_text(color = "black", face = "bold"),  
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))

################################################################################
#--------------------------------- Fit KARMA -----------------------------------
################################################################################

# min AIC  
fit_karma_best <- karma( 
  y,
  ar = c(1, 2, 3, 5, 6),
  ma = c(1, 2, 4, 6),
  h = h,
  diag = 1,
  resid = 1,
  link = "logit",
  prec_start = 5 
)
fit_karma_best_order <- list(
  ar = c(1, 2, 3, 5, 6),
  ma = c(1, 2, 4, 6)
)

# min BIC --------------------- KARMA model selected
fit_karma_best <- karma( 
  y,
  ar = c(4, 5),
  ma = c(1, 4),
  h = h,
  diag = 1,
  resid = 1,
  link = "logit",
  prec_start = 5 
)
fit_karma_best_order <- list(
  ar = c(4, 5),
  ma = c(1, 4)
)

karma_forecast <- fit_karma_best$forecast
karma_fitted <- fit_karma_best$fitted
res_karma <- fit_karma_best$resid3 

#----------------------- Undo transformations and differences ------------------
if (extremes){
  karma_forecast <- ((karma_forecast*n_data)-0.5)/(n_data-1)
  karma_fitted <- ((karma_fitted*n_data)-0.5)/(n_data-1)
}

if (d > 0){
  karma_forecast <- karma_forecast*(b+c-a+c)+(a-c)#karma_forecast*(b-a)+a
  karma_forecast <- diffinv(karma_forecast, 
                            differences = d,
                            xi = df_train[length(df_train)])
  karma_forecast <- karma_forecast[-1]
  karma_fitted <- karma_fitted*(b+c-a+c)+(a-c)#karma_fitted*(b-a)+a
}

# Fitted has max(ar,ma) NAs at the beginning
karma_nas = max(max(fit_karma_best_order$ar), max(fit_karma_best_order$ma))
if (any(is.na(fit_karma_best_order$ar))) karma_nas = max(fit_karma_best_order$ma)
if (any(is.na(fit_karma_best_order$ma))) karma_nas = max(fit_karma_best_order$ar)

# Ignore the NA's as the beginning
karma_fitted <- karma_fitted[(karma_nas+1):length(karma_fitted)]
if (d > 0) {
  karma_fitted <- df_train[
    (karma_nas+1):(karma_nas+length(karma_fitted))
  ]+karma_fitted
}

date_start_forecast = end(df_ts)
date_start_forecast[2] = date_start_forecast[2] - (h-1)

ts_karma_forecast <- ts(karma_forecast, 
                        start = date_start_forecast,
                        frequency = 365)
ts_karma_fitted <- ts(karma_fitted,
                      start = c(2021, karma_nas+2), 
                      frequency = 365)

#---------------------------- Diagnostics and Graphs ---------------------------
# KARMA Residuals
df_resid_karma <- data.frame(Index = 1:length(res_karma), Residuals = res_karma)

# Jarque–Bera test
# H0: joint hypothesis of skewness=0 and excess kurtosis=0
tseries::jarque.bera.test(df_resid_karma$Residuals)
# H0: data came from a normally distributed population
shapiro.test(df_resid_barma$Residuals) # Shapiro-Wilk

# Residuals x Index
resid_index_karma <- 
  ggplot(df_resid_karma, aes(x = Index, y = Residuals)) +
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
       title = "Resíduo quantílico (KARMA)") +
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

# Q-Q plot
karma_qqplot <- 
  ggplot(df_resid_karma, aes(sample = Residuals)) +
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
       title = "Gráfico Q-Q dos resíduos (KARMA)") +
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

# Density
karma_densidade <- 
  ggplot(df_resid_karma, aes(x=Residuals)) +
  geom_density(aes(color = "Densidade dos resíduos"), 
               alpha=.2, fill="#0066cc",
               size = 1.1,
               linetype = "solid") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(df_resid_karma$Residuals), 
                            sd = sd(df_resid_karma$Residuals)), 
                aes(color = "Aproximação da distribuição normal"),
                size = 1.1,
                linetype = "dashed") +
  labs(title='Distribuição dos resíduos (KARMA)', 
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

# ACF
acf_func_karma <- 
  forecast::ggAcf(df_resid_karma$Residuals, lag.max = 15) +
  labs(title = "Função de autocorrelação dos resíduos (KARMA)") +
  theme_bw() +  
  theme(
    plot.title = element_text(color = "black", face = "bold"),  
    axis.title = element_text(color = "black", face = "bold"),  
    axis.text = element_text(color = "black", face = "bold"),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) 

# PACF
acf_func_karma <- 
  forecast::ggPacf(df_resid_karma$Residuals, lag.max = 15) +
  labs(title = "Função de autocorrelação dos resíduos (KARMA)") +
  theme_bw() +  
  theme(
    plot.title = element_text(color = "black", face = "bold"),  
    axis.title = element_text(color = "black", face = "bold"),  
    axis.text = element_text(color = "black", face = "bold"),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) 

# Cumulative periodogram
acum_periodogram_barma <- 
  ggfortify::ggcpgram(df_resid_karma$Residuals) +
  labs(x = "Lag", y = "Valores p", title = "Periodograma acumulado dos resíduos (KARMA)") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),
        axis.title = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

cpgram(df_resid_karma$Residuals, main = "Cumulative Periodogram of Residuals (KARMA)")

#-------------------------------- Portmanteau ----------------------------------  
# Ruey Tsay (2005) Analysis of Financial Time Series, 2nd ed. (Wiley, ch. 2)
# (My suggestion) calculate to m = p' + q' + log(T)... | (m > p + q)

# Initializing response vector
KwanS_P <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
Q4_P <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
LB_P <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
DUF_P <- matrix(rep(NA,30), nrow = 1, ncol = 30) 

KwanS_stat <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
Q4_stat <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
LB_stat <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
DUF_stat <- matrix(rep(NA,30), nrow = 1, ncol = 30) 

# df
karma_fitdf <-  2 + 2 # p'+q'

# Lag m vector
# min(10, T/5) robjhyndman
m <- karma_fitdf + round(log(length(df_resid_karma$Residuals))) # p'+q'+log(T)
vm <- c(seq(from = m, to = 30, by = 1))

for(m in vm)
{
  # Kwan and Sim test(4) 
  KwanS_P[,m] <- Kwan.sim.chest(df_resid_karma$Residuals, 
                                lag = m, 
                                fitdf = karma_fitdf, 
                                type = "correlation",
                                test = 4)$p.value 
  # New portmanteau test statistic (Goodness-of-fit article)
  Q4_P[,m] <- Kwan.sim.chest(df_resid_karma$Residuals, 
                             lag = m, 
                             fitdf = karma_fitdf,
                             type = "partial", 
                             test = 4)$p.value 
  # Ljung-Box
  LB_P[,m] <- Box.test(df_resid_karma$Residuals, 
                       lag=m, 
                       fitdf=karma_fitdf, 
                       type="Ljung-Box")$p.value 
  
  # Robust Dufour and Roy nonparametric test (Residuals ranks)
  DUF_P[,m] <- Dufour.test(df_resid_karma$Residuals, 
                           lag=m, 
                           fitdf=karma_fitdf)$p.value 
  
  # Kwan and Sim test(4) 
  KwanS_stat[,m] <- Kwan.sim.chest(df_resid_karma$Residuals, 
                                   lag = m, 
                                   fitdf = karma_fitdf, 
                                   type = "correlation",
                                   test = 4)$statistic
  # New portmanteau test statistic (Goodness-of-fit article)
  Q4_stat[,m] <- Kwan.sim.chest(df_resid_karma$Residuals, 
                                lag = m, 
                                fitdf = karma_fitdf,
                                type = "partial", 
                                test = 4)$statistic 
  # Ljung-Box
  LB_stat[,m] <- Box.test(df_resid_karma$Residuals, 
                          lag=m, 
                          fitdf=karma_fitdf, 
                          type="Ljung-Box")$statistic 
  
  # Robust Dufour and Roy nonparametric test (Residuals ranks)
  DUF_stat[,m] <- Dufour.test(df_resid_karma$Residuals, 
                              lag=m, 
                              fitdf=karma_fitdf)$statistic
}

nlag <- 30

df_KwanS_P <- data.frame(Lag = 1:nlag, PValue = as.numeric(KwanS_P))
df_q4_P <- data.frame(Lag = 1:nlag, PValue = as.numeric(Q4_P))
df_LB_P <- data.frame(Lag = 1:nlag, PValue = as.numeric(LB_P))
df_DUF_P <- data.frame(Lag = 1:nlag, PValue = as.numeric(DUF_P))

df_KwanS_stat <- data.frame(Lag = 1:nlag, stat = as.numeric(KwanS_stat))
df_q4_stat <- data.frame(Lag = 1:nlag, stat = as.numeric(Q4_stat))
df_LB_stat <- data.frame(Lag = 1:nlag, stat = as.numeric(LB_stat))
df_DUF_stat <- data.frame(Lag = 1:nlag, stat = as.numeric(DUF_stat))

# Kwan and Sim test
pormt_KwanS_karma <- 
  ggplot(df_KwanS_P[!is.na(df_KwanS_P$PValue),], aes(x = Lag, y = PValue)) +
  geom_point(shape = 18, color = "darkblue", size = 3) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
  coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
  labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Kwan and Sim 4 (KARMA)") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),  
        axis.title = element_text(color = "black", face = "bold"),  
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))

# Q4 (BARMA specific)
# pormt_q4_karma <- 
#   ggplot(df_q4_P[!is.na(df_q4_P$PValue),], aes(x = Lag, y = PValue)) +
#   geom_point(shape = 18, color = "darkblue", size = 3) +  
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
#   coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
#   labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Q4 (KARMA)") +
#   theme_bw() +  # Use a theme with a boxed layout
#   theme(plot.title = element_text(color = "black", face = "bold"),  
#         axis.title = element_text(color = "black", face = "bold"),  
#         axis.text = element_text(color = "black"), 
#         panel.border = element_rect(color = "black", fill = NA, size = 1))

# Ljung-Box
pormt_LB_karma <- 
  ggplot(df_LB_P[!is.na(df_LB_P$PValue),], aes(x = Lag, y = PValue)) +
  geom_point(shape = 18, color = "darkblue", size = 3) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
  coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
  labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Ljung-Box (KARMA)") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),  
        axis.title = element_text(color = "black", face = "bold"),  
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))

# Robust Dufour and Roy 
#Dufour, J. M., & Roy, R. (1986). Generalized portmanteau statistics and tests 
#of randomness. Communications in Statistics - Theory and Methods
pormt_DUF_karma <- 
  ggplot(df_DUF_P[!is.na(df_DUF_P$PValue),], aes(x = Lag, y = PValue)) +
  geom_point(shape = 18, color = "darkblue", size = 3) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
  coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
  labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Dufour e Roy com ranks (KARMA)") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),  
        axis.title = element_text(color = "black", face = "bold"),  
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))

################################################################################
#--------------------------------- Fit ARIMA -----------------------------------
################################################################################

# BEST MODEL BY AIC: -1576.45
fit_arima <- forecast::Arima( 
  df_train,
  order = c(6, 1, 6),
  fixed = c(NA, NA, NA, 0, 0, NA, # ar
            0, NA, 0, NA, 0, NA) # ma
)

best_fit_arima <- fit_arima

# #BEST MODEL BY auto.arima
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

# Forecasts
arima_forecast <- forecast::forecast(best_fit_arima, h = h)$mean
arima_forecast <- as.numeric(arima_forecast)

#---------------------------- Diagnostics and Graphs ---------------------------
# ARMA residuals
df_resid_arma <- data.frame(
  Index = 1:length(fit_arima$residuals), 
  Residuals = as.numeric(fit_arima$residuals)
)

lmtest::coeftest(best_fit_arima) 
forecast::checkresiduals(best_fit_arima, test = "LB")
Box.test(residuals(best_fit_arima), type = "Box-Pierce")
Box.test(residuals(best_fit_arima), type = "Ljung-Box")

# Jarque–Bera test
# H0: joint hypothesis of skewness=0 and excess kurtosis=0
tseries::jarque.bera.test(df_resid_arma$Residuals)
# H0: data came from a normally distributed population
shapiro.test(df_resid_barma$Residuals) # Shapiro-Wilk

plotc(df_train,best_fit_arima$fitted)
forecast::accuracy(best_fit_arima$fitted, df_train)
tsdiag(best_fit_arima, gof.lag = 25)
test(best_fit_arima$residuals)
car::qqPlot(best_fit_arima$res)
tsdiag(best_fit_arima)

# Residuals x Index
resid_index_arma <- 
  ggplot(df_resid_arma, aes(x = Index, y = Residuals)) +
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
  coord_cartesian(ylim = c(-0.0025, 0.0025)) +  
  labs(x = "Índice", 
       y = "Resíduo", 
       title = "Resíduos (ARIMA)") +
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

# Q-Q plot
arma_qqplot <- 
  ggplot(df_resid_arma, aes(sample = Residuals)) +
  ggplot2::stat_qq(shape = 4, 
                   color = "black") +  
  qqplotr::stat_qq_line(color = "darkblue", 
                        size = 1) +  
  qqplotr::stat_qq_band(bandType = "ks", 
                        mapping = aes(fill = "Bootstrap"), 
                        conf = 0.95, 
                        alpha = 0.3) +  
  scale_fill_manual(values = c("Bootstrap" = "#0066cc")) +
  guides(fill = FALSE) +
  labs(x = "Quantis teóricos (normal)", 
       y = "Quantis amostrais", 
       title = "Gráfico Q-Q dos resíduos (ARIMA)") +
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
# Create QQ plot
car::qqPlot(df_resid_arma$Residuals, 
            main = "Gráfico Q-Q dos resíduos (ARIMA)", 
            xlab = "Quantis teóricos (normal)", 
            ylab = "Quantis amostrais",
            pch = 4, 
            col = "black", 
            cex = 1, 
            las = 1, 
            grid = TRUE)

# Density
arma_densidade <- 
  ggplot(df_resid_arma, aes(x=Residuals)) +
  geom_density(aes(color = "Densidade dos resíduos"), 
               alpha=.2, fill="#0066cc",
               size = 1.1,
               linetype = "solid") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(df_resid_arma$Residuals), 
                            sd = sd(df_resid_arma$Residuals)), 
                aes(color = "Aproximação da distribuição normal"),
                size = 1.1,
                linetype = "dashed") +
  labs(title='Distribuição dos resíduos (ARIMA)', 
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

# ACF
acf_func_arma <- 
  forecast::ggAcf(df_resid_arma$Residuals, lag.max = 15) +
  labs(title = "Função de autocorrelação dos resíduos (ARIMA)") +
  theme_bw() +  
  theme(
    plot.title = element_text(color = "black", face = "bold"),  
    axis.title = element_text(color = "black", face = "bold"),  
    axis.text = element_text(color = "black", face = "bold"),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) 

# PACF
acf_func_arma <- 
  forecast::ggPacf(df_resid_arma$Residuals, lag.max = 15) +
  labs(title = "Função de autocorrelação dos resíduos (ARIMA)") +
  theme_bw() +  
  theme(
    plot.title = element_text(color = "black", face = "bold"),  
    axis.title = element_text(color = "black", face = "bold"),  
    axis.text = element_text(color = "black", face = "bold"),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) 

# Cumulative periodogram
acum_periodogram_barma <- 
  ggfortify::ggcpgram(df_resid_arma$Residuals) +
  labs(x = "Lag", y = "Valores p", title = "Periodograma acumulado dos resíduos (ARIMA)") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),
        axis.title = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

cpgram(df_resid_arma$Residuals, main = "Cumulative Periodogram of Residuals (ARIMA)")

#-------------------------------- Portmanteau ----------------------------------  
# Ruey Tsay (2005) Analysis of Financial Time Series, 2nd ed. (Wiley, ch. 2)
# (My suggestion) calculate to m = p' + q' + log(T)... | (m > p + q)

# Initializing response vector
KwanS_P <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
Q4_P <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
LB_P <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
DUF_P <- matrix(rep(NA,30), nrow = 1, ncol = 30) 

KwanS_stat <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
Q4_stat <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
LB_stat <- matrix(rep(NA,30), nrow = 1, ncol = 30) 
DUF_stat <- matrix(rep(NA,30), nrow = 1, ncol = 30) 

# df
arma_fitdf <-  1 + 2 # p'+q'

# Lag m vector
# min(10, T/5) robjhyndman
m <- arma_fitdf + round(log(length(df_resid_arma$Residuals))) # p'+q'+log(T)
vm <- c(seq(from = m, to = 30, by = 1))

for(m in vm)
{
  # Kwan and Sim test(4) 
  KwanS_P[,m] <- Kwan.sim.chest(df_resid_arma$Residuals, 
                                lag = m, 
                                fitdf = arma_fitdf, 
                                type = "correlation",
                                test = 4)$p.value 
  # New portmanteau test statistic (Goodness-of-fit article)
  Q4_P[,m] <- Kwan.sim.chest(df_resid_arma$Residuals, 
                             lag = m, 
                             fitdf = arma_fitdf,
                             type = "partial", 
                             test = 4)$p.value 
  # Ljung-Box
  LB_P[,m] <- Box.test(df_resid_arma$Residuals, 
                       lag=m, 
                       fitdf=arma_fitdf, 
                       type="Ljung-Box")$p.value 
  
  # Robust Dufour and Roy nonparametric test (Residuals ranks)
  DUF_P[,m] <- Dufour.test(df_resid_arma$Residuals, 
                           lag=m, 
                           fitdf=arma_fitdf)$p.value 
  
  # Kwan and Sim test(4) 
  KwanS_stat[,m] <- Kwan.sim.chest(df_resid_arma$Residuals, 
                                   lag = m, 
                                   fitdf = arma_fitdf, 
                                   type = "correlation",
                                   test = 4)$statistic
  # New portmanteau test statistic (Goodness-of-fit article)
  Q4_stat[,m] <- Kwan.sim.chest(df_resid_arma$Residuals, 
                                lag = m, 
                                fitdf = arma_fitdf,
                                type = "partial", 
                                test = 4)$statistic 
  # Ljung-Box
  LB_stat[,m] <- Box.test(df_resid_arma$Residuals, 
                          lag=m, 
                          fitdf=arma_fitdf, 
                          type="Ljung-Box")$statistic 
  
  # Robust Dufour and Roy nonparametric test (Residuals ranks)
  DUF_stat[,m] <- Dufour.test(df_resid_arma$Residuals, 
                              lag=m, 
                              fitdf=arma_fitdf)$statistic
}

nlag <- 30

df_KwanS_P <- data.frame(Lag = 1:nlag, PValue = as.numeric(KwanS_P))
df_q4_P <- data.frame(Lag = 1:nlag, PValue = as.numeric(Q4_P))
df_LB_P <- data.frame(Lag = 1:nlag, PValue = as.numeric(LB_P))
df_DUF_P <- data.frame(Lag = 1:nlag, PValue = as.numeric(DUF_P))

df_KwanS_stat <- data.frame(Lag = 1:nlag, stat = as.numeric(KwanS_stat))
df_q4_stat <- data.frame(Lag = 1:nlag, stat = as.numeric(Q4_stat))
df_LB_stat <- data.frame(Lag = 1:nlag, stat = as.numeric(LB_stat))
df_DUF_stat <- data.frame(Lag = 1:nlag, stat = as.numeric(DUF_stat))

# Kwan and Sim test
pormt_KwanS_arma <- 
  ggplot(df_KwanS_P[!is.na(df_KwanS_P$PValue),], aes(x = Lag, y = PValue)) +
  geom_point(shape = 18, color = "darkblue", size = 3) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
  coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
  labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Kwan and Sim 4 (ARIMA)") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),  
        axis.title = element_text(color = "black", face = "bold"),  
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))

# Q4 (BARMA specific)
# pormt_q4_arma <- 
#   ggplot(df_q4_P[!is.na(df_q4_P$PValue),], aes(x = Lag, y = PValue)) +
#   geom_point(shape = 18, color = "darkblue", size = 3) +  
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
#   coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
#   labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Q4 (ARIMA)") +
#   theme_bw() +  # Use a theme with a boxed layout
#   theme(plot.title = element_text(color = "black", face = "bold"),  
#         axis.title = element_text(color = "black", face = "bold"),  
#         axis.text = element_text(color = "black"), 
#         panel.border = element_rect(color = "black", fill = NA, size = 1))

# Ljung-Box
pormt_LB_arma <- 
  ggplot(df_LB_P[!is.na(df_LB_P$PValue),], aes(x = Lag, y = PValue)) +
  geom_point(shape = 18, color = "darkblue", size = 3) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
  coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
  labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Ljung-Box (ARIMA)") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),  
        axis.title = element_text(color = "black", face = "bold"),  
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))

# Robust Dufour and Roy 
#Dufour, J. M., & Roy, R. (1986). Generalized portmanteau statistics and tests 
#of randomness. Communications in Statistics - Theory and Methods
pormt_DUF_arma <- 
  ggplot(df_DUF_P[!is.na(df_DUF_P$PValue),], aes(x = Lag, y = PValue)) +
  geom_point(shape = 18, color = "darkblue", size = 3) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) + 
  coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
  labs(x = "Lag", y = "Valores p", title = "Teste portmanteau Dufour e Roy com ranks (ARIMA)") +
  theme_bw() +  # Use a theme with a boxed layout
  theme(plot.title = element_text(color = "black", face = "bold"),  
        axis.title = element_text(color = "black", face = "bold"),  
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1))

################################################################################
####################### Final results and visualizations #######################
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

# KARMA
df_plot_karma <- data.frame(
  date = time(df_train[1:length(karma_fitted)]),
  observed = as.numeric(df_train[1:length(karma_fitted)])*100000,
  predicted = as.numeric(karma_fitted)*100000
)

df_plot$predicted3 <- df_plot_karma[1:nrow(df_plot), "predicted"]

#-------------- Save high resolution png plot ---------------
#png("new_In_sample_plot_2021.png", units = "in", width = 8, height = 4, res = 300)

# In-sample plot
insample_plot_barma <- 
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

insample_plot_arma <- 
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

insample_plot_karma <- 
  ggplot(df_plot, aes(x = date_formatted)) +
  geom_line(aes(y = observed), 
            color = "#646464",
            size = 0.8) +
  geom_point(aes(y = observed), 
             color = "#646464") +
  geom_line(aes(y = predicted3), 
            color = "#29CE60", 
            linetype = "solid",
            size = 0.9) +
  geom_point(aes(y = predicted3), 
             color = "#29CE60") +
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
           label = "Previsão KARMA", 
           hjust = 0, 
           vjust = 1,
           color = "#29CE60")
#dev.off()

#-------------------------------------------------------------------------------
# Add ARMA, BARMA ans KARMA predictions to data frame
df_plot_test <- data.frame(
  date = time(df_test),
  observed = as.numeric(df_test)*100000,
  predicted = as.numeric(barma_forecast)*100000
)
df_plot_test$arma_predicted <- as.numeric(arima_forecast)*100000
df_plot_test$karma_predicted <- as.numeric(karma_forecast)*100000

# Add a column with the date in Date format
# Make sure the year is correct
date_formatted <- df[(nrow(df) - h + 1):nrow(df), "date"]
df_plot_test$date_formatted <- date_formatted$date

# Define the start and end dates for May
start_date <- as.Date(min(df_plot_test$date_formatted))
end_date <- as.Date(max(df_plot_test$date_formatted))

# Out-of-sample plot
predi_all <- 
ggplot(df_plot_test, aes(x = date_formatted)) +
  geom_line(aes(y = observed, color = "Observado", linetype = "Observado"), 
            size = 0.8) +
  geom_line(aes(y = predicted, color = "βARMA", linetype = "βARMA"), 
            size = 0.8) +
  geom_line(aes(y = karma_predicted, color = "KARMA", linetype = "KARMA"), 
            size = 0.8) +
  geom_line(aes(y = arma_predicted, color = "ARMA", linetype = "ARMA"),
            size = 0.8) +  # Add this line
  labs(x = "\nData (dia/mês)",
       y = "Razão de mortalidade materna em 2021\n (mortes por 100.000 nascidos vivos)\n",
       title = "Previsão fora da amostra treino",
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
                                "KARMA" = "#29CE60",
                                "ARMA" = "blue"),
                     breaks = c("Observado", "βARMA", "KARMA", "ARMA")) +
  scale_linetype_manual(values=c("Observado" = "solid", 
                                 "βARMA" = "dashed",
                                 "KARMA" = "twodash",
                                 "ARMA" = "dotted"),
                        breaks = c("Observado", "βARMA", "KARMA", "ARMA"), 
                        name = "")

#-------------------------------------------------------------------------------
# Join figures

library("patchwork")
#png("pred_in_plot_2021.png", units = "in", width = 14, height = 5, res = 300)
# out and in sample plots
(insample_plot1 +
    insample_plot2 +
    #predi_both +
    plot_layout(ncol = 2))
#dev.off()

#png("pred_out_plot_2021.png", units = "in", width = 8, height = 5, res = 300)
# out and in sample plots
predi_both  
#dev.off()

#png("barma_diag_2021.png", units = "in", width = 13, height = 7, res = 300)
(resid_pondera1 +
    acf_funct_barma +
    p_qqplot +
    p_densidade +
    plot_layout(ncol = 2))
#dev.off()

# portmanteau 


# acf


# normality



################################################################################
#------------------------------ Model comparison -------------------------------
################################################################################

# Initializing the result matrix
table_metrics <- matrix(rep(NA, 9*h), nrow = 9, ncol = h)                           

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
  table_metrics[3,i] <- round(
    mean(abs((df_test[1:i] - karma_forecast[1:i])))*100,
    4
  )
  
  # RMSE
  table_metrics[4,i] <- round(
    sqrt(mean((df_test[1:i] - barma_forecast[1:i])^2))*100,
    4
  )
  table_metrics[5,i] <- round(
    sqrt(mean((df_test[1:i] - arima_forecast[1:i])^2))*100,
    4
  )
  table_metrics[6,i] <- round(
    sqrt(mean((df_test[1:i] - karma_forecast[1:i])^2))*100,
    4
  )
  
  # sMAPE
  table_metrics[7,i] <- round(
    mean(2 * abs(df_test[1:i] - barma_forecast[1:i])/
           (abs(df_test[1:i]) + abs(barma_forecast[1:i]))),
    4
  )
  table_metrics[8,i] <- round(
    mean(2 * abs(df_test[1:i] - arima_forecast[1:i])/
           (abs(df_test[1:i]) + abs(arima_forecast[1:i]))),
    4
  )
  table_metrics[9,i] <- round(
    mean(2 * abs(df_test[1:i] - karma_forecast[1:i])/
           (abs(df_test[1:i]) + abs(karma_forecast[1:i]))),
    4
  )
}

# Set the row names to indicate the metric and method
rownames(table_metrics) <- c("MAE(BARMA)x100", 
                             "MAE(ARIMA)x100",
                             "MAE(KARMA)x100",
                             "RMSE(BARMA)",
                             "RMSE(ARIMA)",
                             "RMSE(KARMA)",
                             "sMAPE(BARMA)", 
                             "sMAPE(ARIMA)",
                             "sMAPE(KARMA)")
colnames(table_metrics) <- c("h=1", 
                             "h=2", 
                             "h=3", 
                             "h=4", 
                             "h=5", 
                             "h=6",
                             "h=7"
                             )

table_metrics
