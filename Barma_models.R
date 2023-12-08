# rm(list = ls()) 

# Valid models up to order 6 can be loaded with: 
# load("models/barma_models_order6.RData")
#-------------------------------------------------------------------------------
# Options for barma():

# resid = 1 : standardized residual
# resid = 2 : standardized residual 2 (predictor scale)
# resid = 3 : standardized weighted residual
# resid = 4 : deviance residual

# diag = 0 : do not plot any graph (useful for simulations)
# diag = 1 : plot graphs
# diag = 2 : save graphs on ps files

# link = "logit"
# link = "probit"
# link = "cloglog"

# h is the prediction window

# Example

# barma(y, ar = p, ma = q, h = h, diag = 0, resid = 1, link = "logit")
#-------------------------------------------------------------------------------
library("tseries")
library("tidyverse")
library("itsmr")
library("tsoutliers")
library("foreach") # parallel loop
library("doParallel") # parallel computation

# Load the data
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
#-------------------------------------------------------------------------------
# Required functions
source("supporting_scripts/barma.r")
source("supporting_scripts/barma.fit.r")
source('supporting_scripts/model_orders.r')

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
#-------------------------- Fit BARMA ------------------------------------------
#Possible combinations including: 0, 1 parameter, 2 parameters, ... 6 parameters
#1 + 6 + 15 + 20 + 15 + 6 + 1 = 64 models 
#64*64 = 4096 models considering, ar 0 up to ar 6 and ma 0 up to ma 6

list_combinations <- list_combinations_BK6 # all combinations up to order 6

# Register the parallel backend
# Detect the phisical number of cores (logical = FALSE detects threads)
n_cores <- detectCores(logical = FALSE) 
cf <- 1 # Number of free CPU cores
num_cores <- n_cores - cf # Reccomended at leat 1 CPU core free  
cl <- makeCluster(num_cores)
registerDoParallel(cl) 

# Use foreach with %:% for nested loops
result <- foreach(k = 1:length(list_combinations),
                  .combine = rbind) %dopar% {
  # Start variables for each outer loop iteration
  valid_models_ar_coef <- list()
  valid_models_ma_coef <- list()
  valid_models_aic <- list()
  valid_models_bic <- list()
  
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
  
  return(data.frame(
    valid_models_ar_coef = I(valid_models_ar_coef),
    valid_models_ma_coef = I(valid_models_ma_coef),
    valid_models_aic = I(valid_models_aic),
    valid_models_bic = I(valid_models_bic)
    )
  )
}
  

# Stop the cluster (REMEMBER TO STOP THE CLUSTER)
stopCluster(cl)

# Combine the results from each iteration
final_result <- do.call(rbind, result)

## list of valid models (54 out of 4096 founded)
final_result[ , ]

# Save valid models
# barma_models_order6 <- final_result
# save(barma_models_order6, file = "models/barma_models_order6.RData")
# load("models/barma_models_order6.RData")

# Get the index of the first, second, third and fourth smallest AIC values
# final_result <- barma_models_order6 -------- HERE AFTER LOADING
#select_index <- final_result["valid_models_aic", ]
#select_index <- final_result["valid_models_bic", ]
index_min <- which.min(select_index) # first
select_index[index_min] <- NA
index_second_min <- which.min(select_index) # second
select_index[index_second_min] <- NA
index_third_min <- which.min(select_index) # third
select_index[index_third_min] <- NA
index_fourth_min <- which.min(select_index) # fourth
#min: AIC -242.229| model selected: AIC -231.837, BIC -209.083

# Position in the list of 1 out of 4 smallest AIC founded
posi <- index_second_min

# Selected model
#final_result["valid_models_aic", posi] # AIC 
#final_result["valid_models_bic", posi] # BIC 
best_barma_combination_ar_aic <- 
  final_result["valid_models_ar_coef", posi][[1]] # best
best_barma_combination_ma_aic <- 
  final_result["valid_models_ma_coef", posi][[1]] # best

# Verify the selected model
fit_barma_best_aic <- barma(
  y,
  ar = best_barma_combination_ar_aic, # 1, 2, 5, 6
  ma = best_barma_combination_ma_aic, # 3, 5
  h = h,
  diag = 1,
  resid = 1,
  link = "logit"
)

## WITHOUT PARALLEL COMPUTATION: SLOW  
# Start variables
# valid_models_ar_coef <- list()
# valid_models_ma_coef <- list()
# valid_models_aic <- list()
# valid_models_bic <- list()
# for (k in 1:length(list_combinations)) {
#   ar <- list_combinations[[k]]
#   for (j in 1:length(list_combinations)) {
#     tryCatch({ # don't stop on errors
#       fit_barma <- barma(
#         y,
#         ar = ar,
#         ma = list_combinations[[j]],
#         h = h,
#         diag = 0,
#         #resid = 1,
#         link = "logit"
#       )
#       coef_df <- as.data.frame(fit_barma$model)
#       coef_pvalues <- coef_df$`Pr(>|z|)` # p-values
#       
#       coef_names <- rownames(coef_df)
#       coef_ar <- coef_df[grepl("phi", coef_names, fixed = TRUE), ]
#       coef_ar_values <- coef_ar$Estimate # ar coef
#       coef_ma <- coef_df[grepl("theta", coef_names, fixed = TRUE), ]
#       coef_ma_values <- coef_ma$Estimate # ma coef
#       
#       # Verify the valid models
#       if (
#         all(coef_pvalues < 0.05) & # 5% significance
#         length(coef_ar_values) != 0 & # At leat 1 ar coef
#         all(coef_ar_values > -0.9 & coef_ar_values < 0.9) & # Causal model
#         all(coef_ma_values > -1.0 & coef_ma_values < 1.0) 
#         ){
#         valid_models_ar_coef <- append(
#           valid_models_ar_coef, 
#           list(ar)
#         )
#         valid_models_ma_coef <- append(
#           valid_models_ma_coef, 
#           list(list_combinations[[j]])
#         )
#         valid_models_aic <- append(
#           valid_models_aic, 
#           list(fit_barma$aic)
#         )
#         valid_models_bic <- append(
#           valid_models_bic, 
#           list(fit_barma$bic)
#         )
#       }
#     }, error = function(e) {
#       
#     }, warning = function(w) {
#       
#     })
#   }
# }