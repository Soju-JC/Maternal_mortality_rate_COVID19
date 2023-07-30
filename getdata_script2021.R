# Load packages
library("tidyverse")
library("data.table")
library("forecast")
library("tseries")

################################################################################

# Load data (SINASC 2021)
# https://opendatasus.saude.gov.br/dataset/sistema-de-informacao-sobre-nascidos-vivos-sinasc/resource/16b72608-a6f0-45f6-8dee-75925ad195ce
df_births <- fread('SINASC_2021.csv', sep = ";")

# Constrain dates to be 8 digits long
df_births$DTNASC <- str_pad(df_births$DTNASC, width = 8, 
                                 side = "left", 
                                 pad = "0")
# Separate the day, month and year information correctly
df_births$DTNASC <- gsub("^(..)(..)(....)$", "\\1/\\2/\\3", 
                              df_births$DTNASC)
# Convert column to date
df_births$DTNASC <- as.Date(df_births$DTNASC,
                                 format = "%d/%m/%Y")

# Groups births for each day of the year 2021
df_births_summary <- df_births %>% 
  group_by(DTNASC) %>%
  summarize(count = n())

# Graph of the number of births per day throughout 2021
ggplot(df_births_summary, 
       aes(x = DTNASC, y = count)) + 
  geom_line()

################################################################################

# Load data (SIM 2021)
# https://opendatasus.saude.gov.br/dataset/sim/resource/2dee2fc9-d19f-41bc-b058-6d4154f126ab
df_deaths <- fread('Mortalidade_Geral_2021.csv', sep = ";")

# Maternal deaths (matches the 3.030 deaths officially reported for 2021)
df_maternal_deaths <- df_deaths %>%
  filter(SEXO == 2 & 
           (
             (CAUSABAS >= "O000" & CAUSABAS <= "O959") | 
             (CAUSABAS >= "O980" & CAUSABAS <= "O999") |
             (CAUSABAS == "A34" & OBITOPUERP != 2) |
             (
               (CAUSABAS >= "B200" & CAUSABAS <= "B249") & 
               (OBITOGRAV == 1 | OBITOPUERP == 1)
               ) |
             (CAUSABAS == "D392" & (OBITOGRAV == 1 | OBITOPUERP == 1)) |
             (CAUSABAS == "E230" & (OBITOGRAV == 1 | OBITOPUERP == 1)) |
             (
               (CAUSABAS >= "F530" & CAUSABAS <= "F539") & 
               (OBITOPUERP != 2 | OBITOPUERP == "")
               ) |
             (CAUSABAS == "M830" & OBITOPUERP != 2)
             )
         )

# Constrain dates to be 8 digits long
df_maternal_deaths$DTOBITO <- str_pad(df_maternal_deaths$DTOBITO, width = 8, 
                            side = "left", 
                            pad = "0")
# Separate the day, month and year information correctly
df_maternal_deaths$DTOBITO <- gsub("^(..)(..)(....)$", "\\1/\\2/\\3", 
                                   df_maternal_deaths$DTOBITO)
# Convert column to date
df_maternal_deaths$DTOBITO <- as.Date(df_maternal_deaths$DTOBITO, 
                                      format = "%d/%m/%Y")

# Groups maternal deaths for each day of the year 2021
df_maternal_deaths_summary <- df_maternal_deaths %>% 
  group_by(DTOBITO) %>%
  summarize(count = n())

# Graph of the number of maternal deaths per day throughout 2021
ggplot(df_maternal_deaths_summary, 
       aes(x = DTOBITO, y = count)) + 
  geom_line()

################################################################################
## df_maternal_deaths_summary is missing one row. Let's fix it!

# Standardize column names
names(df_births_summary) <- c("date", "count")
names(df_maternal_deaths_summary) <- c("date", "count")

# Find row in df_births_summary that is not in df_maternal_deaths_summary
missing_row <- anti_join(df_births_summary, 
                         df_maternal_deaths_summary, 
                         by = "date")

# Find position to insert missing row
position <- which(df_births_summary$date == missing_row$date)

# Insert missing row at specified position
df_maternal_deaths_summary <- df_maternal_deaths_summary %>%
  slice(1:(position - 1)) %>%
  bind_rows(missing_row) %>%
  bind_rows(slice(df_maternal_deaths_summary, (position:n())))

# NA value to the specified position in the df_maternal_deaths_summary
df_maternal_deaths_summary[position, "count"] <- NA

############## Add a value to the specified position using ARIMA ###############

# Slice df_maternal_deaths_summary up to the missing value
df_maternal_deaths_imput <- df_maternal_deaths_summary %>%
  slice(1:(position - 1))

# Transform the data to time series object (day frequency)
df_maternal_deaths_imput_ts <- ts(df_maternal_deaths_imput$count, 
                                  frequency = 365)

# Apply difference transformation to make data stationary
df_maternal_deaths_imput_diff <- diff(df_maternal_deaths_imput_ts)

# Check out autocorrelation function
acf(df_maternal_deaths_imput_ts) # Before difference transformation
acf(df_maternal_deaths_imput_diff) # After difference transformation

# Perform ADF test (indicates stationary data)
adf.test(df_maternal_deaths_imput_diff)

# Perform KPSS test (also indicates stationary data)
kpss.test(df_maternal_deaths_imput_diff)

# Fit ARIMA model (difference order = 1)
model <- auto.arima(df_maternal_deaths_imput_ts, d = 1)

# Predict count for missing date
predicted_count <- forecast(model, h = 1)$mean

# Add the predicted value to the specified position
df_maternal_deaths_summary[position, "count"] <- as.numeric(
  round(predicted_count)
  )

################################################################################
###### Create the maternal mortality rate time series using (NOM/NNV) x 1 ######
################################################################################

df_mr <- df_maternal_deaths_summary[,"date"]
df_mr$rate <- df_maternal_deaths_summary$count/df_births_summary$count

# Graph of the maternal rate per day throughout 2021
ggplot(df_mr, 
       aes(x = date, y = rate)) + 
  geom_line()

# Save the time series as a RDS file
saveRDS(df_mr, file = "maternal_mortality_rate_2021.rds")