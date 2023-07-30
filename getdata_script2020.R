# Load packages
library("tidyverse")
library("data.table")
library("forecast")
library("tseries")
library("imputeTS")

################################################################################

# Load data (SINASC 2020)
# https://opendatasus.saude.gov.br/dataset/sistema-de-informacao-sobre-nascidos-vivos-sinasc/resource/7a139c1b-a4b0-4967-a94c-33a5dfae667b
df_births <- fread('SINASC_2020.csv', sep = ";")

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

# Groups births for each day of the year 2020
df_births_summary <- df_births %>% 
  group_by(DTNASC) %>%
  summarize(count = n())

# Graph of the number of births per day throughout 2020
ggplot(df_births_summary, 
       aes(x = DTNASC, y = count)) + 
  geom_line()

################################################################################

# Load data (SIM 2020)
# https://opendatasus.saude.gov.br/dataset/sim/resource/c622b337-a522-4243-bf19-6c971e809cff
df_deaths <- fread('Mortalidade_Geral_2020.csv', sep = ";")

# Maternal deaths (matches the 1965 deaths officially reported for 2020)
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

# Groups maternal deaths for each day of the year 2020
df_maternal_deaths_summary <- df_maternal_deaths %>% 
  group_by(DTOBITO) %>%
  summarize(count = n())

# Graph of the number of maternal deaths per day throughout 2020
ggplot(df_maternal_deaths_summary, 
       aes(x = DTOBITO, y = count)) + 
  geom_line()

################################################################################
## df_maternal_deaths_summary is missing 3 rows. Let's fix it!

# Standardize column names
names(df_births_summary) <- c("date", "count")
names(df_maternal_deaths_summary) <- c("date", "count")

# Find rows in df_births_summary that are not in df_maternal_deaths_summary
missing_row <- anti_join(df_births_summary, 
                         df_maternal_deaths_summary, 
                         by = "date")

# Find positions to insert missing rows
position1 <- which(df_births_summary$date == missing_row$date[1])
position2 <- which(df_births_summary$date == missing_row$date[2])
position3 <- which(df_births_summary$date == missing_row$date[3])

# Insert missing rows at specified positions
df_maternal_deaths_summary <- df_maternal_deaths_summary %>%
  slice(1:(position1 - 1)) %>%
  bind_rows(missing_row[1,]) %>%
  bind_rows(slice(df_maternal_deaths_summary, (position1:n())))

df_maternal_deaths_summary <- df_maternal_deaths_summary %>%
  slice(1:(position2 - 1)) %>%
  bind_rows(missing_row[2,]) %>%
  bind_rows(slice(df_maternal_deaths_summary, (position2:n())))

df_maternal_deaths_summary <- df_maternal_deaths_summary %>%
  slice(1:(position3 - 1)) %>%
  bind_rows(missing_row[3,]) %>%
  bind_rows(slice(df_maternal_deaths_summary, (position3:n())))

# NA value to the specified positions in the df_maternal_deaths_summary
df_maternal_deaths_summary[position1, "count"] <- NA
df_maternal_deaths_summary[position2, "count"] <- NA
df_maternal_deaths_summary[position3, "count"] <- NA

######### Add values to the specified positions using moving average #########
aux <- na_ma(df_maternal_deaths_summary$count, weighting = "exponential")
df_maternal_deaths_summary$count <- round(aux)

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
saveRDS(df_mr, file = "maternal_mortality_rate_2020.rds")