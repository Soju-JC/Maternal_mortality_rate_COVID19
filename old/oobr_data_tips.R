library(tidyverse)
library(httr)
library(janitor)
library(getPass)
library(repr)
library(data.table)
library(readr)

token = getPass()  #Token de acesso à API da PCDaS (todos os arquivos gerados se encontram na pasta "Databases", no Google Drive)

url_base = "https://bigdata-api.fiocruz.br"

convertRequestToDF <- function(request){
  variables = unlist(content(request)$columns)
  variables = variables[names(variables) == "name"]
  column_names <- unname(variables)
  values = content(request)$rows
  df <- as.data.frame(do.call(rbind,lapply(values,function(r) {
    row <- r
    row[sapply(row, is.null)] <- NA
    rbind(unlist(row))
  } )))
  names(df) <- column_names
  return(df)
}

estados <- c('RO','AC','AM','RR','PA','AP','TO','MA','PI','CE','RN','PB','PE','AL','SE','BA','MG','ES','RJ','SP','PR','SC','RS','MS','MT','GO','DF')
endpoint <- paste0(url_base,"/","sql_query")

# Número de nascidos vivos
df_nascidos_vivos <- df_nascidos_vivos_aux <- data.frame()

for (estado in estados){
  
  params = paste0('{
      "token": {
        "token": "',token,'"
      },
      "sql": {
        "sql": {"query": "SELECT res_REGIAO, res_SIGLA_UF, ano_nasc, COUNT(1)',
                  ' FROM \\"datasus-sinasc_final_1996-2020_preliminar_2021_2022\\"',
                  ' WHERE (res_SIGLA_UF = \'',estado,'\' AND ano_nasc BETWEEN 2012 AND 2020)',
                  ' GROUP BY res_REGIAO, res_SIGLA_UF, ano_nasc",
                        "fetch_size": 65000}
      }
    }')
  
  request <- POST(url = endpoint, body = params, encode = "form")
  df_nascidos_vivos_aux <- convertRequestToDF(request)
  names(df_nascidos_vivos_aux) <- c('regiao', 'uf', 'ano', 'nascidos')
  df_nascidos_vivos <- rbind(df_nascidos_vivos, df_nascidos_vivos_aux)
  
  repeat {
    
    cursor <- content(request)$cursor
    
    params = paste0('{
          "token": {
            "token": "',token,'"
          },
          "sql": {
            "sql": {"query":"SELECT res_REGIAO, res_SIGLA_UF, ano_nasc, COUNT(1)',
                    ' FROM \\"datasus-sinasc_final_1996-2020_preliminar_2021_2022\\"',
                    ' WHERE (res_SIGLA_UF = \'',estado,'\' AND ano_nasc BETWEEN 2012 AND 2020)',
                    ' GROUP BY res_REGIAO, res_SIGLA_UF, ano_nasc",
                           "fetch_size": 65000, "cursor": "',cursor,'"}
          }
        }')
    
    
    request <- POST(url = endpoint, body = params, encode = "form")
    
    if (length(content(request)$rows) == 0)
      break
    else print("oi")
    
    df_nascidos_vivos_aux <- convertRequestToDF(request)
    names(df_nascidos_vivos_aux) <- c('regiao', 'uf', 'ano', 'nascidos')
    df_nascidos_vivos <- rbind(df_nascidos_vivos, df_nascidos_vivos_aux)
  }
}

head(df_nascidos_vivos)

##Óbitos maternos oficiais dos anos de 1996 a 2022
df_obitos_maternos_aux <- dataframe <- data.frame()

for (estado in estados){

  params = paste0('{
      "token": {
        "token": "',token,'"
      },
      "sql": {
        "sql": {"query":" SELECT res_REGIAO, res_SIGLA_UF, res_MUNNOME, res_codigo_adotado, ano_obito, CAUSABAS, causabas_capitulo, causabas_categoria, OBITOGRAV, OBITOPUERP, def_raca_cor, idade_obito_anos, FONTEINV, COUNT(1)',
                        ' FROM \\"datasus-sim_final_1996-2020_preliminar_2021_2022\\"',
                        ' WHERE (res_SIGLA_UF = \'',estado,'\' AND SEXO = 2 AND',
                              ' ((CAUSABAS >= \'O000\'  AND  CAUSABAS <= \'O959\') OR',
                              ' (CAUSABAS >= \'O980\'  AND  CAUSABAS <= \'O999\') OR',
                              ' (CAUSABAS = \'A34\' AND OBITOPUERP != 2) OR',
                              ' ((CAUSABAS >= \'B200\'  AND  CAUSABAS <= \'B249\') AND (OBITOGRAV = 1 OR OBITOPUERP = 1)) OR',
                              ' (CAUSABAS = \'D392\' AND (OBITOGRAV = 1 OR OBITOPUERP = 1)) OR',
                              ' (CAUSABAS = \'E230\' AND (OBITOGRAV = 1 OR OBITOPUERP = 1)) OR',
                              ' ((CAUSABAS >= \'F530\'  AND  CAUSABAS <= \'F539\') AND (OBITOPUERP != 2 OR OBITOPUERP = \' \')) OR',
                              ' (CAUSABAS = \'M830\' AND OBITOPUERP != 2)))',
                        ' GROUP BY res_REGIAO, res_SIGLA_UF, res_MUNNOME, res_codigo_adotado, ano_obito, CAUSABAS, causabas_capitulo, causabas_categoria, OBITOGRAV, OBITOPUERP, def_raca_cor, idade_obito_anos, FONTEINV",
                        "fetch_size": 65000}
      }
    }')

  request <- POST(url = endpoint, body = params, encode = "form")
  dataframe <- convertRequestToDF(request)
  names(dataframe) <- c("regiao", "uf", "municipio", "codigo", "ano", "causabas", "capitulo_cid10", "causabas_categoria", "obitograv", "obitopuerp", "racacor", "idade", "fonteinv", "obitos")
  df_obitos_maternos_aux <- rbind(df_obitos_maternos_aux, dataframe)

  repeat {

    cursor <- content(request)$cursor

    params = paste0('{
          "token": {
            "token": "',token,'"
          },
          "sql": {
            "sql": {"query":" SELECT res_REGIAO, res_SIGLA_UF, res_MUNNOME, res_codigo_adotado, ano_obito, CAUSABAS, causabas_capitulo, causabas_categoria, OBITOGRAV, OBITOPUERP, def_raca_cor, idade_obito_anos, FONTEINV, COUNT(1)',
                            ' FROM \\"datasus-sim_final_1996-2020_preliminar_2021_2022\\"',
                            ' WHERE (res_SIGLA_UF = \'',estado,'\' AND SEXO = 2 AND',
                                  ' ((CAUSABAS >= \'O000\'  AND  CAUSABAS <= \'O959\') OR',
                                  ' (CAUSABAS >= \'O980\'  AND  CAUSABAS <= \'O999\') OR',
                                  ' (CAUSABAS = \'A34\' AND OBITOPUERP != 2) OR',
                                  ' ((CAUSABAS >= \'B200\'  AND  CAUSABAS <= \'B249\') AND (OBITOGRAV = 1 OR OBITOPUERP = 1)) OR',
                                  ' (CAUSABAS = \'D392\' AND (OBITOGRAV = 1 OR OBITOPUERP = 1)) OR',
                                  ' (CAUSABAS = \'E230\' AND (OBITOGRAV = 1 OR OBITOPUERP = 1)) OR',
                                  ' ((CAUSABAS >= \'F530\'  AND  CAUSABAS <= \'F539\') AND (OBITOPUERP != 2 OR OBITOPUERP = \' \')) OR',
                                  ' (CAUSABAS = \'M830\' AND OBITOPUERP != 2)))',
                            ' GROUP BY res_REGIAO, res_SIGLA_UF, res_MUNNOME, res_codigo_adotado, ano_obito, CAUSABAS, causabas_capitulo, causabas_categoria, OBITOGRAV, OBITOPUERP, def_raca_cor, idade_obito_anos, FONTEINV",
                            "fetch_size": 65000, "cursor": "',cursor,'"}
                            }
                    }')


    request <- POST(url = endpoint, body = params, encode = "form")

    if (length(content(request)$rows) == 0)
      break
    else print("oi")

    request <- POST(url = endpoint, body = params, encode = "form")
    dataframe <- convertRequestToDF(request)
    names(dataframe) <- c("regiao", "uf", "municipio", "codigo", "ano", "causabas", "capitulo_cid10", "causabas_categoria", "obitograv", "obitopuerp", "racacor", "idade", "fonteinv", "obitos")
    df_obitos_maternos_aux <- rbind(df_obitos_maternos_aux, dataframe)
  }
}

head(df_obitos_maternos_aux)

df_obitos_maternos <- df_obitos_maternos_aux |>
  mutate(
    tipo_de_morte_materna = if_else(
      condition = (causabas >= "B200" & causabas <= "B249") |
        (causabas >= "O100" & causabas <= "O109") |
        ((causabas >= "O240" & causabas != "O244") & causabas <= "O259") |
        (causabas == "O94") |
        (causabas >= "O980" & causabas <= "O999"),
      true = "Indireta",
      false = if_else(causabas == "O95", true = "Não especificada", false = "Direta")
    ),
    periodo_do_obito = case_when(
      obitograv == "1" & obitopuerp != "1" & obitopuerp != "2" ~ "Durante a gravidez, parto ou aborto",
      obitograv != "1" & obitopuerp == "1" ~ "Durante o puerpério, até 42 dias",
      obitograv != "1" & obitopuerp == "2" ~ "Durante o puerpério, de 43 dias a menos de 1 ano",
      (obitograv == "2" & obitopuerp == "3") | (obitograv == "2" & obitopuerp == "9") | (obitograv == "9" & obitopuerp == "3")  ~ "Não na gravidez ou no puerpério",
      #obitograv == "2" & obitopuerp == "9" ~ "Durante o puerpério, até 1 ano, período não discriminado",
      obitograv == "9" & obitopuerp == "9" ~ "Não informado ou ignorado",
      (obitograv == "1" & obitopuerp == "1") | (obitograv == "1" & obitopuerp == "2") ~ "Período inconsistente"),
    .after = causabas_categoria,
    investigacao_cmm = if_else(fonteinv == "1", true = "Sim", false = if_else(fonteinv == "9", true = "Sem informação", false = "Não",  missing = "Sem informação"), missing = "Sem informação")
  ) |>
  select(!c(codigo, causabas, obitograv, obitopuerp, fonteinv)) |>
  group_by(across(!obitos)) |>
  summarise(obitos = sum(as.numeric(obitos))) |>
  ungroup()

df_obitos_maternos$idade <- as.numeric(df_obitos_maternos$idade)












