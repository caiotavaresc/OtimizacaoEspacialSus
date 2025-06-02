### Step 1: Read the 2022 Hospitalization Data
### Objective: Reconstruct Pereira (2022) Performance Indicators

setwd("D:/Github/OtimizacaoEspacialSus/definition_of_objective_functions")

# Package to read DBF files
install.packages("sf")
library(sf)

# Package to summarize procedures
install.packages("dplyr")
library(dplyr)

# Package to calculate mutual information
install.packages("mpmi")
library(mpmi)

arq_aih_sum <- data.frame()
directories = list("hosp_file(1)","hosp_file(2)","hosp_file(3)","hosp_file(4)","hosp_file(5)","hosp_file(6)","hosp_file(7)")

for (dir in diretorios){
  diretorio = paste0("./input_data/hospitalization_data/", dir)
  
  # List all .dbf files in the directory
  arquivos_dbf <- list.files(path = diretorio, pattern = "\\.dbf$", full.names = TRUE)
  
  arq_aih <- data.frame()
  for (arquivo in arquivos_dbf){
    temp_arq = st_read(arquivo)
    
    # Select only relevant columns
    temp_arq <- temp_arq[, c("MUNIC_RES", "MUNIC_MOV", "PROC_REA", "VAL_TOT")]
    
    # Append data from each file to the full dataset
    arq_aih <- rbind(arq_aih, temp_arq)
  }
  
  # Group and aggregate data
  arq_aih_sum_temp <- data.frame()
  arq_aih_sum_temp <- arq_aih %>%
    group_by(MUNIC_RES, MUNIC_MOV, PROC_REA) %>%
    summarise(
      count = n(),  # Count the number of records (equivalent to COUNT(1))
      sum_val_tot = sum(VAL_TOT, na.rm = TRUE)  # Sum of VAL_TOT values (SUM(VAL_TOT))
    )
  
  # Append the summary of this folder to the full summary dataset
  arq_aih_sum <- rbind(arq_aih_sum, arq_aih_sum_temp)
}

# Read the municipality x health region mapping table
mun_regsau <- read.csv2("./input_data/territorial_data/municipalities_and_health_regions.csv", sep=",", header=TRUE, colClasses="character")

# Fill in the origin and destination health regions
# Perform left join and rename the resulting column to CO_REGSAUD_RES
arq_aih_sum <- left_join(arq_aih_sum, 
                         mun_regsau[, c("CO_MUNICIP", "CO_REGSAUD")], 
                         by = c("MUNIC_RES" = "CO_MUNICIP")) %>%
  rename(CO_REGSAUD_RES = CO_REGSAUD)

# Perform left join and rename the resulting column to CO_REGSAUD_MOV
arq_aih_sum <- left_join(arq_aih_sum, 
                         mun_regsau[, c("CO_MUNICIP", "CO_REGSAUD")], 
                         by = c("MUNIC_MOV" = "CO_MUNICIP")) %>%
  rename(CO_REGSAUD_MOV = CO_REGSAUD)

# Save the arq_aih_sum results to a file for future reference
write.csv(arq_aih_sum, file="arq_aih_sum.csv", row.names=FALSE)

# Dear reader, the lines above demonstrated how to reproduce the calculation of hospitalization indicators by health region.
# Since the full dataset exceeds 8 gigabytes, we simplify the process here by directly reading the already consolidated file: "arq_aih_sum.csv"

arq_aih_sum <- read.csv2("./output_data/arq_aih_sum.csv", sep=",", header=TRUE, colClasses="character")
arq_aih_sum$count <- as.numeric(arq_aih_sum$count)
arq_aih_sum$sum_val_tot <- as.numeric(arq_aih_sum$sum_val_tot)

# Group and aggregate the data
arq_aih_sum_regsau <- arq_aih_sum %>%
  group_by(CO_REGSAUD_RES, CO_REGSAUD_MOV, PROC_REA) %>%
  summarise(
    sum_count = sum(count, na.rm = TRUE),          # Sum of the 'count' column
    sum_sum_val_tot = sum(sum_val_tot, na.rm = TRUE)  # Sum of the 'sum_val_tot' column
  )

# Count the number of distinct procedures per health region
arq_aih_proc_por_regsau <- arq_aih_sum_regsau %>%
  group_by(CO_REGSAUD_MOV) %>%
  summarise(
    num_procedimentos_distintos = n_distinct(PROC_REA),
    total_procedimentos = sum(sum_count, na.rm = TRUE),
    total_valor = sum(sum_sum_val_tot, na.rm = TRUE)
  )

# Calculate the LIFO and LOFI coefficients by health region
arq_atendimentos_sum_regsau <- arq_aih_sum %>%
  group_by(CO_REGSAUD_RES, CO_REGSAUD_MOV) %>%
  summarise(
    sum_count = sum(count, na.rm = TRUE)          # Sum of the 'count' column
  )

# Create the new dataframe arq_atendimentos_internos
arq_atendimentos_internos <- arq_atendimentos_sum_regsau %>%
  filter(CO_REGSAUD_RES == CO_REGSAUD_MOV) %>%  # Filter where CO_REGSAUD_RES == CO_REGSAUD_MOV
  select(COD_REGSAUD = CO_REGSAUD_RES, COUNT_ATEND_MESMA_REG = sum_count)  # Rename columns

# Create the new dataframe arq_atendimentos_vindos_de_fora
arq_atendimentos_vindos_de_fora <- arq_atendimentos_sum_regsau %>%
  filter(CO_REGSAUD_RES != CO_REGSAUD_MOV) %>%  # Filter where CO_REGSAUD_RES != CO_REGSAUD_MOV
  group_by(CO_REGSAUD_MOV) %>%  # Group by CO_REGSAUD_MOV
  summarise(COUNT_ATEND_VIERAM_DE_FORA = sum(sum_count, na.rm = TRUE)) %>%  # Sum the sum_count column
  rename(COD_REGSAUD = CO_REGSAUD_MOV)  # Rename column to COD_REGSAUD

# Create the new dataframe arq_atendimentos_foram_pra_fora
arq_atendimentos_foram_pra_fora <- arq_atendimentos_sum_regsau %>%
  filter(CO_REGSAUD_RES != CO_REGSAUD_MOV) %>%  # Filter where CO_REGSAUD_RES != CO_REGSAUD_MOV
  group_by(CO_REGSAUD_RES) %>%  # Group by CO_REGSAUD_RES
  summarise(COUNT_ATEND_FORAM_PRA_FORA = sum(sum_count, na.rm = TRUE)) %>%  # Sum the sum_count column
  rename(COD_REGSAUD = CO_REGSAUD_RES)  # Rename column to COD_REGSAUD

# Create the dataframe arq_indicadores_regsaud from the three dataframes
arq_indicadores_regsaud <- arq_atendimentos_internos %>%
  left_join(arq_atendimentos_vindos_de_fora, by = "COD_REGSAUD") %>%
  left_join(arq_atendimentos_foram_pra_fora, by = "COD_REGSAUD")

# Calculate the LIFO and LOFI coefficients
arq_indicadores_regsaud <- arq_indicadores_regsaud %>%
  mutate(
    LIFO = COUNT_ATEND_MESMA_REG / (COUNT_ATEND_MESMA_REG + COUNT_ATEND_VIERAM_DE_FORA),
    LOFI = COUNT_ATEND_MESMA_REG / (COUNT_ATEND_MESMA_REG + COUNT_ATEND_FORAM_PRA_FORA)
  )

# Create the columns TX_PERMANENCIA and TX_ATRACAO in the dataframe arq_indicadores_regsaud
arq_indicadores_regsaud <- arq_indicadores_regsaud %>%
  mutate(
    TX_PERMANENCIA = LOFI,              # TX_PERMANENCIA is equal to LOFI
    TX_ATRACAO = 1 - LIFO               # TX_ATRACAO is equal to 1 minus LIFO
  )

# Create the column DESEMPENHO in the dataframe arq_indicadores_regsaud
arq_indicadores_regsaud <- arq_indicadores_regsaud %>%
  mutate(
    DESEMPENHO = ((0.2 + 1) * TX_ATRACAO * TX_PERMANENCIA) / ((0.2) * TX_PERMANENCIA + 1 * TX_ATRACAO)
  )

### Step 2: Work with the 2022 Proadess indicators
### Objective: Define variables of interest

indicadores_proadess <- read.csv2("./input_data/indicadores_estrutura_proadess.csv", sep=";", header=TRUE, row.names = 1)

# 2.1 Removing very sparse columns

# Function to calculate the number of zeros and the percentage of zeros in each variable
calc_zero_info <- function(x) {
  zeros_count <- sum(x == 0)  # Count of zeros
  percent_zeros <- (zeros_count / length(x)) * 100  # Percentage of zeros
  return(c(zeros_count, percent_zeros))
}

# Applying the function to each column of the dataframe
zero_info <- t(apply(indicadores_proadess, 2, calc_zero_info))

# Converting to a dataframe for easier visualization
zero_info_df <- as.data.frame(zero_info)
colnames(zero_info_df) <- c("Zero_Count", "Percent_Zero")

# Adding variable names to the dataframe
zero_info_df$Variable <- rownames(zero_info_df)

# Sorting by number of zeros (from most sparse to least sparse)
zero_info_df <- zero_info_df[order(-zero_info_df$Zero_Count), ]

# Filtering variables with Percent_Zero > 75%
filtered_zero_info_df <- zero_info_df[zero_info_df$Percent_Zero > 75, ]

# Remove columns with high sparsity
columns_to_remove <- filtered_zero_info_df$Variable
indicadores_atu <- indicadores_proadess[, !(colnames(indicadores_proadess) %in% columns_to_remove)]

# 2.2 Obtain the performance index for each health region

# Add the DESEMPENHO column to the dataframe indicadores_atu
indicadores_atu <- indicadores_atu %>%
  mutate(
    DESEMPENHO = arq_indicadores_regsaud$DESEMPENHO[match(rownames(indicadores_atu), arq_indicadores_regsaud$COD_REGSAUD)]
  )
write.csv(indicadores_atu,"output_data/indicadores_regsaud_enriquecidos_desempenho.csv")

### Step 3: Calculate MRMR variables (Maximum Relevance, Minimum Redundancy)

mi_matrix <- cmi(indicadores_atu)
mi_result <- mi_matrix$mi
# Assigns the variable names of the dataframe to the rows and columns of the MI matrix
colnames(mi_result) <- rownames(mi_result) <- colnames(indicadores_atu)

mi_result <- mi_result[!rownames(mi_result) %in% "DESEMPENHO", ]

# Calculate relevance and redundancy for all pairwise combinations
calcular_relevancia <- function(mi_result, conjunto_s, target) {
  # Checks if all names are present in the matrix
  if (!(target %in% colnames(mi_result))) {
    stop("The target variable is not present in the mi_result matrix.")
  }
  if (!all(conjunto_s %in% colnames(mi_result))) {
    stop("One or more variables from conjunto_s are not present in the mi_result matrix.")
  }
  
  # Calculate relevance: MI between each variable in conjunto_s and the target
  inf_mutuas <- sapply(conjunto_s, function(var) {
    mi_result[var, target]
  })
  
  relevancia <- sum(inf_mutuas) / length(conjunto_s)
  return(relevancia)
}

calcular_redundancia <- function(mi_result, conjunto_s) {
  n <- length(conjunto_s)
  
  if (n == 0) return(0)
  
  if (!all(conjunto_s %in% colnames(mi_result))) {
    stop("One or more variables from conjunto_s are not present in the mi_result matrix.")
  }
  
  # Sum of mutual information between all pairs (including i == j)
  soma_redundancia <- sum(sapply(conjunto_s, function(i) {
    sum(sapply(conjunto_s, function(j) {
      mi_result[i, j]
    }))
  }))
  
  # Average with n^2 in the denominator
  redundancia_media <- soma_redundancia / (n^2)
  
  return(redundancia_media)
}

# 4.1 Identify candidate variables (excluding the target)
todas_vars <- colnames(mi_result)
target <- "DESEMPENHO"
variaveis_candidatas <- setdiff(todas_vars, target)

# 4.2 Generate all unique pairs (combinations of 2)
pares <- combn(variaveis_candidatas, 2, simplify = FALSE)

# 4.3 Initialize results list
resultados <- lapply(pares, function(par) {
  # Relevance: average of MI between each variable in the pair and the target
  relevancia <- calcular_relevancia(mi_result, par, target)
  
  # Redundancy: average of MI between the two variables in the pair
  redundancia <- calcular_redundancia(mi_result, par)
  
  # MID: difference between relevance and redundancy
  mid <- relevancia - redundancia
  
  # MIQ: quotient
  miq <- ifelse(redundancia == 0, NA, relevancia / redundancia)
  
  # Return a row of data
  data.frame(
    var1 = par[1],
    var2 = par[2],
    relevancia = relevancia,
    redundancia = redundancia,
    MID = mid,
    MIQ = miq
  )
})

# 4. Combine everything into a single data frame
df_resultados <- do.call(rbind, resultados)
df_resultados <- df_resultados %>%
  arrange(desc(MIQ))

write.csv(df_resultados,"output_data/resultado_analise_MRMR.csv", row.names=FALSE)
