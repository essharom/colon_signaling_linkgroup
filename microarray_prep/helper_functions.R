#helper_functions

#libraries
library("tidyverse")

library("affy")
library("biomaRt")

# sigmoid <- function(x, m, b){
#   1 * (exp(m * x + b)/(1 + exp(m * x + b)))
# }
# 
# breakdown <- function(x){
#   x <- sigmoid(log(x), log(2, base = 10), log(1/2))
#   return(x)
# }

ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl)

run_condition <- function(condition_name, df){
  print(condition_name)
  condition_df <- df %>% 
    filter(condition == condition_name)
  
  platforms <- unique(condition_df$platform)
  
  result <- map_dfr(platforms, run_platform, condition_df) %>% 
    group_by(uniprotswissprot) %>% 
    summarise(expression = median(expression)) %>%  # TODO aggregating probabilities, possibly some better function exists
    mutate(condition = condition_name)
  
  return(result)
}

run_platform <- function(platform_name, df){
  print(platform_name)
  
  platform_samples_df <- df %>% 
    filter(platform == platform_name)
  
  platform_fun <- get_platform_function(platform_name)
  
  result <- platform_fun(platform_samples_df$file_path) %>%
    mutate(expression = scales::rescale(expression, to = c(0, 1)))
  
  return(result)
}

get_platform_function <- function(platform_name){
  if (platform_name == 'GPL570'){
    return (GPL570_read)
  }
  
  if (platform_name == 'GPL571'){
    return (GPL571_read)
  }
  
}


# TODO function to summarize by sample (median) and by protein (sum) instead of doing it at every platform
sample_agg <- function(df){
  df <- df %>% 
    group_by(sample_id, uniprotswissprot) %>% 
    summarise(expression = median(expression)) %>% 
    ungroup()
  
}


#GPL570 
#library("hgu133plus2.db") #GPL570 platform

#file_list = sample_df$file_path[1:2]

GPL570_read <- function(file_list){
  raw_data = ReadAffy(verbose=TRUE, filenames=file_list, cdfname="HG-U133_Plus_2") # read the row data to raw.data
  rma_df = rma(raw_data) # mas5(); normalise the raw data 
  expression_df = exprs(rma_df)	%>% 
    as.data.frame() # get the normalised data in a readable matrix
  
  probes = expression_df %>% 
    mutate(probe_id = row.names(expression_df))

  annot <- getBM(attributes=c('affy_hg_u133_plus_2','uniprotswissprot'), mart = ensembl) %>% 
    dplyr::rename(probe_id = affy_hg_u133_plus_2)
  
  affy <- annot %>% 
    filter(uniprotswissprot!="",
           probe_id!="") %>% 
    left_join(probes) %>%
    dplyr::select(!probe_id) %>% 
    group_by(uniprotswissprot) %>% 
    summarise_all(sum) %>%
    ungroup() %>% 
    pivot_longer(cols = !uniprotswissprot, names_to = "sample_id", values_to = "expression") %>% 
    group_by(uniprotswissprot) %>% 
    summarise(expression = median(expression))
}


#GPL571 
library("hgu133a2.db") #GPL571 platform

GPL571_read <- function(file_list){
  raw_data = ReadAffy(verbose=TRUE, filenames=file_list, cdfname="HG-U133A_2", compress = TRUE) # read the row data to raw.data
  rma_df = rma(raw_data) # mas5(); normalise the raw data 
  expression_df = exprs(rma_df)	%>% 
    as.data.frame() # get the normalised data in a readable matrix
  
  probes = expression_df %>% 
    mutate(probe_id = row.names(expression_df))
  
  annot <- getBM(attributes=c('affy_hg_u133a_2','uniprotswissprot'), mart = ensembl) %>% 
    dplyr::rename(probe_id = affy_hg_u133a_2)
  
  affy <- annot %>% 
    filter(uniprotswissprot!="",
           probe_id!="") %>% 
    left_join(probes) %>%
    dplyr::select(!probe_id) %>% 
    group_by(uniprotswissprot) %>% 
    summarise_all(sum) %>%
    ungroup() %>% 
    pivot_longer(cols = !uniprotswissprot, names_to = "sample_id", values_to = "expression") %>% 
    group_by(uniprotswissprot) %>% 
    summarise(expression = median(expression))
  
}


# <platform_name>_read <- function(file_list){ 
#
#
# output: dataframe of <uniprotswissprot> and <expression> values
# }