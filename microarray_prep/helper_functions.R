#helper_functions

#libraries
library("tidyverse")

library("biomaRt")
library(filesstrings)
library(stringr)

# pre-processing packages
library("oligo")
library("affy")

# sigmoid <- function(x, m, b){
#   1 * (exp(m * x + b)/(1 + exp(m * x + b)))
# }
# 
# breakdown <- function(x){
#   x <- sigmoid(log(x), log(2, base = 10), log(1/2))
#   return(x)
# }

data_download <- function(cancer_type, GSE_ID, GPL_ID, phenotype, pheno_data, pheno_term, cell_files){
  if(length(pheno_term[[1]]) != 0){
    files = ""
    for(j in 1:length(pheno_term[1])){
      files = append(files, cell_files[grep(pheno_term[[j]], pheno_data[['title']])])
    }
    dir.create(paste("./Data/", cancer_type, phenotype,  GPL_ID, sep = ""))
    dir.create(paste("./Data/", cancer_type, phenotype,  GPL_ID, "/", GSE_ID, sep = ""))
    file.copy(paste("./myData/", files[-1], sep =""), 
              paste("./Data/", cancer_type, phenotype, GPL_ID, "/", GSE_ID, "/", files[-1], sep = ""))
  }
}

ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl)


# calculate abundances for a specific condition
run_condition <- function(condition_name, df){
  print(condition_name)
  condition_df <- df %>% 
    filter(condition == condition_name)
  
  platforms <- unique(condition_df$platform)
  
  # TODO median aggregation seems like a good option, but for aggregating probabilities, possibly some better function exists
  result <- map_dfr(platforms, run_platform, condition_df) %>% 
    group_by(uniprotswissprot) %>% 
    summarise(expression = median(expression)) %>%  
    mutate(condition = condition_name)
  
  return(result)
}


# calculate abundances for a platform
run_platform <- function(platform_name, df){
  print(platform_name)
  
  platform_samples_df <- df %>% 
    filter(platform == platform_name)
  
  platform_fun <- get_platform_function(platform_name)
  
  # scale all values on (0 to 1) scale by scales::rescale function
  result <- platform_fun(platform_samples_df$file_path) %>%
    mutate(expression = scales::rescale(expression, to = c(0, 1)))
  
  return(result)
}


# map the correct platform specific function to the CEL file analysis
get_platform_function <- function(platform_name){
  if (platform_name == 'GPL570'){
    return (GPL570_read)
  }
  
  if (platform_name == 'GPL571'){
    return (GPL571_read)
  }
  
  if (platform_name == 'GPL17692'){
    return (GPL17692_read)
  }
  
  if (platform_name == 'GPL6244'){
    return (GPL6244_read)
  }
  
  if (platform_name == 'GPL96'){
    return (GPL96_read)
  }
  
  if (platform_name == 'GPL97'){
    return (GPL97_read)
  }
  
}



#
# FUNCTION blueprint
# <platform_name>_read <- function(file_list){ 
#
#
# output: dataframe of <uniprotswissprot> and <expression> values
# }
#
# add function name to get_platform_function
#


annotation_chip <- function(expression_df, biomart_attribute_name){
  probes = expression_df %>% 
    mutate(probe_id = row.names(expression_df))
  
  annot <- getBM(attributes=c(biomart_attribute_name,'uniprotswissprot'), mart = ensembl) %>% 
    dplyr::rename(probe_id = biomart_attribute_name)
  
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


affy_workflow <- function(file_list, affy_cdfname){
  raw_data = affy::ReadAffy(verbose=TRUE, filenames = file_list, cdfname = affy_cdfname) # read the row data to raw.data
  rma_df = affy::rma(raw_data) # mas5(); normalise the raw data 
  expression_df = exprs(rma_df)	%>% 
    as.data.frame() # get the normalised data in a readable matrix
}

oligo_workflow <- function(file_list){
  raw_data <- oligo::read.celfiles(verbose=TRUE, filenames=file_list)
  rma_df <- oligo::rma(raw_data) # mas5(); normalise the raw data
  expression_df = exprs(rma_df)	%>%
    as.data.frame() # get the normalised data in a readable matrix
}

# GPL96
# [HG-U133A] Affymetrix Human Genome U133A Array

GPL96_read <- function(file_list){
  annotation_chip(affy_workflow(file_list = file_list,
                                affy_cdfname = "HG-U133A"),
                  biomart_attribute_name = 'affy_hg_u133a')
}

# GPL97
# [HG-U133B] Affymetrix Human Genome U133B Array

GPL97_read <- function(file_list){
  annotation_chip(affy_workflow(file_list = file_list,
                                affy_cdfname = "HG-U133B"),
                  biomart_attribute_name = 'affy_hg_u133b')
}

# GPL570 
# [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
#library("hgu133plus2.db") #GPL570 platform

#file_list = sample_df$file_path[3:5]

GPL570_read <- function(file_list){
  annotation_chip(affy_workflow(file_list = file_list,
                                affy_cdfname = "HG-U133_Plus_2"),
                  biomart_attribute_name = 'affy_hg_u133_plus_2')
}

# GPL571 
# [HG-U133A_2] Affymetrix Human Genome U133A 2.0 Array
#library("hgu133a2.db") #GPL571 platform

GPL571_read <- function(file_list){
  annotation_chip(affy_workflow(file_list = file_list,
                                affy_cdfname = "HG-U133A_2"),
                  biomart_attribute_name = 'affy_hg_u133a_2')
}

# GPL5175
# [HuEx-1_0-st] Affymetrix Human Exon 1.0 ST Array [transcript (gene) version]



# GPL6244	[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]

GPL6244_read <- function(file_list){
  annotation_chip(oligo_workflow(file_list = file_list),
                  biomart_attribute_name = 'affy_hugene_1_0_st_v1')
}


# GPL17692	
#	[HuGene-2_1-st] Affymetrix Human Gene 2.1 ST Array [transcript (gene) version]

GPL17692_read <- function(file_list){
  annotation_chip(oligo_workflow(file_list = file_list),
                  biomart_attribute_name = 'affy_hugene_2_0_st_v1')
  
}



