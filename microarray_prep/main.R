source('helper_functions.R')

library(ggfortify)

current_tissue_type = 'eosophagus'

# add the new platform's name, which has analysis function in helper functions
platform_functions_list <- c("GPL570",
                             "GPL571",
                             "GPL17692",
                             "GPL96",
                             "GPL97",
                             "GPL5175",
                             "GPL6244")

sample_root = "C:/Users/kunsi/OneDrive/Data/"

sample_cels = list.files(path=sample_root, pattern = "(CEL|cel)(.gz)?$", recursive = TRUE)

#create design matrix for all samples used in the whole project
sample_df <- sample_cels %>% 
  str_split_fixed(.,pattern = "/", n = 5) %>% 
  as.data.frame() %>% 
  `colnames<-`(c("tissue_type", "condition", "platform", "project", "sample_id")) %>% 
  mutate(file_path = paste0(sample_root, sample_cels))

#GPL17692_read function must be finished but oligo exploration package first

tissue_sample_df <- sample_df %>% 
  dplyr::filter(tissue_type == current_tissue_type)

conditions <- unique(tissue_sample_df$condition)

# run abundance calculation for a tissue type, call necessary functions
result_df <- map(conditions, run_condition, tissue_sample_df) %>%
  purrr::reduce(left_join, by = "uniprotswissprot")

expression = data.frame(result_df, row.names ="uniprotswissprot")
expression = data.frame(t(expression))
rownames(expression) = unlist(lapply(strsplit(rownames(expression), ".", fixed =TRUE), `[[`, 1))
rownames(expression) = unlist(lapply(strsplit(rownames(expression), "_", fixed =TRUE), `[[`, 1))
expression=expression[ , colSums(is.na(expression)) == 0]

sample_data = data.frame(tissue_sample_df, row.names = "sample_id")

rownames(sample_data) = unlist(lapply(strsplit(rownames(sample_data), ".", fixed =TRUE), `[[`, 1))
rownames(sample_data) = unlist(lapply(strsplit(rownames(sample_data), "_", fixed =TRUE), `[[`, 1))


data = merge(x=expression, y= sample_data[,1:4], by ="row.names", all = TRUE)
data = data.frame(data, row.names = "Row.names")

pca_res <- stats::prcomp(data[,1:1496], scale. = TRUE)

autoplot(pca_res, data = data, colour = "condition")
autoplot(pca_res, data = data, colour = "project")
autoplot(pca_res, data = data, colour = "platform")
autoplot(pca_res, data = data, colour = "project", shape = "platform")

# make table human readable
result_df_wide <- result_df %>% 
  pivot_wider(names_from = condition, values_from = expression)

GSM970010_BE0102C1-RE.CEL.gz

# add HUGO GENE SYMBOLS to the `wide` dataframe possibly with biomart
annot_hgnc <- getBM(attributes=c('hgnc_symbol','uniprotswissprot'), mart = ensembl)

result_df_wide <- result_df_wide %>% 
  left_join(annot_hgnc)

# SUPPL_RESULT_Table
# write out result matrix with abundances
write_tsv(result_df_wide, paste0(tissue_type, "_abundances.tsv"))

# SUPPL_Table
# write out design matrix with all used samples for given tissue type
write_tsv(tissue_sample_df, paste0(tissue_type, "_files.tsv"))
