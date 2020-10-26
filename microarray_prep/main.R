source('helper_functions.R')

tissue_type = 'oesophagus'

# add the new platform's name, which has analysis function in helper functions
platform_functions_list <- c("GPL570",
                             "GPL571",
                             "GPL17692")

sample_root = "C://Users/sebes/Dropbox/linkgroup/signaling2020_sample/"

sample_cels = list.files(path=sample_root, pattern = "(CEL|cel)(.gz)?$", recursive = TRUE)

#create design matrix for all samples used in the whole project
sample_df <- sample_cels %>% 
  str_split_fixed(.,pattern = "/", n = 5) %>% 
  as.data.frame() %>% 
  `colnames<-`(c("tissue_type", "condition", "platform", "project", "sample_id")) %>% 
  mutate(file_path = paste0(sample_root, sample_cels))

#GPL17692_read function must be finished but oligo exploration package first

tissue_sample_df <- sample_df %>% 
  filter(tissue_type == tissue_type) %>% 
  filter(platform %in% platform_functions_list)

conditions <- unique(tissue_sample_df$condition)

# run abundance calculation for a tissue type, call necessary functions
result_df <- map_dfr(conditions, run_condition, tissue_sample_df)

# make table human readable
result_df_wide <- result_df %>% 
  pivot_wider(names_from = condition, values_from = expression)

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
