source('helper_functions.R')

tissue_type = 'oesophagus'

sample_root = "C://Users/sebes/Dropbox/linkgroup/signaling2020_sample/"

sample_cels = list.files(path=sample_root, pattern = "(CEL|cel)(.gz)?$", recursive = TRUE)


sample_df <- sample_cels %>% 
  str_split_fixed(.,pattern = "/", n = 5) %>% 
  as.data.frame() %>% 
  `colnames<-`(c("tissue_type", "condition", "platform", "project", "sample_id")) %>% 
  mutate(file_path = paste0(sample_root, sample_cels))


tissue_sample_df <- sample_df %>% 
  filter(tissue_type == tissue_type)

conditions <- unique(tissue_sample_df$condition)

result_df <- map_dfr(conditions, run_condition, tissue_sample_df)

result_df_wide <- result_df %>% 
  pivot_wider(names_from = condition, values_from = expression)

write_tsv(result_df_wide, paste0(tissue_type, "_abundances.tsv"))
write_tsv(tissue_sample_df, paste0(tissue_type, "_files.tsv"))
