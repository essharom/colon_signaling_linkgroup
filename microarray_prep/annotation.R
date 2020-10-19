#sample_script_for GPL571


listMarts(host="http://uswest.ensembl.org")
ensembl = useMart("ENSEMBL_MART_ENSEMBL", host="http://uswest.ensembl.org")
datasets = listDatasets(ensembl)
head(datasets)

ensembl = useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl")
head(searchAttributes(ensembl, pattern="affy"))
u133a2_swissprot_ids <- getBM(attributes=c('affy_hg_u133a_2','uniprotswissprot'), mart = ensembl)

probes <- as.data.frame(ControlC)

affy_GPL571 <- u133a2_swissprot_ids %>% 
  left_join(probes, by = c('affy_hg_u133a_2'='probes')) %>% 
  filter(uniprotswissprot!="") %>% 
  filter(`GSM969967_BE0002-1-A1.CEL` != "") %>% 
  mutate(expression=as.numeric(`GSM969967_BE0002-1-A1.CEL` )) %>% 
  group_by(uniprotswissprot) %>% 
  summarise(mean_expression = mean(expression))


#BiocManager::install("pd.hugene.1.0.st.v1")
library("pd.hugene.1.0.st.v1")
dat = toTable(pmSequence())



library(hgu95av2.db)
dat = toTable(hgu95av2ENSEMBL)
dat[dat[,1]=="32337_at",]
probe_id      ensembl_id
5562 32337_at ENSG00000122026

dim(dat)