library("biomaRt")
library(Biobase)

PlatformDB = function(GPL){
  if(GPL == "GPL96"){
    DB = "affy_hg_u133a"
  } else if(GPL == "GPL97"){
    DB = "affy_hg_u133b"
  } else if(GPL == "GPL570"){
    DB = "affy_hg_u133_plus_2"
  } else if(GPL == "GPL571"){
    DB = "affy_hg_u133a_2"
  } else if(GPL == "GPL5175"){
    DB = "affy_huex_1_0_st_v2"
  } else if(GPL == "GPL6244"){
    DB = "affy_hugene_1_0_st_v1"
  } else if(GPL == "GPL17692"){
    DB = "affy_hugene_2_0_st_v1"
  }
  
  DB
}


ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl)

Probe2UniProt <- function(Data, Info){
  annotationDB = PlatformDB(Info$GEOPlatformID)
  ProbExp = exprs(Data)
  
  annot <- getBM(attributes=c(annotationDB,'uniprotswissprot'), mart = ensembl) %>% 
    dplyr::rename(probe_id = annotationDB) %>% 
    dplyr::mutate(probe_id = as.character(probe_id)) %>% 
    filter(uniprotswissprot!="", probe_id!="")
  
  ProbExp = merge(annot, ProbExp, by.x = "probe_id", by.y = "row.names")
  
  ProtExp = aggregate(.~uniprotswissprot, ProbExp[,-1], FUN = median, na.action = na.omit)
  rownames(ProtExp) <- ProtExp$uniprotswissprot
  ProtExp <- as.matrix(ProtExp[, -1])
  
  ProtExp
}

Probe2Entrez <- function(Data, Info){
  annotationDB = PlatformDB(Info$GEOPlatformID)
  ProbExp = exprs(Data)
  
  probes = row.names(ProbExp)
  
  annot <- getBM(attributes=c(annotationDB,'entrezgene_id'), mart = ensembl) %>% 
    dplyr::rename(probe_id = annotationDB) %>% 
    dplyr::mutate(probe_id = as.character(probe_id)) %>% 
    filter(entrezgene_id!="", probe_id!="")
  
  ProbExp = merge(annot, ProbExp, by.x = "probe_id", by.y = "row.names")
  
  GeneExp = aggregate(.~entrezgene_id, ProbExp[,-1], FUN = median, na.action = na.omit)
  rownames(GeneExp) <- GeneExp$entrezgene_id
  GeneExp <- as.matrix(GeneExp[, -1])
  
  GeneExp
}


Probe2Symbol <- function(Data, Info){
  annotationDB = PlatformDB(Info$GEOPlatformID)
  ProbExp = exprs(Data)
  
  probes = row.names(ProbExp)
  
  annot <- getBM(attributes=c(annotationDB,'external_gene_name'), mart = ensembl) %>% 
    dplyr::rename(probe_id = annotationDB) %>% 
    dplyr::mutate(probe_id = as.character(probe_id)) %>% 
    filter(external_gene_name!="", probe_id!="")
  
  ProbExp = merge(annot, ProbExp, by.x = "probe_id", by.y = "row.names")
  
  GeneExp = aggregate(.~external_gene_name, ProbExp[,-1], FUN = median, na.action = na.omit)
  rownames(GeneExp) <- GeneExp$external_gene_name
  GeneExp <- as.matrix(GeneExp[, -1])
  
  GeneExp
}


Symbol2UniprotID = function(Data){
  annot = getBM(attributes=c('external_gene_name', 'uniprotswissprot'),  mart = ensembl) %>%
  filter(external_gene_name!="", uniprotswissprot != "")
  annot = aggregate(uniprotswissprot~external_gene_name, annot, paste, collapse="|")

  Data = merge(annot, Data, by.x = "external_gene_name", by.y = "row.names", all.y = TRUE)
  Data = data.frame(Data, row.names = "external_gene_name")
  Data = dplyr::rename(Data, "Swissprot ID" = uniprotswissprot)
}

