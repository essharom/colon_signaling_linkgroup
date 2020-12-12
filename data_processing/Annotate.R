library("biomaRt")
library(Biobase)
library(matrixStats)
library(jetset)

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
  } else if(GPL == "GPL8300"){
    DB = "affy_hg_u95av2"
  }
  
  DB
}


PlatformName = function(GPL){
  if(GPL == "GPL96"){
    DB = "hgu133a"
  } else if(GPL == "GPL570"){
    DB = "hgu133plus2"
  }else if(GPL == "GPL8300"){
    DB = "hgu95av2"
  } else if(GPL %in% c("GPL17692", "GPL6244", "GPL5175", "GPL571", "GPL97")){
    DB = NULL
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

Probe2Entrez <- function(Data, Info, method = "median"){
  annotationDB = PlatformDB(Info$GEOPlatformID)
  if(class(Data) == "ExpressionSet")
  {
    ProbExp = exprs(Data)
  } else {
    ProbExp = Data
  }

  annot = getBM(attributes=c(annotationDB, 'entrezgene_id'), 
                filters = annotationDB, 
                values = row.names(ProbExp), 
                mart = ensembl)
  annot = annot[-(which(annot[,1] %in% unique(annot[duplicated(annot[,1]),1]))),]
  ProbExp = ProbExp[which(rownames(ProbExp) %in% annot[,1]),]
  
  if(method == "median"){
    ProbExp = merge(annot, ProbExp, by.x = annotationDB, by.y = "row.names")
    
    GeneExp = aggregate(.~entrezgene_id, ProbExp[,-1], FUN = median, na.action = na.omit)
    
    rownames(GeneExp) <- GeneExp$entrezgene_id
    GeneExp <- as.matrix(GeneExp[, -1])
  } else if (method == "maxIQR"){
    testStat = rowIQRs(ProbExp)
    names(testStat) = rownames(ProbExp)
    
    map = split.default(testStat, annot[which(annot[,1] %in% names(testStat)), 2])
    uniqGenes = sapply(map, function(x) names(which.max(x))) %>% stack()
    colnames(uniqGenes) = c("ProbeID", "EntrezID")
    
    GeneExp = ProbExp[which(rownames(ProbExp) %in% uniqGenes$ProbeID),]
    rownames(GeneExp) = uniqGenes$EntrezID
  }
  
  GeneExp
}


Probe2Symbol <- function(Data, Info, method = "median"){
  annotationDB = PlatformDB(Info$GEOPlatformID)
  if(class(Data) == "ExpressionSet")
  {
    ProbExp = exprs(Data)
  } else {
    ProbExp = Data
  }
  
  annot = getBM(attributes=c(annotationDB, 'hgnc_symbol'), 
                filters = annotationDB, 
                values = row.names(ProbExp), 
                mart = ensembl)
  annot = annot[-(which(annot[,1] %in% unique(annot[duplicated(annot[,1]),1]))),]
  ProbExp = ProbExp[which(rownames(ProbExp) %in% annot[,1]),]
  
  if(method == "median"){
    ProbExp = merge(annot[annot[,2]!="",], ProbExp, by.x = annotationDB, by.y = "row.names")
    
    GeneExp = aggregate(.~hgnc_symbol, ProbExp[,-1], FUN = median, na.action = na.omit)
    
    rownames(GeneExp) <- GeneExp$hgnc_symbol
    GeneExp <- as.matrix(GeneExp[, -1])
  } else if (method == "maxIQR"){
    testStat = rowIQRs(ProbExp)
    names(testStat) = rownames(ProbExp)
    
    map = split.default(testStat, annot[which(annot[,1] %in% names(testStat)), 2])
    uniqGenes = sapply(map, function(x) names(which.max(x))) %>% stack()
    colnames(uniqGenes) = c("ProbeID", "GeneSymbol")
    
    GeneExp = ProbExp[which(rownames(ProbExp) %in% uniqGenes$ProbeID),]
    rownames(GeneExp) = uniqGenes$GeneSymbol
  } else if (method == "jetset"){
    Platform = PlatformName(Info$GEOPlatformID)
    GeneSym = as.vector(annot$hgnc_symbol)
    GeneSym = GeneSym[GeneSym != ""]
    uniqGenes = jmap(Platform, symbol = GeneSym)
    uniqGenes = data.frame(ProbeID = unname(uniqGenes), Symbol = names(uniqGenes))
    uniqGenes = unique(uniqGenes[!is.na(uniqGenes$ProbeID),])
    
    GeneExp = ProbExp[which(rownames(ProbExp) %in% uniqGenes$ProbeID),]
    rownames(GeneExp) = as.vector(uniqGenes$Symbol[which(uniqGenes$ProbeID %in% rownames(GeneExp))])
    
  }
  
  GeneExp
}


Symbol2UniprotID = function(Data){
  annot = getBM(attributes=c('hgnc_symbol', 'uniprotswissprot'),  mart = ensembl) %>%
  filter(hgnc_symbol!="", uniprotswissprot != "")
  annot = aggregate(uniprotswissprot~hgnc_symbol, annot, paste, collapse="|")

  Data = merge(annot, Data, by.x = "hgnc_symbol", by.y = "row.names", all.y = TRUE)
  Data = data.frame(Data, row.names = "hgnc_symbol")
  Data = dplyr::rename(Data, "Swissprot ID" = uniprotswissprot)
}

