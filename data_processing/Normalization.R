library(frma)
library(affy)
library(oligo)
library(arrayQualityMetrics)


#############################################################################
#                                                                           #
#                             Read raw data                                 #
#                                                                           #
#############################################################################

# Read raw data (.CEL) files
ReadCEL = function(Data, Pheno=NULL, QC = FALSE){
  if (is.list(Data)) {
    DataID = as.character(unlist(Data$GEOStudyID))
    GPL = Data$GEOPlatformID
    Path = Data$Path
  } else {
    Path = Data
  }
  
  ##List .CEL file names in current file path
  files = list.files(path=Path,
                      pattern = "(CEL|cel)(.gz)?$",
                      full.names = T,
                      recursive = T)
  
  
  if(QC == TRUE){
    # Remove files with poor quality
    removeFiles = rownames(Pheno)[which(Pheno$Remove == 1)]
    if(length(removeFiles) != 0){
      files = files[!Reduce(`|`, lapply(removeFiles, function(y) endsWith(files, y)))]
      Pheno = Pheno[-which(Pheno$Remove == 1), ]
    }
  }
  
  library(affyio)
  Pheno = cbind(Pheno, ScanDate = get.celfile.dates(files))
  
  ## Read .CEL files
  if(GPL %in% c("GPL6244", "GPL5175", "GPL17692")){
    RawData = oligo::read.celfiles(verbose=TRUE, 
                                   filenames=files, 
                                   phenoData = AnnotatedDataFrame(Pheno))
  } else {
    RawData = ReadAffy(filenames=files, 
                        compress = TRUE, 
                        phenoData = AnnotatedDataFrame(Pheno))
    
    
    if(QC == FALSE){
      QAPath = file.path(Path, "RawData_QA")
    } else {
      QAPath = file.path(Path, "RawData_QA_afterQC")
    }
    
    
    if(!dir.exists(QAPath) | length(list.files(QAPath)) == 0){
      arrayQualityMetrics(RawData, outdir = QAPath)
      gc()
    }
    
  }
  sampleNames(RawData) = unlist(lapply(strsplit(sampleNames(RawData), ".", fixed =TRUE), `[[`, 1))
  sampleNames(RawData) = unlist(lapply(strsplit(sampleNames(RawData), "_", fixed =TRUE), `[[`, 1))
  
  ##Output
  RawData
  
}




#############################################################################
#                                                                           #
#                               Normalize data                              #
#                                                                           #
#############################################################################

## Normalize data
Norm = function(Data, normMethod, backgroundCorrection = T, normalization = T){ #, bioVar = NULL, adjVar = NULL, intVar = NULL, RemoveAdj = T){
  PhenoData = pData(Data)
  GPL = unique(PhenoData$GEOPlatformID)
  
  if (normMethod == 'frma'){
    library(frma)
    ExpSet = frma(Data)
  } else if (normMethod == 'rma') {
    if(GPL %in% c("GPL6244", "GPL5175", "GPL17692")){
      ExpSet = oligo::rma(Data, normalize = normalization, background = backgroundCorrection)
    } else {
      ExpSet = affy::rma(Data, normalize = normalization, background = backgroundCorrection)
    }
  }  else if (normMethod == 'gcrma') {
      ExpSet = affy::gcrma(Data, normalize = normalization, optical.correct = backgroundCorrection)
  } else if (normMethod == 'mas5') {
    ExpSet <- mas5(Data, normalize = normalization)
# } else if (normMethod == 'snm') {
#   bio.var = as.data.frame(PhenoData[, bioVar])
#   adj.var = as.data.frame(PhenoData[, adjVar])
#   int.var = as.data.frame(PhenoData[, intVar])
#   
#   rownames(adj.var) = PhenoData$ArrayDataFile
#   colnames(adj.var) = adjVar
#   colnames(int.var) = intVar
#   
#   int.var$sample = as.factor(int.var$sample)
#   
#   adj.var = model.matrix(~.,adj.var)
#   bio.var = model.matrix(~., bio.var)
#   raw.data = as.matrix(exprs(Data))
#   
#   ExpSet = snm(raw.data, bio.var, adj.var, int.var, rm.adj = RemoveAdj)
#   
#   cat(ExpSet$iter.pi0s)
#   ExpSet = ExpSet$norm.dat
#   colnames(ExpSet) = colnames(raw.data)
  }
  
  ExpSet
}




