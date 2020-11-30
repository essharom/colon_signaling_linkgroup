source("Visualize.R")

library(affyPLM)
library(oligo)

QualityAnalysis = function(Data, DataInfo, PhenoData){
  if(class(Data)[1] %in% c("AffyBatch", "GeneFeatureSet")){
    Path = paste(DataInfo$Path,"/RawData_QC", sep = "")
    if(!dir.exists(Path)){dir.create(Path)}
    
    PSet = PLM(Data, Path)
    logData = logPM(Data)
    
    #rawIntensityImages(Data, Path)
    weightImages(PSet, Path)
    residualImages(PSet, Path)
  } else {
    Path = paste(DataInfo$Path,"/NormData_QC", sep = "")
    if(!dir.exists(Path)){dir.create(Path)}
    
    logData = logNorm(Data)
  }
  PhenoData = pData(Data)
  ExpData = exprs(Data)
  
  DistanceHeatmap(Data, PhenoData, Path, "Phenotype")
  PCAPlot(ExpData, PhenoData, colorCol = "Phenotype", Title = unique(PhenoData$GEOStudyID))
  ggsave(paste(Path, "/PCA.jpg", sep = ""))
  histograms(logData, Path)
  boxplots(logData, Path)
  MAPlot(Data, Path)
}

PLM = function(Data, Path){
  if("AffyBatch" %in% class(Data)){
    Pset = affyPLM::fitPLM(Data, output.param=list(varcov="none"))
    
    jpeg(paste(Path, "/RLE.jpg",sep=""))
    affyPLM::RLE(Pset)
    dev.off()
    jpeg(paste(Path, "/NUSE.jpg",sep=""))
    affyPLM::NUSE(Pset)
    dev.off()
  } else if ("GeneFeatureSet" %in% class(Data)) {
    Pset = oligo::fitProbeLevelModel(Data)
    
    jpeg(paste(Path, "/RLE.jpg",sep=""))
    oligo::RLE(Pset)
    dev.off()
    jpeg(paste(Path, "/NUSE.jpg",sep=""))
    oligo::NUSE(Pset)
    dev.off()
  }
  
  Pset
}


logPM = function(Data){
  if("AffyBatch" %in% class(Data)){
    pmexp = affy::pm(Data)
  } else {
    pmexp = pm(Data)
  }
  
  sampleNames = vector()
  logs = vector()
  
  for (i in 1:length(colnames(Data))){
    sampleNames = c(sampleNames,rep(Data$GEOSampleID[i],dim(pmexp)[1]))
    logs = c(logs,log2(pmexp[,i]))
  }
  
  logData = data.frame(logInt=logs,GEOSampleID=sampleNames)
  
  logData
}

logNorm = function(Data){
  if("ExpressionSet" %in% class(Data)){
    Data = exprs(Data)
  } 
  
  for (i in 1:length(colnames(Data))){
    sampleNames = c(sampleNames,rep(pData(Data)$GEOSampeID[i],dim(Data)[1]))
    logs = c(logs,log2(Data[,i]))
  }
  
  logData = data.frame(logInt=logs,GEOSampleID=sampleNames)
  
  logData
}


