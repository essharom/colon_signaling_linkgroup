library(GEOquery)


## Read in info file about data sets
Info <- function(Path, infoFile, PlatformUsed = NULL, Tissue = NULL, QC = NULL){
  DataSet = read.csv(infoFile, header = TRUE, sep = "\t", na.strings = c("NA", "-", "", "N/A"))
  
  if(!is.null(PlatformUsed)){
    DataSet = DataSet[DataSet$GEOPlatformID == PlatformUsed & DataSet$Tissue == Tissue, ]
  }
  if(!is.null(Tissue)){
    DataSet = DataSet[DataSet$Tissue == Tissue, ]
  }
  if(!is.null(QC)){
    DataSet = DataSet[which(DataSet[, QC] == 0),]
  }
  DataSet = DataSet[, c("GEOStudyID","GEOPlatformID", "PlatformType", "Organism", "Tissue")]
  
  Path = paste(Path, DataSet$GEOStudyID, sep="/")
  DataSet = cbind.data.frame(DataSet, Path)
  DataSet = split(DataSet, seq(nrow(DataSet)))
}


## Read in data set specific phenotype information
Pheno <- function(Data){
  DataID <- as.character(unlist(Data$GEOStudyID))
  GPL = Data$GEOPlatformID
  filePath <- Data$Path
  PhenoFile <- paste(filePath, paste("PhenoData_", GPL, ".txt", sep = ""), sep = "/")
  
  if(file.exists(PhenoFile)){
    PhenoData <- read.csv(PhenoFile, header = TRUE, sep = "\t", na.strings = c("NA", "-", ""), row.names = 1)
  } else {
    cat("No PhenoData.txt file in depository \n Downloading data from NCBI")
    
    temp = getGEO(DataID)
    
    for(j in 1:length(temp)){
      PhenoData = pData(temp[[j]])
      
      GPL = unique(PhenoData[['platform_id']])
      if(any('Homo sapiens'!= PhenoData[['organism_ch1']])){next}
      if(GPL == "GPL6801" || GPL == "GPL7723"){next}
      
      rownames(PhenoData) = unlist(mapply('[[', strsplit(PhenoData$supplementary_file, "/", fixed =TRUE), lengths(strsplit(PhenoData$supplementary_file, "/", fixed =TRUE))))
    
      celFiles = list.files(path=filePath, pattern = "(CEL|cel)(.gz)?$", recursive = TRUE)
      PhenoData = PhenoData[celFiles,]
    
      write.table(PhenoData, file = paste(filePath, "/PhenoData_", GPL, ".txt", sep = ""), 
                  sep = "\t", row.names = TRUE, col.names = TRUE)
    }
  }
  
  
  PhenoData
}
