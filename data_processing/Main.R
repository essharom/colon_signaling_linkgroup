memory.limit(9999999999)

library(purrr)
library(dplyr)

library(preprocessCore)
library(sva)
library(YuGene)
library(ggfortify)
library(heatmaply)



PlatformID   = NULL # Use NULL if you wish to include all of the platforms, or specify GEO platform ID (e.g. "GPL570")
DataSetInfo = "B:/OneDrive/Documents/Microarray data/Data.txt"
TissueType   = "colon" 
DataPath = "B:/OneDrive/Documents/Microarray data/01_Raw data"

source("Info.R")
## Read general info about data sets & extract GEO study ID for the data used
DataInfo = Info(
  Path         = DataPath, 
  infoFile     = DataSetInfo, 
  PlatformUsed = PlatformID, 
  Tissue       = TissueType,
  QC           = "Remove0"
)

## Read phenotype information about data sets used
Description = purrr::map(DataInfo, Pheno)


#############################################################################
#                                                                           #
#                                 Raw data                                  #
#                                                                           #
#############################################################################
source("Normalization.R")
RawExp = map2(DataInfo, Description, ReadCEL, QC="Remove0")
PhenoData = bind_rows(purrr::map(RawExp, pData))


save(list = c("RawExp"), 
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/RawData.RData")
save(list = c("DataInfo", "Description", "PhenoData"), 
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/DataInfo.RData")

#############################################################################
#                                                                           #
#                               Normalized data                             #
#                                                                           #
#############################################################################

## RMA background correction
RMA_BG = purrr::map(RawExp, 
                    Norm, 
                    normMethod = "rma", 
                    normalization = F)

## RMA background correction and normalization
RMA = purrr::map(RawExp, 
                 Norm, 
                 normMethod = "rma", 
                 normalization = T)

## fRMA background correction and normalization
fRMA = purrr::map(RawExp, 
                  Norm, 
                  normMethod = "frma", 
                  fRMATarget = "core")

save(list = c("RMA", "RMA_BG", "fRMA"), 
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/NormData.RData")

#############################################################################
#                                                                           #
#                             Probe to gene mapping                         #
#                                                                           #
#############################################################################
source("Annotate.R")
RMA_BG_GeneExp = map2(RMA_BG, 
                      DataInfo,
                      Probe2Entrez,
                      method = "median")
RMA_GeneExp = map2(RMA, 
                   DataInfo,
                   Probe2Entrez,
                   method = "median")
fRMA_GeneExp = map2(fRMA, 
                    DataInfo,
                    Probe2Entrez,
                    method = "median")
save(list = c("RMA_GeneExp", "RMA_BG_GeneExp", "fRMA_GeneExp"), 
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/GeneExpression_median.RData")

RMA_BG_GeneExp = map2(RMA_BG, 
                      DataInfo,
                      Probe2Entrez,
                      method = "maxIQR")
RMA_GeneExp = map2(RMA, 
                   DataInfo,
                   Probe2Entrez,
                   method = "maxIQR")
fRMA_GeneExp = map2(fRMA, 
                    DataInfo,
                    Probe2Entrez,
                    method = "maxIQR")
save(list = c("RMA_GeneExp", "RMA_BG_GeneExp", "fRMA_GeneExp"), 
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/GeneExpression_maxIQR.RData")


#############################################################################
#                                                                           #
#                             Batch corrected data                          #
#                                                                           #
#############################################################################
## Compare different normalization methods
source("ShapeData.R")

RMA_BG_GeneDF = na.omit(List2DF(RMA_BG_GeneExp))
RMA_GeneDF = na.omit(List2DF(RMA_GeneExp))
fRMA_GeneDF = na.omit(List2DF(fRMA_GeneExp))


## Batch correction
library(preprocessCore)
library(sva)
library(YuGene)

PhenoData = PhenoData[which(rownames(PhenoData) %in% colnames(RMA_GeneDF)),]
fRMAPheno = PhenoData[which(rownames(PhenoData) %in% colnames(fRMA_GeneDF)),]

source("BatchCorrection.R")
RMA_QN= normalize.quantiles(as.matrix(RMA_GeneDF))
RMA_Combat = ComBat(as.matrix(RMA_GeneDF), batch=PhenoData$GEOStudyID, mod=model.matrix(~as.factor(Phenotype), data = PhenoData))
RMA_SVA = SVA(RMA_GeneDF, PhenoData, "Phenotype")
RMA_YuGene = as.matrix(YuGene(as.matrix(RMA_GeneDF)))
colnames(RMA_QN) = colnames(RMA_GeneDF)

fRMA_QN= normalize.quantiles(as.matrix(fRMA_GeneDF))
fRMA_Combat = ComBat(as.matrix(fRMA_GeneDF), batch=fRMAPheno$GEOStudyID, mod=model.matrix(~as.factor(Phenotype), data = fRMAPheno))
fRMA_SVA = SVA(fRMA_GeneDF, fRMAPheno, "Phenotype")
fRMA_YuGene = as.matrix(YuGene(as.matrix(fRMA_GeneDF)))
colnames(fRMA_QN) = colnames(fRMA_GeneDF)

save(list = c("RMA_QN", "RMA_Combat", "RMA_YuGene","fRMA_QN", "fRMA_Combat", "fRMA_YuGene", "RMA_SVA", "fRMA_SVA"), 
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/BatchCorData.RData")

#############################################################################
#                                                                           #
#                                 Gene barcode                              #
#                                                                           #
#############################################################################

source("Annotate.R")

library(frma)

bc = purrr::map(fRMA, barcode)
#bc_DF =  na.omit(List2DF(bc))
#bc_GeneDF =  na.omit(List2DF(map2(bc, DataInfo[-6], Probe2Entrez)))

save(list = c("bc"), 
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/GeneBarcode.RData")
