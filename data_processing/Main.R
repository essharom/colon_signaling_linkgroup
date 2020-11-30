memory.limit(9999999999)

library(purrr)
library(dplyr)

library(preprocessCore)
library(sva)
library(YuGene)
library(ggfortify)
library(heatmaply)



PlatformID   = NULL # Use NULL if you wish to include all of the platforms, or specify GEO platform ID (e.g. "GPL570")
DataSetInfo = "../../01_Data/Data.txt"
TissueType   = "esophagus" 
DataPath = "../../01_Data/01_Raw data"

source("Info.R")
## Read general info about data sets & extract GEO study ID for the data used
DataInfo = Info(
  Path         = DataPath, 
  infoFile     = DataSetInfo, 
  PlatformUsed = PlatformID, 
  Tissue       = TissueType
)

## Read phenotype information about data sets used
Description = purrr::map(DataInfo, Pheno)


#############################################################################
#                                                                           #
#                             Read raw data                                 #
#                                                                           #
#############################################################################
source("Normalization.R")
RawExp = map2(DataInfo, Description, ReadCEL, QC=TRUE)
PhenoData = bind_rows(purrr::map(RawExp, pData))


save(list = c("RawExp"), 
     file = "../../01_Data/02_Normalized data/02_Esophagus/RawData_Corrected.RData")
save(list = c("DataInfo", "Description", "PhenoData"), 
     file = "../../01_Data/02_Normalized data/02_Esophagus/DataInfo_Corrected.RData")


## RMA background correction
RMA_BG = purrr::map(RawExp, 
                    Norm, 
                    normMethod = "rma", 
                    normalization = F)

source("Annotate.R")
RMA_BG_GeneExp = map2(RMA_BG, 
                      DataInfo,
                      Probe2Symbol)
#RMA_BG_ProtExp = map2(RMA_BG, DataInfo, Probe2UniProt)


## RMA background correction and normalization
RMA = purrr::map(RawExp, 
                 Norm, 
                 normMethod = "rma", 
                 normalization = T)
RMA_GeneExp = map2(RMA, 
                   DataInfo,
                   Probe2Symbol)
#RMA_ProtExp = map2(RMA, DataInfo, Probe2UniProt)

for(i in 1:length(DataInfo)){
  QAPath = file.path(DataInfo[[i]]$Path, "RMANormData_QA_afterQC")
  if(!dir.exists(QAPath) | length(list.files(QAPath)) == 0){
    arrayQualityMetrics(RMA[[i]], outdir = QAPath)
    gc()
  }
}

## fRMA background correction and normalization
fRMA = purrr::map(RawExp[-6], 
                  Norm, 
                  normMethod = "frma")
fRMA_GeneExp = map2(fRMA, 
                    DataInfo[-6],
                    Probe2Symbol)
#fRMA_ProtExp = map2(fRMA, DataInfo, Probe2UniProt)

save(list = c("RMA", "RMA_BG", "fRMA"), 
     file = "../../01_Data/02_Normalized data/02_Esophagus/NormData_Corrected.RData")
save(list = c("RMA_GeneExp", "RMA_BG_GeneExp", "fRMA_GeneExp"), 
     file = "../../01_Data/02_Normalized data/02_Esophagus/GeneExpression_Corrected.RData")


#save(list = c("RMA", "RMA_BG", "fRMA"), 
#     file = "../../01_Data/02_Normalized data/02_Esophagus/NormData_Corrected.RData")
#save(list = c("RMA_GeneExp", "RMA_BG_GeneExp", "fRMA_GeneExp"), 
#     file = "../../01_Data/02_Normalized data/02_Esophagus/GeneExpression_Corrected.RData")






## Compare different normalization methods
### Boxplots
source("Visualize.R")
par(mfrow=c(1,3))
BoxPlot(RMA_BG, 
        Logaritmic = "y", 
        Title="RMA background correction")
BoxPlot(RMA, 
        Logaritmic = "y", 
        Title="RMA normalization")
BoxPlot(fRMA, 
        Logaritmic = "y", 
        Title="fRMA normalization")

### PCA plots
source("ShapeData.R")

RMA_BG_GeneDF = na.omit(List2DF(RMA_BG_GeneExp))
RMA_GeneDF = na.omit(List2DF(RMA_GeneExp))
fRMA_GeneDF = na.omit(List2DF(fRMA_GeneExp))


plot1 = PCAPlot(RMA_BG_GeneDF, 
                PhenoData,
                colorCol = "Phenotype",
                shapeCol = "GEOStudyID",
                Title = "RMA background correction")
plot2 = PCAPlot(RMA_GeneDF, 
                PhenoData,
                colorCol = "Phenotype",
                shapeCol = "GEOStudyID",
                Title = "RMA normalization")
plot3 = PCAPlot(fRMA_GeneDF, 
                PhenoData,
                colorCol = "Phenotype",
                shapeCol = "GEOStudyID",
                Title = "fRMA normalization")
plot_grid(plot1, plot2, plot3, ncol=2, nrow = 2, title = "Title")


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
     file = "../../01_Data/02_Normalized data/02_Esophagus/BatchCorData_Corrected.RData")


PCA_RMA_QN = PCAPlot(RMA_QN,
                     PhenoData,
                     colorCol = "Phenotype",
                     shapeCol = "GEOStudyID") + 
  ggtitle("RMA normalized samples after quantile normalization")
PCA_RMA_ComBat = PCAPlot(RMA_Combat, 
                         PhenoData,
                         colorCol = "Phenotype",
                         shapeCol = "GEOStudyID") +
  ggtitle("RMA normalized samples after ComBat batch correction")
PCA_RMA_YuGene = autoplot(pca(t(RMA_YuGene)),
                          data = PhenoData, 
                          colour = "Phenotype", 
                          shape = "GEOStudyID") +
  ggtitle("RMA normalized samples after YuGene batch correction")

plot_grid(PCA_RMA, PCA_RMA_QN, PCA_RMA_ComBat, PCA_RMA_YuGene, 
          ncol=2, 
          nrow = 2, 
          title = "Title")


PCA_fRMA_QN = PCAPlot(fRMA_QN, 
                      PhenoData,
                      colorCol = "Phenotype",
                      shapeCol = "GEOStudyID") +
  ggtitle("fRMA normalized samples after quantile normalization")
PCA_fRMA_ComBat = PCAPlot(fRMA_Combat, 
                          PhenoData,
                          colorCol = "Phenotype",
                          shapeCol = "GEOStudyID") +
  ggtitle("fRMA normalized samples after ComBat batch correction")
PCA_fRMA_YuGene = autoplot(pca(t(fRMA_YuGene)),
                           data = fRMAPheno, 
                           colour = "Phenotype", 
                           shape = "GEOStudyID") + 
  ggtitle("fRMA normalized samples after YuGene batch correction")
plot_grid(PCA_fRMA, PCA_fRMA_QN, PCA_fRMA_ComBat, PCA_fRMA_YuGene,
          ncol=2, 
          nrow = 2, 
          title = "Title")

