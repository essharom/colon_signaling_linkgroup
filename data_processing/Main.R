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
  QC           = "Remove"
)

## Read phenotype information about data sets used
Description = purrr::map(DataInfo, Pheno)

#############################################################################
#                                                                           #
#                                 Raw data                                  #
#                                                                           #
#############################################################################
source("Normalization.R")
RawExp = map2(DataInfo, Description, ReadCEL, QC="Remove")
PhenoData = bind_rows(purrr::map(RawExp, pData))


save(list = c("RawExp"), 
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/RawData_afterQC.RData")
save(list = c("DataInfo", "Description", "PhenoData"), 
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/DataInfo_afterQC.RData")

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
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/NormData_afterQC.RData")

source("Visualize.R")
par(mfrow=c(1,2))
BP_RMA =BoxPlot(RMA, 
                Logaritmic = "y", 
                Title="RMA normalization")
BP_fRMA =  BoxPlot(fRMA, 
                   Logaritmic = "y", 
                   Title="fRMA normalization")

### PCA plots
RMA_DF = na.omit(List2DF(RMA))
fRMA_DF = na.omit(List2DF(fRMA))

PCA_RMA = PCAPlot(RMA_DF,
                  PhenoData,
                  colorCol = "Phenotype",
                  shapeCol = "GEOStudyID") +
  ggtitle("RMA normalization")
PCA_fRMA = PCAPlot(fRMA_DF, 
                   PhenoData,
                   colorCol = "Phenotype",
                   shapeCol = "GEOStudyID") + 
  ggtitle("fRMA normalization")
grid_arrange_shared_legend(PCA_RMA, PCA_fRMA, ncol=2, nrow = 1)

DistanceHeatmap(RMA_DF, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/Normalization",
                c("GEOStudyID", "Phenotype"), 
                "RMA normalization")
DistanceHeatmap(fRMA_DF, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/Normalization",
                c("GEOStudyID", "Phenotype"), 
                "fRMA normalization")


#############################################################################
#                                                                           #
#                             Probe to gene mapping                         #
#                                                                           #
#############################################################################
library(lemon)

source("Annotate.R")

## Mean probe to gene mapping
RMA_BG_GeneExp = map2(RMA_BG, 
                      DataInfo,
                      Probe2Symbol,
                      method = "median")
RMA_GeneExp = map2(RMA, 
                   DataInfo,
                   Probe2Symbol,
                   method = "median")
fRMA_GeneExp = map2(fRMA, 
                    DataInfo,
                    Probe2Symbol,
                    method = "median")
save(list = c("RMA_GeneExp", "RMA_BG_GeneExp", "fRMA_GeneExp"), 
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/GeneExpression_median_afterQC.RData")


source("Visualize.R")
par(mfrow=c(1,2))
BP_RMA =BoxPlot(RMA_GeneExp, 
                Logaritmic = "y", 
                Title="RMA normalization")
BP_fRMA =  BoxPlot(fRMA_GeneExp, 
                   Logaritmic = "y", 
                   Title="fRMA normalization")

### PCA plots
RMA_DF = na.omit(List2DF(RMA_GeneExp))
fRMA_DF = na.omit(List2DF(fRMA_GeneExp))

PCA_RMA = PCAPlot(RMA_DF,
                  PhenoData,
                  colorCol = "Phenotype",
                  shapeCol = "GEOStudyID") +
  ggtitle("RMA normalization")
PCA_fRMA = PCAPlot(fRMA_DF, 
                   PhenoData,
                   colorCol = "Phenotype",
                   shapeCol = "GEOStudyID") + 
  ggtitle("fRMA normalization")
grid_arrange_shared_legend(PCA_RMA, PCA_fRMA, ncol=2, nrow = 1)

DistanceHeatmap(RMA_DF, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/Probe2GeneMapping",
                c("GEOStudyID", "Phenotype"), 
                "Median probe-to-gene mapping of RMA normalized data")
DistanceHeatmap(fRMA_DF, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/Probe2GeneMapping",
                c("GEOStudyID", "Phenotype"), 
                "Medain probe-to-gene mapping of fRMA normalized data")

## Maximal variance probe to gene mapping
RMA_BG_GeneExp = map2(RMA_BG, 
                      DataInfo,
                      Probe2Symbol,
                      method = "maxIQR")
RMA_GeneExp = map2(RMA, 
                   DataInfo,
                   Probe2Symbol,
                   method = "maxIQR")
fRMA_GeneExp = map2(fRMA, 
                    DataInfo,
                    Probe2Symbol,
                    method = "maxIQR")
save(list = c("RMA_GeneExp", "RMA_BG_GeneExp", "fRMA_GeneExp"), 
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/GeneExpression_maxIQR_afterQC.RData")

source("Visualize.R")
par(mfrow=c(1,2))
BP_RMA =BoxPlot(RMA_GeneExp, 
                Logaritmic = "y", 
                Title="RMA normalization")
BP_fRMA =  BoxPlot(fRMA_GeneExp, 
                   Logaritmic = "y", 
                   Title="fRMA normalization")

### PCA plots
RMA_DF = na.omit(List2DF(RMA_GeneExp))
fRMA_DF = na.omit(List2DF(fRMA_GeneExp))

PCA_RMA = PCAPlot(RMA_DF,
                  PhenoData,
                  colorCol = "Phenotype",
                  shapeCol = "GEOStudyID") +
  ggtitle("RMA normalization")
PCA_fRMA = PCAPlot(fRMA_DF, 
                   PhenoData,
                   colorCol = "Phenotype",
                   shapeCol = "GEOStudyID") + 
  ggtitle("fRMA normalization")
grid_arrange_shared_legend(PCA_RMA, PCA_fRMA, ncol=2, nrow = 1)

DistanceHeatmap(RMA_DF, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/Probe2GeneMapping",
                c("GEOStudyID", "Phenotype"), 
                "Maximal variance probe-to-gene mapping of RMA normalized data")
DistanceHeatmap(fRMA_DF, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/Probe2GeneMapping",
                c("GEOStudyID", "Phenotype"), 
                "Maximal variance probe-to-gene mapping of fRMA normalized data")

## Jetset probe to gene mapping
RMA_BG_GeneExp = map2(RMA_BG, 
                      DataInfo,
                      Probe2Symbol,
                      method = "jetset")
RMA_GeneExp = map2(RMA, 
                   DataInfo,
                   Probe2Symbol,
                   method = "jetset")
fRMA_GeneExp = map2(fRMA, 
                    DataInfo,
                    Probe2Symbol,
                    method = "jetset")
save(list = c("RMA_GeneExp", "RMA_BG_GeneExp", "fRMA_GeneExp"), 
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/GeneExpression_jetset_afterQC.RData")

source("Visualize.R")
par(mfrow=c(1,2))
BP_RMA =BoxPlot(RMA_GeneExp, 
                Logaritmic = "y", 
                Title="RMA normalization")
BP_fRMA =  BoxPlot(fRMA_GeneExp, 
                   Logaritmic = "y", 
                   Title="fRMA normalization")

### PCA plots
RMA_DF = na.omit(List2DF(RMA_GeneExp))
fRMA_DF = na.omit(List2DF(fRMA_GeneExp))

PCA_RMA = PCAPlot(RMA_DF,
                  PhenoData,
                  colorCol = "Phenotype",
                  shapeCol = "GEOStudyID") +
  ggtitle("RMA normalization")
PCA_fRMA = PCAPlot(fRMA_DF, 
                   PhenoData,
                   colorCol = "Phenotype",
                   shapeCol = "GEOStudyID") + 
  ggtitle("fRMA normalization")
grid_arrange_shared_legend(PCA_RMA, PCA_fRMA, ncol=2, nrow = 1)

DistanceHeatmap(RMA_DF, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/Probe2GeneMapping",
                c("GEOStudyID", "Phenotype"), 
                "Jetset single probe-to-gene mapping of RMA normalized data")
DistanceHeatmap(fRMA_DF, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/Probe2GeneMapping",
                c("GEOStudyID", "Phenotype"), 
                "Jetset single probe-to-gene mapping of fRMA normalized data")

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
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/BatchCorData_jetset_afterQC.RData")

PCA_RMA_QN = PCAPlot(RMA_QN,
                     PhenoData,
                     colorCol = "Phenotype",
                     shapeCol = "GEOStudyID") + 
  ggtitle("RMA normalization & quantile normalization")
PCA_RMA_ComBat = PCAPlot(RMA_Combat, 
                         PhenoData,
                         colorCol = "Phenotype",
                         shapeCol = "GEOStudyID") +
  ggtitle("RMA normalization & ComBat batch correction")
PCA_RMA_SVA = PCAPlot(RMA_SVA, 
                      PhenoData,
                      colorCol = "Phenotype",
                      shapeCol = "GEOStudyID") +
  ggtitle("RMA normalization & SVA batch correction")
PCA_RMA_YuGene = autoplot(pca(t(RMA_YuGene)),
                          data = PhenoData, 
                          colour = "Phenotype", 
                          shape = "GEOStudyID") +
  ggtitle("RMA normalization & YuGene batch correction")
grid_arrange_shared_legend(PCA_RMA_QN,
                           PCA_RMA_ComBat, 
                           PCA_RMA_SVA,
                           PCA_RMA_YuGene,
                           ncol=2,
                           nrow = 2)

PCA_fRMA_QN = PCAPlot(fRMA_QN, 
                      fRMAPheno,
                      colorCol = "Phenotype",
                      shapeCol = "GEOStudyID") +
  ggtitle("fRMA normalization & quantile normalization")
PCA_fRMA_ComBat = PCAPlot(fRMA_Combat, 
                          fRMAPheno,
                          colorCol = "Phenotype",
                          shapeCol = "GEOStudyID") +
  ggtitle("fRMA normalization & ComBat batch correction")
PCA_fRMA_SVA = PCAPlot(fRMA_SVA, 
                       fRMAPheno,
                       colorCol = "Phenotype",
                       shapeCol = "GEOStudyID") +
  ggtitle("fRMA normalization & SVA batch correction")
PCA_fRMA_YuGene = autoplot(pca(t(fRMA_YuGene)),
                           data = fRMAPheno, 
                           colour = "Phenotype", 
                           shape = "GEOStudyID") + 
  ggtitle("fRMA normalization & YuGene batch correction")
grid_arrange_shared_legend(PCA_fRMA_QN, 
                           PCA_fRMA_ComBat, 
                           PCA_fRMA_SVA,
                           PCA_fRMA_YuGene,
                           ncol=2,
                           nrow = 2)

DistanceHeatmap(RMA_Combat, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/BatchCorrection/afterQC", 
                c("GEOStudyID", "Phenotype"), 
                "ComBat batch correction of Jetset single probe-to-gene mapped RMA normalized data")
DistanceHeatmap(RMA_SVA, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/BatchCorrection/afterQC", 
                c("GEOStudyID", "Phenotype"), 
                "SVA based batch correction of Jetset single probe-to-gene mapped RMA normalized data")
DistanceHeatmap(RMA_YuGene, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/BatchCorrection/afterQC", 
                c("GEOStudyID", "Phenotype"), 
                "YuGene batch correction of Jetset single probe-to-gene mapped RMA normalized data")

DistanceHeatmap(fRMA_Combat, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/BatchCorrection/afterQC", 
                c("GEOStudyID", "Phenotype"), 
                "ComBat batch correction of Jetset single probe-to-gene mapped fRMA normalized data")
DistanceHeatmap(fRMA_SVA, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/BatchCorrection/afterQC", 
                c("GEOStudyID", "Phenotype"), 
                "SVA based batch correction of Jetset single probe-to-gene mapped fRMA normalized data")
DistanceHeatmap(fRMA_YuGene, 
                PhenoData, 
                "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/BatchCorrection/afterQC", 
                c("GEOStudyID", "Phenotype"), 
                "YuGene batch correction of Jetset single probe-to-gene mapped fRMA normalized data")
Normal = PhenoData$GEOSampleID[PhenoData$Phenotype == "Normal"]
Interm = PhenoData$GEOSampleID[PhenoData$Phenotype == "Intermediate"]
Cancer = PhenoData$GEOSampleID[PhenoData$Phenotype == "Cancer"]

NormalExp = apply(RMA_Combat[, which(colnames(RMA_Combat) %in% Normal)],1,median)
IntermExp = apply(RMA_Combat[, which(colnames(RMA_Combat) %in% Interm)],1,median)
CancerExp = apply(RMA_Combat[, which(colnames(RMA_Combat) %in% Cancer)],1,median)

RMACombat_PhenoExp = data.frame(Normal = NormalExp, 
                         Intermediate = IntermExp[which(names(IntermExp) %in% names(NormalExp))], 
                         Cancer = CancerExp[which(names(CancerExp) %in% names(NormalExp))])

RMACombat_PhenoExp = Symbol2UniprotID(RMACombat_PhenoExp)

write.csv(RMACombat_PhenoExp,"B:/OneDrive/Documents/Microarray data/03_Output/01_Colon/JetsetGeneExp.csv", row.names = TRUE)


#############################################################################
#                                                                           #
#                                 Gene barcode                              #
#                                                                           #
#############################################################################

source("Annotate.R")

library(frma)

bc = purrr::map(fRMA, barcode)
bc_DF =  na.omit(List2DF(bc))
bc_GeneDF =  na.omit(List2DF(map2(bc, DataInfo[-6], Probe2Symbol)))

Normal = PhenoData$GEOSampleID[PhenoData$Phenotype == "Normal"]
Interm = PhenoData$GEOSampleID[PhenoData$Phenotype == "Intermediate"]
Cancer = PhenoData$GEOSampleID[PhenoData$Phenotype == "Cancer"]

NormalExp = apply(bc_GeneDF[, which(colnames(bc_GeneDF) %in% Normal)],1,median)
IntermExp = apply(bc_GeneDF[, which(colnames(bc_GeneDF) %in% Interm)],1,median)
CancerExp = apply(bc_GeneDF[, which(colnames(bc_GeneDF) %in% Cancer)],1,median)

bc_PhenoExp = data.frame(Normal = NormalExp, 
                         Intermediate = IntermExp[which(names(IntermExp) %in% names(NormalExp))], 
                         Cancer = CancerExp[which(names(CancerExp) %in% names(NormalExp))])
bc_PhenoExp = Symbol2UniprotID(bc_PhenoExp)

write.csv(bc_PhenoExp,"B:/OneDrive/Documents/Microarray data/03_Output/01_Colon/GeneBarcode.csv", row.names = TRUE)

save(list = c("bc", "bc_DF", "bc_GeneDF", "bc_PhenoExp"), 
     file = "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/GeneBarcode_afterQC.RData")

DistanceHeatmap(bc_GeneDF, PhenoData, "B:/OneDrive/Documents/Microarray data/02_Normalized data/01_Colon/BatchCorrection/afterQC", c("GEOStudyID", "Phenotype"), "Gene barcode")

PhenoData = PhenoData[which(rownames(PhenoData) %in% colnames(bc_GeneDF)),]

heatmaply(bc_GeneDF[apply(bc_GeneDF, 1, var)!=0, ],
          col_side_colors = PhenoData[, c("Phenotype", "GEOStudyID")],
          seriate = "mean", 
          row_dend_left = TRUE,
          plot_method = "plotly",
          main = "Gene barcode heatmap")

DEG = unique(c(rownames(NormIntermDE), rownames(IntermCancerDE), rownames(IntermCancerDE)))
heatmaply(bc_GeneDF[which(rownames(bc_GeneDF) %in% DEG), ],
          col_side_colors = PhenoData[, c("Phenotype", "GEOStudyID")],
          seriate = "mean", 
          row_dend_left = TRUE,
          plot_method = "plotly",
          main = "Gene barcode heatmap")


