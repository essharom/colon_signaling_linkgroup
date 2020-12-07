library(RColorBrewer)
library(htmlwidgets)
library(heatmaply)
library(ggfortify)
library(ggplot2)
library(Biobase)
library(plyr)
#############################################################################
#                                                                           #
#                             Create boxplots                               #
#                                                                           #
#############################################################################

BoxPlot = function(Data, Logaritmic = "y", Title){
  if(class(Data) == "list"){
    if(class(Data[[1]]) %in% c("AffyBatch", "ExpressionSet")){
      Data = purrr::map(Data, exprs)
    }
    ListLength = length(Data)
    BoxLoc = 1:3
    temp = max(unlist(purrr::map(Data, max)))
    xDim = c(0.5,ListLength*3 + 0.5)
    yDim = round_any(temp, 10^ceiling(log10(temp)), f = ceiling)
    cols = brewer.pal(n = ListLength, name = "RdBu")
    
  }
  boxplot(Data[[1]][,1:3],
          at = BoxLoc,
          xlim = xDim,
          ylim = c(1,yDim),
          log = Logaritmic,
          add = F,
          las = 3,
          col = cols[1])
  for(i in 2:ListLength){
    BoxLoc = BoxLoc + 3
    boxplot(Data[[i]][,1:3],
            at=BoxLoc,
            xlim=xDim,
            ylim=c(1,yDim),
            log=Logaritmic,
            add=T,
            las=3,
            col=cols[i])
  }
  title(ylab="expression value",main=Title)
}


#############################################################################
#                                                                           #
#                               PCA analysis                                #
#                                                                           #
#############################################################################

PCAPlot = function(Data, PhenoData = NULL, colorCol, shapeCol = NULL, Title = NULL){
  
  pca = stats::prcomp(t(Data), scale = TRUE)
  if(is.null(shapeCol)){shapeCol=10}
  pca = autoplot(pca, 
            data = PhenoData[which(rownames(PhenoData) %in% colnames(Data)),], 
            colour = colorCol, 
            shape = shapeCol) 
  if(!is.null(Title)){
    pca + ggtitle(Title)
  } else {
    pca
  }
  
}


#############################################################################
#                                                                           #
#                             Raw intensity images                          #
#                                                                           #
#############################################################################

rawIntensityImages = function(Data, Path){
  QAPath = paste(Path,"/RawIntensityImages", sep = "")
  if(!dir.exists(QAPath)){dir.create(QAPath)}
  
  if(length(list.files(QAPath)) == 0){
    for(i in 1:length(colnames(Data))){
      name = paste(QAPath, "/raw_image",i,".jpg",sep="")
      jpeg(name)
      image(Data[,i])
      dev.off()
    }
  }
}


#############################################################################
#                                                                           #
#                   Array pseudo-image based on weights                     #
#                                                                           #
#############################################################################
# weights in probe-level model represent how much the original data contribute 
# to the model,with outlines being strongly downregulated

weightImages = function(Pset, Path){
  QAPath = paste(Path,"/PseudoImages_Weights", sep = "")
  if(!dir.exists(QAPath)){dir.create(QAPath)}
  
  if(length(list.files(QAPath)) == 0){
    for(i in 1:Pset@narrays){
      name = paste(QAPath, "/weight_image",i,".jpg",sep="")
      jpeg(name)
      oligo::image(Pset,which=i)
      dev.off()
    }
  }
}


#############################################################################
#                                                                           #
#                   Array pseudo-image based on weights                     #
#                                                                           #
#############################################################################
# residuals represent the difference between original and ideal data according 
# to probe-level model, with positive residuals indicating larger and negative 
# residuals indicating smaller intensities than ideal value according to the 
# model

residualImages = function(Pset, Path){
  QAPath = paste(Path,"/PseudoImages_Residuals", sep = "")
  if(!dir.exists(QAPath)){dir.create(QAPath)}
  
  if(length(list.files(QAPath)) == 0){
    for(i in 1:Pset@narrays){
      name = paste(QAPath, "/residuals_image",i,".jpg",sep="")
      jpeg(name)
      
      if(class(PSet)[1] == "oligoPLM"){
        oligo::image(Pset,which=i, type = "residuals")
      }else{
        image(Pset,which=i, type = "resids")
      }
      
      dev.off()
    }
  }
}


#############################################################################
#                                                                           #
#                               Density plots                               #
#                                                                           #
#############################################################################

histograms = function(logData, Path){
  dataHist = ggplot(logData, aes(logInt, colour = GEOSampleID)) + geom_density() 
  saveWidget(ggplotly(dataHist), file = file.path(normalizePath(Path, winslash ="/"), "DensityPlot.hthm"))
}


#############################################################################
#                                                                           #
#                                 Box plots                                 #
#                                                                           #
#############################################################################
# boxplots alow assesment of scale and distribuation of the data, with differences
# in scale and center of the boxes indicating that normalization is required

boxplots = function(logData, Path){
  dataBox = ggplot(logData,aes(GEOSampleID,logInt)) + 
    theme(axis.text.x = element_text(angle = 90))
  name = paste(Path, "/Boxplots.jpg",sep="")
  jpeg(name)
  dataBox + geom_boxplot()
  dev.off()
}


#############################################################################
#                                                                           #
#                                 Box plots                                 #
#                                                                           #
#############################################################################

DistanceHeatmap = function(Data, PhenoData, Path, colorCol, Title = NULL) {
  if(class(Data)[1] %in% c("AffyBatch", "GeneFeatureSet", "ExpressionSet")){
    Data = exprs(Data)
  } 
  if(is.null(Title)){Title = ""}
  Dist = dist(t(Data), diag = TRUE)
  PhenoData = PhenoData[which(rownames(PhenoData) %in% colnames(Data)),]
  
  heatmaply(as.matrix(Dist),
            col_side_colors = PhenoData[, colorCol],
            seriate = "mean", 
            row_dend_left = TRUE,
            plot_method = "plotly",
            main = Title, 
            file = paste(Path, "DistanceHeatmap.html", sep="/"))
  
}

#############################################################################
#                                                                           #
#                                 MA plots                                  #
#                                                                           #
#############################################################################
# MA plot for single channel arrays compares array probe intensities to 
# intensities on the median array and idealy the cloud of data points should 
# be centered around M=0, since we assume that majority of genes are not 
# differentially expressed and the number of up- and down-regulates genes is 
# assumed to be the same
# variablity of M should be similar for different A values (average intensities),
# but data normalization may remove some of the dependency of average intensites 
# on M values

MAPlot = function(Data, Path){
  QAPath = paste(Path,"/MAPlots", sep = "")
  if(!dir.exists(QAPath)){dir.create(QAPath)}
  
  if(length(list.files(QAPath)) == 0){
    for (i in 1:length(colnames(Data))){
      name = paste(QAPath, "/MAplot",i,".jpg",sep="")
      jpeg(name)
      if(class(Data)[1] %in% c("AffyBatch", "ExpressionSet")){
        affy::MAplot(Data, which = i)
      } else if (class(Data) %in% "GeneFeatureSet") {
        oligo::MAplot(Data, which = i)
      }
      MAplot(Data,which=i)
      dev.off()
    }
  }
}