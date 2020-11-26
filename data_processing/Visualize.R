library(RColorBrewer)
library(ggfortify)
library(cowplot)
library(Biobase)
library(pvca)
library(plyr)
library(snm)
library(sva)

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

PCAPlot = function(Data, PhenoData = NULL, colorCol, shapeCol = 10, Title){
  if(class(Data) == "list")
  {
    for(i in 1:length(Data)){
      if (class(Data[[i]]) %in% c("AffyBatch", "ExpressionSet")){
        Expr = exprs(Data[[i]])
        Pheno = pData(Data[[i]])
      } else {
        Expr = Data[[i]]
        Pheno = PhenoData[[i]]
      }
      pca = stats::prcomp(t(Expr), 
                          scale. = TRUE)
      autoplot(pca,
               data = Pheno[which(rownames(Pheno) %in% colnames(Expr)),], 
               colour = colorCol,
               shape = shapeCol, 
               title = unique(Pheno[, "GEOStudyID"]))
    } 
    
  } else {
    pca = stats::prcomp(t(Data), 
                        scale = TRUE)
    autoplot(pca, 
             data = PhenoData[which(rownames(PhenoData) %in% colnames(Data)),], 
             colour = colorCol, 
             shape = shapeCol, 
             title = Title)
  }
}



#############################################################################
#                                                                           #
#           Run principal component analysis & plot results                 #
#                                                                           #
#############################################################################

PVCA <- function (Data, PhenoData, Covariates, Threshold = 0.75, title) {
  if (class(Data)[1] == "data.frame") {
    expSet <- Norm(Data, PhenoData, normMethod = "rma", backgroundCorrection  = F, normalization = F)
    
  } else if (class(Data)[1] == "ExpressionSet") {
    expSet <- Data
  }
  PVCAObj <- pvcaBatchAssess(abatch = expSet, threshold = Threshold, batch.factors = Covariates)
  
  bp <- barplot(PVCAObj$dat, ylab = expression(bold("Weighted Average Proportion Variance")), 
                ylim = c(0, 1.1), col = c("navy"), las = 2, font = 2, 
                main = paste("PVCA of", title, sep = " "))
  axis(1, at = bp, labels = PVCAObj$label, cex.axis = 0.6, 
       las = 2, font = 2)
  values <- PVCAObj$dat
  new_values <- round(values, 3)
  text(bp, PVCAObj$dat, labels = new_values, pos = 3, cex = 0.8, font = 2)
  
  PVCAObj
}