if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affyPLM")

memory.limit(9999999999)

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

weightImages = function(Pset, Path){
  QAPath = paste(Path,"/PseudoImages_Weights", sep = "")
  if(!dir.exists(QAPath)){dir.create(QAPath)}
  
  if(length(list.files(QAPath)) == 0){
    for(i in 1:Pset@narrays){
      name = paste(QAPath, "/weight_image",i,".jpg",sep="")
      jpeg(name)
      image(Pset,which=i)
      dev.off()
    }
  }
}


load("../../01_Data/02_Normalized data/04_Gastric/DataInfo.RData")
load("../../01_Data/02_Normalized data/04_Gastric/RawData.RData")

i=4

Data=RawExp[[i]]
class(Data)

Path = paste(DataInfo[[i]]$Path,"/RawData_QC", sep = "")
if(!dir.exists(Path)){dir.create(Path)}

PSet = affyPLM::fitPLM(Data, output.param=list(varcov="none"))

jpeg(paste(Path, "/RLE.jpg",sep=""))
affyPLM::RLE(PSet)
dev.off()
jpeg(paste(Path, "/NUSE.jpg",sep=""))
affyPLM::NUSE(PSet)
dev.off()




weightImages(PSet, Path)
residualImages(PSet, Path)
