source('helper_functions.R')

library(GEOquery)


#file download
input_files = read.delim(file = "Data.txt")

for(i in 1:nrow(input_files)){
  GSE = input_files[i,1]
  type = input_files[i,2]
  normal = strsplit(input_files[i,3], "\\|")
  intermediate = strsplit(input_files[i,4], "\\|")
  cancer = strsplit(input_files[i,5], "\\|")
  
  # download GEO files
  df = getGEO(GSE, getGPL = FALSE, parseCharacteristics = FALSE)
  pheno = df[[1]]@phenoData@data
  getGEOSuppFiles(GSE)
  
  # file unzipping
  untar(paste(GSE, "/", GSE, "_RAW.tar", sep = ""), exdir = "myData")
  cells = list.files(path = "myData", pattern = "(CEL|cel)(.gz)?$")
  
  data_download(type, GSE, "/Normal/", pheno, normal, cells)
  data_download(type, GSE, "/Intermediate/", pheno, intermediate, cells)
  data_download(type, GSE, "/Cancer/", pheno, cancer, cells)
  
  unlink("./myData", recursive = TRUE)
  unlink(paste("./", GSE, sep = ""), recursive = TRUE)
}
