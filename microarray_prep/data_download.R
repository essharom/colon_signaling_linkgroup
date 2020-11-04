memory.limit(9999999999)

source('helper_functions.R')


library(GEOquery)


#file download
input_files = read.delim(file = "Data.txt")
input_files[input_files=="N/A"] = ""
input_files[is.na(input_files)] = ""

for(i in 1:nrow(input_files)){
  GSE = input_files[i,1]
  type = input_files[i,2]
  field = input_files[i,3]
  normal = strsplit(input_files[i,4], "\\|")
  intermediate = strsplit(input_files[i,5], "\\|")
  cancer = strsplit(input_files[i,6], "\\|")

  # download GEO files
  df = getGEO(GSE, getGPL = FALSE, parseCharacteristics = FALSE)
  getGEOSuppFiles(GSE)
  
  # file unzipping
  untar(paste(GSE, "/", GSE, "_RAW.tar", sep = ""), exdir = "temp/")
  cells = list.files(path = "temp/", pattern = "(CEL|cel)(.gz)?$")
  
  for(j in 1:length(df)){
    pheno = df[[j]]@phenoData@data
    GPL = unique(pheno[['platform_id']])
    if(any('Homo sapiens'!=pheno[['organism_ch1']])){next}
    if(GPL == "GPL6801" || GPL == "GPL7723"){next}
    
    if(!dir.exists(file.path("./myData/", type, "/", GSE))){
      dir.create(paste("./myData/", type, "/", GSE, sep = ""))
    }
    
    data_download(type, GSE, GPL, "/Normal/", pheno, field, normal, cells)
    data_download(type, GSE, GPL, "/Intermediate/", pheno, field, intermediate, cells)
    data_download(type, GSE, GPL, "/Cancer/", pheno, field, cancer, cells)
  }
  remove(df)
  unlink("./temp", recursive = TRUE)
  unlink(paste("./", GSE, sep = ""), recursive = TRUE)
}
