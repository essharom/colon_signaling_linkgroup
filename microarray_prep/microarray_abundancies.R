library("affy") # bioconductor package for affy arrays and normalisation 
library("hgu133plus2.db") #GPL570 platform

library("hgu133a2.db") #GPL571 platform
library("biomaRt")
# biomart annotation (creating uniprotac)
ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl)
annot_fullB = getBM(
  attributes = c("affy_hg_u133a_2", "uniprotswissprot"),
  mart = ensembl)


line_dir = paste0("C://Users/sebes/Downloads/GSE39491_RAW","/GSM969967_BE0002-1-A1.CEL/")
line_cels = paste0(line_dir,list.files(path=line_dir, pattern = "CEL$|cel$")) # cels vector for the selected cell line

cels = c(ctrl_cels,line_cels)

cels=line_cels


raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname="HG-U133A_2") # read the row data to raw.data
data=rma(raw.data) # mas5(); normalise the raw data 
data = exprs(data)	# get the normalised data in a readable matrix
data=format(data, digits=5) # keep only the firs five digit after the decimal point
probes=row.names(data) # get the probe names
#EntID = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA)) # get the ENTID-is from the probe names- only 90% is found
ControlC=cbind(probes, data) # bind the probes-entID-expression_data to a matrix

mapping_bm <- getBM(attributes = c("HG-U133A_2",
                                  "ensembl_transcript_id"),
                    filters = "affy_hg_u95av2",
                    values = affyids, mart = ensembl)