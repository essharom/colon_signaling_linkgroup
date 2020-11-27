#############################################################################
#                                                                           #
#                         Phenotypic data summary                           #
#                                                                           #
#############################################################################

PhenoSummary = function(Pheno){
  GSE = unique(Pheno[, "GEOStudyID"])
  GPL = SampleNo = NormalNo = IntermNo = CancerNo = c()
  SampleNo = c()
  for(i in 1:length(GSE)){
    temp = subset(Pheno, eval(expression(GEOStudyID == GSE[i])))
    Platform = unique(temp[, "GEOPlatformID"])
    if(length(Platform)>1){Platform = paste(Platform, collapse = "|")}
    GPL = append(GPL, Platform)
    SampleNo = append(SampleNo, length(temp[, "GEOStudyID"]))
    NormalNo = append(NormalNo, sum(temp[, "Phenotype"] == "Normal"))
    IntermNo = append(IntermNo, sum(temp[, "Phenotype"] == "Intermediate"))
    CancerNo = append(CancerNo, sum(temp[, "Phenotype"] == "Cancer"))
  }
  PhenoDF = data.frame(GSE, GPL, SampleNo, NormalNo, IntermNo, CancerNo)
  
  sumrow = as.data.frame(lapply(PhenoDF, function(z) if (is.numeric(z)) sum(z) else ''))
  sumrow[1] <- "Total"
  
  PhenoDF = rbind(PhenoDF, sumrow)
  names(PhenoDF) = c("GEO study ID", "GEO Platform ID", "# of samples",
                     "# of normal samples", "# of intermediate samples", 
                     "# of cancer samples")
  
  PhenoDF
}

#############################################################################
#                                                                           #
#                         Expression data summary                           #
#                                                                           #
#############################################################################


#############################################################################
#                                                                           #
#                         Combine expression data                           #
#                                                                           #
#############################################################################

List2DF = function(Data){
  if("matrix" %in% class(Data[[1]])){
    expr = Data
  } else {
    expr = purrr::map(Data, exprs)
  }
  
  expr = map2(purrr::map(expr, rownames), expr, cbind) %>% 
    imap(~  {colnames(.x)[1] <- "RowNames"; .x}) %>% 
    purrr::map(data.frame) %>% 
    purrr::reduce(left_join, by = "RowNames")
  rnames = expr$RowNames
  expr = data.frame(lapply(expr[,-1], function(x) as.numeric(as.character(x))))
  rownames(expr) = rnames
  
  expr
}
