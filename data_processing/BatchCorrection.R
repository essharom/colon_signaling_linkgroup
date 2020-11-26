library(sva)

#############################################################################
#                                                                           #
#              Surrogate variable analysis batch correction                 #
#                                                                           #
#############################################################################

SVA <- function(Data, PhenoData, bioVar) {
  if (class(Data)[1] == "data.frame") {
    expMat = as.matrix(Data)
  } else if (class(Data)[1] == "ExpressionSet") {
    expMat = as.matrix(exprs(expSet))
  }
  
  PhenoData = PhenoData[which(rownames(PhenoData) %in% colnames(expMat)),]
  
  ## Preform SVA analysis
  mod <- model.matrix(~as.factor(PhenoData[, bioVar]), data = PhenoData)
  mod0 <- model.matrix(~1, data = PhenoData)

  SVAObject <- sva(expMat, mod, mod0)
  
  W <- as.matrix(SVAObject$sv)
  alpha <- solve(t(W) %*% W) %*% t(W) %*% t(expMat)
  corrected <- t(t(expMat) - W %*% alpha)
  
  corrected
}