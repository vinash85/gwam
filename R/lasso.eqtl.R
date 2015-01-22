lasso.EQTL <- function (x, y, gene.map=NULL, nthreads = 1)
{   
  require(glmnet)
  require(doMC)
  require(foreach)
  registerDoMC(cores = nthreads)
  n <- nrow(x)
  if(is.null(gene.map)){
    genes <- seq(ncol(y))
    coefs <- foreach(gene = genes, .inorder = T) %dopar% {
      fit <- cv.glmnet(x = x, y = y[, gene], alpha = 0.5,
                       parallel = T)
      coef = coef(fit, s = "lambda.min")
      coef = coef[coef[,1]!=0,]
      coef
    }
#     coefs <- list()
#   for(gene in genes) {
#       fit <- cv.glmnet(x = x, y = y[, gene], alpha = 1,
#                        parallel = F)
#       coef = coef(fit, s = "lambda.min")
#       coef = coef[coef[,1]!=0,]
#       coefs[[gene]] = coef
#     } 
  }else{
  genes <- unique(gene.map)
  coefs <- foreach(gene = genes, .inorder = T) %dopar% {
    curr.inx <- which(gene.map == gene)
    fit <- cv.glmnet(x = x[, curr.inx],  y = y[, gene], alpha = 0.5, nfolds=20,
                     parallel = T)
    coef = coef(fit, s = "lambda.min")
    coef <- coef[-1]
    coef
  }
 
  }
  coefs
}
