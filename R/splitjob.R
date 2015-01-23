args=(commandArgs(trailingOnly=TRUE))
curinx=eval(parse(text=args[1]))
library(glmnet)
load("expressions.RData")
exp.diff.donor = exp.diff.donor[,curinx]
exp.diff.disease = exp.diff.disease[,curinx]
load("genotype.RData") # contains genotype.donor and genotype.disease
# genotype.donor=genotype[donor.inx,]
# genotype.disease=genotype[disease.inx,]
gc()
source("/cbcbhomes/vinash85/project/gwam/R/lasso.eqtl.R")
donor.eqtl = lasso.EQTL(x=genotype.donor, y=exp.diff.donor, nthreads=10)
disease.eqtl = lasso.EQTL(x=genotype.disease, y=exp.diff.disease, nthreads=10)
save(file=paste0("eqtlsDir/eqtls",curinx[1],curinx[length(curinx)], ".RData"), donor.eqtl, disease.eqtl)