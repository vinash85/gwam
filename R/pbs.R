# 5000 differetional genes
library(data.table)
library(avinash)
load("~/shortcuts/Data/MikeData/MAGnet_combat_eset.Rdata")
require(Biobase)
require(BiocGenerics)
require(parallel)
expressions <- exprs(eset.combat)
expression.sub <- as.data.table(t(expressions))
genes.col <- paste0("X", colnames(expression.sub))
setnames(expression.sub,seq(ncol(expression.sub)), genes.col)
expression.sub$sample_name <- colnames(expressions)
setkey(expression.sub, sample_name)
physiological.tab <- fread("~/project/gtps/gtps/doc/magnet_redcap.txt")
setkey(physiological.tab, sample_name)
expression.annot <- merge(x=expression.sub, y=physiological.tab, by="sample_name")
expression.annot[,disease.p:=ifelse(etiology=="donor",1,2)]
expression.annot[,history_of_afib_aflutter.p:=sapply(history_of_afib_aflutter, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
expression.annot[,history_of_diabetes.p:=sapply(history_of_diabetes, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
expression.annot[,history_of_hypertension.p:=sapply(history_of_hypertension, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]

library(limma)
exp.mat  <- as.matrix(expressions)
design <- cbind(grp1=1, disease=expression.annot$disease.p) 
fit  <-  lmFit(exp.mat, design)
# Fold-change thresholding
fit2 <- treat(fit,lfc=0.1)
topTreat(fit2,coef=2)
fit2 <- eBayes(fit2) 
topGene <- topTreat(fit2,coef=2,10000)
topGene <- topGene[topGene$adj.P.Val < 1e-5,] # 6670 genes with adj p.value  < 1e-5



##################################################################
#find in donors (and heart failures) SNP explaining variance by lasso
##################################################################
expression.diff = t(expressions[rownames(topGene),])
exp.diff.donor = expression.diff[expression.annot$disease.p==1,]
load("/fs/sh-project/Projects/Cardio/Data/MikeData/MAGnet_eQTL_impute/genotype1.RData")
genotype.sample = rownames(genotype)
genotype = genotype[colnames(expressions),]
source("/cbcbhomes/vinash85/project/gwam/R/lasso.eqtl.R")
exp.diff.donor = expression.diff[expression.annot$disease.p==1,]
genotype.donor = genotype[expression.annot$disease.p==1,]
donor.eqtl = lasso.EQTL(x=genotype.donor, y=t(exp.diff.donor), threads=64)

exp.diff.disease = expression.diff[expression.annot$disease.p==2,]
genotype.disease = genotype[expression.annot$disease.p==2,]
disease.eqtl = lasso.EQTL(x=genotype.disease, y=t(exp.diff.disease), threads=64)
save(file="eqtls.RData", donor.eqtl, disease.eqtl)
