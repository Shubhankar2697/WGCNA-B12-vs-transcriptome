library(metan)
getwd()
setwd("/B12_WGCNA/nutshell/")
traitData = read.csv("clin_var.csv");
dim(traitData)
names(traitData)
Wpheno<-column_to_rownames(traitData, var = "child_no")
Wpheno
All <- corr_coef(Wpheno)
All
plot(All)
