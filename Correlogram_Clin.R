library(metan)
getwd()
setwd("C:/Priya")
traitData = read.csv("Clinical.csv");
dim(traitData)
names(traitData)
Wpheno<-column_to_rownames(traitData, var = "child_no")
Wpheno
All <- corr_coef(Wpheno)
All
plot(All)
