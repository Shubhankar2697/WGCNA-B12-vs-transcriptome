library(tibble)
library(WGCNA)
options(stringsAsFactors = FALSE)
getwd()
setwd("C:/Priya")
gene1 <-read.csv("neonatal_expr_data_yy.csv")
#read the clinical feature data
#take a quick look on what is the dataset.
dim(gene1)
names(gene1)
#fix(gene1)
#fix(gene2)
gene1<-column_to_rownames(gene1, var = "gene_symbol")

t_gene1<-t(as.data.frame(gene1))
rownames(t_gene1)<-sub("*\\X", "", rownames(t_gene1))
#Checking data for excessive missing values and identification of outlier microarray samples

gsg = goodSamplesGenes(t_gene1, verbose = 3);
gsg$allOK
#Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers
sampleTree = hclust(dist(t_gene1), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 2e+07, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight =2e+07 , minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = t_gene1[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#loading clinical cluster data
traitData = read.csv("B12_79.csv");
dim(traitData)
names(traitData)
Wpheno<-column_to_rownames(traitData, var = "child_no")
Wpheno

# Form a data frame analogous to expression data that will hold the clinical traits.
Samples<-rownames(datExpr)
view(Samples)
rownames(datExpr)
traitRows = match(Samples,rownames(Wpheno))
view(traitRows)
datTraits = Wpheno[traitRows,]
collectGarbage();


#regrouping samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(Wpheno, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(Wpheno),
                    main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "cordEXP.RData")
lnames = load(file ="cordEXP.RData")
lnames
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
net = blockwiseModules(datExpr, power = 7, networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Cord B12 expr",
                       verbose = 3)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Signed_7.RData")
# Load network data saved in the second part.
vnames = load(file = "Signed_7.RData");
vnames
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 15,6, 6));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#Define variable weight contaning the weight column of datTrait
B12 = as.data.frame(datTraits$cord.blood..B12.levels);
names(B12) = "B12"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
geneTraitSignificance = as.data.frame(cor(datExpr, B12, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(B12), sep="");
names(GSPvalue) = paste("p.GS.", names(B12), sep="");


########MM-GS for yellow  module which is positively associated with B12 levels########
module = "yellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for B12 levels",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

########MM-GS for turquoise module which is positively associated with B12 levels#######
module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for B12 levels",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


###########MM-GS for brown module which is negatively associated with B12 levels###
module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for B12 levels",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

############# To extract modules associated with cord blood B12 levels##########
yellow = colnames(datExpr1)[moduleColors=="yellow"]
write.csv(yellow, file = "Yellow_b12_signed.csv")

turquoise = colnames(datExpr1)[moduleColors=="turquoise"]
write.csv(turquoise, file = "Turquoise_b12_signed.csv")

brown = colnames(datExpr1)[moduleColors=="brown"]
write.csv(brown, file = "Brown_B12_signed.csv")
