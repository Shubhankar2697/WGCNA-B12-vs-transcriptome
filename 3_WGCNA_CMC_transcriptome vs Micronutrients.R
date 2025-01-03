## Author: Shubhankar A Pawar. Contact: shubhankarpawar3@gmail.com
getwd()
setwd("/B12_WGCNA/nutshell/")
gene1 <-read.csv("vst_no__batch_correction.csv")

library(tibble)
library(WGCNA)
options(stringsAsFactors = FALSE)

dim(gene1) #read the clinical feature data
names(gene1) #take a quick look on what is the dataset.


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
abline(h = 140, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight =140 , minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = t_gene1[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)
#write.csv(datExpr, "expr_out_vst.csv")
traitData = read.csv("clin_var.csv");#loading clinical cluster data
dim(traitData)
names(traitData)
Wpheno<-column_to_rownames(traitData, var = "child_no")
Wpheno
# Form a data frame analogous to expression data that will hold the clinical traits.
Samples<-rownames(datExpr)
#view(Samples)
rownames(datExpr)
traitRows = match(Samples,rownames(Wpheno))
#view(traitRows)
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

save(datExpr, datTraits, file = "cordEXP_vst_nb_new.RData")
#--------------


lnames = load(file ="cordEXP_vst_nb_new.RData")
lnames
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, 
                        networkType = "signed", 
                        enableWGCNAThreads(nThreads = 15), 
                        verbose = 5)
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
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



cor<-WGCNA::cor
net = blockwiseModules(datExpr, power = 7, networkType = "signed",
                       TOMType = "signed", minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Cord_varriant_power7_expr_vst",
                       verbose = 3,
                       nThreads = 15)

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
     file = "Signed_7_vst_nb.RData")

# Load network data saved in the second part.
vnames = load(file = "Signed_7_vst_nb.RData");
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

# Adjust margins and text size
par(mar = c(8, 10,3, 5));
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
##################################################################

B12 = as.data.frame(datTraits$Cord.blood..B12)
names(B12) = "B12"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
geneTraitSignificance = as.data.frame(cor(datExpr, B12, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(B12), sep="");
names(GSPvalue) = paste("p.GS.", names(B12), sep="");
#########   MM-GS ####### B12 associated clusters
# for power 14)
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



