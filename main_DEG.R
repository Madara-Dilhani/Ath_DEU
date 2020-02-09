########################################################################################################
## title: "Differential Gene Expression Analysis"
## author: "Madara Hetti-Arachchilage; mhettiar 'at' illinois.edu"
## date: "February 3, 2020"
########################################################################################################

##### Importing libraries and custom scripts 
# Loading R source code which has defined custom functions needed for later steps
# Loading R/Biconductor packages (checks if the requested package is already installed in 
# the environment. If not, the package is installed from the appropriate source and then 
# loaded into the workspace

library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(VennDiagram)
library(HTSFilter)

########################################################################################################
##### Step 1: Exploratory analysis on raw read counts 
## Import Data 
rawCountTable <- read.csv("Output/Ath_rawCount.csv", header=TRUE, row.names=1)
sampleInfo <- read.table("Input/At_timeSeries_metadata.csv", header=TRUE, sep=",")

## Check read Data
head(rawCountTable)
nrow(rawCountTable)

## Check Experimental Design
head(sampleInfo)
sampleInfo$Sample <- as.factor(sampleInfo$Sample)

## Create EdgeR DGE Object
dgeFull <- DGEList(rawCountTable, remove.zeros = TRUE, group = sampleInfo$Sample)
dgeFull

## Visualize exploratory analysis
## pseudoCounts
pseudoCounts <- log2(dgeFull$counts + 1)
head(pseudoCounts)

##  histogram for PseudoCounts
hist(pseudoCounts[ ,"Cond.T1.Rep.1"], main = "", xlab = "counts")

## Box Plot for PseudoCounts
par(mar = c(8,4,1,2))
boxplot(pseudoCounts, col = "gray", las = 3, cex.names = 1)

## MA Plot for PseudoCounts
limma::plotMA(pseudoCounts[ ,1:2], xlab = "M", ylab = "A", main = "")
abline(h = 0, col = "red")

## MDS from PseudoCounts
colConditions <- brewer.pal(10, "Set3")
colConditions <- colConditions[match(sampleInfo$time,
                                     unique(sampleInfo$time))]
pchGenotypes <- c(8, 15, 16)[match(sampleInfo$organ,
                                   levels(sampleInfo$organ))]
plotMDS(pseudoCounts, pch = pchGenotypes, col = colConditions) + 
  legend("bottomleft", lwd = 2, col = brewer.pal(10, "Set3"), legend = unique(sampleInfo$time))+  
  legend("bottomright", pch = c(8, 15, 16), legend = levels(sampleInfo$organ))

## CIM from PseudoCounts
sampleDists <- as.matrix(dist(t(pseudoCounts)))
sampleDists
cimColor <- colorRampPalette(rev(brewer.pal(9, "Reds")))(16)
cim(sampleDists, color = cimColor, symkey = FALSE, row.cex = 0.7, col.cex = 0.7)
########################################################################################################

##### Step 2: Exploratory analysis on normalized read counts
## estimate Normalization Factors
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull
head(dgeFull$counts)

## estimate Normalized Counts
normCounts <- cpm(dgeFull, log = TRUE)
write.csv(normCounts, "featureCounts/Ath_Normalized_log-CPM.csv", quote = F, row.names = T)
pseudoNormCounts <- cpm(dgeFull, log = TRUE, prior.count = 1)
par(mar = c(8,4,1,2))
boxplot(pseudoNormCounts, col = "gray", las = 3, cex.names = 1)

## MDS for normalized counts
plotMDS(pseudoNormCounts, pch = pchGenotypes, col = colConditions)
legend("topright", lwd = 2, col = brewer.pal(3, "Set2")[1:2], 
       legend = levels(sampleInfo$condition))
legend("bottomright", pch = c(8, 15, 16), 
       legend = levels(sampleInfo$genotype))

##############################################################################################

##### Step 3: Find DEGs
## Filter lowly expressed genes
exp_log <- log(normCount)
melted_expr <- melt(exp_log)
p <- qplot(value, geom = "density", data = melted_expr) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(TPM)") +
  ylab("Density")

activeGenes <- rowSums(normCount > 0.5) >= nrow(sampleInfo)/2
normCount <- normCount[activeGenes,]

## DGEList Creation
dgeFull.group <- DGEList(normCount, remove.zeros = TRUE, 
                         group = dgeFull$samples$group)
dgeFull$samples$genotype <- sampleInfo$genotype
dgeFull
dgeFull.group$samples$norm.factors <- dgeFull$samples$norm.factors
dgeFull.group

## estimate Dispersion
dgeFull.group <- estimateCommonDisp(dgeFull.group)
dgeFull.group <- estimateTagwiseDisp(dgeFull.group)
dgeFull.group

## plot BCV
plotBCV(dgeFull.group)

## fisher Exact
dgeExactTest <- exactTest(dgeFull.group, pair=c("R0","R5"))
dgeExactTest

## get top Tag
resExactTest <- topTags(dgeExactTest, n = nrow(dgeExactTest$table))
head(resExactTest$table)

## histogram of PVal from Fisher Exact test
par(mfrow = c(1,2))
hist(resExactTest$table$PValue, xlab = "p-value", main = "raw p-values")
hist(resExactTest$table$FDR, xlab = "p-value", main = "adjusted p-values")

## DEGs from Fisher Exact
selectedET <- resExactTest$table$FDR < 0.05 & abs(resExactTest$table$logFC) > 2
selectedET <- resExactTest$table[selectedET, ]
nrow(selectedET)
head(selectedET)

## Up- and Down- regulated genes
selectedET$updown <- factor(ifelse(selectedET$logFC > 0, "up", "down"))
head(selectedET)

## Save DEG list
write.table(selectedET, file = "AthDEG.csv", sep = ",")

## volcanoPlot
volcanoData <- cbind(selectedET$logFC, -log10(selectedET$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
DEGs <- selectedET$FDR < 0.05 & abs(selectedET$logFC) > 2
point.col <- ifelse(DEGs, "red", "black")
plot(volcanoData, pch = 16, col = point.col, cex = 0.5)
###########################################################################


