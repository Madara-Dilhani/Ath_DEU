########################################################################################################
## title: "RNA-seq data preprocessing"
## author: "Madara Hetti-Arachchilage; mhettiar 'at' illinois.edu"
## date: "February 3, 2020"
########################################################################################################

##### Importing libraries and custom scripts 
# Loading R source code which has defined custom functions needed for later steps
# Loading R/Biconductor packages (checks if the requested package is already installed in 
# the environment. If not, the package is installed from the appropriate source and then 
# loaded into the workspace

source("get_CPM.R")
load_pkg("edgeR", bioconductor = TRUE)
load_pkg("ggplot2", bioconductor = TRUE)
load_pkg("factoextra", bioconductor = FALSE)
load_pkg("matrixStats", bioconductor = FALSE)
########################################################################################################

##### Read all featureCount files in the given directory and save count matrix into a .csv file
# Calling "CountMat" function which reads multiple read count text files generated from featureCounts tool 
# and save in to count matrix, where rows are genes and columns are samples. This function takes 2 arguments, 
# folder name where featurecount files are saved and given pattern for filenames

FeatureC_File <- "Output/Ath_rawCount.csv"
folder = "Input/featureCounts/"
pattern = "*counts.txt"
files = list.files(folder, pattern, full.names = T)
out.file <- CountMat(files)

print("Renaming columns....")
df <- colnames(out.file)[-1]
cc <- sapply(strsplit(sapply(strsplit(df, "_paired_prinseq_good_hisat2"), "[[",1), "scratch.MH_Data.Sorted_Bam."), "[[",2)
meta <- read.csv("/At_timeSeries_metadata.csv", header = T, as.is = T)
ff <- meta$Sample[match(cc, meta$Run, nomatch = 0)]
colnames(out.file) <- c("Geneid", cc)
print("Renaming columns....Done")
write.csv(x = out.file, file = FeatureC_File, row.names = FALSE, quote = F)
########################################################################################################

## Normalizing featureCount raw read matrix for experimental samples and saving CPM values to a.CSV file
# The "calc_cpm" function uses "edgeR' package to normalize featureCounts raw reads to TMM values and 
# then transform into CPM values. This function takes 2 arguments, raw count matrix and a table of 
# sample information (Control and Treatment group should be listed in "Sample" Column).

CPM_Count_File <- "Output/Ath_NormCount.csv"
out.file <- read.csv("Output/Ath_rawCount.csv", header = T, as.is = T)

row.names(out.file) <- out.file[,1]
out.file <- out.file[,-1]
expr_norm <- calc_cpm(out.file, meta$Sample)
write.csv(x = expr_norm, file = CPM_Count_File, quote=F, row.names = T)
########################################################################################################

##### This step perform hierachical clustering and plot the dendogram to visualize sample similarity
# Filtering not expressed genes
df <- expr_norm[rowSums(expr_norm > 5) >= 1, ]
# Log transformation
df[df == 0] <- 0.0000000001
df <- log2(df)
df[df < 0] <- 0
# Quantile normalization for across tissue comparison
df_q <- preprocessCore::normalize.quantiles(as.matrix(df))
row.names(df_q) <- row.names(df)
colnames(df_q) <- colnames(df)

# Optional; Filtering top highly variable genes 
rv <- rowVars(x = as.matrix(df_q))
select <- order(rv, decreasing=T)[seq_len(min(5000,length(rv)))]
df_q <- df_q[select,]

# Hierachial clusterin on tanspose of scaled expression values to see sample clustering
res.hc <- eclust(t(df_q), "hclust")
# Create dendogram
fviz_dend(res.hc, rect = F, type = "rectangle", horiz = F, color_labels_by_k = TRUE, rect_border = "jco", 
          palette = "jco", labels_track_height = 1e+3, rect_lty = 2.5, cex = 1, rect_fill = F, 
          ggtheme = theme_grey(), lower_rect = -300) 
# Create Scatter plot
fviz_cluster(res.hc, labelsize = 7) 
########################################################################################################

##### Perform PCA
pca <- prcomp(t(df_q))
scores <- as.data.frame(pca$x)
summary(pca)
# Hierachical clustering on PCA 
res.hc <-  eclust(scores, "hclust")
# Create dendogram
fviz_dend(res.hc, rect = F, type = "rectangle", horiz = F, color_labels_by_k = TRUE, rect_border = "jco", 
          palette = "jco", labels_track_height = 1e+3, rect_lty = 2.5, cex = 1, rect_fill = F, 
          ggtheme = theme_grey(), lower_rect = -300) 
# Create Scatter plot
fviz_cluster(res.hc, labelsize = 7) 
########################################################################################################

