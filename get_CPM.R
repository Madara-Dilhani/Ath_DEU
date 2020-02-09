########################################################################################################
## Title: R modules for Generating CPM values from featureCount text files 
## Author: Madara Hetti-Arachchilage
## date: "February 3, 2020"
########################################################################################################
# Defining 'load_pkg' function here
load_pkg <- function (pkg_name, bioconductor = FALSE) {
  # This function checks if the requested package is already installed in the
  # environment. If not, the package is installed from the appropriate source
  # and then loaded into the workspace
  
  if (!pkg_name %in% installed.packages()) {
    if (bioconductor == TRUE) {
      source("https://bioconductor.org/biocLite.R")
      biocLite(pkg_name)
    } else {
      install.packages(pkg_name, quiet = T, verbose = F)
    }
  }
  if (pkg_name %in% installed.packages()) {
    library(pkg_name, character.only = TRUE, quietly = T, verbose = F)
  } else {
    stop(paste("Cannot install/load package", pkg_name, sep = ": "))
  }
}
  
# Function to read count files in a given folder with given pattern for filenames
CountMat <- function(files){
  
  print("Preparing count matrix.....")
  #Merge the first to files and store
  file1 = read.table(files[1], header = T)[,c(1,7)]
  file2 = read.table(files[2], header = T)[,c(1,7)]
  out.file <- merge(file1, file2, by= "Geneid")
  #For loop to merge contents of remaining files
  for(i in 3:length(files)){
    file = read.table(files[i], header = T)[,c(1,7)]
    out.file <- merge(out.file, file, by=c("Geneid"))
  }
  print("Preparing count matrix.....Done")
  return(out.file)
}


#  Defining the Function to get normalized CPM values. This function takes count matrix, meta file
calc_cpm <- function(count_matrix, Sample){
  
  y <- DGEList(counts = count_matrix, group= Sample)
  dge <- calcNormFactors(y, method= "TMM")
  cpm <- cpm(dge, log = FALSE)
  return(cpm)
  
}

# Defining the function to Draw CPM value distribution to decide a CPM threshold for expressed/not expressed gene 
plotCPM <- function(expr_norm) {
  # Convert CPM values to log-CPM for visualization
  exp_log <- log(expr_norm)
  # Convert count matrix to two column table for plotting
  melted_expr <- melt(exp_log)
  # Density plot for distribution of CPM values. Generally, this plot gives a bimodal distibution which can be
  # used to decide a CPM threshold for expressed gene 
  p <- qplot(value, geom = "density", data = melted_expr) +
    stat_function(fun = dnorm, size = 0.5, color = 'red') +
    xlab("Standardized log(FPKM)") +
    ylab("Density")
  return(p)
}

#############################################################################################