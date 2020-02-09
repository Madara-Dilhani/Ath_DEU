###############################################################################################################

# Defining the function to calculate differential exon usage
# This function takes a list of count files, dataframe containing meta information, gene annotation file, 
# no. of cores to be used, false discovery rate threshold
CalcDEU <- function(countFiles, sampleTable, flattenedFile, workers, FDR){
  suppressPackageStartupMessages( library( "DEXSeq" ) )
  n_dxd1 = DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData=sampleTable,
    design= ~ sample + exon + Condition:exon,
    flattenedfile=flattenedFile )
  
  n_dxd1 = na.omit(n_dxd1)
  snow <- SnowParam(workers)
  options(bphost="localhost")
  
  n_dxd1 = estimateSizeFactors( n_dxd1 )
  n_dxd = estimateDispersions( n_dxd1, BPPARAM=snow)
  n_deu = testForDEU( n_dxd, BPPARAM=snow)
  n_exf = estimateExonFoldChanges( n_deu, fitExpToVar="Condition", BPPARAM=snow)
  n_dxr1 = DEXSeqResults(n_exf)
  return(n_dxr1)
}
