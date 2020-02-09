######################################################################################################
library(DEXSeq)
library(BiocParallel)
library(htmltools)
######################################################################################################

##### Load Input files
# list exon count files
countFiles_all = list.files(path = "Input/Exon_Count_Files/", pattern=glob2rx("*.txt"), full.names=TRUE)
basename(countFiles_all)
# list gene annotation file
flattenedFile = list.files(path = "Input/", pattern="gff$", full.names=TRUE)
basename(flattenedFile)
# read meta information file
sampleTable_all = read.table(file = "Input/Ath_timeSeries_metadata.csv", header = T, sep = ",")

######################################################################################################

##### Calculate DEU for each condition (treatment*time)
# Define time points and organ tested
con <- unique(sampleTable_all$time)
org <- c("R", "S")

# Run CalcDEU for each condition
for (u in 1:length(org)){
  for (i in 1:length(con)){
    sampleTable <- sampleTable_all[sampleTable_all$time == con[i] & sampleTable_all$org == org[u],]
    countFiles <- subset(countFiles_all, countFiles_all %in% paste("./", sampleTable$Run, ".txt", sep=""))
    sampleTable$Condition <- factor(sampleTable$Condition)
    
    # Creating output file names
    name1= paste("N_", org[u], "_", con[i], "min_condition_DXResult_11.26.19.csv", sep="")
    name2= paste("N_", org[u], "_", con[i], "min_HTML_11.26.19.html", sep="")
    folder= paste(getwd(), org[u], con[i], sep="/")
    # creating folder to save output
    dir.create(path = paste(getwd(), org[u], con[i], sep="/"))
    # Calculating DEUs and saving output
    CalcDEU(countFiles=countFiles, sampleTable=sampleTable, flattenedFile=flattenedFile, workers=12, FDR=0.1)
    write.table(nitrogen_dxr1, file=name1, quote=F, sep=",",row.names=T, col.names=T)
    nitrogen_DEXSeqHTML = DEXSeqHTML(object = nitrogen_dxr1, fitExpToVar = "Condition",
                                     path=folder , file=name2, FDR=FDR, color=c("#FF000080", "#0000FF80"),BPPARAM=snow)
  }
}