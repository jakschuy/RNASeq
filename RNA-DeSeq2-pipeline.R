Sys.setenv(LANG="en")
rm(list=ls())
options(scipen=999)

#input section
NameOfControl = "ctrl"
ProjectName <- "RNA-DeSeq2-pipeline_all-runs_HG38"

#Working directory
setwd(paste(rstudioapi::getSourceEditorContext()$path, "/../..", sep = ""))

#create folder for extracted data like lists and dataframes
if (!dir.exists(paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/", sep=""))) {
  dir.create(paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/", sep=""))     }

#The differential expression analysis tools
suppressMessages(library("DESeq2"))

# #to install the packages use Bioconductor Manager
# BiocManager::install("tximport")

#Other functions for reading the salmon input
library("readr")
library("tximport")

#generate sample list
ListWithSpecifications <- list.files(path = "Input_files/Specifications", pattern = ".txt")

#load all genes-ENSG ID-connection
warning("Using short annotated tx2gene file!")
tx2gene <- read.table(file = "Input_files/tx2gene_hg38_short.csv", sep = ",", header = TRUE)

#loop for input of all samples
i <- 1
for(i in 1:length(ListWithSpecifications)){
  subdat <- read.table(file = paste("Input_files/Specifications/", ListWithSpecifications[i], sep=""), sep = " ", header = TRUE)
  rownames(subdat) <- subdat$run
  
  #reorder factors because of DeSeq2 uses factors ordered by alphabet
  PositionCTRL <- which(levels(factor(subdat$sample)) == NameOfControl)
  PositionNoCTRL <- which(levels(factor(subdat$sample)) != NameOfControl)
  subdat$sample <- factor(subdat$sample, levels(factor(subdat$sample))[c(PositionCTRL, PositionNoCTRL)])
  
  
  #variables defining specifications
  files <- file.path(subdat$path)
  names(files) <- subdat$run
  
  #prepare data structure before feeding into DeSeq2
  #= combining ENSG with data
  txi <- tximport(files, type="salmon", tx2gene=tx2gene)
  
  #run DE analysis
  condition = factor(subdat$sample)
  ExpDesign <- data.frame(row.names=colnames(txi$counts), condition = condition)
  ddsTxi <- DESeqDataSetFromTximport(txi, ExpDesign, ~condition)
  dds <- DESeq(ddsTxi)
  res <- results(dds)
  df_export <- as.data.frame(res)
  
  #label and export
  NameNoDim <- sub('\\.txt$', '', ListWithSpecifications[i])
  write.table(df_export, 
              file = paste("R stuff/ExtractedData_", ProjectName, "/", NameNoDim, ".tsv" , sep = ""), row.names = TRUE)
  # write.table(df_export,
  #           file = paste("R/ExtractedData/", NameNoDim, ".txt" , sep = ""))

}









