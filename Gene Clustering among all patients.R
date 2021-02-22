Sys.setenv(LANG="en")
rm(list=ls())
options(scipen=999)

# library(ggplot2)
# library(ggrepel)
 library(dplyr)
# library("tximport")
# library(reshape)
library(gplots)

#input settings
ProjectName <- "GeneClusterDelAli9"
RefOfSeqDat <- "CTRL"
BaseMean_threshold <- 25
colorsForLegend <- c("grey", "black", "dodgerblue", "black", "red", "green", "palegreen1", "palegreen4", "darkgoldenrod1", "darkorchid1")
thresholdForPval <- 0.05
thresholdForExpr <- 1.25
thresholdForBaseMean <- 25



#-------start----

#Working directory
setwd(paste(rstudioapi::getSourceEditorContext()$path, "/../..", sep = ""))

#create folder for output graphs
if (!dir.exists(paste(getwd(), "/R stuff/Graphs_", ProjectName, "/", sep="")))  {
  dir.create(paste(getwd(), "/R stuff/Graphs_", ProjectName, "/", sep=""))      }
#create folder for extracted data like lists and dataframes
if (!dir.exists(paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/", sep=""))) {
  dir.create(paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/", sep=""))     }

#load all samples
if (file.exists(paste("R Stuff/ExtractedData_", ProjectName, "/saveEnvironment_", ProjectName, ".RData", sep = ""))) {
  load(paste("R Stuff/ExtractedData_", ProjectName, "/saveEnvironment_", ProjectName, ".RData", sep = ""))} else {


#create a list with input files, sort for filename
ListWithDeSeq2Files <- list.files(path = "R Stuff/ExtractedData_RNA-DeSeq2-pipeline_all runs/", pattern = ".tsv")

#import data, RNA seq from tsv. files
i <- 1
RNASeqData <- vector()
ListPatient <- vector()

for (i in 1:length(ListWithDeSeq2Files)) {

  subdat <- read.delim(file = paste("R Stuff/ExtractedData_RNA-DeSeq2-pipeline_all runs/", ListWithDeSeq2Files[i], sep = ""), sep = " ", header = TRUE)
  subdat$Patient <- strsplit(ListWithDeSeq2Files[i], "vs")
  subdat$Patient <- sapply(subdat$Patient, function(x) x[1])
  ListPatient <- c(ListPatient, subdat$Patient[1])  #this is a variable for later in correlation calculations
  RNASeqData <- rbind(RNASeqData, subdat)
}

#rearrange and rename
RNASeqData$ENSG <- rownames(RNASeqData)
RNASeqData <- RNASeqData[,c(8,1:7)]
row.names(RNASeqData) <- NULL
RNASeqData$ENSG <- strsplit(RNASeqData$ENSG, ".", fixed=TRUE)
RNASeqData$ENSG <- sapply(RNASeqData$ENSG, function(x) x[1])

#import data, genes for current Project (here: allChr)
geneNames <- read.table(file = paste("R stuff/mart_export_", ProjectName, ".txt", sep = ""), sep = "\t", header = TRUE)

#combine geneNames(contains all genes for project, here: all chr)  with RNA SeqData, exclude Data from other chr
RNASeqData_annot <- merge(geneNames, RNASeqData, by.x = "Gene.stable.ID", by.y = "ENSG", all = TRUE)

#drop list: a) RNA seq hits not part of loaded ENSG library (found in all sampleTripl.), b) library hits not found in RNA seq hits
RNASeqData_annot_drop <- RNASeqData_annot[is.na(RNASeqData_annot$Gene.name) | is.na(RNASeqData_annot$log2FoldChange),]  
RNASeqData_annot <- RNASeqData_annot[!is.na(RNASeqData_annot$Gene.name) & !is.na(RNASeqData_annot$log2FoldChange),]

#convert log2 to fold change
RNASeqData_annot$FoldChange<- 2^RNASeqData_annot$log2FoldChange

#some values in pVal adjusted and pValue are NAs. to overcome this problem we generate an artifical p val of 1
RNASeqData_annot$pvalue[is.na(RNASeqData_annot$pvalue) & RNASeqData_annot$baseMean != 0] <- 1
RNASeqData_annot$padj[is.na(RNASeqData_annot$padj) & RNASeqData_annot$baseMean != 0] <- 1

#add pValue threshold for colouring, creating two columns, for both pvalues(adjusted and unadjusted)
RNASeqData_annot$PValAdjusted <- "non-significant"
RNASeqData_annot$PValAdjusted[RNASeqData_annot$padj <= thresholdForPval] <- "significant"
RNASeqData_annot$PValNONAdjusted <- "non-significant"
RNASeqData_annot$PValNONAdjusted[RNASeqData_annot$pvalue <= thresholdForPval] <- "significant"

#add low Expression indicator, abs() make the "Betrag" as boundaries
RNASeqData_annot$DeltaLowExpr <- NA
RNASeqData_annot$DeltaLowExpr[abs(RNASeqData_annot$log2FoldChange) >= log(thresholdForExpr, 2)] <- FALSE
RNASeqData_annot$DeltaLowExpr[abs(RNASeqData_annot$log2FoldChange) <  log(thresholdForExpr, 2)] <- TRUE

#Identifier for Expression Pattern
RNASeqData_annot$ExpressionPatternVSCtrl <- NA
RNASeqData_annot$ExpressionPatternVSCtrl[RNASeqData_annot$DeltaLowExpr == TRUE] <- "unchanged"
RNASeqData_annot$ExpressionPatternVSCtrl[RNASeqData_annot$log2FoldChange <=  -log(thresholdForExpr, 2)] <- "downregulated"
RNASeqData_annot$ExpressionPatternVSCtrl[RNASeqData_annot$log2FoldChange >=  log(thresholdForExpr, 2)] <- "upregulated"

#add conversion of bases into kilo bases
conversionFactorKB <- 1000
RNASeqData_annot$CoordStart_KB <- RNASeqData_annot$Gene.start..bp. / conversionFactorKB
conversionFactorMB <- 1000000
RNASeqData_annot$CoordStart_MB <- RNASeqData_annot$Gene.start..bp. / conversionFactorMB

#add protein gene type identifier for plot
RNASeqData_annot$logicalProteinType <- FALSE
RNASeqData_annot$logicalProteinType[RNASeqData_annot$Gene.type == "protein_coding"] <- TRUE

#clean environment and save progress
RNASeqData_clean <- RNASeqData_annot
remove(geneNames, RNASeqData, RNASeqData_annot,subdat, ListWithDeSeq2Files, i)
save.image(file = paste("R Stuff/ExtractedData_", ProjectName, "/saveEnvironment_", ProjectName, ".RData", sep = ""))

}

#backup from clean data
RNASeqData_clean_B <- RNASeqData_clean
#RNASeqData_clean <-RNASeqData_clean_B


#get gene identifier for filterung
firstDeletion_1 <- c(43414908,44797076)
firstDeletion_2 <- c(44797226,44866986)
firstDeletion <- c(firstDeletion_1[1], firstDeletion_2[2])
secondDeletion <- c(45781450,48110000)
surroundingsDeletion <- c(42414908,47110000)

#add deletion section
RNASeqData_clean$Position <- "Unaffected"
RNASeqData_clean$Position[RNASeqData_clean$Chromosome.scaffold.name == 21 & RNASeqData_clean$Gene.start..bp. > surroundingsDeletion[1] & RNASeqData_clean$Gene.end..bp. < surroundingsDeletion[2]] <- "surroundings_1MB"
RNASeqData_clean$Position[RNASeqData_clean$Chromosome.scaffold.name == 21 & RNASeqData_clean$Gene.start..bp. > firstDeletion[1] & RNASeqData_clean$Gene.end..bp. < firstDeletion[2]]  <- "First/Second Deletion"
RNASeqData_clean$Position[RNASeqData_clean$Chromosome.scaffold.name == 21 & RNASeqData_clean$Gene.start..bp. > secondDeletion[1] & RNASeqData_clean$Gene.end..bp. < secondDeletion[2]]  <- "Third Deletion"

#Colour Code for Gene in and outside of deletions
RNASeqData_clean$ColGenes <- "blue"
RNASeqData_clean$ColGenes[RNASeqData_clean$Position %in% c("First/Second Deletion", "Third Deletion")] <- "green"

#remove unsignificant entries for better overview
specificsFilter <- vector()
FilterNonSignificance <- FALSE
if(FilterNonSignificance == TRUE){
  RNASeqData_clean <- RNASeqData_clean[RNASeqData_clean$PValAdjusted == "significant",]
  specificsFilter <- c(specificsFilter, "onlySignif")
}

#focus on protein-coding
FilterProteinGenes <- TRUE
if(FilterProteinGenes == TRUE){
  RNASeqData_clean <- RNASeqData_clean[RNASeqData_clean$Gene.type == "protein_coding",]
  specificsFilter <- c(specificsFilter, "ProteinCoding")
}

#focus on chr21
FilterChr21 <- TRUE
if(FilterChr21 == TRUE){
  RNASeqData_clean <- RNASeqData_clean[RNASeqData_clean$Chromosome.scaffold.name == 21,]
  specificsFilter <- c(specificsFilter, "Chr21")
}

#focus on ALI9 and DSB
FilterALI9 <- TRUE
if(FilterALI9 == TRUE){
RNASeqData_clean <- RNASeqData_clean[RNASeqData_clean$Patient %in% c("ALI09", "DSB2", "DSB3", "DSBmerge"),]
specificsFilter <- c(specificsFilter, "focusDSB-ALI9")
}

#focus on Deletions
FilterDeletion <- FALSE
if(FilterDeletion == TRUE){
RNASeqData_clean <- RNASeqData_clean[RNASeqData_clean$Gene.name %in% c("First/Second Deletion", "Third Deletion"),]
specificsFilter <- c(specificsFilter, "only Dels")
}

specificsFiltercollapse = paste(specificsFilter, collapse = "_")

#----matrix building----
PlotName <- "ClusterAll"
Clean_ClusterAll <- RNASeqData_clean

#exclude all genes with no data, baseMean=0, only for plot to avoid unnecessary messages
Clean_ClusterAll <- Clean_ClusterAll[Clean_ClusterAll$baseMean >thresholdForBaseMean,]  
 
#ordering for correct data implementation
Clean_ClusterAll <- Clean_ClusterAll[order(Clean_ClusterAll$Patient),]
Clean_ClusterAll <- Clean_ClusterAll[order(Clean_ClusterAll$Gene.name),]
  
#filter against genes which are not present in all samples
Clean_ClusterAll_sub <- Clean_ClusterAll %>% group_by(Gene.name) %>% tally()
GenesPresentInAllSamples <- Clean_ClusterAll_sub$Gene.name[Clean_ClusterAll_sub$n == length(levels(factor(Clean_ClusterAll$Patient)))]
Clean_ClusterAll <- Clean_ClusterAll[Clean_ClusterAll$Gene.name %in% GenesPresentInAllSamples,]

#get order for gene colouring
ColDelGenes <- Clean_ClusterAll$ColGenes[!duplicated(Clean_ClusterAll$Gene.name)]

#create matrix for all genes
CCA_M<- matrix(data=Clean_ClusterAll$log2FoldChange, 
                   nrow = nrow(Clean_ClusterAll) / length(levels(factor(Clean_ClusterAll$Patient))),
                   ncol = length(levels(factor(Clean_ClusterAll$Patient))),
                   byrow=TRUE)
colnames(CCA_M) <- levels(factor(Clean_ClusterAll$Patient))
rownames(CCA_M) <- Clean_ClusterAll$Gene.name[seq(1,nrow(Clean_ClusterAll), by = length(levels(factor(Clean_ClusterAll$Patient))))]
#CCA_M <- log(CCA_M, 10)
#CCA_M[is.infinite(CCA_M)] <- NA
#CCA_M <- CCA_M[apply(CCA_M, MARGIN = 1, FUN = function(x) length(x[!is.na(x)]) > (ncol(CCA_M)/2)),]  #removes rows with more than half entries of NAs rows


#----Plotting----

#VectorResolution <- c(2500,25000)
#VectorResolution <- c(2500,5000)
VectorResolution <- c(5000)
m <- 1
for (m in 1:length(VectorResolution)) {

  png(paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_", VectorResolution[m], specificsFiltercollapse, "_NOscaleRow", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""),    # create PNG for the heat map        
      width = 2500, height = VectorResolution[m], res = 300, pointsize = 12)
  

  test2 <- heatmap.2(CCA_M, srtCol = 90, offsetRow = -0.3, offsetCol = -0.3, scale = "none", na.rm=TRUE,
          col=bluered, na.color = "grey",
          #trace = "none", 
          tracecol="black", 
          RowSideColors = ColDelGenes,
          density.info=c("density"), densadj = 0.5, key.title = "RNA Seq all Genes", keysize = 1)
  test2
  dev.off()
  
  png(paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_", VectorResolution[m], "_", specificsFiltercollapse, "_WITHscaleRow", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""),    # create PNG for the heat map
      width = 2500, height = VectorResolution[m], res = 300, pointsize = 12)


  test2 <- heatmap.2(CCA_M, srtCol = 90, offsetRow = -0.3, offsetCol = -0.3, scale = "row", na.rm=TRUE,
                     col=bluered, na.color = "grey",
                     #trace = "none", 
                     tracecol="black",
                     RowSideColors = ColDelGenes,
                     density.info=c("density"), densadj = 0.5, key.title = "RNA Seq all Genes", keysize = 1)
  test2
  dev.off()
  
}

#-----------------testing----------
## Some input sample matrix

y <- matrix(rnorm(50), 10, 5, dimnames=list(paste("g", 1:10, sep=""),
                                            paste("t", 1:5, sep="")))
CCA_M

## Run heatmap.2 on this matrix
test <- heatmap.2(y)
test3 <- heatmap.2(CCA_M)

y[rev(test$rowInd), test$colInd]
CCA_M[rev(test3$rowInd), test3$colInd]
CCA_M_Cluster <- CCA_M[rev(test2$rowInd), test2$colInd]

#------getting cluster tree----
#example
hm <- heatmap.2( y )
hc <- as.hclust( hm$rowDendrogram )
cutree( hc, h=5 )

for (i in seq(nrow(hc),0,-1)) {
  print(i)
  print(cutree( hc, h=i ))
  print("")
}
#-for own Matrix
test2
test2_hm <- as.hclust( test2$rowDendrogram )

#for (i in seq(nrow(test2_hm),0,-1)) {
for (i in seq(nrow(CCA_M),0,-1000)) {
  print(i)
  print(cutree( test2_hm, h=i ))
  print("")
}

max(cutree(test2_hm, h=10 ))























https://support.bioconductor.org/p/48624/
