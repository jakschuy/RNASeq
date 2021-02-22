Sys.setenv(LANG="en")
rm(list=ls())
options(scipen=999)

library(ggplot2)
library(ggrepel)
#library(xlsx)
library(writexl)
library(readxl)
library(dplyr)
library(cowplot)

#input
ProjectName <- "all-runs_HG38"
RefOfSeqDat <- "CTRL"
#NameInPutData <- "Sequencing data/ALI9_vs_ctrl.results.tsv"
firstDeletion <- c(1,8204967)
secondDeletion <- c(41994799, 43447106)
thirdDeletion <- c(44361567, 46690088)
Inversion <- c(43447107, 44361566)
surroundingsDeletion <- c(40994799, 47690088)
surroundingsDeletionP <- c(-1000000,9204967)
Region_PAR1 <- c(10000, 2781479)
Region_PAR2 <- c(56887902, 57217415)
colorsForLegend <- c("grey", "black", "dodgerblue", "black", "red", "green", "palegreen1", "palegreen4", "darkgoldenrod1", "darkorchid1")
thresholdForPval <- 0.05
thresholdForExpr <- 1.25
thresholdForBaseMean <- 25
useBothPVal_Adj_NonAdj <- FALSE
ListRiskGenesDS <- c("APP","ADAMTS1","COL6A1","DSCAM","DSCR1","DYRK1A","ERG","ETS2","ITSN1","JAM2","OLIG1","OLIG2","PTTG1IP","SIM2","SLC19A1","SYNJ1")
ListNESCells <- c("PAX6","NES")
ListNeuronalGenes <- c("LMX1B", "NGFR",  "MAP2",  "NTRK1", "NTRK2", "NTRK3",  "TUBB3")
ListGliaCells <- c("GFAP", "VIM", "S100B", "TAGLN2", "FABP7")
ListPARGenes <- c("AKAP17A", "ASMT", "ASMTL", "ASMTL-AS1", "CD99", "CD99P1", "CRLF2", "CSF2RA", "DHRSX", "DHRSX-IT1", "FABP5P13", "GTPBP6", "IL3RA", "LINC00102", "LINC00106", "LINC00685", "MIR3690", "MIR6089", "PLCXD1", "PPP2R3B", "P2RY8", "SHOX", "SLC25A6", "XG", "ZBED1", "AMD1P2", "DDX11L16", "DPH3P2", "IL9R", "SPRY3", "ELOCP24",  "TRPC6P", "VAMP7", "WASH6P", "WASIR1")
filterForRiskGenesDS <- FALSE
filterForRiskGenesDS_inPlot <- TRUE
SkipExportSection <- TRUE
ChosenPatient <- "ALI9"
RemoveUnSigForPlots <- TRUE
ListGOtermsNeuron <- c("GO:0022008", "GO:0021872", "GO:0048699",  "GO:0007407", "GO:0014016", "GO:0007405", "GO:0030182", "GO:0001764", "GO:0050767") 
ListGOtermsGlia <- c("GO:0042063", "GO:0010001", "GO:0042065", "GO:0008347", "GO:0014009", "GO:0070444")
ListGOterms_all <- c(ListGOtermsNeuron, ListGOtermsGlia)

#required functions
tophits <- function(x,n=5,decreasing = TRUE, na.last = TRUE){
  # Function tophits: goes through a numeric vector(!), sorts it and returns the tio entries with a list length of n,
  # can deal with NAs (na.last=F/T)
  # can deal with reverse order (decreasing = T/F)
  
  if(!is.numeric(x)) {
    stop('\nI am so sorry, but this function only works for numeric input!\n',
         'You have provided an object of class: ', class(x)[1])
  }
  
  listSort <- sort(x, decreasing = decreasing, na.last = na.last)
  result <- head(listSort, n=n)
  
  return(result)
}

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
if (file.exists(paste(getwd(), "/R Stuff/ExtractedData_", ProjectName, "/saveEnvironment.RData", sep=""))) {
  load(paste(getwd(), "/R Stuff/ExtractedData_", ProjectName, "/saveEnvironment.RData", sep=""))} else {

#create a list with input files, sort for .tsv
ListWithSamples <- list.files(path = paste("Sequencing data/", "all-runs", sep = ""), pattern = ".tsv")

#import data, genes for current Project (here: allChr)
geneNames <- read.delim(file = "Sequencing data/Ensemble-hg38_noFilter.txt", sep = "\t", header = TRUE)
geneNamesCellTypes_manual <- read.delim(file = "Sequencing data/marker_neuron_glia_hg38.txt", sep = "\t", header = TRUE)
ListMarkersCellTypes_manual <- c(ListNESCells, ListNeuronalGenes, ListGliaCells)
geneNamesCellTypes_GOterms <- read.delim(file = "Sequencing data/mart_export_GO-0022008 neurogenesis_hg38.txt", sep = "\t", header = TRUE)


#remove duplicates in geneNames
geneNamesWithGO_NCBI <- geneNames  #backup with GO terms
geneNames <- geneNames[,-which(names(geneNames) %in% c("NCBI.gene..formerly.Entrezgene..accession","GO.term.name"))]  #remove GO terms
geneNames <- geneNames[!duplicated(geneNames),]
geneNames2 <- unique(geneNames$Gene.stable.ID.version[duplicated(geneNames$Gene.stable.ID.version)])  #list of only duplicates
geneNames3 <- vector()
geneNames_new <- geneNames[!geneNames$Gene.stable.ID.version %in% geneNames2,]  #take not duplicated

i <- 1
for (i in 1:length(geneNames2)) {
  subdat <- geneNames[as.character(geneNames$Gene.stable.ID.version) == geneNames2[i],]   #isolate duplicate
  
  if (any(levels(factor(subdat$Gene.name)) == levels(factor(subdat$HGNC.symbol))) == TRUE) {    #choose more official name
    keep <- which(levels(factor(subdat$Gene.name)) == levels(factor(subdat$HGNC.symbol)))
    subdat <- subdat[keep,]
    if (nrow(subdat) > 1) {
      warning(message = paste("check duplicates with", geneNames2[i], sep = " "))
    }   #quality control
    
  }else{
    subdat <- subdat[!duplicated(subdat),] #delete duplicated rows
  }
  if (nrow(subdat) > 1) {
    warning(message = paste("check manually with", geneNames2[i], sep = " "))
  } #quality control
  
  geneNames3 <- rbind(geneNames3, subdat)
  
}
geneNames_new <- rbind(geneNames_new, geneNames3)
geneNames <- geneNames_new[order(geneNames_new$Gene.stable.ID.version),]
rownames(geneNames) <- NULL
geneNames$Gene.stable.ID <- strsplit(as.character(geneNames$Gene.stable.ID.version), ".", fixed=TRUE)
geneNames$Gene.stable.ID <- sapply(geneNames$Gene.stable.ID, function(x) x[1])
geneNames <- geneNames[, -which(names(geneNames) == "Gene.stable.ID.version")]


#import data, RNA seq from tsv. files
i <- 1
RNASeqData <- vector()
ListPatient <- vector()
for(i in 1:length(ListWithSamples)){
  
  subdat <- read.delim(file = paste("Sequencing data/", "all-runs", "/", ListWithSamples[i], sep=""), sep = " ", header = TRUE)
  subdat$Patient <- strsplit(ListWithSamples[i], "vs")
  subdat$Patient <- sapply(subdat$Patient, function(x) x[1])
  ListPatient <- c(ListPatient, subdat$Patient[1])  #this is a variable for later in correlation calculations
  RNASeqData <- rbind(RNASeqData, subdat)
}

#rearrange and rename
RNASeqData$ENSG <- rownames(RNASeqData)
RNASeqData <- RNASeqData[,c(8,1:7)]
rownames(RNASeqData) <- NULL
RNASeqData$ENSG <- strsplit(RNASeqData$ENSG, ".", fixed=TRUE)
RNASeqData$ENSG <- sapply(RNASeqData$ENSG, function(x) x[1])

#combine geneNames(contains all genes for project, here: all chr)  with RNA SeqData, exclude Data from other chr
Chr_all <- merge(RNASeqData, geneNames, by.x = "ENSG", by.y = "Gene.stable.ID", all.x = TRUE, all.y = FALSE)

#convert log2 to fold change
Chr_all$FoldChange<- 2^Chr_all$log2FoldChange

#some values in pVal adjusted and pValue are NAs. to overcome this problem we generate an artifical p val of 1
Chr_all$pvalue[is.na(Chr_all$pvalue) & Chr_all$baseMean != 0] <- 1
Chr_all$padj[is.na(Chr_all$padj) & Chr_all$baseMean != 0] <- 1

#add pValue threshold for colouring, creating two columns, for both pvalues(adjusted and unadjusted)
Chr_all$PValAdjusted <- "non-significant"
Chr_all$PValAdjusted[Chr_all$padj <= thresholdForPval] <- "significant"
Chr_all$PValNONAdjusted <- "non-significant"
Chr_all$PValNONAdjusted[Chr_all$pvalue <= thresholdForPval] <- "significant"

#add low Expression indicator, abs() make the "Betrag" as boundaries
Chr_all$DeltaLowExpr <- TRUE
Chr_all$DeltaLowExpr[Chr_all$log2FoldChange >=  log(thresholdForExpr, 2)] <- FALSE
Chr_all$DeltaLowExpr[Chr_all$log2FoldChange <=  log(1-(thresholdForExpr-1), 2)] <- FALSE

#Identifier for Expression Pattern
Chr_all$ExpressionPatternVSCtrl <- NA
Chr_all$ExpressionPatternVSCtrl[Chr_all$DeltaLowExpr == TRUE] <- "unchanged"
Chr_all$ExpressionPatternVSCtrl[Chr_all$log2FoldChange <=  log(1-(thresholdForExpr-1), 2)] <- "downregulated"
Chr_all$ExpressionPatternVSCtrl[Chr_all$log2FoldChange >=  log(thresholdForExpr, 2)] <- "upregulated"


#add deletion section
Chr_all$Position <- "Unaffected"
Chr_all$Position[Chr_all$Gene.start..bp. > Inversion[1] & Chr_all$Gene.start..bp. < Inversion[2] & Chr_all$Chromosome.scaffold.name == "21"]  <- "Inversion"
Chr_all$Position[Chr_all$Gene.start..bp. > surroundingsDeletion[1] & Chr_all$Gene.start..bp. < surroundingsDeletion[2] & Chr_all$Chromosome.scaffold.name  == "21"] <- "surroundings_1MB"
Chr_all$Position[Chr_all$Gene.start..bp. > surroundingsDeletionP[1] & Chr_all$Gene.start..bp. < surroundingsDeletionP[2] & Chr_all$Chromosome.scaffold.name  == "21"] <- "surroundings_1MB"
Chr_all$Position[Chr_all$Gene.start..bp. > firstDeletion[1] & Chr_all$Gene.start..bp. < firstDeletion[2] & Chr_all$Chromosome.scaffold.name  == "21"]  <- "First Deletion"
Chr_all$Position[Chr_all$Gene.start..bp. > secondDeletion[1] & Chr_all$Gene.start..bp. < secondDeletion[2] & Chr_all$Chromosome.scaffold.name  == "21"]  <- "Second Deletion"
Chr_all$Position[Chr_all$Gene.start..bp. > thirdDeletion[1] & Chr_all$Gene.start..bp. < thirdDeletion[2] & Chr_all$Chromosome.scaffold.name  == "21"]  <- "Third Deletion"

#reorder Position
Chr_all$Position <- factor(Chr_all$Position, levels(factor(Chr_all$Position))[c(1,2,4,3,5)])

#add protein gene type identifier for plot
Chr_all$logicalProteinType <- FALSE
Chr_all$logicalProteinType[Chr_all$Gene.type == "protein_coding"] <- TRUE

#add conversion of bases into kilo bases
conversionFactorKB <- 1000
Chr_all$CoordStart_KB <- Chr_all$Gene.start..bp. / conversionFactorKB
conversionFactorMB <- 1000000
Chr_all$CoordStart_MB <- Chr_all$Gene.start..bp. / conversionFactorMB

##


#stuff for plot
#section for pValNON adjusted input
vectpValue <- c("PValAdjusted", "PValNONAdjusted")
vectpValueNumber <- c("padj", "pvalue")
if(useBothPVal_Adj_NonAdj == TRUE)+
{usevectpValue <- c(1:2)}else+
{usevectpValue <- c(1)}

#reorder Chromosome factors
if(length(levels(factor(Chr_all$Chromosome.scaffold.name))) != length(c(1,12, 16:22, 2:11, 13:15, 24:25, 23))){
  warning("Levels in CHR wrong")
  }else{
  message("Levels in CHR okay")
  }
Chr_all$Chromosome.scaffold.name <- factor(Chr_all$Chromosome.scaffold.name, levels(factor(Chr_all$Chromosome.scaffold.name))[c(1,12, 16:22, 2:11, 13:15, 24:25, 23)])

#build DSgenes List
temp <- Chr_all
temp <- temp[temp$Patient %in% c("ALI09", "DSBmerge"),]
temp <- temp[temp$Gene.name %in% c(
  "JAM2", "APP", "ADAMTS1", "SYNJ1", "OLIG2","OLIG1","ITSN1", "MRPS6", "RCAN1", 
  "RUNX1","CBR1", "SIM2","DYRK1A", "ERG","ETS2","HMGN1","DSCAM","PFKL", "PTTG1IP",
  "SLC19A1", "COL6A1","S100B","CRLF2"),]
write_xlsx(temp, path = "R stuff/ExtractedData_all-runs_HG38/DSgenesForALI09-DSBmerge.xlsx")

##for export, all samples = patients__ BEFORE BASEMEAN CUT
if(SkipExportSection == FALSE){
  ListExport <- Chr_all
  #export the gene list with all computed values
  write_xlsx(ListExport, paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/", "ListExport_all_AllSamples-vs", RefOfSeqDat, "NO_BaseMEAN CUTOFF.xlsx", sep=""))
}

#exclude all genes with no data, baseMean=0, only for plot to avoid unnecessary messages
Chr_all <- Chr_all[Chr_all$baseMean >= 5,]

##for export, all samples = patients
if(SkipExportSection == FALSE){
  ListExport <- Chr_all
  #export the gene list with all computed values
  write_xlsx(ListExport, paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/", "ListExport_all_AllSamples-vs", RefOfSeqDat, ".xlsx", sep=""))
}

#Backup
Chr_all.red <- Chr_all

#remove unsignificant entries for better overview
Chr_all.redSig <- Chr_all.red
if(RemoveUnSigForPlots ==TRUE){
  Chr_all.redSig <- Chr_all.red[Chr_all.red$PValAdjusted == "significant",]
}



#Sobol, Chr21_all, not reduced, is not used so far
#NPC: "OLIG1", "OLIG2", "APP", "CBR3", "NCAM2" 
#DiffNPCs/Chr21/trisomy21: "ATP5O", "PDXK", "C21orf91", "C21orf33", "C2CD2", "CYYR1"
#both stages, DEGs: "PTTG1IP", "COL8A1", "RUNX1"
df_From_Sobol2019 <- Chr_all[Chr_all$Patient != "ALI9" & Chr_all$Gene.name %in% c("OLIG1", "OLIG2", "APP", "CBR3", "NCAM2", "ATP5O", "PDXK", "C21orf91", "C21orf33", "C2CD2", "CYYR1", "PTTG1IP", "COL8A1", "RUNX1"),]


#save environment for processed data
rm(geneNames2, geneNames3, geneNames_new, subdat,temp,RNASeqData,i, keep, ListExport)
save.image(file = paste(getwd(), "/R Stuff/ExtractedData_", ProjectName, "/saveEnvironment.RData", sep=""))

} #close load command



#----all Chromosomes----

#plot and loop
i <- 1
# for(i in 1:length(usevectpValue)){
for(i in 1:length(levels(factor(Chr_all.redSig$Patient)))){
  #exclude all Patients but chosen one
  Chr_all.redSig_sub <- Chr_all.redSig[Chr_all.redSig$Patient == levels(factor(Chr_all.redSig$Patient))[i],]
  
  #Plot
  plot.data <- ggplot(Chr_all.redSig_sub[Chr_all.redSig_sub$Chromosome.scaffold.name %in% c(1:22, "X"),], aes(x=CoordStart_MB, y=FoldChange))+
    geom_hline(yintercept = 2^log(thresholdForExpr, 2), alpha = 0.5,size=0.2)+
    geom_hline(yintercept = 2^log(1-(thresholdForExpr-1), 2), alpha = 0.5,size=0.2)+
    geom_hline(linetype=2, yintercept = 1, alpha = 0.2,size=0.2)+
    facet_wrap(~Chromosome.scaffold.name, scales = "free_x", ncol = 6)+
    geom_point(aes(shape=Patient, color = logicalProteinType))+
    geom_smooth(method = "loess", show.legend = FALSE, col = colorsForLegend[3], se = FALSE, span = 0.8)+
    scale_color_manual(values=colorsForLegend[c(2,5)], labels = c("other", "protein coding"))+
    theme_classic(base_size=20)+
   scale_y_continuous(trans='log10')+
    labs(x = "Chromosome location [mb]", y = "Fold Change", color = "Gene type") +
    #labs(x = "Chromosomes, y = "Fold Change") +
    guides(
      fill=FALSE, 
      color=guide_legend(order=1), 
      shape=FALSE 
      #alpha=guide_legend(order=3)
    )+
    theme(
      element_line(size = 1.0),
      plot.subtitle = element_text(colour="black",face="italic", size = 15),
      #axis.line = element_line(size = 2),
      axis.ticks = element_line(colour="black"),
      #axis.ticks = element_line(size = 1, colour="black"),
      #axis.text.x = element_text(colour="black",face="bold", hjust=0.5, size = 25),
      axis.text.x = element_text(colour="black", hjust=0.5),
      #axis.text.y = element_text(colour="black",face="bold", size = 25),
      axis.text.y = element_text(colour="black", hjust=0.5),
      #legend.title = element_text(colour="black", size = 25),
      #strip.background = element_rect(size=2, linetype="blank", colour ="red"),
      legend.position = c(0.9, 0.1)
      #legend.title = element_blank(),
    )+
    #xlim(c(-5,15))+
    #ylim(c(-5,15))+
    # geom_text_repel(label = ifelse(Chr_all.redSig$FoldChange %in% c(tophits(Chr_all.redSig$FoldChange[Chr_all.redSig[,vectpValue[i]] == "significant"]),
    #                                                            tophits(Chr_all.redSig$FoldChange[Chr_all.redSig[,vectpValue[i]] == "significant"], decreasing = FALSE)), 
    #                               paste(Chr_all.redSig$Gene.name), 
    #                                ""),
    # # geom_text_repel(label = ifelse(Chr_all.redSig[,vectpValue[i]] == "significant" & Chr_all.redSig$Position != "Unaffected", paste(Chr_all.redSig$Gene.name), ""),
    #                 col = "black", size = 3.5, #for label settings
    #                 #ylim = c(1.8, NA), #for y-axis movements
    #                 #xlim=c(NA, 57000000),
    #                 xlim=c(NA, 48000000),
  #                 #force=3, #no change
  #                 #vjust = -1, #no change
  #                 hjust = 0.3, #no change
  #                 segment.alpha = 0.2, #line transparency
  #                 # nudge_y = 2.1 - Chr_all.redSig$FoldChange[!is.na(Chr_all.redSig$FoldChange)], #no need for this
  #                 nudge_y = 0.1, 
  #                 # box.padding = 4,
  #                 #direction = "x", 
  #                 angle = 0)+ # to which direction should the text be aligend to,turn text
  
  ggtitle("RNA Seq Data, NES, Chr_all", subtitle = paste(levels(factor(Chr_all.redSig$Patient))[i], "vs", RefOfSeqDat, ", threshold low deltaExpr: |",thresholdForExpr, "|x; threshold pVal: ", thresholdForPval, "; ", "p-valAdjusted", sep =""))
  #plot.data
  
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", levels(factor(Chr_all.redSig$Patient))[i], "vs", RefOfSeqDat, "_", "p-valAdjusted", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
  #ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/temp.png", sep=""),
  #      height = 18, width = 30, unit="cm")
}

#plot.data

ggsave(filename = paste("R stuff/test_all.png", sep=""),
       height = 100, width = 100, unit="cm",limitsize = FALSE)





###----check for global effects----
PlotName <- "GlobalEffect"

Chr_globalEf <- Chr_all

#change order for aes in plot
Chr_globalEf$PValAdjusted <- relevel(as.factor(Chr_globalEf$PValAdjusted), ref = "significant")
Chr_globalEf$logicalProteinType <- relevel(as.factor(Chr_globalEf$logicalProteinType), ref = "TRUE")

# Chr_globalEf$Chromosome.scaffold.name <- droplevels(Chr_globalEf$Chromosome.scaffold.name)
# Chr_globalEf$Chromosome.scaffold.name <- factor(Chr_globalEf$Chromosome.scaffold.name, levels(factor(Chr_globalEf$Chromosome.scaffold.name))[c(1,12, 16:22, 2:11, 13:15, 82:83, 81, 23:80)])

#remove unsignificant entries for better overview
if(RemoveUnSigForPlots ==TRUE){
  Chr_globalEf <- Chr_globalEf[Chr_globalEf$PValAdjusted == "significant",]
}

#filter against unchanged expression
Chr_globalEf <- Chr_globalEf[Chr_globalEf$ExpressionPatternVSCtrl != "unchanged",]

#export data
ListExport_globalEf <- Chr_globalEf

if(SkipExportSection == FALSE) {
    write_xlsx(ListExport_globalEf, paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/", "ListExport_", PlotName, ".xlsx", sep=""))
}


#prepare all vs all loop
j <- 1
MatPatients <- combn(ListPatient[!ListPatient %in% c("CTRL11", "CTRL7", "CTRL9", "CTRLmerge", "DSB2", "DSB3")], 2)
Chr_glEf_c_m_numbers <- as.data.frame(matrix(data=NA, nrow=3, ncol=ncol(MatPatients)))

SkipExportSection <- FALSE
SkipPlotSection <- FALSE

for (j in 1:ncol(MatPatients)){
  
  #define both variables for the whole function (could be skipped and the code be implemented in the corresponding lines but the for loop was added afterwards and doing so it saves time and reduce complexity)
  m <- MatPatients[1,j]
  n <- MatPatients[2,j]

##filter for common genes: DSmerge and loop[i]
#create df with selected two patients
Chr_globalEf_common <- Chr_globalEf[Chr_globalEf$Patient %in% c(m,n),]

keep <- vector()  #robust against double entries with two dif annotations for same patient!
keep <- c(keep, unique(as.character(Chr_globalEf_common$ENSG[Chr_globalEf_common$Patient == m])))
keep <- c(keep, unique(as.character(Chr_globalEf_common$ENSG[Chr_globalEf_common$Patient == n])))
DifExpGenesInclNotShared <- unique(keep) #for ratio alculations: mirrored genes/detected genes
keep <- keep[duplicated(keep)]

#Chr_globalEf_common <- Chr_globalEf_common[Chr_globalEf_common$ENSG %in% keep,] 
Chr_globalEf_common$SharedG <- FALSE
Chr_globalEf_common$SharedG[Chr_globalEf_common$ENSG %in% keep] <- TRUE

#PlotSection common genes
if (SkipPlotSection == FALSE) {
  

#Plot
plot.data <- ggplot(Chr_globalEf_common[Chr_globalEf_common$Chromosome.scaffold.name %in% c(1:22, "X"),], aes(x=CoordStart_MB, y=FoldChange))+
  geom_hline(yintercept = 2^log(thresholdForExpr, 2), alpha = 0.5,size=0.2)+
  geom_hline(yintercept = 2^log(1-(thresholdForExpr-1), 2), alpha = 0.5,size=0.2)+
  geom_hline(linetype=2, yintercept = 1, alpha = 0.2,size=0.2)+
  facet_wrap(~Chromosome.scaffold.name, scales = "free_x", ncol = 6)+
  geom_point(aes(color=Patient,  shape = logicalProteinType, alpha = PValAdjusted))+
  # geom_smooth(method = "loess", show.legend = FALSE, col = colorsForLegend[3], se = FALSE, span = 0.8)+
  scale_color_manual(values=colorsForLegend[c(2,5)])+
  scale_alpha_manual(values = c(1,0.1))+
  scale_shape_manual(values = c(16, 1))+
  theme_classic(base_size=20)+
  scale_y_continuous(trans='log10')+
  labs(x = "Chromosome location [Mb]", y = "Fold Change")+
  guides(
    fill=FALSE,
    color=guide_legend(order=1),
    shape=guide_legend(order=2),
    alpha=guide_legend(order=3)
   )+
  theme(
    element_line(size = 1.0),
    plot.subtitle = element_text(colour="black",face="italic", size = 15),
    axis.ticks = element_line(colour="black"),
    axis.text.x = element_text(colour="black", hjust=0.5),
    axis.text.y = element_text(colour="black", hjust=0.5),
    #legend.position = c(0.9, 0.1)
  )+

ggtitle("RNA Seq Data, NES, GlobalEffect", subtitle = paste(m,"vs", n, ", threshold low deltaExpr: |",thresholdForExpr, "|x; threshold pVal: ", thresholdForPval, "; ", "p-valAdjusted", sep =""))


# ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_", m,"vs", n, "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_", m,"vs", n, ".png", sep=""), 
       height = 30, width = 50, unit="cm", limitsize = FALSE)
}

#filter for DSmerge up/down + ALI09 down/up

  #going up, only shared genes
keep1 <- vector()
keep1 <- c(keep1, as.character(Chr_globalEf_common$ENSG[
                      Chr_globalEf_common$SharedG == TRUE &
                      (Chr_globalEf_common$Patient == m & Chr_globalEf_common$ExpressionPatternVSCtrl == "upregulated")
                                                                  ]
                               )
           )

keep1 <- c(keep1, as.character(Chr_globalEf_common$ENSG[
                      Chr_globalEf_common$SharedG == TRUE & 
                      (Chr_globalEf_common$Patient == n & Chr_globalEf_common$ExpressionPatternVSCtrl == "downregulated")
                                                                  ]
                               )
           )
keep1 <- keep1[duplicated(keep1)]

  #going down, only shared genes
keep2 <- vector()
keep2 <- c(keep2, as.character(Chr_globalEf_common$ENSG[
                      Chr_globalEf_common$SharedG == TRUE &
                      (Chr_globalEf_common$Patient == m & Chr_globalEf_common$ExpressionPatternVSCtrl == "downregulated")
                                                                  ]
                                )
           )

keep2 <- c(keep2, as.character(Chr_globalEf_common$ENSG[
                      Chr_globalEf_common$SharedG == TRUE &
                      (Chr_globalEf_common$Patient == n & Chr_globalEf_common$ExpressionPatternVSCtrl == "upregulated")
                                                                  ]
                                )
           )
keep2 <- keep2[duplicated(keep2)]

  #unite
keep <- c(keep1, keep2)

  #clear df for keep
Chr_globalEf_c_mirror <- Chr_globalEf_common[Chr_globalEf_common$ENSG %in% keep,]
Chr_globalEf_common$mirrored <- Chr_globalEf_common$ENSG %in% keep
Chr_globalEf_common$mirrored <- relevel(factor(Chr_globalEf_common$mirrored), ref = "TRUE")


#save df
if (SkipExportSection == FALSE) {
    
    # write_xlsx(Chr_globalEf_c_mirror, paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/", "ListExport_", PlotName, "_", m, "vs", n, "_MirrorEf.xlsx", sep=""))
    write_xlsx(Chr_globalEf_common, paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/", "ListExport_", PlotName, "_", m, "vs", n, ".xlsx", sep=""))

}

#PlotSection common genes+mirrored marked
if (SkipPlotSection == FALSE) {
  

  #Plot _ GlobalEffect_mmarked
  plot.data <- ggplot(Chr_globalEf_common[Chr_globalEf_common$Chromosome.scaffold.name %in% c(1:22, "X"),], aes(x=CoordStart_MB, y=FoldChange))+
    geom_hline(yintercept = 2^log(thresholdForExpr, 2), alpha = 0.5,size=0.2)+
    geom_hline(yintercept = 2^log(1-(thresholdForExpr-1), 2), alpha = 0.5,size=0.2)+
    geom_hline(linetype=2, yintercept = 1, alpha = 0.2,size=0.2)+
    facet_wrap(~Chromosome.scaffold.name, scales = "free_x", ncol = 6)+
    geom_point(aes(color=Patient,  shape = logicalProteinType, alpha = mirrored))+
    # geom_smooth(method = "loess", show.legend = FALSE, col = colorsForLegend[3], se = FALSE, span = 0.8)+
    scale_color_manual(values=colorsForLegend[c(2,5)])+
    scale_alpha_manual(values = c(1,0.1))+
    scale_shape_manual(values = c(16, 1))+
    theme_classic(base_size=20)+
    scale_y_continuous(trans='log10')+
    labs(x = "Chromosome location [mb]", y = "Fold Change")+
    guides(
      fill=FALSE,
      color=guide_legend(order=1),
      shape=guide_legend(order=2),
      alpha=guide_legend(order=3)
    )+
    theme(
      element_line(size = 1.0),
      plot.subtitle = element_text(colour="black",face="italic", size = 15),
      axis.ticks = element_line(colour="black"),
      axis.text.x = element_text(colour="black", hjust=0.5),
      axis.text.y = element_text(colour="black", hjust=0.5),
      #legend.position = c(0.9, 0.1)
    )+
    
    ggtitle("RNA Seq Data, NES, GlobalEffect_mmarked", subtitle = paste(m,"vs", n, ", threshold low deltaExpr: |",thresholdForExpr, "|x; threshold pVal: ", thresholdForPval, "; ", "p-valAdjusted; mirrored genes: ", sum(as.logical(Chr_globalEf_common$mirrored))/2, sep =""))
  
  
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_mmarked_", m,"vs", n, "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         height = 30, width = 50, unit="cm", limitsize = FALSE)
}#close plot


#PlotSection mirrored genes
if (SkipPlotSection == FALSE) {
  
  
#plot - only mirrored
#Plot
plot.data <- ggplot(Chr_globalEf_c_mirror[Chr_globalEf_c_mirror$Chromosome.scaffold.name %in% c(1:22, "X"),], aes(x=CoordStart_MB, y=FoldChange))+
  geom_hline(yintercept = 2^log(thresholdForExpr, 2), alpha = 0.5,size=0.2)+
  geom_hline(yintercept = 2^log(1-(thresholdForExpr-1), 2), alpha = 0.5,size=0.2)+
  geom_hline(linetype=2, yintercept = 1, alpha = 0.2,size=0.2)+
  facet_wrap(~Chromosome.scaffold.name, scales = "free_x", ncol = 6)+
  geom_point(aes(color=Patient,  shape = logicalProteinType, alpha = PValAdjusted))+
  # geom_smooth(method = "loess", show.legend = FALSE, col = colorsForLegend[3], se = FALSE, span = 0.8)+
  scale_color_manual(values=colorsForLegend[c(2,5)])+
  scale_alpha_manual(values = c(1,0.1))+
  scale_shape_manual(values = c(16, 1))+
  theme_classic(base_size=20)+
  scale_y_continuous(trans='log10')+
  labs(x = "Chromosome location [mb]", y = "Fold Change") +
  guides(
    fill=FALSE, 
    color=guide_legend(order=1), 
    shape=guide_legend(order=2),
    alpha=guide_legend(order=3)
  )+
  theme(
    element_line(size = 1.0),
    plot.subtitle = element_text(colour="black",face="italic", size = 15),
    axis.ticks = element_line(colour="black"),
    axis.text.x = element_text(colour="black", hjust=0.5),
    axis.text.y = element_text(colour="black", hjust=0.5),
    #legend.position = c(0.9, 0.1)
  )+
  
  ggtitle("RNA Seq Data, NES, GlobalEffect_MirrorExpression", subtitle = paste(m,"vs", n, RefOfSeqDat, ", threshold low deltaExpr: |",thresholdForExpr, "|x; threshold pVal: ", thresholdForPval, "; ", "p-valAdjusted; uniq.genes: ", nrow(Chr_globalEf_c_mirror)/2, sep =""))

ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_mirror", m,"vs", n, "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       height = 30, width = 50, unit="cm", limitsize = FALSE)
}#close plot

#add number of mirrored genes
colnames(Chr_glEf_c_m_numbers)[j] <- paste(m, "_", n, sep="")
Chr_glEf_c_m_numbers[1,j] <- nrow(Chr_globalEf_c_mirror)/2
Chr_glEf_c_m_numbers[2,j] <- length(DifExpGenesInclNotShared)
Chr_glEf_c_m_numbers[3,j] <- Chr_glEf_c_m_numbers[1,j] / Chr_glEf_c_m_numbers[2,j]

# test <- as.data.frame(t(Chr_glEf_c_m_numbers))
# test$comp <- rownames(test)
# test$ali9 <- FALSE
# test$ali9[grepl(test$comp, pattern = "ALI05")] <- TRUE
# ggplot(test, aes(x= comp, y=V1, col=ali9)) + geom_point()

} #close all_runs loop









#----Mirror effect_SVs----
#!! must be run after gloabl effect chapter
#this only generates a plot!

PlotName <- "Mirror_SV"

#prepare list with correct patients 
Chr_Mirror_SV <- Chr_globalEf_common

#Filter for SV
Chr_Mirror_SV <- Chr_Mirror_SV[Chr_Mirror_SV$Position != "Unaffected",]


#Filter for protein coding
Chr_Mirror_SV <- Chr_Mirror_SV[Chr_Mirror_SV$logicalProteinType ==TRUE,]

#add colour code
Chr_Mirror_SV$Position <- as.character(Chr_Mirror_SV$Position)
Chr_Mirror_SV$Position[Chr_Mirror_SV$mirrored == FALSE] <- "not mirrored"

#reorder levels
if(length(levels(factor(Chr_Mirror_SV$Position))) != length(c(2,4,3,1))){
  warning("Levels in CHR wrong")
}else{
  message("Levels in CHR okay")
}
Chr_Mirror_SV$Position <- factor(Chr_Mirror_SV$Position, levels(factor(Chr_Mirror_SV$Position))[c(2,4,3,1)])



# #Filter for selected patients
# SelectedPatientsMirrorSV <- c("ALI09", "DSBmerge")
# Chr_Mirror_SV <- Chr_Mirror_SV[Chr_Mirror_SV$Patient %in% SelectedPatientsMirrorSV,]


#Plot
Chr_Mirror_SVQ <- Chr_Mirror_SV[Chr_Mirror_SV$CoordStart_KB > 40000,]
# plot.data.Q <- ggplot(Chr_Mirror_SVQ, aes(x=CoordStart_KB, y=FoldChange))+
plot.data.Q <- ggplot(Chr_Mirror_SVQ, aes(x=CoordStart_KB, y=log2FoldChange))+
  # geom_hline(yintercept = 2^log(thresholdForExpr, 2), alpha = 0.5,size=0.2)+
  # geom_hline(yintercept = 2^log(1-(thresholdForExpr-1), 2), alpha = 0.5,size=0.2)+
  # geom_hline(yintercept = 1, linetype=2, alpha = 0.2,size=0.2)+
  geom_hline(yintercept = log(thresholdForExpr, 2), alpha = 0.5,size=0.2)+
  geom_hline(yintercept = log(1-(thresholdForExpr-1), 2), alpha = 0.5,size=0.2)+
  geom_hline(yintercept = 0, linetype=2, alpha = 0.2,size=0.2)+
  geom_point(aes(color= Position,  shape = Patient), size=6, alpha = 0.4)+
  scale_color_manual(labels = c("Second Deletion", "Third Deletion", "SV proximity","not mirrored"), values=colorsForLegend[c(5,6,3,1)])+
  scale_alpha_manual(values = c(1,0.1))+
  scale_shape_manual(values = c(16, 17))+
  theme_classic(base_size=20)+
  scale_x_continuous(limits = c(40900,47000), breaks = seq(42000, 46000, 1000))+
  # scale_y_continuous(trans='log2', limits= c(0.08, 8.5), breaks = c(0.1, 0.2, 0.3, 0.5, 0.7,1,1.5, 2,3,4))+
  scale_y_continuous(limits= c(-3, 3), breaks = seq(-2.5,2, 0.5))+
  labs(x = "Chromosome 21 [kb]", y = "Fold Change", col = "Colour code") +
  labs(x = "Chromosome 21 [kb]", y = "Fold Change [log2]", col = "Colour code") +
  guides(
    color=guide_legend(order=1),
    shape=guide_legend(order=2),
    fill=guide_legend(order=3),
    alpha=guide_legend(order=4)
  )+
  theme(
    element_line(size = 1.0),
    plot.subtitle = element_text(colour="black",face="italic", size = 15),
    axis.ticks = element_line(colour="black"),
    axis.text.x = element_text(colour="black", hjust=0.5),
    axis.text.y = element_text(colour="black", hjust=0.5),
  )+
  geom_text_repel(label = ifelse((Chr_Mirror_SVQ$FoldChange > 1) & (Chr_Mirror_SVQ$Position != "not mirrored"),
    as.character(Chr_Mirror_SVQ$Gene.name), ""),
    col = "black", size = 6, #for label settings
    xlim=c(40900, 47000),
    ylim=c(2.1,NA),
    vjust = 0, #no change
    hjust =0, #no change
    segment.alpha = 0.2, #line transparency
    box.padding = 0.1,
    direction = "x",
    max.overlaps = 100,
    angle = 35)+ # to which direction should the text be aligend to,turn text

  # annotate(geom = "rect", xmin = secondDeletion[1]/conversionFactorKB, xmax = secondDeletion[2]/conversionFactorKB, ymin = 0, ymax = 4, alpha=0.08, fill = colorsForLegend[5])+
  # annotate(geom = "rect", xmin = thirdDeletion[1]/conversionFactorKB, xmax = thirdDeletion[2]/conversionFactorKB, ymin = 0, ymax = 4, alpha=0.08, fill = colorsForLegend[6])+
  # annotate(geom = "rect", xmin = surroundingsDeletion[1]/conversionFactorKB, xmax = secondDeletion[1]/conversionFactorKB, ymin = 0, ymax = 4, alpha=0.08, fill = colorsForLegend[3])+
  # annotate(geom = "rect", xmin = secondDeletion[2]/conversionFactorKB, xmax = thirdDeletion[1]/conversionFactorKB, ymin = 0, ymax = 4, alpha=0.08, fill = colorsForLegend[3])+
  annotate(geom = "rect", xmin = secondDeletion[1]/conversionFactorKB, xmax = secondDeletion[2]/conversionFactorKB, ymin = -3, ymax = 2, alpha=0.08, fill = colorsForLegend[5])+
  annotate(geom = "rect", xmin = thirdDeletion[1]/conversionFactorKB, xmax = thirdDeletion[2]/conversionFactorKB, ymin = -3, ymax = 2, alpha=0.08, fill = colorsForLegend[6])+
  annotate(geom = "rect", xmin = surroundingsDeletion[1]/conversionFactorKB, xmax = secondDeletion[1]/conversionFactorKB, ymin = -3, ymax = 2, alpha=0.08, fill = colorsForLegend[3])+
  annotate(geom = "rect", xmin = secondDeletion[2]/conversionFactorKB, xmax = thirdDeletion[1]/conversionFactorKB, ymin = -3, ymax = 2, alpha=0.08, fill = colorsForLegend[3])+

  
  ggtitle("RNA Seq Data, NES, SV-Mirrored-marked, q-arm", subtitle = paste(m,"vs", n, ", threshold low deltaExpr: |",thresholdForExpr, "|x; threshold pVal: ", thresholdForPval, "; ", "p-valAdjusted; mirrored genes: ", sum(as.logical(Chr_Mirror_SV$mirrored[Chr_Mirror_SV$CoordStart_KB > 40000]))/2, sep =""))
plot.data.Q

ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "SV-Mirrored-marked_", m,"vs", n, "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       height = 18, width = 30, unit="cm", limitsize = FALSE, dpi = 600)

#P-arm only
Chr_Mirror_SVP <- Chr_Mirror_SV[Chr_Mirror_SV$CoordStart_KB < 20000,]
# plot.data.P <- ggplot(Chr_Mirror_SVP, aes(x=CoordStart_KB, y=FoldChange))+
plot.data.P <- ggplot(Chr_Mirror_SVP, aes(x=CoordStart_KB, y=log2FoldChange))+
  # geom_hline(yintercept = 2^log(thresholdForExpr, 2), alpha = 0.5,size=0.2)+
  # geom_hline(yintercept = 2^log(1-(thresholdForExpr-1), 2), alpha = 0.5,size=0.2)+
  # geom_hline(yintercept = 1, linetype=2, alpha = 0.2,size=0.2)+
  geom_hline(yintercept = log(thresholdForExpr, 2), alpha = 0.5,size=0.2)+
  geom_hline(yintercept = log(1-(thresholdForExpr-1), 2), alpha = 0.5,size=0.2)+
  geom_hline(yintercept = 0, linetype=2, alpha = 0.2,size=0.2)+
  geom_point(aes(color= Position,  shape = Patient), size=6, alpha = 0.4)+
  scale_color_manual(labels = c("First Deletion","not mirrored", "SV proximity"), values = colorsForLegend[c(9,1,3)])+
  scale_alpha_manual(values = c(1,0.1))+
  scale_shape_manual(values = c(16, 17))+
  theme_classic(base_size=20)+
  scale_x_continuous(limits = c(0,9300), breaks = seq(0, 10000, 2000), n.breaks = 10)+
  # scale_y_continuous(trans='log2', limits= c(0.08, 6), breaks = c(0.1, 0.2, 0.3, 0.5, 0.7,1,1.5, 2,3,4))+
  scale_y_continuous(limits = c(-11, 5), breaks = seq(-11,3,1))+
  
  # labs(x = "Chromosome 21 [kb]", y = "Fold Change", col = "Colour code") +
  labs(x = "Chromosome 21 [kb]", y = "Fold Change [log2]", col = "Colour code") +
  guides(
    color=guide_legend(order=1),
    shape=guide_legend(order=2),
    fill=guide_legend(order=3),
    alpha=guide_legend(order=4)
  )+
  theme(
    element_line(size = 1.0),
    plot.subtitle = element_text(colour="black",face="italic", size = 15),
    axis.ticks = element_line(colour="black"),
    axis.text.x = element_text(colour="black", hjust=0.5),
    axis.text.y = element_text(colour="black", hjust=0.5),
  )+
  geom_text_repel(label = ifelse((Chr_Mirror_SVP$FoldChange > 1) & (Chr_Mirror_SVP$Position != "not mirrored"),
                                 as.character(Chr_Mirror_SVP$Gene.name), ""),
                  col = "black", size = 6, #for label settings
                  xlim=c(NA, 9000),
                  ylim=c(3.3,NA),
                  vjust = 0, #no change
                  hjust =0, #no change
                  segment.alpha = 0.2, #line transparency
                  box.padding = 0.1,
                  direction = "x",
                  angle = 35)+ # to which direction should the text be aligend to,turn text
  
  # annotate(geom = "rect", xmin = firstDeletion[1]/conversionFactorKB, xmax = firstDeletion[2]/conversionFactorKB, ymin = 0, ymax = 4, alpha=0.08, fill = colorsForLegend[9])+
  # annotate(geom = "rect", xmin = firstDeletion[2]/conversionFactorKB, xmax = surroundingsDeletionP[2]/conversionFactorKB, ymin = 0, ymax = 4, alpha=0.08, fill = colorsForLegend[3])+
  annotate(geom = "rect", xmin = firstDeletion[1]/conversionFactorKB, xmax = firstDeletion[2]/conversionFactorKB, ymin = -11, ymax = 3.3, alpha=0.08, fill = colorsForLegend[9])+
  annotate(geom = "rect", xmin = firstDeletion[2]/conversionFactorKB, xmax = surroundingsDeletionP[2]/conversionFactorKB, ymin = -11, ymax = 3.3, alpha=0.08, fill = colorsForLegend[3])+
  
  
  ggtitle("RNA Seq Data, NES, SV-Mirrored-marked, p-arm", subtitle = paste(m,"vs", n, ", threshold low deltaExpr: |",thresholdForExpr, "|x; threshold pVal: ", thresholdForPval, "; ", "p-valAdjusted; mirrored genes: ", sum(as.logical(Chr_Mirror_SV$mirrored[Chr_Mirror_SV$CoordStart_KB < 20000]))/2, sep =""))

plot.data.P

ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "SV-Mirrored-marked_", m,"vs", n, "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       height = 18, width = 30, unit="cm", limitsize = FALSE, dpi = 600)
 








#----Permuations for global effects----

PlotName <- "GlobalEffect_Permutations"

#Chosen Comparison/Patient for permutation
PermutationsPatient <- "ALI09"
PermutationsReference <- "DSBmerge"

#grep for selected mirrored genes comparison, important for statistics
ListSamplesGlobal <- list.files(path = paste("R stuff/ExtractedData_", ProjectName, "/", sep = ""), pattern = c("GlobalEffect"))
ListSamplesGlobal <- ListSamplesGlobal[grepl(x=ListSamplesGlobal, pattern = ".xlsx")]
ListSamplesGlobal <- ListSamplesGlobal[grepl(x=ListSamplesGlobal, pattern = paste0(PermutationsPatient,".*",PermutationsReference))]

#load mirrored genes data calculations, important for statistics
MirroredTableRef <-read_excel(path = paste("R stuff/ExtractedData_", ProjectName, "/", ListSamplesGlobal, sep=""))
MirroredTableRef <- MirroredTableRef[MirroredTableRef$logicalProteinType == TRUE,]


#filter for Chosen Patient
Chr_Permutations <- Chr_all[Chr_all$Patient %in% c(PermutationsPatient,PermutationsReference),]

#trim geneNames
geneNamesTrimmed <- geneNames[!duplicated(geneNames$Gene.stable.ID),c(8,5,7)]
colnames(geneNamesTrimmed)[1] <- "ENSG"
geneNamesTrimmed$ENSG <- as.character(geneNamesTrimmed$ENSG)


#backup
Chr_PermutationsNew <- Chr_Permutations

#apply filter (if filter is before or after loop)
#only protein-coding
#Chr_PermutationsNew <- Chr_PermutationsNew[Chr_PermutationsNew$logicalProteinType == TRUE,]

#remove unsignificant entries for better overview
#Chr_PermutationsNew <- Chr_PermutationsNew[Chr_PermutationsNew$PValAdjusted == "significant",]

#filter against unchanged expression
#Chr_PermutationsNew <- Chr_PermutationsNew[!Chr_PermutationsNew$ExpressionPatternVSCtrl == "unchanged",]




#filter for essential columns to enter loop, geneID, patient, pval, expression
Chr_PermutationsNew <- Chr_PermutationsNew[,c(1,8,17,20)]


i <- 1
no <- 1
nPseudoPatients <- 1000 
set.seed(no)
Chr_glEf_P_c_m_numbers <- as.data.frame(matrix(data=NA, nrow=nPseudoPatients, ncol=4))
SkipPermutating <- TRUE

if (SkipPermutating == FALSE) {

for (i in 1:nPseudoPatients) {
  
  #subdat
  Chr_Permutations_subdat <- Chr_PermutationsNew[Chr_PermutationsNew$Patient == PermutationsPatient,]
  
  #create pseudo patient
  PseudoPatient <- paste(PermutationsPatient, "_", i, sep = "")
  Chr_Permutations_subdat$Patient <- PseudoPatient

  #generate pseudo data
  #Chr_Permutations_subdat[(1+(i-1)*nrow(Chr_Permutations)): (i*nrow(Chr_Permutations))] <- sample(as.numeric(Chr_Permutations$log2FoldChange, size = nrow(Chr_Permutations)))
  Chr_Permutations_subdat$ENSG <- sample(Chr_Permutations_subdat$ENSG, replace = FALSE)

  #merge with reference
  Chr_Permutations_subdat <- rbind(Chr_Permutations_subdat, Chr_PermutationsNew[Chr_PermutationsNew$Patient == PermutationsReference,])
  
  #add old annotation to permutated data
  Chr_Permutations_subdat <- inner_join(Chr_Permutations_subdat, geneNamesTrimmed, by = "ENSG", all=FALSE) 
  
  #filter for protein-coding genes
  Chr_Permutations_subdat <- Chr_Permutations_subdat[Chr_Permutations_subdat$Gene.type == "protein_coding",]
  
  #remove unsignificant entries for better overview
  Chr_Permutations_subdat <- Chr_Permutations_subdat[Chr_Permutations_subdat$PValAdjusted == "significant",]

  #filter against unchanged expression
  Chr_Permutations_subdat <- Chr_Permutations_subdat[!Chr_Permutations_subdat$ExpressionPatternVSCtrl == "unchanged",]

  #filter for common genes: DSmerge and loop[i]
  Chr_glEf_Permutations_common <- Chr_Permutations_subdat
  
  #robust against double entries with two dif annotations for same patient!
  nGenesPseudoPatient_keep <- c(unique(as.character(Chr_glEf_Permutations_common$ENSG[Chr_glEf_Permutations_common$Patient == PseudoPatient])))
  nGenesReference_keep <- c(unique(as.character(Chr_glEf_Permutations_common$ENSG[Chr_glEf_Permutations_common$Patient == PermutationsReference])))
  keep <- c(nGenesPseudoPatient_keep,nGenesReference_keep)
  DifExpGenesInclNotShared <- unique(keep) #number of dif. expr. genes; needed for ratio calculations
  DifExpGenesOnlyShared <- keep[duplicated(keep)]
  Chr_glEf_Permutations_common <- Chr_glEf_Permutations_common[Chr_glEf_Permutations_common$ENSG %in% DifExpGenesOnlyShared,] 
  
#idea 1: filter for DSmerge up/down + ALI09 down/up
  #going up
  keep1 <- vector()
  keep1 <- c(keep1, as.character(Chr_glEf_Permutations_common$ENSG[(Chr_glEf_Permutations_common$Patient == PseudoPatient & 
                                                                       Chr_glEf_Permutations_common$ExpressionPatternVSCtrl == "upregulated")]))
  keep1 <- c(keep1, as.character(Chr_glEf_Permutations_common$ENSG[(Chr_glEf_Permutations_common$Patient == PermutationsReference & 
                                                                       Chr_glEf_Permutations_common$ExpressionPatternVSCtrl == "downregulated")]))
  keep1 <- keep1[duplicated(keep1)]
  
  #going down
  keep2 <- vector()
  keep2 <- c(keep2, as.character(Chr_glEf_Permutations_common$ENSG[(Chr_glEf_Permutations_common$Patient == PseudoPatient & 
                                                                       Chr_glEf_Permutations_common$ExpressionPatternVSCtrl == "downregulated")]))
  keep2 <- c(keep2, as.character(Chr_glEf_Permutations_common$ENSG[(Chr_glEf_Permutations_common$Patient == PermutationsReference & 
                                                                       Chr_glEf_Permutations_common$ExpressionPatternVSCtrl == "upregulated")]))
  keep2 <- keep2[duplicated(keep2)]
  
  #unite
  keep <- c(keep1, keep2)
  
  #clear df for keep
  Chr_glEf_Per_c_mirror <- Chr_glEf_Permutations_common[Chr_glEf_Permutations_common$ENSG %in% keep,]




  #add number of mirrored genes
  row.names(Chr_glEf_P_c_m_numbers)[i] <- paste(PseudoPatient, "_", PermutationsReference, sep="")
  Chr_glEf_P_c_m_numbers[i,1] <-  nrow(Chr_glEf_Per_c_mirror)/2     #sum of mirrored genes per patient, V1
  Chr_glEf_P_c_m_numbers[i,2] <- length(DifExpGenesInclNotShared)   #sum of unshared DE genes, V2
  Chr_glEf_P_c_m_numbers[i,3] <- length(DifExpGenesOnlyShared)      #sum of shared DE genes, V3
  Chr_glEf_P_c_m_numbers[i,4] <- length(nGenesPseudoPatient_keep)   #
  
  

} #close all_runs loop

#calculate ratio for random mirrored genes
Chr_glEf_P_c_m_numbers$ratioAll <- Chr_glEf_P_c_m_numbers$V1 / Chr_glEf_P_c_m_numbers$V2
Chr_glEf_P_c_m_numbers$ratioShared <- Chr_glEf_P_c_m_numbers$V1 / Chr_glEf_P_c_m_numbers$V3
Chr_glEf_P_c_m_numbers$ratioPseudoPatient <- Chr_glEf_P_c_m_numbers$V1 / Chr_glEf_P_c_m_numbers$V4
Chr_glEf_P_c_m_numbers$SampleName <- rownames(Chr_glEf_P_c_m_numbers)

#save df
write_xlsx(Chr_glEf_P_c_m_numbers, paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/","Permutation for all genes", "/", "ListExport_", PlotName, "_MirrorEf_", no, ".xlsx", sep=""))
write.table(Chr_glEf_P_c_m_numbers, paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/","Permutation for all genes", "/", "ListExport_", PlotName, "_MirrorEf_", no, ".tsv", sep=""), quote = FALSE, row.names = FALSE)
}
# 
# #test plot
# test <- as.data.frame(Chr_glEf_P_c_m_numbers)
# test$comp <- rownames(test)
# test$ali9 <- FALSE
# test$ali9[grepl(test$comp, pattern = "ALI05")] <- TRUE
# ggplot(test, aes(x= comp, y=V1, col=ali9)) + geom_point()
# 


#load all files in case the permutation was outsourced(this case: yes)
ListSamplesPer <- list.files(path = paste("R stuff/ExtractedData_", ProjectName, "/","Permutation for all genes", "/", sep = ""), pattern = c("Permutations_MirrorEf"))
ListSamplesPer <- ListSamplesPer[grepl(x=ListSamplesPer, pattern = ".tsv")]
ListSamplesPer <- ListSamplesPer[!grepl(x=ListSamplesPer, pattern = "statsOutput")]

#import data, Permutation files
i <- 1
Chr_glEf_P_c_m_numbers <- as.data.frame(matrix(data=NA, nrow=length(ListSamplesPer)*nPseudoPatients, ncol=8))
colnames(Chr_glEf_P_c_m_numbers) <- colnames(x = read.delim(file = paste("R stuff/ExtractedData_", ProjectName, "/","Permutation for all genes", "/", ListSamplesPer[i], sep=""), sep = " ", header = TRUE))
for(i in 1:length(ListSamplesPer)){
  
  Chr_glEf_P_c_m_numbers[((i-1)*nPseudoPatients+1):(i*nPseudoPatients),] <- read.delim(file = paste("R stuff/ExtractedData_", ProjectName, "/","Permutation for all genes", "/", ListSamplesPer[i], sep=""), sep = " ", header = FALSE, skip = 1, stringsAsFactors = FALSE)

}


#mir: number of mirrored genes per patient[i] for V1
Length <- length(Chr_glEf_P_c_m_numbers)
Per_binomial <- data.frame(Length)
Per_binomial$Min_mir <- as.vector(summary(Chr_glEf_P_c_m_numbers$V1))[1]
Per_binomial$Max_mir <- as.vector(summary(Chr_glEf_P_c_m_numbers$V1))[6]
Per_binomial$Range_mir <- Per_binomial$Max_mir - Per_binomial$Min_mir
Per_binomial$Median_mir <- as.vector(summary(Chr_glEf_P_c_m_numbers$V1))[3]
Per_binomial$Mean_mir <- as.vector(summary(Chr_glEf_P_c_m_numbers$V1))[4]

#all: number of total dif.expressed(sign.) genes in patient[i] for V2
Per_binomial$Min_all <- as.vector(summary(Chr_glEf_P_c_m_numbers$V2))[1]
Per_binomial$Max_all <- as.vector(summary(Chr_glEf_P_c_m_numbers$V2))[6]
Per_binomial$Range_all <- Per_binomial$Max_all - Per_binomial$Min_all
Per_binomial$Median_all <- as.vector(summary(Chr_glEf_P_c_m_numbers$V2))[3]
Per_binomial$Mean_all <- as.vector(summary(Chr_glEf_P_c_m_numbers$V2))[4]

#ratio: ratio of mirrored genes/total genes; =FDR per patient [i]
Per_binomial$Min_ratioAll <- as.vector(summary(Chr_glEf_P_c_m_numbers$ratioAll))[1]
Per_binomial$Max_ratioAll <- as.vector(summary(Chr_glEf_P_c_m_numbers$ratioAll))[6]
Per_binomial$Range_ratioAll <- Per_binomial$Max_ratioAll - Per_binomial$Min_ratioAll
Per_binomial$Median_ratioAll <- as.vector(summary(Chr_glEf_P_c_m_numbers$ratioAll))[3]
Per_binomial$Mean_ratioAll <- as.vector(summary(Chr_glEf_P_c_m_numbers$ratioAll))[4]

#ratio: ratio of ratioShared
Per_binomial$Median_ratioShared <- as.vector(summary(Chr_glEf_P_c_m_numbers$ratioShared))[3]
Per_binomial$Mean_ratioShared <- as.vector(summary(Chr_glEf_P_c_m_numbers$ratioShared))[4]

#ratio: ratio of ratioPseudoPatient
Per_binomial$Median_ratioPseudoPatient <- as.vector(summary(Chr_glEf_P_c_m_numbers$ratioPseudoPatient))[3]
Per_binomial$Mean_ratioPseudoPatient <- as.vector(summary(Chr_glEf_P_c_m_numbers$ratioPseudoPatient))[4]

#Standard error for global FDR
se <- function(x) sqrt(var(x)/length(x))
Per_binomial$SE_ratio <- se(Chr_glEf_P_c_m_numbers$ratioAll)

#save table
# write.table(x=Per_binomial, file = paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/", "ListExport_", PlotName, "_MirrorEf_statsOutput.tsv", sep=""), quote = FALSE, row.names = FALSE)
write_xlsx(Per_binomial, paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/", "Permutation for all genes", "/", "ListExport_", PlotName, "_MirrorEf_statsOutput_filterAfterPer.xlsx", sep=""))




#Binomial test


Test_mirroredOnRatioAll <- binom.test(
  sum(as.logical(MirroredTableRef$mirrored))/2,
  sum(MirroredTableRef$Patient == PermutationsPatient),
  Per_binomial$Mean_ratioAll,
  alternative="greater")
Test_mirroredOnRatioAll
Test_mirroredOnRatioAll$p.value


Test_mirroredOnRatioShared <- binom.test(
  sum(as.logical(MirroredTableRef$mirrored))/2,
  sum(MirroredTableRef$Patient == PermutationsPatient),
  Per_binomial$Mean_ratioPseudoPatient,
  alternative="greater")
Test_mirroredOnRatioShared
Test_mirroredOnRatioShared$p.value
#with hg38: 0.0011










#----only Chr21 - all Patients----
#2nd plot: only chr21, with loop for each patient
# sig and nonsig hits

PlotName <- "Chr21Del"

#take only Chr21, decide between significant or all points
#Chr_all.21 <- Chr_all.red[Chr_all.red$Chromosome.scaffold.name == "21",]
Chr_all.21 <- Chr_all.redSig[!is.na(Chr_all.redSig$Chromosome.scaffold.name),]
Chr_all.21 <- Chr_all.21[Chr_all.21$Chromosome.scaffold.name == "21",]


#add colour code
Chr_all.21$Position <- as.character(Chr_all.21$Position)
Chr_all.21$Position[Chr_all.21$logicalProteinType == FALSE] <- "not protein-coding"

#Filter for DE
Chr_all.21 <- Chr_all.21[Chr_all.21$ExpressionPatternVSCtrl != "unchanged",]

#reorder levels
Chr_all.21$Position <- factor(Chr_all.21$Position, levels(factor(Chr_all.21$Position))[c(1,3,5,4,6,2)])

#rename one level
if ("surroundings_1MB" %in% levels(factor(Chr_all.21$Position))) {
  levels(Chr_all.21$Position)[levels(Chr_all.21$Position)=="surroundings_1MB"] <- "SV proximity"
}else{
  warning("level was not renamed due absence and therefore not renamed!")
}

#----chr21 only_new plot
i <- 1

# i <- 8 #ALI09
for(i in 1:length(levels(factor(Chr_all.red$Patient)))){
  
  #exclude all Patients but chosen one
  Chr_all.21_sub <- Chr_all.21[Chr_all.21$Patient == levels(factor(Chr_all.red$Patient))[i],]

  #Plot
  plot.21only <- ggplot(Chr_all.21_sub, aes(x=CoordStart_KB, y=log2FoldChange))+
    geom_hline(yintercept = log(thresholdForExpr, 2), alpha = 0.5,size=0.2)+
    geom_hline(yintercept = log(1-(thresholdForExpr-1), 2), alpha = 0.5,size=0.2)+
    geom_hline(yintercept = 0, linetype=2, alpha = 0.2,size=0.2)+
    geom_point(aes(color = factor(Position)), size=6, alpha=0.3, shape=16)+
    scale_color_manual(values=colorsForLegend[c(9,5,6,3,2,1)])+
    theme_classic(base_size=20)+
    scale_x_continuous(limits = c(0,46691), breaks = c(0,10000, 20000, 30000, 40000, 46000))+
    scale_y_continuous(limits= c(-11,6), breaks = seq(-11,6,1))+
    labs(x = "Chromosome 21 [kb]", y = "Fold Change [log2]", color = "Colour code", alpha = "Gene type") +
    guides(
     color=guide_legend(order=1),
     shape=guide_legend(order=2),
     alpha=guide_legend(order=3)
    )+
    theme(
     element_line(size = 1.0),
     plot.subtitle = element_text(colour="black",face="italic", size = 15),
     axis.ticks = element_line(colour="black"),
     axis.text.x = element_text(colour="black", hjust=0.5),
     axis.text.y = element_text(colour="black", hjust=0.5)
    )+
    geom_text_repel(label = ifelse(Chr_all.21_sub$FoldChange %in% c(tophits(Chr_all.21_sub$FoldChange[Chr_all.21_sub$PValAdjusted == "significant" & Chr_all.21_sub$logicalProteinType == TRUE]),
                                                                    tophits(Chr_all.21_sub$FoldChange[Chr_all.21_sub$PValAdjusted == "significant" & Chr_all.21_sub$logicalProteinType == TRUE], decreasing = FALSE)),
                                   paste(Chr_all.21_sub$Gene.name),
                                   ""),
                    col = "black", size = 3.5, #for label settings
                    xlim=c(NA, 48000000),
                    # ylim=c(NA, 30),  
                    hjust = 0.3, #no change
                    segment.alpha = 0.13, #line transparency
                    nudge_y = 0.3,
                    box.padding = 1,
                    max.overlaps = 100,
                    angle = 0)+ # to which direction should the text be aligend to,turn text
    
    annotate(geom = "rect", xmin = firstDeletion[1]/conversionFactorKB, xmax = firstDeletion[2]/conversionFactorKB, ymin = -11, ymax = 6, alpha=0.08, fill = colorsForLegend[9])+
    annotate(geom = "rect", xmin = firstDeletion[2]/conversionFactorKB, xmax = surroundingsDeletionP[2]/conversionFactorKB, ymin = -11, ymax = 6, alpha=0.08, fill = colorsForLegend[3])+
    annotate(geom = "rect", xmin = secondDeletion[1]/conversionFactorKB, xmax = secondDeletion[2]/conversionFactorKB, ymin = -11, ymax = 6, alpha=0.08, fill = colorsForLegend[5])+
    annotate(geom = "rect", xmin = thirdDeletion[1]/conversionFactorKB, xmax = thirdDeletion[2]/conversionFactorKB, ymin = -11, ymax = 6, alpha=0.08, fill = colorsForLegend[6])+
    annotate(geom = "rect", xmin = surroundingsDeletion[1]/conversionFactorKB, xmax = secondDeletion[1]/conversionFactorKB, ymin = -11, ymax = 6, alpha=0.08, fill = colorsForLegend[3])+
    annotate(geom = "rect", xmin = secondDeletion[2]/conversionFactorKB, xmax = thirdDeletion[1]/conversionFactorKB, ymin = -11, ymax = 6, alpha=0.08, fill = colorsForLegend[3])+
    
    ggtitle(paste("RNA Seq Data, NES, ", PlotName, sep=""), 
            subtitle = paste(levels(factor(Chr_all.red$Patient))[i], "vs", RefOfSeqDat, ", threshold low deltaExpr: |",thresholdForExpr, "|x; threshold pVal: ", thresholdForPval, "; ", "PValAdjusted", sep =""))
  
  # plot.21only
  
   ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_", levels(factor(Chr_all.red$Patient))[i], "vs", RefOfSeqDat, "_", "PValAdjusted", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         height = 18, width = 30, unit="cm", limitsize = FALSE, dpi = 600)
  
}


# ggsave(filename = paste("R stuff/test_all.png", sep=""),
#        height = 100, width = 100, unit="cm",limitsize = FALSE)



#----Deleted part on chr21, all patients-------------
PlotName <- "Chr21DeletedRegion"

#take only Chr21, decide between significant or all points
Chr_all.21 <- Chr_all.redSig[Chr_all.redSig$Chromosome.scaffold.name == "21",]

#delete genes upstream of deletions
Chr_all.21 <- Chr_all.21[Chr_all.21$Gene.start..bp. > surroundingsDeletion[1]-1,]

#add colour code
Chr_all.21$Position <- as.character(Chr_all.21$Position)
Chr_all.21$Position[Chr_all.21$logicalProteinType == FALSE] <- "not protein-coding"

#Filter for DE
Chr_all.21 <- Chr_all.21[Chr_all.21$ExpressionPatternVSCtrl != "unchanged",]

#rename one level
if ("surroundings_1MB" %in% levels(factor(Chr_all.21$Position))) {
  Chr_all.21$Position[Chr_all.21$Position == "surroundings_1MB"] <- "SV proximity"
}else{
  warning("level was not renamed due absence and therefore not renamed!")
}

#reorder levels
if(length(levels(factor(Chr_all.21$Position))) != length(c(2,4,3,1))){
  warning("Levels in CHR wrong")
}else{
  message("Levels in CHR okay")
}
Chr_all.21$Position <- factor(Chr_all.21$Position, levels(factor(Chr_all.21$Position))[c(2,4,3,1)])


i <- 1

# for(i in 1:length(levels(factor(Chr_all.red$Patient)))){
for(i in which(levels(factor(Chr_all.21$Patient)) %in% c("ALI09", "DSB2", "DSB3", "DSBmerge"))){
  
  #exclude all Patients but chosen one
  Chr_all.21_sub <- Chr_all.21[Chr_all.21$Patient == levels(factor(Chr_all.21$Patient))[i],]
  
    #Plot
    plot.21only <- ggplot(Chr_all.21_sub, aes(x=CoordStart_KB, y=log2FoldChange))+
      geom_hline(yintercept = log(thresholdForExpr, 2), alpha = 0.5,size=0.2)+
      geom_hline(yintercept = log(1-(thresholdForExpr-1), 2), alpha = 0.5,size=0.2)+
      geom_hline(yintercept = 0, linetype=2, alpha = 0.2,size=0.2)+
      # geom_smooth(method = "loess", show.legend = FALSE, col = colorsForLegend[9], se = FALSE, size=0.6, span =1)+
      geom_point(aes(color = factor(Position)), size=6, alpha=0.3, shape=16)+
      # scale_alpha_discrete(range=c(0.1,0.9), labels = c("other", "protein coding"))+
      scale_color_manual(values=colorsForLegend[c(5,6,3,1)])+
      theme_classic(base_size=20)+
      scale_x_continuous(limits = c(40900,46691), breaks = seq(40000, 46000, 1000))+
      scale_y_continuous(limits= c(-4, 4), breaks = seq(-3,3,1))+
      labs(x = "Chromosome 21 [kb]", y = "Fold Change [log2]", color = "Colour Code", alpha = "Gene type") +
      guides(
        color=guide_legend(order=1),
        shape=guide_legend(order=2),
        alpha=guide_legend(order=3)
      )+
      theme(
        element_line(size = 1.0),
        plot.subtitle = element_text(colour="black",face="italic", size = 15),
        axis.ticks = element_line(colour="black"),
        axis.text.x = element_text(colour="black", hjust=0.5),
        axis.text.y = element_text(colour="black", hjust=0.5)
      )+
      geom_text_repel(label = ifelse(Chr_all.21_sub$Gene.name %in% Chr_all.21_sub$Gene.name[
                        order(Chr_all.21_sub$FoldChange)[(length(Chr_all.21_sub$FoldChange)-round((length(Chr_all.21_sub$FoldChange)-length(Chr_all.21_sub$FoldChange)/2),0)):length(Chr_all.21_sub$FoldChange)]],
                        as.character(Chr_all.21_sub$Gene.name), ""),
                      col = "black", size = 6, #for label settings
                      xlim=c(NA, 46000),
                      ylim=c(3,NA),
                      # force=5, #no change
                      vjust = 0, #no change
                      hjust =0, #no change
                      segment.alpha = 0.2, #line transparency
                      #nudge_x = 1,
                      box.padding = 0.1,
                      #nudge_y = 0.1,
                      direction = "x",
                      max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                      angle = 35)+ # to which direction should the text be aligend to,turn text
      geom_text_repel(label = ifelse(Chr_all.21_sub$Gene.name %in% Chr_all.21_sub$Gene.name[
                        order(Chr_all.21_sub$FoldChange)[1:round((length(Chr_all.21_sub$FoldChange)-length(Chr_all.21_sub$FoldChange)/2),0)]],
                        as.character(Chr_all.21_sub$Gene.name), ""),
                      col = "black", size = 6, #for label settings
                      xlim=c(NA, 46000),
                      ylim=c(NA,-3.2),
                      #force=0.4, #no change
                      vjust = 1, #no change
                      hjust =0.5, #no change
                      segment.alpha = 0.2, #line transparency
                      #nudge_y = 1,
                      box.padding = 0.2,
                      #nudge_y = 0.1,
                      direction = "x",
                      max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                      angle = 35)+ # to which direction should the text be aligend to,turn text

    annotate(geom = "rect", xmin = secondDeletion[1]/conversionFactorKB, xmax = secondDeletion[2]/conversionFactorKB, ymin = -3, ymax = 3, alpha=0.08, fill = colorsForLegend[5])+
    annotate(geom = "rect", xmin = thirdDeletion[1]/conversionFactorKB, xmax = thirdDeletion[2]/conversionFactorKB, ymin = -3, ymax = 3, alpha=0.08, fill = colorsForLegend[6])+
    annotate(geom = "rect", xmin = surroundingsDeletion[1]/conversionFactorKB, xmax = secondDeletion[1]/conversionFactorKB, ymin = -3, ymax = 3, alpha=0.08, fill = colorsForLegend[3])+
    annotate(geom = "rect", xmin = secondDeletion[2]/conversionFactorKB, xmax = thirdDeletion[1]/conversionFactorKB, ymin = -3, ymax = 3, alpha=0.08, fill = colorsForLegend[3])+
      
    ggtitle(paste("RNA Seq Data, NES, ", PlotName, sep=""), 
              subtitle = paste(levels(factor(Chr_all.red$Patient))[i], "vs", RefOfSeqDat, ", threshold low deltaExpr: |",thresholdForExpr, "|x; threshold pVal: ", thresholdForPval, "; ", "PValAdjusted", sep =""))
    
  # plot.21only

    ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_", levels(factor(Chr_all.red$Patient))[i], "vs", RefOfSeqDat, "_", "PValAdjusted", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
           height = 18, width = 35, unit="cm", limitsize = FALSE, dpi = 600)
    
}





#------volcano-patients separated/no circles-------------
PlotName <- "VolcanoSep"

#volcano plot
Chr_all.volcano <- Chr_all.redSig
Chr_all.volcano$DifExpr <- FALSE
Chr_all.volcano$DifExpr[Chr_all.volcano$ExpressionPatternVSCtrl != "unchanged"] <- TRUE
Chr_all.volcano$groups <- ifelse(Chr_all.volcano$log2FoldChange > log(thresholdForExpr, 2) & Chr_all.volcano$padj < thresholdForPval, "upregulated",
                                 ifelse(Chr_all.volcano$log2FoldChange < log(1-(thresholdForExpr-1),2) & Chr_all.volcano$padj < thresholdForPval, "downregulated",
                                        "unchanged/ns"))

#for loop for each patient alone
h <- 1

for(h in 1:length(ListPatient)){
  Chr_all.volcano.subset <- Chr_all.volcano[Chr_all.volcano$Patient == ListPatient[h],]
  
  #ggplot(Chr_all.redSig, aes(x=log2FoldChange, y=-log10(padj), col= ExpressionPatternVSCtrl, shape=Position))+
  ggplot(Chr_all.volcano.subset, aes(x = log2FoldChange, y = -log10(padj)))+
    theme_classic(base_size=20)+
    geom_point(aes(color= groups))+
    #geom_point()+
    scale_color_manual(values=c("red", "grey", "blue", "black"))+
    labs(x="log2 foldChange", y="-log10 p-value", shape = "Position")+
    guides(color=FALSE)+
    scale_x_continuous(limits = c(-25.5,12), breaks = seq(-25,10,5))+
    scale_y_continuous(limits= c(0, 300), breaks = seq(0,300,50))+
    geom_text_repel(aes(x = log2FoldChange,y = -log10(padj)),
                    label= ifelse((Chr_all.volcano.subset$padj %in% tophits(Chr_all.volcano.subset$padj, n = 7, decreasing = FALSE)) |
                                    (abs(Chr_all.volcano.subset$log2FoldChange) %in% tophits(abs(Chr_all.volcano.subset$log2FoldChange[Chr_all.volcano.subset$PValAdjusted == "significant"]), n = 7, decreasing = TRUE)),
                                  paste(Chr_all.volcano.subset$Gene.name),
                                  ""),
                    col = "black", size = 3.5, #for label settings
                    #ylim = c(1.8, NA), #for y-axis movements
                    #xlim=c(NA, 57000000),
                    #xlim=c(NA, 48000000),
                    #force=3, #no change
                    #vjust = -1, #no change
                    hjust = .2, #horizontal offset, similar to nudge_x
                    segment.alpha = 0.2, #line transparency
                    # nudge_y = 2.1 - Chr_all.redSig$FoldChange[!is.na(Chr_all.redSig$FoldChange)], #no need for this
                    nudge_y = .1,
                    #box.padding = .1,
                    #direction = "x",
                    max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                    angle = 0)+ # to which direction should the text be aligend to,turn text
    # --
    theme(
      plot.subtitle = element_text(colour="black",face="italic", size = 14),
      axis.ticks = element_line(colour="black"),
      axis.text = element_text(colour="black", hjust=0.5),
      legend.position = c(0.15, 0.6),
      legend.background = element_rect(linetype="solid", colour ="black"))+
    ggtitle("RNA Seq Data, NES, Chr_all", subtitle = paste(ListPatient[h], ", threshold low deltaExpr: |",thresholdForExpr, "|x; threshold pVal: ", thresholdForPval, "; PValadjusted", sep =""))
  
  
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, ListPatient[h],"_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         height = 18, width = 30, unit="cm")
}

ggsave(filename = paste("R stuff/test_all.png", sep=""),
       height = 18, width = 30, unit="cm",limitsize = FALSE)



#####continue here with the liftover, check for the next thing to do!!!!!

#------skewness_single_manual
#skewness analysis, AI09, DS2, DS3
PlotName <- "Skewness"

#get full list with significant and unsignificant hits, all chromosomes, all patients, no BaseMean=0
Chr_all.skewness.manual <- Chr_all.red[Chr_all.red$ENSG %in% geneNamesCellTypes_manual$ENSG,]

#create empty columns to fill later
Chr_all.skewness.manual <- cbind(Chr_all.skewness.manual, data.frame(NES=F, Neuron=F, Glia=F, CellType=NA))

#fill groupings
Chr_all.skewness.manual$NES[Chr_all.skewness.manual$Gene.name %in% ListNESCells] <- TRUE
Chr_all.skewness.manual$Neuron[Chr_all.skewness.manual$Gene.name %in% ListNeuronalGenes] <- TRUE
Chr_all.skewness.manual$Glia[Chr_all.skewness.manual$Gene.name %in% ListGliaCells] <- TRUE

Chr_all.skewness.manual$CellType[Chr_all.skewness.manual$NES == TRUE] <- "NES"
Chr_all.skewness.manual$CellType[Chr_all.skewness.manual$Neuron == TRUE] <- "Neuron"
Chr_all.skewness.manual$CellType[Chr_all.skewness.manual$Glia == TRUE] <- "Glia"

#reorder factors for labelling
Chr_all.skewness.manual$Gene.name <- factor(Chr_all.skewness.manual$Gene.name, levels = levels(factor(Chr_all.skewness.manual$Gene.name))[c(5,10,  3:4,6:9,13,  1:2,11:12,14)])
Chr_all.skewness.manual$CellType <- factor(Chr_all.skewness.manual$CellType, levels = levels(factor(Chr_all.skewness.manual$CellType))[c(2,3,1)])

#plot_detail

#------skewness_summarized_manual
#plot_summarized
plot_summarized <- ggplot(Chr_all.skewness.manual, aes(x=Chr_all.skewness.manual$Patient, y=Chr_all.skewness.manual$log2FoldChange, fill = Chr_all.skewness.manual$Patient)) +
  geom_hline(yintercept = 0, linetype =3, alpha = 0.5,size=0.2)+
  geom_hline(yintercept = log2(thresholdForExpr), linetype =2, alpha = 0.5, size=0.2)+
  geom_hline(yintercept = log2(1-(thresholdForExpr-1)), linetype =2, alpha = 0.5, size=0.2)+
  theme_classic(base_size=20) +
  facet_grid(~CellType, scales = "free_x", space = "free_x") +
  geom_boxplot(outlier.shape = NA, size = 0.8)+
  geom_jitter(width = 0.1, show.legend = NA)+
  stat_summary(fun.y=mean, colour="red3", geom="point", shape=18, size=4, show.legend = NA)+
  labs(x = "", y = "log2Foldchange", fill = element_blank()) +
  scale_fill_manual(values=colorsForLegend[c(3,7,8,6)])+
  theme(
    plot.subtitle = element_text(colour="black",face="italic", size = 15),
    axis.ticks = element_line(colour="black"),
    axis.text.x = element_text(colour="black", hjust=0.5),
    axis.text.y = element_text(colour="black", hjust=0.5),
    #legend.title = element_text(colour = "black", size = 15),
    strip.background = element_rect(size=1.5, linetype="solid", colour ="black"),
    legend.position = "bottom") +
  guides(fill=FALSE)+
  ggtitle(paste("RNA Seq Data, NES, Chr_all, ", PlotName, sep=""), subtitle = paste( "AllSamples vs. ", RefOfSeqDat, ", threshold pVal: ", thresholdForPval, "; ", "p-valAdjusted", sep =""))

plot_summarized

ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_", "allSamples_summarized", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       height = 18, width = 35, unit="cm", limitsize = FALSE)

# ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/temp.png", sep=""),
#      height = 18, width = 30, unit="cm")




#------skewness_GO terms------------
PlotName <- "SkewnessGO"

#reduce Lists to work with it better, exclude "GOA", list is implemented in "load biomart lists
geneNamesCellTypes_GOterms_red <- subset(geneNamesCellTypes_GOterms, select = -c(GOSlim.GOA.Accession.s., GOSlim.GOA.Description, GO.domain))
geneNamesCellTypes_GOterms_red <- geneNamesCellTypes_GOterms_red[!duplicated(geneNamesCellTypes_GOterms_red),]
#!duplicated() is equal to unique(), but duplicated returns a logical vector, unique returns a vector

#get full list with significant and unsignificant hits, all chromosomes, all patients, no BaseMean=0, GO:0022008 neurogenesis
Chr_all.skewness.GO <- Chr_all.red[Chr_all.red$ENSG %in% geneNamesCellTypes_GOterms_red$ENSG,]

#merge columns, add GO terms to genes
Chr_all.skewness.GO <- merge(Chr_all.red, geneNamesCellTypes_GOterms_red, 
              by.x = c("ENSG", "Gene.start..bp.", "Gene.end..bp.", "Chromosome.scaffold.name", "Gene.name", "HGNC.symbol"), 
              by.y = c("ENSG", "Gene.start..bp.", "Gene.end..bp.", "Chromosome.scaffold.name", "Gene.name", "HGNC.symbol"))

##reorder factors for labelling
#get the GO terms from GO.list specified in the input section: GO:xx -> GO-Name[all from input]
ListOrderNames <- Chr_all.skewness.GO$GO.term.name[Chr_all.skewness.GO$GO.term.accession %in% ListGOterms_all]
ListOrderNames <- unique(ListOrderNames)

# #Get Index of GO terms in current df[GO:xx]
#[GO-Name]
ListOrderPosition  <- which(levels(Chr_all.skewness.GO$GO.term.name) %in% ListOrderNames)

#define new order, [GO-Name]
levels(Chr_all.skewness.GO$GO.term.name)[ListOrderPosition]
ListOrderPosition <- ListOrderPosition[c(9,2,1,7,8,10,11,12,6,3,4,5)]
ListOrderOrderNOT <- which(!levels(Chr_all.skewness.GO$GO.term.name) %in% ListOrderNames)
Chr_all.skewness.GO$GO.term.name <- factor(Chr_all.skewness.GO$GO.term.name, 
                                                levels = levels(factor(Chr_all.skewness.GO$GO.term.name))[c(ListOrderPosition,ListOrderOrderNOT)])

#filter for significant hits
if(RemoveUnSigForPlots == TRUE){
  Chr_all.skewness.GO <- Chr_all.skewness.GO[Chr_all.skewness.GO$PValAdjusted =="significant",]
}

#checking that all chosen GO terms have data for each , and every patient should have at least one data point
Chr_all.skewness.GO.subdat <- Chr_all.skewness.GO[Chr_all.skewness.GO$GO.term.accession %in% ListGOterms_all,]

Chr_all.skewness.GO.subdat <- Chr_all.skewness.GO.subdat %>% 
  group_by(GO.term.name, add = TRUE) %>%
  group_by(Patient, add = TRUE) %>%
  summarise(counted = length(GO.term.name))

ListSkewness.GO.drop <- names(summary(as.factor(Chr_all.skewness.GO.subdat$GO.term.name)))[as.vector(summary(as.factor(Chr_all.skewness.GO.subdat$GO.term.name)) != length(levels(factor(Chr_all.skewness.GO.subdat$Patient))))]
Chr_all.skewness.GO <- Chr_all.skewness.GO[!Chr_all.skewness.GO$GO.term.name %in% ListSkewness.GO.drop,]


#------X-Y plots-----------
#Y plots
PlotName <- "YChr"

#get dataframe, including non-significant, chr-y, PAR genes(Biomart excludes them from Y), inclduing basemean > 5, all other chr, no x
Chr_all.ChrXY <- Chr_all.red[Chr_all.red$Chromosome.scaffold.name == "Y" | Chr_all.red$Gene.name %in% ListPARGenes,]

#label PAR Genes
Chr_all.ChrXY$ParRegion <- FALSE
Chr_all.ChrXY$ParRegion[Chr_all.ChrXY$Gene.name %in% ListPARGenes] <- TRUE

#manual labelling PAR2 Genes, assumption chromosome:GRCh38:Y:56887903 - 57217415 is shared with X: 155701383 - 156030895 (PAR2)
ChrY_max <- 57217415
Chr_all.ChrXY$Gene.start..bp.[Chr_all.ChrXY$Gene.start..bp. > ChrY_max] <- ChrY_max

#Group for low Basemean
Chr_all.ChrXY$GroupBaseMean <- "low"
Chr_all.ChrXY$GroupBaseMean[Chr_all.ChrXY$baseMean > thresholdForBaseMean] <- paste(">", thresholdForBaseMean, sep="")
Chr_all.ChrXY$GroupBaseMean <- as.factor(Chr_all.ChrXY$GroupBaseMean)

#filter for Patient: ALI9, DSB2, DSB3, DSm
Chr_all.ChrXY <- Chr_all.ChrXY[Chr_all.ChrXY$Patient %in% c("ALI09", "DSB2", "DSB3", "DSBmerge"),]

#plot
plot.ChrXY <-ggplot(Chr_all.ChrXY, aes(x=`Gene.start..bp.`, y=log2(FoldChange)))+
  geom_point(aes(shape=Patient, color = ParRegion, alpha=PValAdjusted, fill=GroupBaseMean), stroke=1.5, size = 3)+
  # geom_point(inherit.aes = FALSE, aes(x=`Gene.start..bp.`, y=FoldChange), shape = 1, size = 4, alpha = 0.5)+
  scale_shape_manual(values=c(21, 24, 25, 22))+
  scale_alpha_discrete(range=c(.1,1))+
  scale_color_manual(labels = c("outside", "inside"), values=colorsForLegend[c(2,5)])+
  scale_fill_manual(values=colorsForLegend[c(9,10)])+
  theme_classic(base_size=20)+
  scale_x_continuous(limits = c(0,ChrY_max), breaks = c(seq(0,ChrY_max, 10000000),ChrY_max))+
  #scale_y_continuous(trans='log10', limits= c(0.01, 700), breaks = c(0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50))+
  labs(x = "Chromosome Y", y = "Fold Change (log2)", color = "PAR-Genes", shape = "Patient", alpha = "p-value", fill = "BaseMean") +
  guides(
    fill=guide_legend(order=1, override.aes=list(shape=21)),
    color=guide_legend(order=3),
    shape=guide_legend(order=2),
    alpha=guide_legend(order=4)
  )+
  theme(
    element_line(size = 1.0),
    plot.subtitle = element_text(colour="black",face="italic", size = 15),
    axis.ticks = element_line(colour="black"),
    axis.text.x = element_text(colour="black", hjust=0.5),
    axis.text.y = element_text(colour="black", hjust=0.5)
  )+
  annotate("text", label = list("PAR2: manually", "adjusted coordinates!"), x = ChrY_max, y = log2(c(50, 30)), size = 4, colour = "black", hjust=1)+
  geom_text_repel(label = ifelse(Chr_all.ChrXY$Gene.name == "SRY", as.character(Chr_all.ChrXY$Gene.name), ""),
                   col = "black", size = 4, direction ="x", xlim=Chr_all.ChrXY$Gene.start..bp.[Chr_all.ChrXY$Gene.name == "SRY"][1],
                   box.padding = 1,
                   segment.alpha = 0.2)+ #line transparency)
  
  #xlim(c(NA,40000000))+
  #ylim(c(-5,15))+
  ggtitle(paste("RNA Seq Data, NES, Chr_all, ", PlotName, sep=""), subtitle = paste("AllSamples", "vs", RefOfSeqDat, ", threshold pVal: ", thresholdForPval, "; ", "PvalAdjusted", sep =""))

plot.ChrXY

ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       height = 18, width = 35, unit="cm", limitsize = FALSE)



