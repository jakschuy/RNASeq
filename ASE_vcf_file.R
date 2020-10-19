

library(ggplot2)
library("xlsx")
library(dplyr)
library(tidyr)
library(purrr)

#input
ProjectName <- "vcf"
thresholdForCounts <-  50
ListDesiredVariants <- c("missense_variant", "synonymous_variant", "protein_altering_variant", "stop_gained")
ListUnDesiredVariants <- c("3_prime", "5_prime", "UTR", "upstream", "downstream", "intergenic", "intron", "mature_miRNA", "miRNA", "non_coding")
firstDeletion <- c(43415001,44867000)
secondDeletion <- c(45781001,48110000)
surroundingsDeletion <- c(42415000,47110001)

colorsForLegend <- c("grey", "black", "dodgerblue", "black", "red", "green", "palegreen1", "palegreen4", "darkgoldenrod1", "darkorchid1")

#required functions
tophits <- function(x,n=5,decreasing = TRUE, na.last = TRUE){
  # Function tophits: goes through a numeric vector(!), sorts it and returns the tio entries with a list length of n,
  # can deal with NAs (na.last=F/T)
  # can deal with reverse order (decreasing = T/F)
  # v.1.0 (15.6.2020)
  
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


#create folder structure for loading data
ListFolder  <- list.dirs(path = "Z:/jaksch-sens2017106/vcf_files")[-1]

#create list with files
i <- 1
ListSamples <- vector()
ListSamples2 <- vector()
for(i in 1:length(ListFolder)){
  subdat <- list.files(path = ListFolder[i], pattern = ".annotated.GATKASE.vcf")
  subdat2 <- paste(ListFolder[i], "/", subdat, sep = "")
  ListSamples <- c(ListSamples, subdat2)
  ListSamples2 <- c(ListSamples2, subdat)
}

#import data, load vcf files)
i <- 1
Merge_Dat <- vector()
for(i in 1:length(ListSamples)){
  subdat <- read.table(file = ListSamples[i], sep = "\t", skip=55, comment.char="", header = TRUE)
  colnames(subdat)[10] <- "FORMAT_CONTENT"
  #subdat$INFO <- strsplit(as.character(subdat$INFO), ";", fixed=TRUE)
  #unnest(test2, cols = "INFO",keep_empty = TRUE, names_sep = "asdf")
  subdat$Sample <- ListSamples2[i]
  assign(ListSamples2[i], subdat)
  Merge_Dat <- rbind(Merge_Dat, subdat)
  
  
}


#extract information from allele specific column
Merge_Dat$FORMAT_CONTENT <- strsplit(as.character(Merge_Dat$FORMAT_CONTENT), ":", perl =TRUE)
Merge_Dat$GT <- sapply(Merge_Dat$FORMAT_CONTENT, function(x) x[1])
Merge_Dat$AD <- sapply(Merge_Dat$FORMAT_CONTENT, function(x) x[2])
Merge_Dat$DP <- sapply(Merge_Dat$FORMAT_CONTENT, function(x) as.numeric(x[3]))
Merge_Dat$GQ <- sapply(Merge_Dat$FORMAT_CONTENT, function(x) as.numeric(x[4]))
Merge_Dat$PL <- sapply(Merge_Dat$FORMAT_CONTENT, function(x) x[5])
Merge_Dat$BT <- sapply(Merge_Dat$FORMAT_CONTENT, function(x) x[6])

#clean environment and save progress
remove(subdat, subdat2, i)
save.image(file = paste(getwd(), "/R Stuff/ExtractedData_", ProjectName, "/saveEnvironment.RData", sep=""))


} #close load-Environment command

#load processed data
if (file.exists(paste(getwd(), "/R Stuff/ExtractedData_", ProjectName, "/saveEnvironment_processed.RData", sep=""))) {
  load(paste(getwd(), "/R Stuff/ExtractedData_", ProjectName, "/saveEnvironment_processed.RData", sep=""))} else {

#clean imported raw data
Merge_Clean <- Merge_Dat[, -c(3,7,9:10)]
Merge_Clean$AD_ref <- strsplit(Merge_Clean$AD, ",")
Merge_Clean$AD_alt <- as.numeric(sapply(Merge_Clean$AD_ref, function(x) x[2]))
Merge_Clean$AD_ref <- as.numeric(sapply(Merge_Clean$AD_ref, function(x) x[1]))

Merge_Clean$BT_altCount <- strsplit(Merge_Clean$BT, ",")
Merge_Clean$BT_totalCount <- as.numeric(sapply(Merge_Clean$BT_altCount, function(x) x[2]))
Merge_Clean$BT_Binominal <- as.numeric(sapply(Merge_Clean$BT_altCount, function(x) x[3]))
Merge_Clean$BT_pValNonparametric <- as.numeric(sapply(Merge_Clean$BT_altCount, function(x) x[4]))
Merge_Clean$BT_altCount <- as.numeric(sapply(Merge_Clean$BT_altCount, function(x) x[1]))

#Filter for non-performed Nonparamtetric test
Merge_Clean <- Merge_Clean[!Merge_Clean$BT_pValNonparametric %in% c(-1, NA),]

#filter for Counts, as quality check to make sure to have a certain rad depth
Merge_Clean <- Merge_Clean[Merge_Clean$BT_totalCount >= thresholdForCounts,]


# quick plotting to see the data
plot(Merge_Clean$DP, Merge_Clean$BT_totalCount)
plot(Merge_Clean$AD_alt, Merge_Clean$BT_altCount)


#design function to extract
readINFO <- function(x, match, separator, separatorV, num = TRUE, save = FALSE, saveDF){
  fetchMatch <- paste(match, separator, sep="")
  # v.1.0 (15.6.2020)
  
  #checks input
  if(class(match) != "character"){
    warning(" Match argument must be of type [character]! \n\n Stopping script now! \n")
    stop()
  } 
  if(class(separator) != "character"){
    warning(" Seperator argument must be of type [character]! \n\n Stopping script now! \n")
    stop()
  }
  if(class(x) != "list"){
    warning(" Input is not of type [list]! \n\n Trying to convert to a list \n", call.=FALSE)
    sub <- strsplit(as.character(x), separatorV)
    
  }else{
    sub <- x
  }

  #extraction
  findings <- lapply(sub, function(sub) strsplit(sub[startsWith(sub, match)], fetchMatch)) %>% unlist() %>% paste(sep="")
  
  #convert output to class numeric if assigned
  if(num == TRUE){
    findings <- as.numeric(findings[findings != ""])
    
  }else{
    findings <- findings[findings != ""]
  }

  #quality control
  if(length(findings) != length(sub)){
    warning(" Extracted data does not match working dataframe! \n Check empty values and length(findings) vs. nrow(df) \n\n Stopping script now! \n")
    stop()
  }
  
  #output
  if(save == TRUE){
    #check for destination df
    if(!exists(as.character(substitute(saveDF)))){
      warning(" Saving option==TRUE but no/invalid destination chosen! \n\n Stopping script now! \n", call.=FALSE)
      stop()
    }
    
    #save findings
    savesub <- cbind(saveDF, findings)
    colnames(savesub)[ncol(savesub)] <- paste("Extr_", match, sep = "")
    assign("saveDF", savesub)
    return(savesub)
    
  }else{
    
    return(findings)
  }
  
  ##debug section
  if(debug == TRUE){
  #define variables
  x <- Merge_Clean$INFO
  match <- "CSQ"
  separator <- "="
  separatorV <- ";"
  save <- TRUE
  saveDF <- Merge_Clean
  
  #clean environment
  remove(x, match, separator, separatorV, save, saveDF, savesub, sub, findings, fetchMatch)
  }
}



#use function to extract information
Merge_Clean <- readINFO(x = Merge_Clean$INFO, match = "AC", separator = "=", separatorV = ";",  save = TRUE, saveDF = Merge_Clean)
Merge_Clean <- readINFO(x = Merge_Clean$INFO, match = "AF", separator = "=", separatorV = ";", save = TRUE, saveDF = Merge_Clean)
Merge_Clean <- readINFO(x = Merge_Clean$INFO, match = "AN", separator = "=", separatorV = ";", save = TRUE, saveDF = Merge_Clean)
Merge_Clean <- readINFO(x = Merge_Clean$INFO, match = "CSQ", separator = "=", separatorV = ";", save = TRUE, saveDF = Merge_Clean, num = FALSE)
#Merge_Clean <- readINFO(x = Merge_Clean$Extr_CSQ, match = "CSQ", separator = "=", separatorV = ";", save = TRUE, saveDF = Merge_Clean, num = FALSE)

#get gene type, gene name, ENSG ID
Merge_Clean2 <- Merge_Clean
Merge_Clean2$Variant_type <- strsplit(as.character(Merge_Clean2$Extr_CSQ), ",") 
Merge_Clean2 <- unnest(data = Merge_Clean2[, -c(6,13,23)], cols = "Variant_type")
Merge_Clean2$Variant_type <- strsplit(Merge_Clean2$Variant_type, "|", fixed = TRUE)
Merge_Clean2$Gene.name <- sapply(Merge_Clean2$Variant_type, function (x) x[4])
Merge_Clean2$Gene.stable.ID <- sapply(Merge_Clean2$Variant_type, function (x) x[5])
Merge_Clean2$Variant_type <- sapply(Merge_Clean2$Variant_type, function (x) x[2])

#add location column
Merge_Plot <- Merge_Clean2
Merge_Plot$Position <- "Unaffected"
Merge_Plot$Position[Merge_Plot$POS > surroundingsDeletion[1] & Merge_Plot$POS < surroundingsDeletion[2] & Merge_Plot$X.CHROM == "chr21"] <- "surroundings_1MB"
Merge_Plot$Position[Merge_Plot$POS > firstDeletion[1] & Merge_Plot$POS < firstDeletion[2] & Merge_Plot$X.CHROM == "chr21"]  <- "First Deletion"
Merge_Plot$Position[Merge_Plot$POS > secondDeletion[1] & Merge_Plot$POS < secondDeletion[2] & Merge_Plot$X.CHROM == "chr21"]  <- "Second Deletion"
Merge_Plot$Position <- as.factor(Merge_Plot$Position)
levels(Merge_Plot$Position) <- levels(Merge_Plot$Position)[1:4]

#filter for variants
SearchParameter <- paste(ListUnDesiredVariants, collapse ="|")
SearchOutput <- grep(pattern = SearchParameter, x = levels(factor(Merge_Plot$Variant_type)), value = TRUE)
Merge_Plot <- Merge_Plot[!Merge_Plot$Variant_type %in% SearchOutput,]

#save environment for processed data
save.image(file = paste(getwd(), "/R Stuff/ExtractedData_", ProjectName, "/saveEnvironment_processed.RData", sep=""))
}


##Chr21_allele specific expression----
#filter for chromosomes
Merge_Plot21 <- Merge_Plot[Merge_Plot$X.CHROM == "chr21",]

PlotName <- "Allele specific expression_Chr21"


skipThisPart <- TRUE
if (skipThisPart != TRUE) {

Plot21 <- ggplot(Merge_Plot21, aes(x=POS, y=QUAL, col = Sample))+
  geom_point()
Plot21
# Plot21 <- ggplot(Merge_Plot21, aes(x=POS, y=QUAL, col = Position))+
#   geom_point()+
#   facet_grid(~Sample)
# Plot21
Plot21 <- ggplot(Merge_Plot21, aes(x=POS, y=QUAL, col = Position))+
  geom_point()+
  scale_color_manual(values=colorsForLegend[c(5,6,3,2)])+
  facet_grid(~Sample)
Plot21

#shape for heterocygocity
Plot21 <- ggplot(Merge_Plot21, aes(x=POS, y=QUAL, col = Position, shape = GT))+
  geom_point()+
  scale_color_manual(values=colorsForLegend[c(5,6,3,2)])+
  facet_grid(~Sample)
Plot21

#only ALI09
Plot21 <- ggplot(Merge_Plot21[Merge_Plot21$Sample == "ALI9.annotated.GATKASE.vcf",], aes(x=POS, y=QUAL, col = Position, shape = GT))+
  geom_point()+
  scale_color_manual(values=colorsForLegend[c(5,6,3,2)])+
  facet_grid(~Sample)
Plot21

}

##AD_BT_correlation coefficient----
#calculate allele ratio
Merge_Plot21$BT_ratio <- Merge_Plot21$BT_altCount / Merge_Plot21$BT_totalCount      #BT comes from ASERead counter
Merge_Plot21$AD_ratio <- Merge_Plot21$AD_alt / Merge_Plot21$DP        #AD comes from GATK haplotype caller

#statistics
options(scipen = 0)
shapiro.test(Merge_Plot21$BT_ratio) #data is not normally distibuted
shapiro.test(Merge_Plot21$AD_ratio) #data is not normally distibuted
Corr.CoefSpear <- cor.test(Merge_Plot21$BT_ratio, Merge_Plot21$AD_ratio, method = "spear")[["estimate"]]
Corr.Coef_pSpear <- cor.test(Merge_Plot21$BT_ratio, Merge_Plot21$AD_ratio, method = "spear")[["p.value"]]    #data is not noramlly distributed therefore spearman
Corr.CoefPear <- cor.test(Merge_Plot21$BT_ratio, Merge_Plot21$AD_ratio, method = "pearson")[["estimate"]]     #data is not normally distributed therefore spearman
Corr.Coef_pPear <- cor.test(Merge_Plot21$BT_ratio, Merge_Plot21$AD_ratio, method = "pearson")[["p.value"]]


#correlation between both ratios
Plot21_correlation <- ggplot(Merge_Plot21, aes(x = BT_ratio, y = AD_ratio))+
  geom_point(size = 8, shape = 1, stroke = 1)+
  ylim(c(0,1))+
  theme_classic()+
  geom_text(x=0, y=1,   label=paste("rho(Pearson):   ", round(Corr.CoefPear, 5), ", p= ", Corr.Coef_pPear, sep=""), hjust = "left", color="black", show.legend = FALSE)+
  geom_text(x=0, y=0.98, label=paste("rho(Spearman): ", round(Corr.CoefSpear, 5), ", p= ", Corr.Coef_pSpear, sep=""), hjust = "left", color="black", show.legend = FALSE)+
  theme_classic(base_size=20)
Plot21_correlation
options(scipen = 999)



#manually convert deletes regeion into heterocygous
Merge_Plot21$GT[Merge_Plot21$Position %in% c("First Deletion", "Second Deletion")] <- "0/1"

#filter for heterocygocity, we are interested in the allele specific expression
Merge_Plot21 <- Merge_Plot21[Merge_Plot21$GT == "0/1",]

#fiter against indels, reason: noise, for allele specific expression we need high quality SNVs
Merge_Plot21 <- Merge_Plot21[nchar(as.character(Merge_Plot21$REF)) == 1,]

#add colour code for ratio grouping
thresholdForRatioGrouping <- 0.5
rangeRatioGrouping <- 0.2
Merge_Plot21$GroupRatio <- "normal"
Merge_Plot21$GroupRatio[Merge_Plot21$BT_ratio > thresholdForRatioGrouping + rangeRatioGrouping | Merge_Plot21$AD_ratio > thresholdForRatioGrouping + rangeRatioGrouping] <- "changed"
Merge_Plot21$GroupRatio[Merge_Plot21$BT_ratio < thresholdForRatioGrouping - rangeRatioGrouping | Merge_Plot21$AD_ratio < thresholdForRatioGrouping - rangeRatioGrouping] <- "changed"

#all patients
df_GroupRatio <- Merge_Plot21 %>% 
  group_by(GroupRatio, add = TRUE) %>%
  group_by(Sample, add = TRUE) %>%
  summarise(counted = length(GroupRatio))
df_GroupRatio

#add better x axis for position, convert to KB
Merge_Plot21$POS_10E3 <- Merge_Plot21$POS / 10^3
Merge_Plot21$POS_10E6 <- Merge_Plot21$POS / 10^6


#plot with BT
Plot21 <- ggplot(Merge_Plot21[Merge_Plot21$Sample == "ALI9.annotated.GATKASE.vcf",], aes(x=POS_10E3, y=BT_ratio, col = Position, shape = GT))+
  theme_classic()+
  annotate(geom = "rect", xmin = firstDeletion[1]/ 10^3, xmax = firstDeletion[2]/ 10^3, ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[5])+
  annotate(geom = "rect", xmin = secondDeletion[1]/ 10^3, xmax = secondDeletion[2]/ 10^3, ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[6])+
  annotate(geom = "rect", xmin = surroundingsDeletion[1]/ 10^3, xmax = firstDeletion[1]/ 10^3, ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[3])+
  annotate(geom = "rect", xmin = firstDeletion[2]/ 10^3, xmax = secondDeletion[1]/ 10^3, ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[3])+
  geom_hline(yintercept = 0.5, linetype =2, alpha = 0.25, size=0.3)+
  geom_hline(yintercept = 0.75, linetype =2, alpha = 0.25, size=0.3)+
  geom_hline(yintercept = 0.25, linetype =2, alpha = 0.25, size=0.3)+
  geom_point(size = 5, alpha = 0.8)+
  scale_color_manual(values=colorsForLegend[c(5,6,3,2)])+
  #facet_grid(~Sample)+
  ylim(0,1)+
  labs(x="Chromosome 21 [kb]", y="BT: nVariant / ntotalVariant")+
  ggtitle(paste("VCF (from RNA Seq Data) - ", PlotName, sep=""), subtitle = paste( "filtered for minReadCounts (", thresholdForCounts, "), ", "only heterocygous, deleted region manually described->0/1, Indels excluded",  sep =""))
Plot21
#save
ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       height = 15, width = 40, unit="cm", limitsize = FALSE)


#plot with AD
Plot21 <- ggplot(Merge_Plot21[Merge_Plot21$Sample == "ALI9.annotated.GATKASE.vcf",], aes(x=POS_10E3, y=AD_ratio, col = Position, shape = GT))+
  theme_classic()+
  annotate(geom = "rect", xmin = firstDeletion[1]/ 10^3, xmax = firstDeletion[2]/ 10^3, ymin = 1, ymax = 0, alpha=0.2, fill = colorsForLegend[5])+
  annotate(geom = "rect", xmin = secondDeletion[1]/ 10^3, xmax = secondDeletion[2]/ 10^3, ymin = 1, ymax = 0, alpha=0.2, fill = colorsForLegend[6])+
  annotate(geom = "rect", xmin = surroundingsDeletion[1]/ 10^3, xmax = firstDeletion[1]/ 10^3, ymin = 1, ymax = 0, alpha=0.2, fill = colorsForLegend[3])+
  annotate(geom = "rect", xmin = firstDeletion[2]/ 10^3, xmax = secondDeletion[1]/ 10^3, ymin = 1, ymax = 0, alpha=0.2, fill = colorsForLegend[3])+
  geom_hline(yintercept = 0.5, linetype =2, alpha = 0.25, size=0.3)+
  geom_hline(yintercept = 0.75, linetype =2, alpha = 0.25, size=0.3)+
  geom_hline(yintercept = 0.25, linetype =2, alpha = 0.25, size=0.3)+
  geom_point(size = 4, alpha = 0.7)+
  scale_color_manual(values=colorsForLegend[c(5,6,3,2)])+
  scale_x_continuous(breaks = seq(15000, 45000, by = 5000), labels = c("", "20000", "", "30000", "", "40000", ""))+
  #facet_grid(~Sample)+
  ylim(0,1)+
  labs(x="Chromosome 21 [kb]", y="AD: nVariant / ntotalVariant")+
  theme(
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(colour="black"),
    axis.text.x = element_text(colour="black", hjust=0.5),
    axis.text.y = element_text(colour="black")
  )+
  ggtitle(paste("VCF (from RNA Seq Data) - ", PlotName, sep=""), subtitle = paste( "filtered for minReadCounts (", thresholdForCounts, "), ", "only heterocygous, deleted region manually described->0/1, Indels excluded",  sep =""))
Plot21
#save
ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       height = 7, width = 15, unit="cm", limitsize = FALSE)


##only Deleletions
#plot with BT
Merge_Plot21Del <- Merge_Plot21[Merge_Plot21$POS > 40000000,]
Plot21 <- ggplot(Merge_Plot21Del[Merge_Plot21Del$Sample == "ALI9.annotated.GATKASE.vcf",], aes(x=POS, y=BT_ratio, col = Position, shape = GT))+
  theme_classic()+
  annotate(geom = "rect", xmin = firstDeletion[1], xmax = firstDeletion[2], ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[5])+
  annotate(geom = "rect", xmin = secondDeletion[1], xmax = secondDeletion[2], ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[6])+
  annotate(geom = "rect", xmin = surroundingsDeletion[1], xmax = firstDeletion[1], ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[3])+
  annotate(geom = "rect", xmin = firstDeletion[2], xmax = secondDeletion[1], ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[3])+
  geom_point()+
  scale_color_manual(values=colorsForLegend[c(5,6,3,2)])+
  facet_grid(~Sample)+
  ylim(0,1)+
  labs(x="Chromosome 21", y="BT: nVariant / ntotalVariant")+
  ggtitle(paste("VCF (from RNA Seq Data) - ", PlotName, "_Del", sep=""), subtitle = paste( "filtered for minReadCounts (", thresholdForCounts, "), ", "only heterocygous, deleted region manually described->0/1, Indels excluded",  sep =""))
Plot21

#plot with AD
Plot21 <- ggplot(Merge_Plot21Del[Merge_Plot21Del$Sample == "ALI9.annotated.GATKASE.vcf",], aes(x=POS, y=AD_ratio, col = Position, shape = GT))+
  theme_classic()+
  annotate(geom = "rect", xmin = firstDeletion[1], xmax = firstDeletion[2], ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[5])+
  annotate(geom = "rect", xmin = secondDeletion[1], xmax = secondDeletion[2], ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[6])+
  annotate(geom = "rect", xmin = surroundingsDeletion[1], xmax = firstDeletion[1], ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[3])+
  annotate(geom = "rect", xmin = firstDeletion[2], xmax = secondDeletion[1], ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[3])+
  geom_point()+
  scale_color_manual(values=colorsForLegend[c(5,6,3,2)])+
  facet_grid(~Sample)+
  ylim(0,1)+
  labs(x="Chromosome 21", y="AD: nVariant / ntotalVariant")+
  ggtitle(paste("VCF (from RNA Seq Data) - ", PlotName, "_Del", sep=""), subtitle = paste( "filtered for minReadCounts (", thresholdForCounts, "), ", "only heterocygous, deleted region manually described->0/1, Indels excluded",  sep =""))
Plot21

#save
ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       height = 10, width = 30, unit="cm", limitsize = FALSE)

ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/temp.png", sep=""),
       height = 10, width = 30, unit="cm", limitsize = FALSE)





##focus on genomewide
Merge_Plot_genomewide <- Merge_Plot
levels(Merge_Plot_genomewide$X.CHROM) <- levels(Merge_Plot_genomewide$X.CHROM)[c(1,12, 16:22, 2:11, 13:15,23,24)]



PlotName <- "Genomewide"

#----genomewide_Variant counting----
#get list with variant types and respective abundance
#all patients
df_Variant_all <- Merge_Plot_genomewide %>% 
  group_by(Variant_type, add = TRUE) %>%
  group_by(X.CHROM, add = TRUE) %>%
  summarise(counted = length(Variant_type))

#plot for variant counting
Plot_df_Variant <- ggplot(df_Variant_all, aes(x=X.CHROM, y=Variant_type))+
  geom_point(aes(size=log10(df_Variant_all$counted)), alpha = .8)+
  ggtitle(paste("VCF (from RNA Seq Data) - ", PlotName, "_nCounts_allPatients", sep=""), subtitle = paste( "filtered for minReadCounts (", thresholdForCounts, ")",  sep =""))
Plot_df_Variant

ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_variants-counted-allPatients", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       height = 15, width = 40, unit="cm", limitsize = FALSE)

#only ALI09
df_Variant_ALI09 <- Merge_Plot_genomewide[Merge_Plot_genomewide$Sample == "ALI9.annotated.GATKASE.vcf",] %>% 
  group_by(Variant_type, add = TRUE) %>%
  group_by(X.CHROM, add = TRUE) %>%
  summarise(counted = length(Variant_type))

#plot ALI09
Plot_df_Variant <- ggplot(df_Variant_ALI09, aes(x=X.CHROM, y=Variant_type))+
  geom_point(aes(size=log10(df_Variant_ALI09$counted)), alpha = .8)+
  ggtitle(paste("VCF (from RNA Seq Data) - ", PlotName, "_nCounts_ALI09", sep=""), subtitle = paste( "filtered for minReadCounts (", thresholdForCounts, ")",  sep =""))
Plot_df_Variant
ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_variants-counted-ALI09", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       height = 15, width = 40, unit="cm", limitsize = FALSE)


#----genomewide_ASE----
PlotName <- "Allele-specific"
options(scipen=2)


Merge_Plot_genomewide <- Merge_Plot
Merge_Plot_genomewide$X.CHROM <- factor(Merge_Plot_genomewide$X.CHROM, levels = levels(factor(Merge_Plot_genomewide$X.CHROM))[c(1,12, 16:22, 2:11, 13:15,23,24)])


##calculations
#add Variant ratio
Merge_Plot_genomewide$BT_ratio <- Merge_Plot_genomewide$BT_altCount / Merge_Plot_genomewide$BT_totalCount
Merge_Plot_genomewide$AD_ratio <- Merge_Plot_genomewide$AD_alt / Merge_Plot_genomewide$DP 

#add colour code for ratio grouping
thresholdForRatioGrouping <- 0.5
rangeRatioGrouping <- 0.1
Merge_Plot_genomewide$GroupRatio <- "normal"
Merge_Plot_genomewide$GroupRatio[Merge_Plot_genomewide$BT_ratio > thresholdForRatioGrouping + rangeRatioGrouping | Merge_Plot_genomewide$AD_ratio > thresholdForRatioGrouping + rangeRatioGrouping] <- "changed"
Merge_Plot_genomewide$GroupRatio[Merge_Plot_genomewide$BT_ratio < thresholdForRatioGrouping - rangeRatioGrouping | Merge_Plot_genomewide$AD_ratio < thresholdForRatioGrouping - rangeRatioGrouping] <- "changed"

#all patients
df_GroupRatio <- Merge_Plot_genomewide %>% 
  group_by(GroupRatio, add = TRUE) %>%
  group_by(X.CHROM, add = TRUE) %>%
  summarise(counted = length(GroupRatio))

#add better x axis for position
Merge_Plot_genomewide$POS_10E6 <- Merge_Plot_genomewide$POS / 10^6

##Quality control
#filter for heterocygocity, we are interested in the allele specific expression, false assigned hits will be lost
Merge_Plot_genomewide <- Merge_Plot_genomewide[Merge_Plot_genomewide$GT == "0/1",]

#fiter against indels, reason: noise, for allele specific expression we need high quality SNVs
Merge_Plot_genomewide <- Merge_Plot_genomewide[nchar(as.character(Merge_Plot_genomewide$REF)) == 1,]

#Extract List of Genes
write.xlsx(Merge_Plot_genomewide[Merge_Plot_genomewide$GroupRatio == "changed", c("X.CHROM", "POS", "Gene.stable.ID", "Gene.name", "Sample", "REF", "ALT", "BT_ratio", "AD_ratio")],
           paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/", "ListGeneWithDifExprAllele_all", ".xls", sep=""))
extractGenes        <- Merge_Plot_genomewide[Merge_Plot_genomewide$GroupRatio == "changed", c("X.CHROM", "POS", "Gene.stable.ID", "Gene.name", "Sample", "REF", "ALT", "BT_ratio", "AD_ratio")]
extractGenes_unique <- extractGenes$Gene.name[!duplicated(extractGenes$Gene.name)]

#Plot
i <- 1
for (i in 1:length(levels(factor(Merge_Plot_genomewide$Sample)))) {
  
Plot21 <- ggplot(Merge_Plot_genomewide[Merge_Plot_genomewide$Sample == levels(factor(Merge_Plot_genomewide$Sample))[i],], aes(x=POS_10E6, y=AD_ratio, col = GroupRatio))+
  geom_point(show.legend = TRUE, size = 1)+
  geom_smooth(method = "loess", show.legend = FALSE, col = colorsForLegend[3], se = TRUE, span = 0.37, size = 0.5, na.rm =TRUE, inherit.aes = TRUE, fill = "grey50")+
  #geom_text(x=0, y=1,   label=paste("rho(Pearson):   ", round(Corr.CoefPear, 5), ", p= ", Corr.Coef_pPear, sep=""), hjust = "left", color="black", show.legend = FALSE)+
  scale_color_manual(values=colorsForLegend[c(5,2,6)])+
  facet_wrap(~X.CHROM, scales = "free_x", ncol = 6)+
  #scale_x_continuous(labels = scales::comma, limits = c(0, NA))+
  scale_x_continuous(labels = scales::comma)+
  ylim(0,1)+
  labs(x="Position [bp x10^6]", y="AD: nVariant / ntotalVariant")+
  theme(
    plot.subtitle = element_text(colour="black",face="italic"),
    axis.ticks = element_line(colour="black"),
    axis.text = element_text(colour="black"),
    legend.position = c(0.9, 0.1),
    legend.background = element_rect(linetype="solid", colour ="black"))+
    ggtitle(paste("VCF (from RNA Seq Data) - ", PlotName, "_Genome_", levels(factor(Merge_Plot_genomewide$Sample))[i], sep=""), subtitle = paste( "filtered for minReadCounts (", thresholdForCounts, "), indels excluded, homocygous variants excluded, RatioRange: ", rangeRatioGrouping*100, "%",  sep =""))

#save
ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_", levels(factor(Merge_Plot_genomewide$Sample))[i], "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       height = 18, width = 30, unit="cm", limitsize = FALSE)

}









