Sys.setenv(LANG="en")
rm(list=ls())
options(scipen=999)

library(Hapi)
library(ggplot2)
library(writexl)
library(dplyr)
library(tidyr)
library(purrr)
library(car)
library(pastecs)
library(gmodels)



#input
ProjectName <- "vcf_trio"
thresholdForCounts <-  50
thresholdForQual <- 40
colorsForLegend <- c("grey", "black", "dodgerblue", "black", "red", "green", "palegreen1", "palegreen4", "darkgoldenrod1", "darkorchid1")
SkipPlotSection <- TRUE
SkipExportSection <- TRUE
firstDeletion <- c(1,8204967)
secondDeletion <- c(41994799, 43447106)
thirdDeletion <- c(44361567, 46699872)
Inversion <- c(43447107, 44361566)
surroundingsDeletion <- c(40994799, 47690088)
surroundingsDeletionP <- c(-1000000,9204967)



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


#design function to extract
readINFO <- function(x, match, separator, separatorV, num = FALSE, save = FALSE, saveDF){
  # v.1.1 (15.6.2020)
  
  #build search mask
  fetchMatch <- paste(match, separator, sep="")
  
  #checks input on data format
  if(class(match) != "character"){
    warning(" Match argument must be of type [character]! \n\n Stopping script now! \n")
    stop()
  }
  if(class(separator) != "character"){
    warning(" Seperator argument must be of type [character]! \n\n Stopping script now! \n")
    stop()
  }
  if(class(x) != "list"){
    message(" INFO: Input is not of type [list]!  Trying to convert to a list")
    sub <- strsplit(as.character(x), separatorV)
    
  }else{
    sub <- x
  }
  
  #checks whether match is part of data
  if(all(!grepl(fetchMatch, sub))){
    stop(" Cannot find search mask in data! (\"", fetchMatch, "\") \n Stopping script now! \n")
  }
  
  #extraction
  findings1 <- lapply(sub, function(sub) sub[startsWith(sub, fetchMatch)])  #lapply or sapply, doesnt matter
  findings1[findings1 == "character(0)"] <- paste(match, NA, sep=separator)
  #findings1[[3]] <- paste(findings1[[3]], separator, "bullshit", sep="")    #for debug
  findings2 <- sapply(findings1, function(x) strsplit(x[startsWith(x, match)], separator))
  findings3 <- vapply(findings2, function(x) paste(x[-1], collapse = "; "), FUN.VALUE = character(1))
  
  
  #convert output to class numeric if assigned
  if(num == TRUE){
    #
    
    #checks if findings3 has numbers/punctuations, without any letters/,/: but with NAs
    if(all(grepl("[[:digit:][:punct:]]", findings3) & !grepl("[[:alpha:]]", findings3) & !grepl(",", findings3) & !grepl(":", findings3)| grepl("NA", findings3))){
      suppressWarnings(
        findings3 <- as.numeric(findings3)
      )
      message(" INFO: Conversion into numeric: PASS!")
    } else if(any(grepl("[[:alpha:]]", findings3))){
      message(" INFO: Cannot convert output into number -> output contains letters. >Output remains as type character")
    } else if(any(grepl(",", findings3) | grepl(":", findings3))){
      message(" INFO: Cannot convert output into number -> output contains  \",\" and/or \":\" >Output remains as type character")
    } else {
      message(" INFO: Unknown error in output conversion into numeric. Please run debug section. >Output remains as type character")
    }
  }
  
  #Checking length, then Quality check
  if(length(findings3) != length(sub)){
    findings <- findings3[findings3 != ""]
    message(" INFO: Data length does not match working dataframe! Trying to remove empty elements")
  }
  
  #quality control for data length
  if(length(findings3) == length(sub)){
    message(" INFO: Quality Control for Data length: PASS!")
  } else {
    
    warning(" Extracted data does not match working dataframe! \n Check empty values and length(findings) vs. nrow(df) \n\n Stopping script now! \n", )
    stop()
  }
  
  #output section
  if(save == TRUE){
    #check for destination df
    if(!exists(as.character(substitute(saveDF)))){
      warning(" Saving option==TRUE but no/invalid destination chosen! \n\n Stopping script now! \n", call.=FALSE)
      stop()
    }
    
    #save findings
    savesub <- cbind(saveDF, findings3)
    colnames(savesub)[ncol(savesub)] <- paste("Extr_", match, sep = "")
    assign("saveDF", savesub)
    return(savesub)
    
  } else {
    
    return(findings3)
  }
  
  
  ##debug section
  #define variables
  # x <- Merge_Clean$INFO
  # match <- "RankScore"
  # num <- FALSE
  # separator <- "="
  # separatorV <- ";"
  # save <- FALSE
  # saveDF <- Merge_Clean
  
  #clean environment
  remove(x, match, separator, separatorV, save, saveDF, savesub, sub, findings, findings1, findings2, findings3, fetchMatch)
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
if (file.exists(paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/saveEnvironment.RData", sep=""))) {
  load(paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/saveEnvironment.RData", sep=""))} else {
    
    #load coordintates for hg19
    data("hg19")
    
    #load coordintates for hg38
    hg38_1 <- read.delim("Sequencing data/Modeled_regions_for_GRCh38_from ncbi.nlm.nih.govgrchuman.tsv", sep = "\t", skip=0, comment.char="", header = TRUE)
    hg38_1 <- hg38_1[!duplicated(hg38_1$X.region_name),]
    colnames(hg38_1) <- c(colnames(hg38_1[-1]), "cen_length")
    
    hg38_2 <- read.delim("Sequencing data/GCA_000001405.15_GRCh38_assembly_report.txt", sep = "\t", skip=0, comment.char="#", header = TRUE)
    hg38_2 <- hg38_2[1:24, c(9:10)]
    
    hg38 <- data.frame(hg38_1, hg38_2)
    rownames(hg38) <- hg38$UCSC.style.name
    hg38 <- hg38[, c(5,2,3)]
    colnames(hg38) <- colnames(hg19)
    
    hg38[24,c(2)] <- hg38[24,c(2)] - 750000 #manual extension of centromer of Y for visualization
    hg38[24,c(3)] <- hg38[24,c(3)] + 750000
    
    
    #import data, load vcf files, chr21)
    vcfsa21 <- read.table(file = "//sshfs/jaksch-sens2017106@bianca-sftp.uppmax.uu.se/jaksch-sens2017106/vcf_files_hg38/phased_possorted_ALI9.LReads.hg38.REF_chr21.bam_marked_duplicated.ReadGroups.bam.vcf", sep = "\t", skip=0, comment.char="#", header = TRUE)
    vcfmo21 <- read.table(file = "//sshfs/jaksch-sens2017106@bianca-sftp.uppmax.uu.se/jaksch-sens2017106/vcf_files_hg38/phased_possorted_mo.REF_chr21.bam_marked_duplicated.ReadGroups.bam.vcf", sep = "\t", skip=0, comment.char="#", header = TRUE)
    vcffa21 <- read.table(file = "//sshfs/jaksch-sens2017106@bianca-sftp.uppmax.uu.se/jaksch-sens2017106/vcf_files_hg38/phased_possorted_fa.REF_chr21.bam_marked_duplicated.ReadGroups.bam.vcf", sep = "\t", skip=0, comment.char="#", header = TRUE)
    vcfsa21PCRf <- read.table(file = "//sshfs/jaksch-sens2017106@bianca-sftp.uppmax.uu.se/jaksch-sens2017106/vcf_files_hg38/P2109_188.ALI9.PCRfree.hg38.REF_chr21.bam_marked_duplicated.ReadGroups.bam.vcf", sep = "\t", skip=0, comment.char="#", header = TRUE)

    #rename columns
    vcfsa21$name <- "ALI09sa_Lr"
    colnames(vcfsa21)[10] <- "FORMAT_CONTENT"
    vcfmo21$name <- "ALI09mo_Lr"
    colnames(vcfmo21)[10] <- "FORMAT_CONTENT"
    vcffa21$name <- "ALI09fa_Lr"
    colnames(vcffa21)[10] <- "FORMAT_CONTENT"
    vcfsa21PCRf$name <- "ALI09sa_PCRf"
    colnames(vcfsa21PCRf)[10] <- "FORMAT_CONTENT"

    # #unite df
    # Merge_Dat21 <- rbind(vcfsa21, vcfmo21, vcffa21)
    
    #info about vcf file
    {
      ##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
      ##FILTER=<ID=LowQual,Description="Low quality">
      ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
      ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
      ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
      ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
      ##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
      ##INFO=<ID=RO,Number=1,Type=Integer,Description="Reference allele observed count">
      ##INFO=<ID=AO,Number=.,Type=Integer,Description="Alternate allele observed count">
      ##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
      ##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
      ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
      ##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
      ##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
      ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
      ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
      ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
      ##INFO=<ID=AS_InbreedingCoeff,Number=A,Type=Float,Description="allele specific heterozygosity as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation; relate to inbreeding coefficient">
      ##INFO=<ID=AS_QD,Number=1,Type=Float,Description="Allele-specific Variant Confidence/Quality by Depth">
      ##INFO=<ID=AS_RAW_BaseQRankSum,Number=1,Type=String,Description="raw data for allele specific rank sum test of base qualities">
      ##INFO=<ID=AS_RAW_MQ,Number=A,Type=Float,Description="Allele-specfic raw data for RMS Mapping Quality">
      ##INFO=<ID=AS_RAW_MQRankSum,Number=1,Type=String,Description="Allele-specific raw data for Mapping Quality Rank Sum">
      ##INFO=<ID=AS_RAW_ReadPosRankSum,Number=1,Type=String,Description="allele specific raw data for rank sum test of read position bias">
      ##INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description="Allele-specific forward/reverse read counts for strand bias tests">
      ##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
      ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
      ##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
      ##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
      ##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
      ##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
      ##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
      ##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
      ##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
      ##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
      ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
      ##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
      ##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
      ##INFO=<ID=RAW_MQ,Number=1,Type=Float,Description="Raw data for RMS Mapping Quality">
      ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
      ##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
    }
    
    # https://cran.r-project.org/web/packages/Hapi/vignettes/Hapi.html
    
    #load SNV data from parents
    PaFa_noMo <-  read.table(file = "//sshfs/jaksch-sens2017106@bianca-sftp.uppmax.uu.se/jaksch-sens2017106/vcf_files_hg38/SaFa_noMo.vcf", sep = "\t", skip=0, comment.char="", header = FALSE)
    PaMo_noFa <-  read.table(file = "//sshfs/jaksch-sens2017106@bianca-sftp.uppmax.uu.se/jaksch-sens2017106/vcf_files_hg38/SaMo_noFa.vcf", sep = "\t", skip=0, comment.char="", header = FALSE)
    FaMo_noPa <-  read.table(file = "//sshfs/jaksch-sens2017106@bianca-sftp.uppmax.uu.se/jaksch-sens2017106/vcf_files_hg38/FaMo_noPa.vcf", sep = "\t", skip=0, comment.char="", header = FALSE)
    
    #clean environment and save progress
    rm(hg38_1, hg38_2)
    save.image(file = paste(getwd(), "/R Stuff/ExtractedData_", ProjectName, "/saveEnvironment.RData", sep=""))
    
    
  } #close load-Environment command



#----analyse SNVs in vicinity of breakpoints(as part of start sequence)
#prepare table for Haplo layout
PlotName <- "breakpointAnalysis"

#build backup
bpjchr21 <- vcfsa21

#focus on junctions
bpjchr21$Position <- "Unaffected"
# bpjchr21$Position[bpjchr21$CoordStart_MB > surroundingsDeletion[1] & bpjchr21$CoordStart_MB < surroundingsDeletion[2] & bpjchr21$CHROM  == "chr21"] <- "surroundings_1MB"
# bpjchr21$Position[bpjchr21$CoordStart_MB > surroundingsDeletionP[1] & bpjchr21$CoordStart_MB < surroundingsDeletionP[2] & bpjchr21$CHROM  == "chr21"] <- "surroundings_1MB"
bpjchr21$Position[bpjchr21$POS > Inversion[1] & bpjchr21$POS < Inversion[2] & bpjchr21$CHROM == "chr21"]  <- "Inversion"
bpjchr21$Position[bpjchr21$POS >= hg38$cenStart[21] & bpjchr21$POS <= hg38$cenEnd[21]] <- "Centromer"
bpjchr21$Position[bpjchr21$POS > firstDeletion[1] & bpjchr21$POS < firstDeletion[2] & bpjchr21$CHROM  == "chr21"]  <- "First Deletion"
bpjchr21$Position[bpjchr21$POS > secondDeletion[1] & bpjchr21$POS < secondDeletion[2] & bpjchr21$CHROM  == "chr21"]  <- "Second Deletion"
bpjchr21$Position[bpjchr21$POS > thirdDeletion[1] & bpjchr21$POS < thirdDeletion[2] & bpjchr21$CHROM  == "chr21"]  <- "Third Deletion"

#reorder Position
bpjchr21$Position <- factor(bpjchr21$Position, levels(factor(bpjchr21$Position))[c(2,1,4,5,3,6)])

#add conversion of bases into kilo bases
conversionFactorKB <- 1000
bpjchr21$CoordStart_KB <- bpjchr21$POS / conversionFactorKB
conversionFactorMB <- 1000000
bpjchr21$CoordStart_MB <- bpjchr21$POS / conversionFactorMB

#quick plot
if (SkipPlotSection == FALSE) {
  plot_position <- ggplot(bpjchr21,aes(x=CoordStart_MB, y=log2(QUAL), col = Position)) +
    geom_point()+
    # scale_x_continuous(limits = c(Deletion_chr21p[1], secondDeletion[2]), breaks = seq(1, 50000000, by=5000000))+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_positions", subtitle = "noQualFilter")
  
  plot_position
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_patient_Qual-positions", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
}



##extract information from allele specific column
i <- 1
bpjchr21_FORMAT <- data.frame(matrix(nrow = 0, ncol = ncol(bpjchr21)))
colnames(bpjchr21_FORMAT) <- colnames(bpjchr21)

for (i in 1:length(levels(factor(bpjchr21$FORMAT)))) {
  #subset for dif levels
  bpjchr21_sub <- bpjchr21[bpjchr21$FORMAT == levels(factor(bpjchr21$FORMAT))[i],]
  #split info
  bpjchr21_sub$FORMAT <- strsplit(as.character(bpjchr21_sub$FORMAT), ":", perl =TRUE)
  bpjchr21_sub$FORMAT_CONTENT <- strsplit(as.character(bpjchr21_sub$FORMAT_CONTENT), ":", perl =TRUE)
  #read no of elements
  numberElements <- length(bpjchr21_sub$FORMAT[[1]])
  Elements <- bpjchr21_sub$FORMAT[[1]] 
  #loop for extraction
  j <-1
  for (j in 1:numberElements){
    bpjchr21_sub[Elements[j]] <- sapply(bpjchr21_sub$FORMAT_CONTENT, "[[", j)
  }
  # bpjchr21_FORMAT <- merge(bpjchr21_FORMAT, bpjchr21_sub, by = colnames(bpjchr21_FORMAT), all.y = TRUE, all.x = TRUE)
  bpjchr21_FORMAT <- merge(bpjchr21_FORMAT, bpjchr21_sub, all.y = TRUE, all.x = TRUE)
}

#generate DP from AD 
bpjchr21_FORMAT$AD <-  strsplit(as.character(bpjchr21_FORMAT$AD), ",", perl =TRUE)
bpjchr21_FORMAT$RO <- as.numeric(sapply(bpjchr21_FORMAT$AD, function (x) x[1]))
bpjchr21_FORMAT$AO <- as.numeric(sapply(bpjchr21_FORMAT$AD, function (x) x[2]))
bpjchr21_FORMAT$DP <- bpjchr21_FORMAT$RO + bpjchr21_FORMAT$AO

#convert values into numbers
bpjchr21_FORMAT$GQ <- as.numeric(bpjchr21_FORMAT$GQ)

#delete extracted columns
bpjchr21 <- bpjchr21_FORMAT[,-c(3,7,9,10,16)]
rm(bpjchr21_FORMAT, bpjchr21_sub)

#filter against poor quality
bpjchr21_clean <- bpjchr21[bpjchr21$QUAL > thresholdForQual,]

#quick plot
if (SkipPlotSection == FALSE) {
  plot_position <- ggplot(bpjchr21_clean,aes(x=CoordStart_MB, y=log2(QUAL), col = Position)) +
    geom_point()+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_positions", subtitle = "QualFilt40")
  
  plot_position
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_positions", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
}




#end of "Start"



#----heterozygosity----

#quick plotting genotype + overall read depth + ref/alt observations
if (SkipPlotSection == FALSE) {
  plot_GT <- ggplot(bpjchr21_clean,aes(x=CoordStart_MB, y=as.factor(GT), col = Position)) +
    geom_point()+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_GT", subtitle = "QualFilt40")
  
  plot_GT
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_GT", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
  
  
  #read depth
  plot_ReadDepth <- ggplot(bpjchr21_clean[bpjchr21_clean$DP < 100 & !is.na(bpjchr21_clean$DP),],aes(x=CoordStart_MB, y=DP, col = Position,)) +
    geom_point()+
    # geom_smooth(method = "loess", show.legend = FALSE, col = colorsForLegend[3], se = FALSE, span = 0.8)+
    geom_smooth(mapping = aes(group=Position), col = colorsForLegend[2], span = 0.5)+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_ReadDepth", subtitle = "QualFilter40, DP<100")
  
  plot_ReadDepth
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_ReadDepth", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
  
  #ref obs
  plot_ReadDepth <- ggplot(bpjchr21_clean[bpjchr21_clean$DP < 100 & !is.na(bpjchr21_clean$DP),],aes(x=CoordStart_MB, y=RO, col = Position,)) +
    geom_point()+
    # geom_smooth(method = "loess", show.legend = FALSE, col = colorsForLegend[3], se = FALSE, span = 0.8)+
    geom_smooth(mapping = aes(group=Position), col = colorsForLegend[2], span = 0.5)+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_refObs", subtitle = "QualFilter40, DP<100")
  
  plot_ReadDepth
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_RefObs", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)  
  
  #alt obs
  plot_ReadDepth <- ggplot(bpjchr21_clean[bpjchr21_clean$DP < 100 & !is.na(bpjchr21_clean$DP),],aes(x=CoordStart_MB, y=AO, col = Position,)) +
    geom_point()+
    # geom_smooth(method = "loess", show.legend = FALSE, col = colorsForLegend[3], se = FALSE, span = 0.8)+
    geom_smooth(mapping = aes(group=Position), col = colorsForLegend[2], span = 0.5)+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_AltObs", subtitle = "QualFilter40, DP<100")
  
  plot_ReadDepth
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_altObs", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
}

#remove na in DP
bpjchr21_filtered <- bpjchr21_clean[!is.na(bpjchr21_clean$DP),]

##define highly repetitive regions
thresholdReadDepth <- median(bpjchr21_filtered$DP) *3
thresholdRepBoundary <- 1000
PositionsRep <- bpjchr21_filtered$POS[bpjchr21_filtered$DP > thresholdReadDepth]


i <- 1
PositionsRepPOS <- c(rep(NA, length(bpjchr21_filtered$POS)))
for (i in 1:length(bpjchr21_filtered$POS)) {
  PositionsRepPOS[i] <- if_else(any(abs(bpjchr21_filtered$POS[i] - PositionsRep) <= thresholdRepBoundary), "repetitive", "normal")
}

#exclude highly repetitve regions
bpjchr21_filtered$rep <- PositionsRepPOS

#quick plotting for rep region filter
if (SkipPlotSection == FALSE) {
  plot_ReadDepth <- ggplot(bpjchr21_filtered[bpjchr21_filtered$DP < 100,],aes(x=CoordStart_MB, y=DP, col = rep,)) +
    geom_point()+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_FilterRepRegion", subtitle = paste("QualFilter40, DP<100, marked for repRegionBoundary of ", thresholdRepBoundary, sep = ""))
  
  plot_ReadDepth
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_FilterRepRegion", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
  
  #remove rep regions
  plot_ReadDepth <- ggplot(bpjchr21_filtered[bpjchr21_filtered$DP < 100 & bpjchr21_filtered$rep != "repetitive",],aes(x=CoordStart_MB, y=DP, col = Position,)) +
    geom_point()+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_FilterRepRegion", subtitle = paste("QualFilter40, DP<100, REMOVED! for repRegionBoundary of ", thresholdRepBoundary, sep = ""))
  
  plot_ReadDepth
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_RemovedRepRegion", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
  
}




#----analyse SNVs in vicinity of breakpoints_mother----
#----prepare table for Haplo layout
PlotName <- "breakpointAnalysis_mother"

#build backup
bpjchr21mo <- vcfmo21

#focus on junctions
bpjchr21mo$Position <- "Unaffected"
bpjchr21mo$Position[bpjchr21mo$POS > Inversion[1] & bpjchr21mo$POS < Inversion[2] & bpjchr21mo$CHROM == "chr21"]  <- "Inversion"
bpjchr21mo$Position[bpjchr21mo$POS >= hg38$cenStart[21] & bpjchr21mo$POS <= hg38$cenEnd[21]] <- "Centromer"
bpjchr21mo$Position[bpjchr21mo$POS > firstDeletion[1] & bpjchr21mo$POS < firstDeletion[2] & bpjchr21mo$CHROM  == "chr21"]  <- "First Deletion"
bpjchr21mo$Position[bpjchr21mo$POS > secondDeletion[1] & bpjchr21mo$POS < secondDeletion[2] & bpjchr21mo$CHROM  == "chr21"]  <- "Second Deletion"
bpjchr21mo$Position[bpjchr21mo$POS > thirdDeletion[1] & bpjchr21mo$POS < thirdDeletion[2] & bpjchr21mo$CHROM  == "chr21"]  <- "Third Deletion"

#reorder Position
bpjchr21mo$Position <- factor(bpjchr21mo$Position, levels(factor(bpjchr21mo$Position))[c(2,1,4,5,3,6)])

#add conversion of bases into kilo bases
conversionFactorKB <- 1000
bpjchr21mo$CoordStart_KB <- bpjchr21mo$POS / conversionFactorKB
conversionFactorMB <- 1000000
bpjchr21mo$CoordStart_MB <- bpjchr21mo$POS / conversionFactorMB

#quick plot
if (SkipPlotSection == FALSE) {
  plot_position <- ggplot(bpjchr21mo,aes(x=CoordStart_MB, y=log2(QUAL), col = Position)) +
    geom_point()+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_positions_mother", subtitle = "noQualFilter")
  
  plot_position
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_mother_Qual-positions", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
}




##extract information from allele specific column
i <- 1
bpjchr21mo_FORMAT <- data.frame(matrix(nrow = 0, ncol = ncol(bpjchr21mo)))
colnames(bpjchr21mo_FORMAT) <- colnames(bpjchr21mo)

for (i in 1:length(levels(factor(bpjchr21mo$FORMAT)))) {
  #subset for dif levels
  bpjchr21mo_sub <- bpjchr21mo[bpjchr21mo$FORMAT == levels(factor(bpjchr21mo$FORMAT))[i],]
  #split info
  bpjchr21mo_sub$FORMAT <- strsplit(as.character(bpjchr21mo_sub$FORMAT), ":", perl =TRUE)
  bpjchr21mo_sub$FORMAT_CONTENT <- strsplit(as.character(bpjchr21mo_sub$FORMAT_CONTENT), ":", perl =TRUE)
  #read no of elements
  numberElements <- length(bpjchr21mo_sub$FORMAT[[1]])
  Elements <- bpjchr21mo_sub$FORMAT[[1]] 
  #loop for extraction
  j <-1
  for (j in 1:numberElements){
    bpjchr21mo_sub[Elements[j]] <- sapply(bpjchr21mo_sub$FORMAT_CONTENT, "[[", j)
  }
  # bpjchr21mo_FORMAT <- merge(bpjchr21mo_FORMAT, bpjchr21mo_sub, by = colnames(bpjchr21mo_FORMAT), all.y = TRUE, all.x = TRUE)
  bpjchr21mo_FORMAT <- merge(bpjchr21mo_FORMAT, bpjchr21mo_sub, all.y = TRUE, all.x = TRUE)
}

#generate DP from AD 
bpjchr21mo_FORMAT$AD <-  strsplit(as.character(bpjchr21mo_FORMAT$AD), ",", perl =TRUE)
bpjchr21mo_FORMAT$RO <- as.numeric(sapply(bpjchr21mo_FORMAT$AD, function (x) x[1]))
bpjchr21mo_FORMAT$AO <- as.numeric(sapply(bpjchr21mo_FORMAT$AD, function (x) x[2]))
bpjchr21mo_FORMAT$DP <- bpjchr21mo_FORMAT$RO + bpjchr21mo_FORMAT$AO

#convert values into numbers
bpjchr21mo_FORMAT$GQ <- as.numeric(bpjchr21mo_FORMAT$GQ)

#delete extracted columns
bpjchr21mo <- bpjchr21mo_FORMAT[,-c(3,7,9,10,16)]
rm(bpjchr21mo_FORMAT, bpjchr21mo_sub)

#filter against poor quality
bpjchr21mo_clean <- bpjchr21mo[bpjchr21mo$QUAL > thresholdForQual,]

#quick plot
if (SkipPlotSection == FALSE) {
  plot_position <- ggplot(bpjchr21mo,aes(x=CoordStart_MB, y=log2(QUAL), col = Position)) +
    geom_point()+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_positions_mother", subtitle = "QualFilt40")
  
  plot_position
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_positions_mother", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
}


#----heterozygosity_mother----


#quick plotting genotype + overall read depth + ref/alt observations
if (SkipPlotSection == FALSE) {
  plot_GT <- ggplot(bpjchr21mo_clean,aes(x=CoordStart_MB, y=as.factor(GT), col = Position)) +
    geom_point()+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_GT_mother", subtitle = "QualFilter40")
  
  plot_GT
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_GT_mother", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
  
  
  #read depth
  plot_ReadDepth <- ggplot(bpjchr21mo_clean[bpjchr21mo_clean$DP < 100 & !is.na(bpjchr21mo_clean$DP),],aes(x=CoordStart_MB, y=DP, col = Position,)) +
    geom_point()+
    # geom_smooth(method = "loess", show.legend = FALSE, col = colorsForLegend[3], se = FALSE, span = 0.8)+
    geom_smooth(mapping = aes(group=Position), col = colorsForLegend[2], span = 0.5)+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_ReadDepth_mother", subtitle = "QualFilter40, DP<100")
  
  plot_ReadDepth
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_ReadDepth_mother", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
  
  #ref obs
  plot_ReadDepth <- ggplot(bpjchr21mo_clean[bpjchr21mo_clean$DP < 100 & !is.na(bpjchr21mo_clean$DP),],aes(x=CoordStart_MB, y=RO, col = Position,)) +
    geom_point()+
    # geom_smooth(method = "loess", show.legend = FALSE, col = colorsForLegend[3], se = FALSE, span = 0.8)+
    geom_smooth(mapping = aes(group=Position), col = colorsForLegend[2], span = 0.5)+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_refObs_mother", subtitle = "QualFilter40, DP<100")
  
  plot_ReadDepth
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_RefObs_mother", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)  
  
  #alt obs
  plot_ReadDepth <- ggplot(bpjchr21mo_clean[bpjchr21mo_clean$DP < 100 & !is.na(bpjchr21mo_clean$DP),],aes(x=CoordStart_MB, y=AO, col = Position,)) +
    geom_point()+
    # geom_smooth(method = "loess", show.legend = FALSE, col = colorsForLegend[3], se = FALSE, span = 0.8)+
    geom_smooth(mapping = aes(group=Position), col = colorsForLegend[2], span = 0.5)+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_AltObs_mother", subtitle = "QualFilter40, DP<100")
  
  plot_ReadDepth
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_altObs_mother", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
}

#remove na in DP
bpjchr21mo_filtered <- bpjchr21mo_clean[!is.na(bpjchr21mo_clean$DP),]

##define highly repetitive regions
thresholdReadDepth <- median(bpjchr21mo_filtered$DP) *3
thresholdRepBoundary <- 1000
PositionsRep <- bpjchr21mo_filtered$POS[bpjchr21mo_filtered$DP > thresholdReadDepth]


i <- 1
PositionsRepPOS <- c(rep(NA, length(bpjchr21mo_filtered$POS)))
for (i in 1:length(bpjchr21mo_filtered$POS)) {
  PositionsRepPOS[i] <- if_else(any(abs(bpjchr21mo_filtered$POS[i] - PositionsRep) <= thresholdRepBoundary), "repetitive", "normal")
}

#exclude highly repetitve regions
bpjchr21mo_filtered$rep <- PositionsRepPOS

#quick plotting for rep region filter
if (SkipPlotSection == FALSE) {
  plot_ReadDepth <- ggplot(bpjchr21mo_filtered[bpjchr21mo_filtered$DP < 100,],aes(x=CoordStart_MB, y=DP, col = rep,)) +
    geom_point()+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_FilterRepRegion_mother", subtitle = paste("QualFilter40, DP<100, marked for repRegionBoundary of ", thresholdRepBoundary, sep = ""))
  
  plot_ReadDepth
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_FilterRepRegion_mother", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
  
  #remove rep regions
  plot_ReadDepth <- ggplot(bpjchr21mo_filtered[bpjchr21mo_filtered$DP < 100 & bpjchr21mo_filtered$rep != "repetitive",],aes(x=CoordStart_MB, y=DP, col = Position,)) +
    geom_point()+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle("SNV_chr21_ALI9_FilterRepRegion_mother", subtitle = paste("QualFilter40, DP<100, REMOVED! for repRegionBoundary of ", thresholdRepBoundary, sep = ""))
  
  plot_ReadDepth
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_RemovedRepRegion_mother", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
  
}







#----binning SVs for het fraction----
PlotName <- "HetFractionBinning"

#get data
binning_sa <- bpjchr21_clean

#remove NAs in data before loop
binning_sa_loop <- binning_sa[(!is.na(binning_sa$RO) & !is.na(binning_sa$AO)),]

#remove zero counts in AO and RO
binning_sa_loop <- binning_sa[!(binning_sa$RO == 0 & binning_sa$AO == 0),]

#loop through df, generate bins defined by binsize and binlength
binlength <- 50
binsize <- ceiling(nrow(binning_sa_loop) / binlength)
i <- 1
binning_sa_bins <- data.frame(matrix(nrow=binsize, ncol=8))
colnames(binning_sa_bins) <- c("CHROM", "CoordStart_MB", "Position", "DP", "fMean", "fMedian", "fMean2", "fMedian2")

for (i in 1:binsize) {
  #generate bins
  if (nrow(binning_sa_loop) >= (i*binlength)) {
    binning_sa_temp <- binning_sa_loop[(1+(i-1)*binlength):(i*binlength), c(1,8,10,16,17,18)] 
  }else{
    binning_sa_temp <- binning_sa_loop[(1+(i-1)*binlength):nrow(binning_sa_loop),c(1,8,10,16,17,18)]
  }
  
  #calculate fraction of heterozygosity
  # fractionMean <- mean(binning_sa_temp$AO/binning_sa_temp$DP)
  # fractionMedian <- median(binning_sa_temp$AO/binning_sa_temp$DP)
  
  fractionMean2 <- mean(binning_sa_temp$AO/(binning_sa_temp$AO+binning_sa_temp$RO))
  fractionMedian2 <- median(binning_sa_temp$AO/(binning_sa_temp$AO+binning_sa_temp$RO))

  #correct for read depth
  binning_sa_temp$AO_corr <- binning_sa_temp$AO / binning_sa_temp$DP
  binning_sa_temp$RO_corr <- binning_sa_temp$RO / binning_sa_temp$DP
  
  fractionMean <- mean(binning_sa_temp$AO_corr/(binning_sa_temp$AO_corr+binning_sa_temp$RO_corr))
  fractionMedian <- median(binning_sa_temp$AO_corr/(binning_sa_temp$AO_corr+binning_sa_temp$RO_corr))
  
  #write output for i-th bin
  output <- c(paste(binning_sa_temp$CHROM[1]), mean(binning_sa_temp$CoordStart_MB), as.character(binning_sa_temp$Position[1]), median(binning_sa_temp$DP), fractionMean, fractionMedian, fractionMean2, fractionMedian2)
  
  #write data
  binning_sa_bins[i,] <- output
  
}

#clean loop
# rm(output, fractionMean, fractionMedian, binning_sa_temp, binlength, binsize, i, binning_sa_loop)
binning_sa_bins$CoordStart_MB <- as.numeric(binning_sa_bins$CoordStart_MB)
binning_sa_bins$DP <- as.numeric(binning_sa_bins$DP)
binning_sa_bins$fMean <- as.numeric(binning_sa_bins$fMean)
binning_sa_bins$fMedian <- as.numeric(binning_sa_bins$fMedian)
binning_sa_bins$fMean2 <- as.numeric(binning_sa_bins$fMean2)
binning_sa_bins$fMedian2 <- as.numeric(binning_sa_bins$fMedian2)

#reorder Position
binning_sa_bins$Position <- factor(binning_sa_bins$Position, levels(factor(binning_sa_bins$Position))[c(2,1,4,5,3,6)])


#plots
if (SkipPlotSection == FALSE) {

plot_fractionVSposMean2 <- ggplot(binning_sa_bins,aes(x=CoordStart_MB, y=fMean2, col = Position,)) +
  geom_point()+
  # scale_x_continuous(breaks = seq(0, 46, by=5))+
  scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
  ggtitle("SNV_chr21_ALI9_Binning", subtitle = paste("QualFilter40, SNV-mean2", sep = ""))
plot_fractionVSposMean2
ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", ProjectName, "_", PlotName, "_", "patient_SNV-mean2", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       #height = 18, width = 350, unit="cm", limitsize = FALSE)
       height = 30, width = 50, unit="cm", limitsize = FALSE)


plot_fractionVSposMedian2 <- ggplot(binning_sa_bins,aes(x=CoordStart_MB, y=fMedian2, col = Position,)) +
  geom_point()+
  # scale_x_continuous(limits = c(Deletion_chr21p[1], secondDeletion[2]), breaks = seq(1, 50000000, by=5000000))+
  scale_color_manual(values=colorsForLegend[c(2,9,6,3,5,1)])+
  ggtitle("SNV_chr21_ALI9_Binning", subtitle = paste("QualFilter40, SNV-median2", sep = ""))

ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", ProjectName, "_", PlotName, "_", "patient_SNV-median2", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       #height = 18, width = 350, unit="cm", limitsize = FALSE)
       height = 30, width = 50, unit="cm", limitsize = FALSE)

}

#statistics
stats <- binning_sa_bins
#test for normal distribution
for (i in 1:length(levels(factor(stats$Position)))) {
  print(paste0(levels(factor(stats$Position))[i],":"))
  print(shapiro.test(stats$fMedian2[stats$Position == levels(factor(stats$Position))[i]]))
  print(shapiro.test(stats$fMean2[stats$Position == levels(factor(stats$Position))[i]]))
} #non-normal distributed!

#test for homogeneity of variance
for (i in 1:length(levels(factor(stats$Position)))) {
  print(paste0(levels(factor(stats$Position))[i],":"))
  print(leveneTest(stats$fMedian2, stats$Position, center=median))
  print(leveneTest(stats$fMean2, stats$Position, center=median))
} #non-heterogeneous variance

#stats desc
for (i in 1:length(levels(factor(stats$Position)))) {
  print(paste0(levels(factor(stats$Position))[i],":"))
  print(stat.desc(stats$fMedian2[stats$Position == levels(factor(stats$Position))[i]]))
  print(stat.desc(stats$fMean2[stats$Position == levels(factor(stats$Position))[i]]))
}


#prepare non-parametreic tests for all vs all loop
p <- 1
GroupsPositions <- combn(levels(factor(stats$Position)), 2)
GroupsPositions_out <- as.data.frame(matrix(data=NA, nrow=4, ncol=ncol(GroupsPositions)))
GroupsPositions_out$Desc <- c("Mean_W-statistics", "Mean_p-value","Median_W-statistics", "Median_p-value")
GroupsPositions_out <- GroupsPositions_out %>% select(Desc, everything())

SkipExportSection <- TRUE
SkipPlotSection <- TRUE

for (p in 1:ncol(GroupsPositions)){
  
  #define both variables for the whole function (could be skipped and the code be implemented in the corresponding lines but the for loop was added afterwards and doing so it saves time and reduce complexity)
  m <- GroupsPositions[1,p]
  n <- GroupsPositions[2,p]
  
  #do the comparison
  newModel<-wilcox.test(fMean2 ~ Position, data = stats[stats$Position %in% c(m,n),], paired = FALSE)
  newModel2<-wilcox.test(fMedian2 ~ Position, data = stats[stats$Position %in% c(m,n),], paired = FALSE)

  #add statistics to output
  colnames(GroupsPositions_out)[p+1] <- paste(m, "_vs_", n, sep="")
  GroupsPositions_out[1,p+1] <- newModel$statistic
  GroupsPositions_out[2,p+1] <- newModel$p.value
  GroupsPositions_out[3,p+1] <- newModel2$statistic
  GroupsPositions_out[4,p+1] <- newModel2$p.value
}

if (SkipExportSection == FALSE) {
write_xlsx(GroupsPositions_out, paste(getwd(), "/R stuff/ExtractedData_", ProjectName, "/", PlotName, "Nonparametric-test_violinplot_patient",  ".xlsx", sep=""))
}
  
# SelectCol <- "Unaffected"
# temp <- as.vector(colnames(GroupsPositions_out))
# StatsForPlot <- as.data.frame(t(GroupsPositions_out[,which(grepl(pattern = SelectCol, colnames(GroupsPositions_out)))]))
# StatsForPlot$Comp <- rownames(StatsForPlot)
# rownames(StatsForPlot) <- NULL
# StatsForPlot <- StatsForPlot[c(2,1,),]
TextforStatsMean <- c("***", "***", "REF", "***", "ns", "***")
TextforStatsMean2 <- c("1.706E-07", "6.597E-10", "REF", "4.218E-24", "0.067", "3.350E-28")
TextforStatsMedian <- c("***", "***", "REF", "***", "ns", "***")
TextforStatsMedian2 <- c("5.828E-06", "1.296E-10", "REF", "1.303E-17", "0.069", "2.633E-21")

#reorder Position
binning_sa_bins$Position <- factor(binning_sa_bins$Position, levels(factor(binning_sa_bins$Position))[c(1,2,6,3,5,4)])

if (SkipPlotSection == FALSE) {
  # plot_fractionVSposMean2 <- ggplot(binning_sa_bins,aes(x=CoordStart_MB, y=fMean2, col = Position,)) +
  plot_fractionVSposMean2 <- ggplot(binning_sa_bins, aes(y=fMean2, x=Position, fill=Position)) +
    geom_violin(draw_quantiles = c(0.1,0.5,0.9), trim = TRUE, scale="width", )+
    theme_classic(base_size=30) +
    scale_fill_manual(values=colorsForLegend[c(9,2,1,6,3,5)])+
    labs(x = "Region on Chromosome 21", y = "Bins: Fraction of Homozygosity (Mean)") +
    ylim(c(0,1))+
    annotate("text", label =paste("p-value:"), x = 0, y = 1, size = 10, colour = "black", vjust=-0.2, hjust=-0.1)+
    annotate("text", label =paste(TextforStatsMean2)[1], x = 1, y = 1, size = 10, colour = "black", vjust=-0.2)+
    annotate("text", label =paste(TextforStatsMean2)[2], x = 2, y = 1, size = 10, colour = "black", vjust=-0.2)+
    annotate("text", label =paste(TextforStatsMean2)[3], x = 3, y = 1, size = 10, colour = "black", vjust=-0.2)+
    annotate("text", label =paste(TextforStatsMean2)[4], x = 4, y = 1, size = 10, colour = "black", vjust=-0.2)+
    annotate("text", label =paste(TextforStatsMean2)[5], x = 5, y = 1, size = 10, colour = "black", vjust=-0.2)+
    annotate("text", label =paste(TextforStatsMean2)[6], x = 6, y = 1, size = 10, colour = "black", vjust=-0.2)+
    theme(
      plot.subtitle = element_text(colour="black",face="italic", size = 15),
      strip.text = element_text(size = 14),
      axis.ticks = element_line(colour="black"),
      axis.text.x = element_text(colour="black"),
      axis.text.y = element_text(colour="black"),
      #legend.title = element_text(colour = "black", size = 15),
      legend.position = "none")+
    ggtitle("SNV_chr21_ALI9_Binning", subtitle = paste("QualFilter40, SNV-mean2_violin", sep = ""))
  plot_fractionVSposMean2
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", ProjectName, "_", PlotName, "_", "patient_SNV-mean2_boxplot", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
  
  
  plot_fractionVSposMedian2 <- ggplot(binning_sa_bins,aes(x=Position, y=fMedian2, fill = Position,)) +
    geom_violin(draw_quantiles = c(0.1,0.5,0.9), trim = TRUE, scale="width", )+
    theme_classic(base_size=30) +
    scale_fill_manual(values=colorsForLegend[c(9,2,1,6,3,5)])+
    labs(x = "Region on Chromosome 21", y = "Bins: Fraction of Homozygosity (Median)") +
    ylim(c(0,1))+
    annotate("text", label =paste("p-value:"), x = 0, y = 1, size = 10, colour = "black", vjust=-0.2, hjust=-0.1)+
    annotate("text", label =paste(TextforStatsMedian2)[1], x = 1, y = 1, size = 10, colour = "black", vjust=-0.2)+
    annotate("text", label =paste(TextforStatsMedian2)[2], x = 2, y = 1, size = 10, colour = "black", vjust=-0.2)+
    annotate("text", label =paste(TextforStatsMedian2)[3], x = 3, y = 1, size = 10, colour = "black", vjust=-0.2)+
    annotate("text", label =paste(TextforStatsMedian2)[4], x = 4, y = 1, size = 10, colour = "black", vjust=-0.2)+
    annotate("text", label =paste(TextforStatsMedian2)[5], x = 5, y = 1, size = 10, colour = "black", vjust=-0.2)+
    annotate("text", label =paste(TextforStatsMedian2)[6], x = 6, y = 1, size = 10, colour = "black", vjust=-0.2)+
    theme(
      plot.subtitle = element_text(colour="black",face="italic", size = 15),
      strip.text = element_text(size = 14),
      axis.ticks = element_line(colour="black"),
      axis.text.x = element_text(colour="black"),
      axis.text.y = element_text(colour="black"),
      #legend.title = element_text(colour = "black", size = 15),
      legend.position = "none")+
    ggtitle("SNV_chr21_ALI9_Binning", subtitle = paste("QualFilter40, SNV-median2_violin", sep = ""))
  plot_fractionVSposMedian2
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", ProjectName, "_", PlotName, "_", "patient_SNV-median2_boxplot", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       #height = 18, width = 350, unit="cm", limitsize = FALSE)
       height = 30, width = 50, unit="cm", limitsize = FALSE)
}




#----phasing with SNPs----
PlotName <- "PhasingSNPs"

#load data, loop over mother and father
m <- 1
ListForLoop <- c("vcffa21", "vcfmo21", "vcfsa21PCRf")

for (m in 1:length(ListForLoop)) {

#build backup
dataForLoop <- as.data.frame(mget(x = ListForLoop[m]))
colnames(dataForLoop) <- gsub(paste0(ListForLoop[m],"."),"",colnames(dataForLoop),)
dataForLoop$QUAL <- as.numeric(dataForLoop$QUAL)

#focus on junctions
dataForLoop$Position <- "Unaffected"
dataForLoop$Position[dataForLoop$POS > Inversion[1] & dataForLoop$POS < Inversion[2] & dataForLoop$CHROM == "chr21"]  <- "Inversion"
dataForLoop$Position[dataForLoop$POS >= hg38$cenStart[21] & dataForLoop$POS <= hg38$cenEnd[21]] <- "Centromer"
dataForLoop$Position[dataForLoop$POS > firstDeletion[1] & dataForLoop$POS < firstDeletion[2] & dataForLoop$CHROM  == "chr21"]  <- "First Deletion"
dataForLoop$Position[dataForLoop$POS > secondDeletion[1] & dataForLoop$POS < secondDeletion[2] & dataForLoop$CHROM  == "chr21"]  <- "Second Deletion"
dataForLoop$Position[dataForLoop$POS > thirdDeletion[1] & dataForLoop$POS < thirdDeletion[2] & dataForLoop$CHROM  == "chr21"]  <- "Third Deletion"

#reorder Position
dataForLoop$Position <- factor(dataForLoop$Position, levels(factor(dataForLoop$Position))[c(2,1,4,5,3,6)])

#add conversion of bases into kilo bases
conversionFactorKB <- 1000
dataForLoop$CoordStart_KB <- dataForLoop$POS / conversionFactorKB
conversionFactorMB <- 1000000
dataForLoop$CoordStart_MB <- dataForLoop$POS / conversionFactorMB

#quick plot
if (SkipPlotSection == FALSE) {
  plot_position <- ggplot(dataForLoop,aes(x=CoordStart_MB, y=log2(QUAL), col = Position)) +
    geom_point()+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle(paste("SNV_chr21_ALI9_positions_", ListForLoop[m]), subtitle = "noQualFilter")
  
  plot_position
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_", ListForLoop[m], "_Qual-positions", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
  
  
  plot_position <- ggplot(dataForLoop[dataForLoop$QUAL>40,],aes(x=CoordStart_MB, y=log2(QUAL), col = Position)) +
    geom_point()+
    scale_x_continuous(breaks = seq(0, 46, by=5))+
    scale_color_manual(values=colorsForLegend[c(9,2,6,5,3,1)])+
    ggtitle(paste("SNV_chr21_ALI9_positions_", ListForLoop[m]), subtitle = "QualFilter40")
  
  plot_position
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_", ListForLoop[m], "_Qual-positions", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
}

##extract information from allele specific column
i <- 1
dataForLoop_FORMAT <- data.frame(matrix(nrow = 0, ncol = ncol(dataForLoop)))
colnames(dataForLoop_FORMAT) <- colnames(dataForLoop)

for (i in 1:length(levels(factor(dataForLoop$FORMAT)))) {
  #subset for dif levels
  dataForLoop_sub <- dataForLoop[dataForLoop$FORMAT == levels(factor(dataForLoop$FORMAT))[i],]
  #split info
  dataForLoop_sub$FORMAT <- strsplit(as.character(dataForLoop_sub$FORMAT), ":", perl =TRUE)
  dataForLoop_sub$FORMAT_CONTENT <- strsplit(as.character(dataForLoop_sub$FORMAT_CONTENT), ":", perl =TRUE)
  #read no of elements
  numberElements <- length(dataForLoop_sub$FORMAT[[1]])
  Elements <- dataForLoop_sub$FORMAT[[1]] 
  #loop for extraction
  j <-1
  for (j in 1:numberElements){
    dataForLoop_sub[Elements[j]] <- sapply(dataForLoop_sub$FORMAT_CONTENT, "[[", j)
  }
  # dataForLoop_FORMAT <- merge(dataForLoop_FORMAT, dataForLoop_sub, by = colnames(dataForLoop_FORMAT), all.y = TRUE, all.x = TRUE)
  dataForLoop_FORMAT <- merge(dataForLoop_FORMAT, dataForLoop_sub, all.y = TRUE, all.x = TRUE)
}

#generate DP from AD 
dataForLoop_FORMAT$AD <-  strsplit(as.character(dataForLoop_FORMAT$AD), ",", perl =TRUE)
dataForLoop_FORMAT$RO <- as.numeric(sapply(dataForLoop_FORMAT$AD, function (x) x[1]))
dataForLoop_FORMAT$AO <- as.numeric(sapply(dataForLoop_FORMAT$AD, function (x) x[2]))
dataForLoop_FORMAT$DP <- dataForLoop_FORMAT$RO + dataForLoop_FORMAT$AO

#convert values into numbers
dataForLoop_FORMAT$GQ <- as.numeric(dataForLoop_FORMAT$GQ)

#delete extracted columns
dataForLoop <- dataForLoop_FORMAT[,-c(3,7,9,10,16)]
rm(dataForLoop_FORMAT, dataForLoop_sub)

#filter against poor quality
dataForLoop_clean <- dataForLoop[dataForLoop$QUAL > thresholdForQual,]

#extract data from loop
assign(paste0(ListForLoop[m], "_extracted"), dataForLoop)
assign(paste0(ListForLoop[m], "_extracted_clean"), dataForLoop_clean)


rm(j, i, m, dataForLoop_clean, dataForLoop)
}



##label snvs in patient for origin(mother/father/denovo)
#merge trio
Phasing_Merge <- rbind(bpjchr21_clean, vcffa21_extracted_clean, vcfmo21_extracted_clean)

#clean data 
Phasing_Merge <- Phasing_Merge[,-c(5,6,9,12,13,14,15)]
Phasing_Merge <- Phasing_Merge[!is.na(Phasing_Merge$AO) & !is.na(Phasing_Merge$RO) &!is.na(Phasing_Merge$DP),]
Phasing_Merge <- Phasing_Merge[!Phasing_Merge$DP == 0,]

#generate levels for SNVs
Phasing_Merge$SNV <- paste(Phasing_Merge$POS, Phasing_Merge$ALT, sep = "-")
Phasing_Merge <- Phasing_Merge[Phasing_Merge$GT %in% c("0/1", "1/1"),]
Phasing_Merge$origin <- "unknown"

#focus on patient
Phasing_Merge_sa <- Phasing_Merge[Phasing_Merge$name == "ALI09sa_Lr",]
Phasing_Merge_mo <- Phasing_Merge[Phasing_Merge$name == "ALI09mo_Lr",]
Phasing_Merge_fa <- Phasing_Merge[Phasing_Merge$name == "ALI09fa_Lr",]
Phasing_Merge_saPCRf <- Phasing_Merge[Phasing_Merge$name == "ALI09sa_Lr",]

i <- 1
for (i in 1:length(levels(factor(Phasing_Merge_sa$SNV)))) {
  sub <- Phasing_Merge_sa[Phasing_Merge_sa$SNV == Phasing_Merge_sa$SNV[i],]
          if ((sub$SNV %in% Phasing_Merge_mo$SNV) & (sub$SNV %in% Phasing_Merge_fa$SNV)) {
      Phasing_Merge_sa$origin[Phasing_Merge_sa$SNV == Phasing_Merge_sa$SNV[i]] <- "parental"
    }else if ((sub$SNV %in% Phasing_Merge_mo$SNV) & (!sub$SNV %in% Phasing_Merge_fa$SNV)) {
      Phasing_Merge_sa$origin[Phasing_Merge_sa$SNV == Phasing_Merge_sa$SNV[i]] <- "maternal"
    }else if ((!sub$SNV %in% Phasing_Merge_mo$SNV) & (sub$SNV %in% Phasing_Merge_fa$SNV)) {
      Phasing_Merge_sa$origin[Phasing_Merge_sa$SNV == Phasing_Merge_sa$SNV[i]] <- "paternal"
    }else if ((!sub$SNV %in% Phasing_Merge_mo$SNV) & (!sub$SNV %in% Phasing_Merge_fa$SNV)) {
      Phasing_Merge_sa$origin[Phasing_Merge_sa$SNV == Phasing_Merge_sa$SNV[i]] <- "denovo"
    }
}

  
 #remove hom SNVs
# Phasing_Merge_sa <- Phasing_Merge_sa[Phasing_Merge_sa$GT == "0|1",] #not feasable because hemizygosity is regarded as homozygous

#reorder Position
Phasing_Merge_sa$Position <- factor(Phasing_Merge_sa$Position, levels(factor(Phasing_Merge_sa$Position))[c(1,2,6,3,5,4)])



##statistics
#prepare table for binominal test
PhaseStats <- Phasing_Merge_sa

v <- 1
picked <- "maternal"
compare <- "paternal"
PhaseStats_out <- as.data.frame(matrix(NA, nrow = 0, ncol = 8))
colnames(PhaseStats_out) <- 	c("Position", "origin", "counted", "picked", "comparedTo", "pvalLower", "CILower", "nullValLower")
for (v in 1:length(levels(factor(PhaseStats$Position)))) {
  #create df for testing
  PhaseStats_count <- as.data.frame(matrix(NA, nrow = 0, ncol = 8))
  colnames(PhaseStats_count) <- 	c("Position", "origin", "counted", "picked", "comparedTo", "pvalLower", "CILower", "nullValLower")
  
  PhaseStats_count <- merge(PhaseStats_count,(PhaseStats[PhaseStats$Position == levels(factor(PhaseStats$Position))[v],] %>% 
  group_by(Position, add = TRUE) %>%
  group_by(origin, add = TRUE) %>%
  summarise(counted = length(Position))), by = c("Position", "origin", "counted"), all = TRUE)

  #claim comparison
  PhaseStats_count[PhaseStats_count$origin == picked,"picked"] <- picked
  PhaseStats_count[PhaseStats_count$origin == picked,"comparedTo"] <- compare
  
  #run stat test, hypothesis, observed fraction is lower than expected
  test.binomial <- binom.test(
                      PhaseStats_count$counted[PhaseStats_count$origin == picked],
                      sum(PhaseStats_count$counted),
                      PhaseStats_count$counted[PhaseStats_count$origin == compare]/sum(PhaseStats_count$counted),
                      alternative="less")
  #extract data
  PhaseStats_count[PhaseStats_count$origin == picked,"pvalLower"] <- test.binomial$p.value
  PhaseStats_count[PhaseStats_count$origin == picked,"CILower"] <- paste(test.binomial$conf.int, collapse = ",")
  PhaseStats_count[PhaseStats_count$origin == picked,"nullValLower"] <- test.binomial$null.value
  

  PhaseStats_out <- rbind(PhaseStats_out, PhaseStats_count)
}


#remove untested 
PhaseStats_out_Export <- PhaseStats_out[!is.na(PhaseStats_out$picked),]
#save table
write.table(x = PhaseStats_out_Export, file = paste("R stuff/ExtractedData_", ProjectName, "/", PlotName, "_binomialtesting_fraction_NOparental.csv", sep=""), sep = ";", row.names = FALSE, quote = FALSE)

#prepare data for plot
PhaseStats_out$pvalPlot <- round(PhaseStats_out$pvalLower, 4)
PhaseStats_out$pvalPlot[PhaseStats_out$pvalLower < 0.0001] <- "<0.0001"

PhaseStats_out$counted2 <- PhaseStats_out$counted
PhaseStats_out$counted2[PhaseStats_out$counted == max(PhaseStats_out$counted)] <- 15000


#quick plot
if (SkipPlotSection == FALSE) {
  
  ggplot(PhaseStats_out,aes(x=origin, y=counted,  fill=origin)) +
    theme_classic(base_size=30) +
    facet_grid(~Position)+
    geom_col()+
    # scale_fill_manual(values=colorsForLegend[c(9,2,1,6,3,5)])+
    scale_fill_manual(values=colorsForLegend[c(9,8,1,3)])+
    geom_text(aes(x=origin, y=counted+500, label=pvalPlot), size=5, na.rm = TRUE)+
    theme(
      plot.subtitle = element_text(colour="black",face="italic", size = 15),
      strip.text = element_text(size = 20),
      strip.background = element_blank(),
      axis.ticks = element_line(colour="black"),
      axis.text.x = element_text(colour="black", angle = 45, hjust=0.9),
      axis.text.y = element_text(colour="black"),
      legend.title = element_blank(),
      legend.position = c(0.75,1.16)
    )+
    guides(fill = guide_legend(nrow = 1))+
    ggtitle("SNV_chr21_ALI9_Phasing", subtitle = "QualFilter40")
  # ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_SNVs", "_", ".png", sep=""), 
  #        #height = 18, width = 350, unit="cm", limitsize = FALSE)
  #        height = 30, width = 50, unit="cm", limitsize = FALSE)
  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_SNVs", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
  
  
  ggplot(PhaseStats_out[!PhaseStats_out$Position=="Unaffected",],aes(x=origin, y=counted,  fill=origin)) +
    theme_classic(base_size=30) +
    facet_grid(~Position)+
    geom_col()+
    scale_fill_manual(values=colorsForLegend[c(9,8,1,3)])+
    geom_text(aes(x=origin, y=counted+100, label=pvalPlot), size=5, na.rm = TRUE)+
    theme(
      plot.subtitle = element_text(colour="black",face="italic", size = 15),
      strip.text = element_text(size = 20),
      strip.background = element_blank(),
      axis.ticks = element_line(colour="black"),
      axis.text.x = element_text(colour="black", angle = 45, hjust=0.9),
      axis.text.y = element_text(colour="black"),
      legend.title = element_blank(),
      legend.position = c(0.75,1.16)
    )+
    guides(fill = guide_legend(nrow = 1))+
    ggtitle("SNV_chr21_ALI9_Phasing", subtitle = "QualFilter40, no Unaffected")

  ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", "SNV_chr21_ALI9_SNVs", "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
         #height = 18, width = 350, unit="cm", limitsize = FALSE)
         height = 30, width = 50, unit="cm", limitsize = FALSE)
}


#----prepareSNVsfrom ALI9 for VEP----
#merge PCRf and Linked-reads
PlotName <- "VEP-List"

#load data, loop over mother and father
m <- 1
ListForLoop <- c("vcfsa21", "vcfsa21PCRf")

for (m in 1:length(ListForLoop)) {
  
  #build backup
  dataForLoop <- as.data.frame(mget(x = ListForLoop[m]))
  colnames(dataForLoop) <- gsub(paste0(ListForLoop[m],"."),"",colnames(dataForLoop),)
  dataForLoop$QUAL <- as.numeric(dataForLoop$QUAL)
  
  #focus on junctions
  dataForLoop$Position <- "Unaffected"
  dataForLoop$Position[dataForLoop$POS > Inversion[1] & dataForLoop$POS < Inversion[2] & dataForLoop$CHROM == "chr21"]  <- "Inversion"
  dataForLoop$Position[dataForLoop$POS >= hg38$cenStart[21] & dataForLoop$POS <= hg38$cenEnd[21]] <- "Centromer"
  dataForLoop$Position[dataForLoop$POS > firstDeletion[1] & dataForLoop$POS < firstDeletion[2] & dataForLoop$CHROM  == "chr21"]  <- "First Deletion"
  dataForLoop$Position[dataForLoop$POS > secondDeletion[1] & dataForLoop$POS < secondDeletion[2] & dataForLoop$CHROM  == "chr21"]  <- "Second Deletion"
  dataForLoop$Position[dataForLoop$POS > thirdDeletion[1] & dataForLoop$POS < thirdDeletion[2] & dataForLoop$CHROM  == "chr21"]  <- "Third Deletion"
  
  #reorder Position
  dataForLoop$Position <- factor(dataForLoop$Position, levels(factor(dataForLoop$Position))[c(2,1,4,5,3,6)])
  
  #add conversion of bases into kilo bases
  conversionFactorKB <- 1000
  dataForLoop$CoordStart_KB <- dataForLoop$POS / conversionFactorKB
  conversionFactorMB <- 1000000
  dataForLoop$CoordStart_MB <- dataForLoop$POS / conversionFactorMB
  
  ##extract information from allele specific column
  i <- 1
  dataForLoop_FORMAT <- data.frame(matrix(nrow = 0, ncol = ncol(dataForLoop)))
  colnames(dataForLoop_FORMAT) <- colnames(dataForLoop)
  
  for (i in 1:length(levels(factor(dataForLoop$FORMAT)))) {
    #subset for dif levels
    dataForLoop_sub <- dataForLoop[dataForLoop$FORMAT == levels(factor(dataForLoop$FORMAT))[i],]
    #split info
    dataForLoop_sub$FORMAT <- strsplit(as.character(dataForLoop_sub$FORMAT), ":", perl =TRUE)
    dataForLoop_sub$FORMAT_CONTENT <- strsplit(as.character(dataForLoop_sub$FORMAT_CONTENT), ":", perl =TRUE)
    #read no of elements
    numberElements <- length(dataForLoop_sub$FORMAT[[1]])
    Elements <- dataForLoop_sub$FORMAT[[1]] 
    #loop for extraction
    j <-1
    for (j in 1:numberElements){
      dataForLoop_sub[Elements[j]] <- sapply(dataForLoop_sub$FORMAT_CONTENT, "[[", j)
    }
    # dataForLoop_FORMAT <- merge(dataForLoop_FORMAT, dataForLoop_sub, by = colnames(dataForLoop_FORMAT), all.y = TRUE, all.x = TRUE)
    dataForLoop_FORMAT <- merge(dataForLoop_FORMAT, dataForLoop_sub, all.y = TRUE, all.x = TRUE)
  }
  
  #generate DP from AD 
  dataForLoop_FORMAT$AD <-  strsplit(as.character(dataForLoop_FORMAT$AD), ",", perl =TRUE)
  dataForLoop_FORMAT$RO <- as.numeric(sapply(dataForLoop_FORMAT$AD, function (x) x[1]))
  dataForLoop_FORMAT$AO <- as.numeric(sapply(dataForLoop_FORMAT$AD, function (x) x[2]))
  dataForLoop_FORMAT$DP <- dataForLoop_FORMAT$RO + dataForLoop_FORMAT$AO
  
  #convert values into numbers
  dataForLoop_FORMAT$GQ <- as.numeric(dataForLoop_FORMAT$GQ)
  
  #delete extracted columns
  dataForLoop <- dataForLoop_FORMAT[,-c(3,7,9,10,16)]
  rm(dataForLoop_FORMAT, dataForLoop_sub)
  
  #filter against poor quality
  dataForLoop_clean <- dataForLoop[dataForLoop$QUAL > thresholdForQual,]
  
  #extract data from loop
  assign(paste0(ListForLoop[m], "_extracted"), dataForLoop)
  assign(paste0(ListForLoop[m], "_extracted_clean"), dataForLoop_clean)
  
  
  rm(j, i, m, dataForLoop_clean, dataForLoop)
}



##label snvs in patient for origin(mother/father/denovo)
#merge trio
Phasing_Merge <- rbind(vcfsa21_extracted_clean, vcfsa21PCRf_extracted_clean)

#clean data 
Phasing_Merge <- Phasing_Merge[,-c(5,6,9,12,13,14,15)]
Phasing_Merge <- Phasing_Merge[!is.na(Phasing_Merge$AO) & !is.na(Phasing_Merge$RO) &!is.na(Phasing_Merge$DP),]
Phasing_Merge <- Phasing_Merge[!Phasing_Merge$DP == 0,]

#generate levels for SNVs
Phasing_Merge$SNV <- paste(Phasing_Merge$POS, Phasing_Merge$ALT, sep = "-")
Phasing_Merge <- Phasing_Merge[Phasing_Merge$GT %in% c("0/1", "1/1"),]
Phasing_Merge$origin <- "notshared"

#focus on patient
Phasing_Merge_sa <- Phasing_Merge[Phasing_Merge$name == "ALI09sa_Lr",]
Phasing_Merge_saPCRf <- Phasing_Merge[Phasing_Merge$name == "ALI09sa_PCRf",]

i <- 1
for (i in 1:length(levels(factor(Phasing_Merge_sa$SNV)))) {
  sub <- Phasing_Merge_sa[Phasing_Merge_sa$SNV == Phasing_Merge_sa$SNV[i],]
  if ((sub$SNV %in% Phasing_Merge_saPCRf$SNV)) {
    Phasing_Merge_sa$origin[Phasing_Merge_sa$SNV == Phasing_Merge_sa$SNV[i]] <- "shared"
  }
}


export <- Phasing_Merge_sa[Phasing_Merge_sa$origin == "shared",]
export <- export[!export$Position %in% c("Unaffected", "Centromer","Inversion"),]
export$Chr <- "21"
export$Start <- export$POS
export$End <- export$POS + (nchar(export$ALT) - nchar(export$REF))
export$Variant <- paste(export$REF, export$ALT, sep = "/")
export$strand <- "+"
export <- export[,-c(1:13)]

#----ASE_Chr21_patient----
#inherited from hg19 script but adapted to hg38 GATK-ASE run
PlotName <- "ASE_Chr21"

stop("Make sure that the input data is expression not wgs data!")

#build backup
ASE_chr21 <- 

#fiter against indels, reason: noise, for allele specific expression we need high quality SNVs
ASE_chr21 <- ASE_chr21[nchar(as.character(ASE_chr21$REF)) == 1,]
ASE_chr21 <- ASE_chr21[nchar(as.character(ASE_chr21$ALT)) == 1,]

#filter further against poor quality, calculated from original data set 
thresholdForQual_ASE <- median(bpjchr21$QUAL)
ASE_chr21 <- ASE_chr21[ASE_chr21$QUAL > thresholdForQual_ASE,]

#filter against ReadDepth, calculated, between 1st and 3rd quartile
thresholdReadDepth <- c(as.vector(quantile(ASE_chr21$DP))[2],as.vector(quantile(ASE_chr21$DP))[3])
ASE_chr21 <- ASE_chr21[ASE_chr21$DP > thresholdReadDepth[1] & ASE_chr21$DP < thresholdReadDepth[2],]

#manually convert deletes region into heterocygous
ASE_chr21$GT[ASE_chr21$Position %in% c("Second Deletion", "Third Deletion")] <- "0/1"

#filter for heterocygocity, we are interested in the allele specific expression
ASE_chr21 <- ASE_chr21[ASE_chr21$GT == "0/1",]

#calculate allele ratio
ASE_chr21$AD_ratio <- ASE_chr21$AO / ASE_chr21$DP        #AD comes from GATK haplotype caller

#reorder Position
ASE_chr21$Position <- factor(ASE_chr21$Position, levels(factor(ASE_chr21$Position))[c(1,2,6,3,5,4)])


#plot with AD
Plot21 <- ggplot(ASE_chr21, aes(x=CoordStart_MB, y=AD_ratio, col = Position, shape = GT))+
  theme_classic()+
  # annotate(geom = "rect", xmin = firstDeletion[1]/ 10^3, xmax = firstDeletion[2]/ 10^3, ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[5])+
  # annotate(geom = "rect", xmin = secondDeletion[1]/ 10^3, xmax = secondDeletion[2]/ 10^3, ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[6])+
  # annotate(geom = "rect", xmin = surroundingsDeletion[1]/ 10^3, xmax = firstDeletion[1]/ 10^3, ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[3])+
  # annotate(geom = "rect", xmin = firstDeletion[2]/ 10^3, xmax = secondDeletion[1]/ 10^3, ymin = 1, ymax = 0, alpha=0.08, fill = colorsForLegend[3])+
  geom_hline(yintercept = 0.5, linetype =2, alpha = 0.25, size=0.3)+
  geom_hline(yintercept = 0.75, linetype =2, alpha = 0.25, size=0.3)+
  geom_hline(yintercept = 0.25, linetype =2, alpha = 0.25, size=0.3)+
  geom_point(size = 5, alpha = 0.8)+
  scale_x_continuous(breaks = seq(5, 45, by = 5), labels = c("", "10", "", "20","", "30","", "40", ""))+
  scale_color_manual(values=colorsForLegend[c(9,2,1,6,3,5)])+
  ylim(0,1)+
  labs(x="Chromosome 21 [Mb]", y="AD: nVariant / ntotalVariant")+
  ggtitle(paste("VCF (from RNA Seq Data) - ", PlotName, sep=""), 
          subtitle = paste( "filter: QUAL(", thresholdForQual_ASE, "), ReadDepth (", 
          paste(thresholdReadDepth, collapse = ","), "), only heterocygous, deleted region manually described->0/1, Indels excluded",  sep =""))
Plot21
#save
ggsave(filename = paste("R stuff/Graphs_", ProjectName, "/", PlotName, "_", gsub(":", "-", format(Sys.time(), "%Y%b%d_%X")), ".png", sep=""), 
       height = 15, width = 40, unit="cm", limitsize = FALSE)

