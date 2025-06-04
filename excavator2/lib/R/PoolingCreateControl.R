#!/usr/bin/env Rscript
library(argparse)
library(yaml)

parser <- ArgumentParser(prog = "PoolingCreateControl.R", description = "Creates controls pool data to be used in subsequent analysis")

parser$add_argument("programFolder", help = "Path to Excavator2 folder")
parser$add_argument("outputFolder", help = "Path to output folder")
parser$add_argument("targetFolder", help = "Path to target folder")
parser$add_argument("samplesFile", help = "Path to samples file")
parser$add_argument("mode", help = "Experimental design mode [\"pooling\" (default) or \"paired\"]", default = "pooling")
parser$add_argument("inputFolder", help = "Path to input folder")

args <- parser$parse_args()

#vars.tmp <- commandArgs()
#vars <- vars.tmp[4:length(vars.tmp)]
#split.vars <- unlist(strsplit(vars,","))

##  Setting input paths for normalized read count and experimental design ###
#ProgramFolder <- split.vars[1]
#DataFolder <- split.vars[2]
#TargetFolder <- split.vars[3]
#ExperimentalFile <- split.vars[4]
#ExperimentalDesign <- split.vars[5]
#TargetName <- split.vars[6]
#InputFolder <- split.vars[7]

### Load and set experimental design ###
#ExperimentalTable <- read.table(ExperimentalFile, sep = " ", quote = "", header = F)
#LabelName <- as.character(ExperimentalTable[,1])
#PathInVec <- as.character(ExperimentalTable[,2])
#ExpName <- as.character(ExperimentalTable[,3])
ExpName <- unlist(yaml.load_file(input = args$samplesFile))
PathInVec <- file.path(args$inputFolder, ExpName)
LabelName <- names(ExpName)

### Loading target chromosomes ###
TargetChrom <- Sys.glob(file.path(args$targetFolder, "*_chromosome.txt"))
CHR <- readLines(con = TargetChrom, n = -1L, ok = TRUE, warn = TRUE, encoding = "unknown")
unique.chrom <- strsplit(CHR, " ")[[1]]

Path2ExomeRC <- file.path(args$programFolder, "lib/R/LibraryExomeRC.R")
source(Path2ExomeRC)

### Create the vector for the experimental design ###
if (args$mode == "pooling") {
  indC <- grep("C", LabelName)
  PathInVecC <- PathInVec[indC]
  ExpNameC <- ExpName[indC]
  ### Create the RC matrix with all the Experiments ### 
  
  PathRC <- file.path(PathInVecC[1], "RCNorm")
  FileIn <- file.path(PathRC, paste0(ExpNameC[1], ".NRC.RData"))
  load(FileIn)
  
  MetaData <- MatrixNorm[, 1:5]
  Class <- MatrixNorm[, 7]
  RCNorm <- as.numeric(MatrixNorm[, 6])
  if (length(PathInVecC) > 1) {
    for (i in 2:length(PathInVecC)) {
      PathRC <- file.path(PathInVecC[i], "RCNorm")
      FileIn <- file.path(PathRC, paste0(ExpNameC[i], ".NRC.RData"))
      load(FileIn)
      RCNorm <- RCNorm + as.numeric(MatrixNorm[, 6])
    }
    RCNorm <- RCNorm / length(PathInVecC)
  }
  MatrixNorm <- cbind(MetaData, RCNorm, Class)
}

dir.create(file.path(args$outputFolder, "Control", "RCNorm"), recursive = T)
FileOut <- file.path(args$outputFolder, "Control", "RCNorm", "Control.NRC.RData")
save(MatrixNorm, file = FileOut)
