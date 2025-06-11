#!/usr/bin/env Rscript
library(argparse)
library(stringr)

parser <- ArgumentParser(prog = "DataAnalysisParallel.R", description = "Splits bash jobs for multi-processor data analysis")

parser$add_argument("programFolder", help = "Path to Excavator2 folder")
parser$add_argument("samplesFile", help = "Path to samples file")
parser$add_argument("outputFolder", help = "Path to output folder")
parser$add_argument("threads", type = "integer", help = "Number of threads to use")
parser$add_argument("experiment", help = "Experimental design [\"pooling\" (default) or \"paired\"] to let me understand", default = "pooling")
parser$add_argument("inputFolder", help = "Path to input folder")
parser$add_argument("targetPath", help = "Path to target folder")
parser$add_argument("parametersFile", help = "Path to parameters file")

args <- parser$parse_args()

#vars.tmp <- commandArgs()
#vars <- vars.tmp[4:length(vars.tmp)]
#split.vars <- gsub(pattern = "\\s", replacement = "\\ ", x = unlist(strsplit(vars, ",")))

###  Setting input paths for normalized read count and experimental design ###
#ProgramFolder <- split.vars[1]
#ExperimentalFile <- split.vars[2]
#OutputFolder <- split.vars[3]
#TargetName <- split.vars[4]
#Assembly <- split.vars[5]
#Processors <- as.numeric(split.vars[6])
#experiment <- split.vars[7]
#InputFolder <- split.vars[8]
#TargetsFolder <- split.vars[9]
#ParametersFile <- split.vars[10]
#CentromeresFile <- split.vars[11]

PathOut <- file.path(args$outputFolder, ".tmp")
SettingLabel <- paste0(
  " -i ", args$inputFolder,
  " -o ", args$outputFolder,
  " -t ", args$targetPath,
  " -e ", args$experiment,
  " -p ", args$parametersFile
)

### Load and set experimental design ###
ExperimentalTable <- readLines(con = args$samplesFile)
ExperimentalTable <- grep(pattern = "^[^#]+", x = ExperimentalTable, value = T)

if (args$experiment == "pooling") {
  indT <- grep("^T", ExperimentalTable)
  ExperimentalTableFinal <- ExperimentalTable[indT]
  
  NExp <- length(ExperimentalTableFinal)
  if (args$threads > NExp) args$threads <- NExp
  Q <- NExp %/% args$threads
  R <- NExp %% args$threads
  ExpPart <- c(rep(Q + 1, times = R), rep(Q, times = args$threads - R))
  
  ShVector <- c()
  StartInd <- 1
  for (i in 1:length(ExpPart)) {
    EndInd <- StartInd + ExpPart[i] - 1
    ExperimentalTableSplit <- ExperimentalTableFinal[StartInd:EndInd]
    FileOut <- file.path(PathOut, sub("\\.([^.]+)$", paste0(".", i, ".\\1"), basename(args$samplesFile)))
    write.table(ExperimentalTableSplit, file = FileOut, col.names = F, row.names = F, quote = F)
    ShVector <- c(ShVector, paste0("perl ", args$programFolder, "/lib/perl/AnalyzePerla.pl ", SettingLabel, " -s ", FileOut))
    StartInd <- EndInd + 1
  }
} else if (args$experiment == "paired") {
  ExperimentalTable <- str_sort(ExperimentalTable, numeric = TRUE)
  L <- length(ExperimentalTable)
  if (L %% 2 == 1) stop("In \"paired\" mode the number of test and control samples must be the same!")
  expT <-data.frame("C" = ExperimentalTable[1:(L/2)] ,
                    "T" = ExperimentalTable[((L/2)+1):L],stringsAsFactors = F) 
  expTnumber <- data.frame("C" = gsub("^C|:.*","",expT$C),
                           "T" = gsub("^T|:.*","",expT$T),stringsAsFactors = F)
  if (!all(expTnumber$C == expTnumber$T)) stop("In \"paired\" mode test and control samples must match!")
  if (any(duplicated(expTnumber))) stop("Test and control samples must be unique!")
  
  NExp <- nrow(expT) 
  if (args$threads > NExp) args$threads <- NExp
  Q <- NExp %/% args$threads
  R <- NExp %% args$threads
  ExpPart <- c(rep(Q + 1, times = R), rep(Q, times = args$threads - R)) * 2
  
  ShVector <- c()
  StartInd <- 1
  for (i in 1:length(ExpPart)) {
    EndInd <- StartInd + ExpPart[i] - 1
    ExperimentalTableSplit <- as.vector(t(expT))[StartInd:EndInd]
    FileOut <- file.path(PathOut, sub("\\.([^.]+)$", paste0(".", i, ".\\1"), basename(args$samplesFile)))
    write.table(ExperimentalTableSplit, file = FileOut, col.names = F, row.names = F, quote = F)
    ShVector <- c(ShVector, paste0("perl ", args$programFolder, "/lib/perl/AnalyzePerla.pl ", SettingLabel, " -s ", FileOut))
    StartInd <- EndInd + 1
  }
}


FileShOut <- file.path(PathOut, "ParallelAnalyzePerla.sh")
write.table(ShVector, FileShOut, col.names = F, row.names = F, sep = "\t", quote = F)
