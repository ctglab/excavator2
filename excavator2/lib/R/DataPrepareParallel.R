#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(prog = "DataPrepareParallel.R", description = "Splits bash jobs for multi-processor data preparation")

parser$add_argument("programFolder", help = "Path to Excavator2 folder")
parser$add_argument("samplesFile", help = "Path to samples list file")
parser$add_argument("outputFolder", help = "Path to output folder")
parser$add_argument("targetFolder", help = "Path to target folder")
parser$add_argument("threads", type = "integer", help = "Number of threads to use")

args <- parser$parse_args()

# vars.tmp <- commandArgs()
# vars <- vars.tmp[4:length(vars.tmp)]
# split.vars <- gsub(pattern = "\\s", replacement = "\\ ", x = unlist(strsplit(vars, ",")))

###  Setting input paths for normalized read count and experimental design ###
#ProgramFolder <- split.vars[1]
#ExperimentalFile <- split.vars[2]
#PathOut <- split.vars[3]
#TargetName <- split.vars[4]
#Assembly <- split.vars[5]
#Threads <- as.numeric(split.vars[6])
#TargetsPath <- split.vars[7]

SettingLabel <- paste0(
  " -o ", args$outputFolder,
  " -t ", args$targetFolder
)

tmp <- file.path(args$outputFolder, ".tmp")
ExperimentalTable <- readLines(con = args$samplesFile)
ExperimentalTable <- grep(pattern = "^[^#]+", x = ExperimentalTable, value = T)

NExp <- length(ExperimentalTable)
if (args$threads > NExp) args$threads <- NExp

Q <- NExp %/% args$threads
R <- NExp %% args$threads
ExpPart <- c(rep(Q + 1, times = R), rep(Q, times = args$threads - R))

ShVector <- c()
StartInd <- 1
for (i in 1:length(ExpPart)) {
  EndInd <- StartInd + ExpPart[i] - 1
  ExperimentalTableSplit <- ExperimentalTable[StartInd:EndInd]
  FileOut <- file.path(tmp, sub("\\.([^.]+)$", paste0(".", i, ".\\1"), basename(args$samplesFile)))
  write.table(ExperimentalTableSplit, file = FileOut, col.names = F, row.names = F, quote = F)
  ShVector <- c(ShVector, paste0("perl ", file.path(args$programFolder, "lib", "perl", "ReadPerla.pl"), " -s ", FileOut, SettingLabel))
  StartInd <- EndInd + 1
}

FileShOut <- file.path(tmp, "ParallelReadPerla.sh")
write.table(ShVector, file = FileShOut, col.names = F, row.names = F, sep = "\t", quote = F)

