#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(prog = "SaveFRB.R", description = "Saves chromosome FRB information to output folder")

parser$add_argument("chromosome", help = "Selected chromosome")
parser$add_argument("outputFolder", help = "Path to output folder")

args <- parser$parse_args()

# vars.tmp <- commandArgs()
# vars <- vars.tmp[length(vars.tmp)]
# split.vars <- unlist(strsplit(vars, ","))
# chrsel <- split.vars[1]
# FRBPath <- split.vars[2]

File_In <- file.path(args$outputFolder, "FRB.txt")
FRBMat <- read.table(File_In, sep = "\n", quote = "\"", fill = T, header = F)
FRBMat <- as.character(FRBMat[, 1])
indCoord <- seq(1, length(FRBMat), by = 2)
indSeq <- seq(2, length(FRBMat), by = 2)

FRBData <- cbind(
  unlist(strsplit(
    x = unlist(strsplit(FRBMat[indCoord], ":"))[seq(2, length(indCoord) * 2, by = 2)],
    split = "-"
  ))[seq(1, length(indCoord) * 2, by = 2)],
  FRBMat[indSeq]
)

save(FRBData, file = file.path(args$outputFolder, paste0("FRB.", args$chromosome, ".RData")))
