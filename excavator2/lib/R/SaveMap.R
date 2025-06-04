#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(prog = "SaveMap.R", description = "Saves chromosome Mappability information to output folder")

parser$add_argument("chromosome", help = "Selected chromosome")
parser$add_argument("outputFolder", help = "Path to output folder")

args <- parser$parse_args()

# vars.tmp <- commandArgs()
# vars <- vars.tmp[length(vars.tmp)]
# split.vars <- unlist(strsplit(vars,","))
# chrsel <- split.vars[1]
# MAPPath <- split.vars[2]

FileIn <- file.path(args$outputFolder, "Mapout.txt")

MappingTable <- read.table(FileIn, sep = "\t", quote = "", header = F)
MapMed <- MappingTable[, 6]

save(MapMed, file = file.path(args$outputFolder, paste0("Map.", args$chromosome, ".RData")))
