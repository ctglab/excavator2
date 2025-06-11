#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(prog = "SaveGCC.R", description = "Saves chromosome GCC information to output folder")

parser$add_argument("chromosome", help = "Selected chromosome")
parser$add_argument("outputFolder", help = "Path to output folder")

args <- parser$parse_args()

# vars.tmp <- commandArgs()
# vars <- vars.tmp[length(vars.tmp)]
# split.vars <- unlist(strsplit(vars, ","))
# chrsel <- split.vars[1]
# GCCPath <- split.vars[2]

FileIn <- file.path(args$outputFolder, "GCC.txt")
GCContent <- scan(FileIn)

save(GCContent, file = file.path(args$outputFolder, paste0("GCC.", args$chromosome, ".RData")))
