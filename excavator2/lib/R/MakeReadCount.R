#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(prog = "MakeReadCount.R", description = "Creates read count data by sample")

parser$add_argument("sampleName", help = "Name of selected sample")
parser$add_argument("outputFolder", help = "Path to output folder")
parser$add_argument("chromosome", help = "Name of selected chromosome")
parser$add_argument("programFolder", help = "Path to Excavator2 folder")
parser$add_argument("targetFolder", help = "Path to target folder")

args <- parser$parse_args()

# vars.tmp <- commandArgs()
# vars <- vars.tmp[4:length(vars.tmp)]
# split.vars <- unlist(strsplit(vars, ","))
#Sample_Name <- split.vars [1]
#Output_Folder <- split.vars[2]
#chr <- split.vars[3]
#Program_Folder <- split.vars[4]
#target.name <- split.vars[5]
#assembly.name <- split.vars[6]
#target.folder <- split.vars[7]

targetData <- Sys.glob(file.path(args$targetFolder, "*.RData"))
load(targetData)
dyn.load(paste0(args$programFolder, "/lib/F77/F4R.so"))
fileIN <- file.path(args$outputFolder, ".tmp", ".FilteredBamtmp.000")
fileconIN <- file(fileIN, open = "r")

chrom <- as.character(MyTarget[, 1])
ix.chr <- which(chrom == args$chromosome)
start <- as.numeric(MyTarget[ix.chr, 2])
end <- as.numeric(MyTarget[ix.chr, 3])
Nexome <- length(start)
step <- 500000
jex <- 1
residual <- 0
RC <- rep(0, Nexome)
ConLogic <- 1
while (ConLogic != 0) {
  tt <- readLines(fileconIN, n = step)
  pari <- seq(2, length(tt) * 2, by = 2)
  ReadVector <- as.integer(unlist(strsplit(tt, "\t"))[pari])
  N <- length(ReadVector)
  control <- (ReadVector[N] - end[jex])
  if (control <= 0) {
    residual <- N
    if (jex > Nexome) RC[jex] <- RC[jex] + residual
  }
  if (control > 0) {
    out <- .Fortran(
	  "EXOMECOUNT",
	  as.integer(ReadVector),
	  as.integer(Nexome),
	  as.integer(start),
	  as.integer(end),
	  as.integer(N),
	  as.integer(residual),
	  as.integer(RC),
	  as.integer(jex)
	)
    residual <- out[[6]]
    RC <- out[[7]]
    jex <- out[[8]]
  }
  if (N < step | jex > Nexome) ConLogic <- 0
}

close.connection(fileconIN)  
save(RC, file = file.path(args$outputFolder, "RC", paste0(args$sampleName, ".RC.", args$chromosome, ".RData")))
