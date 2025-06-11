#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(prog = "FilterTarget.R", description = "Filters target regions from input BED file")

parser$add_argument("inputBedFile", help = "Path to input BED file")
parser$add_argument("outputFolder", help = "Path to output folder")
parser$add_argument("targetName", help = "Name of targeted resequencing panel")
parser$add_argument("assembly", help = "Reference genome assembly name (\"hg19\", \"hg38\", ...)")
parser$add_argument("windowSize", type = "integer", help = "Size of the window in base pairs")
parser$add_argument("chromosomeCoordinates", help = "Path to chromosome coordinates file of reference genome")
parser$add_argument("gapFile", help = "Path to reference genome gap coordinates file")

args <- parser$parse_args()

#vars <- commandArgs(trailingOnly = TRUE)
#split.vars <- unlist(x = strsplit(x = vars, split = ","))
#BedIn <- split.vars[1]
#ProgramFolder <- split.vars[2] ### REMOVED
#OutputFolder <- split.vars[3]
#target.name <- split.vars[4]
#assembly <- split.vars[5]
step <- args$windowSize #as.numeric(split.vars[6]) ### Shorter
#CoordIn <- split.vars[7]
#FileGap <- split.vars[8]
flank <- 200

options("scipen" = 20)

CoordTable <- read.table(file = args$chromosomeCoordinates, sep = "\t", quote = "\"", fill = T, header = F)
ChrCoord <- as.character(x = CoordTable[, 1])
StartCoord <- as.numeric(CoordTable[, 2])
EndCoord <- as.numeric(CoordTable[, 3])

BedTable <- read.table(file = args$inputBedFile, sep = "\t", quote = "\"", fill = T, header = F)

Chr <- as.character(x = BedTable[, 1])
Start <- as.numeric(x = BedTable[, 2])
End <- as.numeric(x = BedTable[, 3])

ChrVec <- ChrCoord <- c(1:22, "X")
if (nchar(x = Chr[1]) > 3) {
  ChrVec <- ChrCoord <- paste0("chr", ChrVec)
}

emptyS <- emptyE <- emptyChr <- c()

stepP <- step + 2 * flank
Totalempty <- c()
TotalTarget <- c()
for (i in 1:length(x = ChrVec)) {
  
  indCoordC <- which(ChrCoord == ChrVec[i])
  StartCoordC <- StartCoord[indCoordC]
  EndCoordC <- EndCoord[indCoordC]
  indC <- which(x = Chr == ChrVec[i]) 
  if (length(indC)!=0){
      StartC <- Start[indC]
      EndC <- End[indC]
      emptyS <- c(StartCoordC, EndC[1:(length(EndC))])
      emptyE <- c(StartC[1:length(StartC)], EndCoordC)
      
      emptySize <- emptyE - emptyS
      indF <- which(emptySize >= stepP)
      StartCF <- emptyS[indF]
      EndCF <- emptyE[indF]
      emptySizeF <- emptySize[indF]
      
      TotalTarget <- rbind(TotalTarget, cbind(Chr[indC], StartC, EndC, rep("IN", length(indC))))
      
      TotalemptyC <- c()
      for (k in 1:length(indF)) {
        numstep <- floor(x = emptySizeF[k] / step)
        seqstep <- seq(from = 0, to = numstep * step, by = step)
        emptystep <- StartCF[k] + seqstep + flank
        emptystepS <- emptystep[1:(length(emptystep) - 1)] + 1
        emptystepE <- emptystep[2:length(emptystep)]
        TotalemptyC <- rbind(TotalemptyC, cbind(rep(ChrVec[i], length(emptystep) - 1), emptystepS, emptystepE, rep("OUT", length(emptystep) - 1)))
      }
    }
  if (length(indC)==0){
    TotalSize<-EndCoordC-StartCoordC
    NumStep<-floor(TotalSize/step)
    IntervalStep<-seq(StartCoordC,EndCoordC,by=step)
    TotalemptyC<-cbind(rep(ChrVec[i],NumStep),IntervalStep[-length(IntervalStep)],IntervalStep[-1]-1,rep("OUT",NumStep))
  }
  Totalempty <- rbind(Totalempty, TotalemptyC)
  
}

FullTarget <- rbind(TotalTarget, Totalempty)
ChrFull <- as.character(FullTarget[, 1])
StartFull <- as.numeric(FullTarget[, 2])
EndFull <- as.numeric(FullTarget[, 3])
StatusFull <- as.character(FullTarget[, 4])

FinalTarget <- c()
for (i in 1:length(ChrVec)) {
  indC <- which(ChrFull == ChrVec[i]) 
  StartFullC <- StartFull[indC]
  EndFullC <- EndFull[indC]
  StatusFullC <- StatusFull[indC]
  indS <- sort(StartFullC, index.return = T)$ix
  FinalTarget <- rbind(FinalTarget, cbind(ChrFull[indC], StartFullC[indS], EndFullC[indS], StatusFullC[indS]))
  
}

### Filtro Gapmeri e Gap##
GapTable <- read.table(file = args$gapFile, sep = "\t", header = F, skip = 1, fill = TRUE, quote = "")

if (nchar(Chr[1]) > 3) {
  GapChrIn <- as.character(GapTable[, 2])
} else {
  GapChrIn <- sub(pattern = "chr", replacement = "", x = as.character(GapTable[, 2]))
} 

indAlt <- grep(pattern = "_", x = GapChrIn) 
GapChr <- GapChrIn[-indAlt]
GapStart <- as.numeric(GapTable[-indAlt, 3])
GapEnd <- as.numeric(GapTable[-indAlt, 4])

FinalTargetF <- c()
for (c in 1:length(ChrVec)) {
  
  cat(paste("Filtering chromosome", ChrVec[c], "..."), file = stderr())
  
  indGap <- which(GapChr == ChrVec[c])
  GapPos <- cbind(GapStart[indGap], GapEnd[indGap])
  indC <- which(FinalTarget[, 1] == ChrVec[c])
  FinalTargetC <- cbind(FinalTarget[indC, ], rep(x = 0, length(indC)))
  for (i in 1:nrow(FinalTargetC)) {
    if (FinalTargetC[i, 5] != 1) { 
      for (j in 1:nrow(GapPos)) {
        if (
          length(
            intersect(
              x = findInterval(
                x = as.vector(FinalTargetC[i, c(2, 3)]),
                vec = as.vector(GapPos[j, c(1, 2)])
              ),
              y = 1
            )
          ) != 0
        ) {
          FinalTargetC[i, 5] <- 1
        }
      }
    }
  }
  indF <- which(FinalTargetC[, 5] == 0)
  FinalTargetF <- rbind(FinalTargetF, FinalTargetC[indF, c(1, 2, 3, 4)])
  
  cat(" done.\n", file = stderr())
}

#new format
a <- paste0("a", seq(from = 1, to = nrow(FinalTargetF)))
FinalTargetF <- cbind(FinalTargetF[, c(1:3)], a, FinalTargetF[, 4])

MyTarget <- FinalTargetF
MyChr <- unique(FinalTargetF[, 1])

dir.create(path = args$outputFolder, recursive = T)
dir.create(path = file.path(args$outputFolder, "GCC"))
dir.create(path = file.path(args$outputFolder, "MAP"))
dir.create(path = file.path(args$outputFolder, "FRB"))

write.table(
  x = data.frame(rbind(MyChr)),
  file = file.path(args$outputFolder, paste0(args$targetName, "_chromosome.txt")),
  col.names = F,
  row.names = F,
  quote = F
)
save(MyTarget, file = file.path(args$outputFolder, paste0(args$targetName, ".RData")))
write.table(
  x = MyTarget,
  file = file.path(args$outputFolder, "Filtered.txt"),
  col.names = F,
  row.names = F,
  sep = "\t",
  quote = F
)

