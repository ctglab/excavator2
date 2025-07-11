#!/usr/bin/env Rscript
library(argparse)
library(yaml)

parser <- ArgumentParser(prog = "EXCAVATORInferenceExome.R", description = "Performs segmentation and calling for CNV analysis.")

parser$add_argument("outputFolder", help = "Path to output folder")
parser$add_argument("targetPath", help = "Path to target folder")
parser$add_argument("samplesFile", help = "Path to samples file")
parser$add_argument("mode", help = "Experimental design mode [\"pooling\" (default) or \"paired\"]", default = "pooling")
parser$add_argument("programFolder", help = "Path to Excavator2 folder")
parser$add_argument("assembly", help = "Reference genome assembly name (\"hg19\", \"hg38\", ...)")
parser$add_argument("inputFolder", help = "Path to input folder")
parser$add_argument("parametersFile", help = "Path to parameters file")
parser$add_argument("centromeresFile", help = "Path to centromeres file")

args <- parser$parse_args()

#vars.tmp <- commandArgs()
#vars <- vars.tmp[4:length(vars.tmp)]
#split.vars <- unlist(strsplit(vars, ","))

##  Setting input paths for normalized read count and experimental design ###
#DataFolder <- split.vars[1]
#TargetFolder <- split.vars[2]
#ExperimentalFile <- split.vars[3]
#ExperimentalDesign <- split.vars[4]
#TargetName <- split.vars[5]
#ProgramFolder <- split.vars[6]
#Assembly <- split.vars[7]
#InputFolder <- split.vars[8]
#ParametersFile <- split.vars[9]
#CentromeresFile <- split.vars[10]

### Load and set experimental design ###
ExpName <- unlist(yaml.load_file(input = args$samplesFile))
PathInVec <- file.path(args$inputFolder, ExpName)
LabelName <- names(ExpName)
ExperimentalTable <- data.frame(label = LabelName, path = PathInVec, sample = ExpName)
#ExperimentalTable <- read.table(args$samplesFile, sep = " ", quote = "", header = F)
#LabelName <- as.character(ExperimentalTable[, 1])
#PathInVec <- as.character(ExperimentalTable[, 2])
#ExpName <- as.character(ExperimentalTable[, 3])

### Create the vector for the experimental design ###
if (args$mode == "pooling") {
  ExpTest <- ExpName
  PathInVecTS <- PathInVec
  FileRef <- file.path(args$outputFolder, "Control", "RCNorm", "Control.NRC.RData")
  load(file = FileRef)
  DataSeqRef <- as.numeric(MatrixNorm[, 6])
}

if (args$mode == "paired") {
  indC <- grep("C", LabelName)
  indT <- grep("T", LabelName)
  LabelNameC <- LabelName[indC]
  LabelNameT <- LabelName[indT]
  ExpNameC <- ExpName[indC]
  ExpNameT <- ExpName[indT]
  PathInVecC <- PathInVec[indC]
  PathInVecT <- PathInVec[indT]
  NumC <- as.numeric(substr(LabelNameC, 2, 100000))
  NumT <- as.numeric(substr(LabelNameT, 2, 100000))
  indCS <- sort(NumC, index.return = T)$ix
  indTS <- sort(NumT, index.return = T)$ix
  ExpTest <- ExpNameT[indTS]
  ExpControl <- ExpNameC[indCS]
  PathInVecTS <- PathInVecT[indTS]
  PathInVecCS <- PathInVecC[indCS]
}

### Loading target chromosomes ###
TargetChrom <- Sys.glob(file.path(args$targetPath, "*_chromosome.txt"))
CHR <- scan(TargetChrom, what = "character") #, quiet = T
unique.chrom <- paste0(CHR, ".")

source(file.path(args$programFolder, "/lib/R/LibraryFastCall.R"))
dyn.load(file.path(args$programFolder, "/lib/F77/FastJointSLMLibraryI.so"))
source(file.path(args$programFolder, "/lib/R/LibraryJSLMIn.R"))

### Loading centromere file ###
CentromereTable <- read.table(file = args$centromeresFile, sep = "\t", quote = "", header = T)
CentroChr <- as.character(CentromereTable[, 1])
if (nchar(CHR[1]) < 4) CentroChr <- substr(CentroChr, 4, 100000)

CentroStart <- as.numeric(CentromereTable[, 2])
CentroEnd <- as.numeric(CentromereTable[, 3])

for (zz in 1:length(ExpTest)) {
  ExpLabelOut <- ExpTest[zz]
  ### Loading normalized read count for reference sample when args$mode = paired or pooling ####
  if (args$mode == "paired") {
    FileRef <- file.path(PathInVecCS[zz], "RCNorm", paste0(ExpControl[zz], ".NRC.RData"))
    load(file = FileRef)
    RefMatrixNorm <- MatrixNorm
    DataSeqRef <- as.numeric(RefMatrixNorm[, 6])
  }
  
  ### Loading normalized read count for test sample ###
  FileTest <- file.path(PathInVecTS[zz], "RCNorm", paste0(ExpTest[zz], ".NRC.RData"))
  
  load(file = FileTest)
  
  TestMatrixNorm <- MatrixNorm
  DataSeqTest <- as.numeric(TestMatrixNorm[, 6])
  Position <- as.integer(TestMatrixNorm[, 2])
  chrom <- as.character(TestMatrixNorm[, 1])
  start <- as.numeric(TestMatrixNorm[, 3])
  end <- as.numeric(TestMatrixNorm[, 4])
  Gene <- as.character(TestMatrixNorm[, 5])
  Class <- as.character(TestMatrixNorm[, 7])
  indInTarget <- which(Class == "IN")
  indOutTarget <- which(Class == "OUT")
  
  ### Calculating Log2-ratio and lowess normalization for In e Out ###
  A <- 0.5 * log2(DataSeqTest * DataSeqRef)
  M <- log2(DataSeqTest / DataSeqRef)
  
  AIn <- A[indInTarget]
  MIn <- M[indInTarget]
  LogDataNormIn <- MIn - median(MIn)
  
  AOut <- A[indOutTarget]
  MOut <- M[indOutTarget]
  LogDataNormOut <- MOut - median(MOut)
  
  LogDataNorm <- rep(NA, length(DataSeqTest))
  LogDataNorm[indInTarget] <- LogDataNormIn
  LogDataNorm[indOutTarget] <- LogDataNormOut
  
  ### Setting starting parameters of HSLM ###
  parameters <- yaml.load_file(input = args$parametersFile)
  omega <- parameters$HSLM$Omega
  eta <- as.numeric(parameters$HSLM$Theta)
  stepeta <- as.numeric(parameters$HSLM$D_norm)
  cell <- parameters$FastCall$Cellularity
  thrd <- parameters$FastCall$d
  thru <- parameters$FastCall$u
  FW <- parameters$FastCall$minExons
  
  ### Calculating parameters of the HSLM ###
  mw <- 1
  ParamList <- ParamEstSeq(rbind(LogDataNorm), omega)
  mi <- ParamList$mi
  smu <- ParamList$smu
  sepsilon <- ParamList$sepsilon
  muk <- MukEst(rbind(LogDataNorm), mw)
  
  ###  Segmentation of the log2-ratio profiles with HSLM ###
  MatrixSeg <- matrix(NA, nrow = length(LogDataNorm), ncol = 8)
  indCountS <- 1
  for (i in 1:length(CHR)) {
    chr <- as.character(unlist(CHR))[i]
    
    indchr <- which(chrom == chr)
    seqChrom <- LogDataNorm[indchr]
    PosChrom <- Position[indchr]
    startChrom <- start[indchr]
    endChrom <- end[indchr]
    GeneChrom <- Gene[indchr]
    ClassChrom <- Class[indchr]
    
    indCountE <- indCountS + length(indchr) - 1
    
    seqChrom <- rbind(seqChrom)
    
    splitchrom1 <- CentroStart[which(CentroChr == chr)]
    splitchrom2 <- CentroEnd[which(CentroChr == chr)]
    
    splitind1 <- tail(which(PosChrom < splitchrom1), 1)
    splitind2 <- which(PosChrom > splitchrom2)[1]
    ChrSeg <- c()
    if (length(splitind1) != 0) {
      ind1 <- c(1:splitind1)
      ind2 <- c(splitind2:length(seqChrom))
      
      DataSeq1 <- rbind(seqChrom[ind1])
      DataSeq2 <- rbind(seqChrom[ind2])
      Pos1 <- PosChrom[1:splitind1]
      Pos2 <- PosChrom[splitind2:length(seqChrom)]
      start1 <- startChrom[1:splitind1]
      end1 <- endChrom[1:splitind1]
      start2 <- startChrom[splitind2:length(seqChrom)]
      end2 <- endChrom[splitind2:length(seqChrom)]
      Gene1 <- GeneChrom[ind1]
      Gene2 <- GeneChrom[ind2]
      Class1 <- ClassChrom[ind1]
      Class2 <- ClassChrom[ind2]
      
      TotalPredBreak1 <- JointSegIn(DataSeq1, muk, mi, smu, sepsilon, Pos1, omega, eta, stepeta)
      if (length(which(is.na(TotalPredBreak1))) == 0) {
        TotalPredBreak1 <- FilterSeg(TotalPredBreak1, FW)
        DataSeg1 <- SegResults(DataSeq1, TotalPredBreak1)
        ChrSeg <- rbind(
          ChrSeg,
          cbind(
            rep(chr, length(ind1)),
            Pos1,
            Gene1,
            start1,
            end1,
            t(DataSeq1),
            t(DataSeg1),
            Class1
          )
        )
      }
      if (length(which(is.na(TotalPredBreak1))) != 0) {
        message(
          "The HSLM analysis of short arm of chromosome ",
          chr,
          " was aborted because the total number of EMRC is too small: ",
          length(ind1)
        )
        ChrSeg <- rbind(
          ChrSeg,
          cbind(
            rep(chr, length(ind1)),
            Pos1,
            Gene1,
            start1,
            end1,
            matrix(NA, nrow = length(ind1), ncol = 2),
            Class1
          )
        )
      }
      
      TotalPredBreak2 <- JointSegIn(DataSeq2, muk, mi, smu, sepsilon, Pos2, omega, eta, stepeta)
      if (length(which(is.na(TotalPredBreak2))) == 0) {
        TotalPredBreak2 <- FilterSeg(TotalPredBreak2, FW)
        DataSeg2 <- SegResults(DataSeq2, TotalPredBreak2)
        ChrSeg <- rbind(
          ChrSeg,
          cbind(
            rep(chr, length(ind2)),
            Pos2,
            Gene2,
            start2,
            end2,
            t(DataSeq2),
            t(DataSeg2),
            Class2
          )
        )
      }
      if (length(which(is.na(TotalPredBreak2))) != 0) {
        message(
          "The HSLM analysis of long arm of chromosome ",
          chr,
          " was aborted because the total number of EMRC is too small: ",
          length(ind1)
        )
        ChrSeg <- rbind(
          ChrSeg,
          cbind(
            rep(chr, length(ind2)),
            Pos2,
            Gene2,
            start2,
            end2,
            matrix(NA, nrow = length(ind2), ncol = 2),
            Class2
          )
        )
      }
    }
    
    if (length(splitind1) == 0) {
      TotalPredBreak <- JointSegIn(seqChrom, muk, mi, smu, sepsilon, PosChrom, omega, eta, stepeta)
      if (length(which(is.na(TotalPredBreak))) == 0) {
        TotalPredBreak <- FilterSeg(TotalPredBreak, FW)
        DataSeg <- SegResults(seqChrom, TotalPredBreak)
        ChrSeg <- cbind(
          rep(chr, length(indchr)),
          PosChrom,
          cbind(GeneChrom),
          startChrom,
          endChrom,
          t(seqChrom),
          t(DataSeg),
          ClassChrom
        )
      }
      if (length(which(is.na(TotalPredBreak))) != 0) {
        message(
          "The HSLM analysis of chromosome ",
          chr,
          " was aborted because the total number of EMRC is too small: ",
          length(indchr)
        )
        ChrSeg <- rbind(
          ChrSeg,
          cbind(
            rep(chr, length(indchr)),
            PosChrom,
            cbind(GeneChrom),
            startChrom,
            endChrom,
            matrix(NA, nrow = length(indchr), ncol = 2),
            ClassChrom
          )
        )
      }
    }
    
    MatrixSeg[c(indCountS:indCountE), ] <- ChrSeg
    indCountS <- indCountE + 1
  }
  
  MatrixSeg1 <- MatrixSeg[, c(1, 2, 4, 5, 6, 7, 8)]
  colnames(MatrixSeg1) <-
    c("Chromosome", "Position", "Start", "End", "Log2R", "SegMean", "Class")
  FileOutSeg <- file.path(args$outputFolder, "Results", ExpTest[zz], paste0("HSLMResults_", ExpTest[zz], ".txt"))
  write.table( MatrixSeg1, FileOutSeg, col.names = T, row.names = F, sep = "\t", quote = F)
  
  #### Filtering Not-Segmented Data ##
  DataFilt <- as.character(MatrixSeg[, 6])
  indFilt <- which(is.na(DataFilt))
  if (length(indFilt) != 0) MatrixSeg <- MatrixSeg[-indFilt, ]
  
  MatrixSegFC <- MatrixSeg[, -c(2, 8)]
  
  ###  FastCall analysis ###
  AnalisiList <- MakeData(MatrixSegFC, infoPos.StartEnd = TRUE)
  
  MetaData <- AnalisiList$MetaTable
  SummaryData <- AnalisiList$SummaryData
  
  mdata <- SummaryData[, 4]
  
  if (cell < 1) {
    datac <- (2 ^ (mdata) / (cell) - (1 - cell) / cell)
    thrc <- 2 ^ (-5)
    datac[which(datac < thrc)] <- thrc
    mdata <- log2(datac)
  }
  
  ### EM algorithm ####
  ResultsEM <- EMFastCall(mdata, thru, thrd)
  muvec <- ResultsEM$muvec
  sdvec <- ResultsEM$sdvec
  prior <- ResultsEM$prior
  bound <- ResultsEM$bound
  
  P0 <- PosteriorP(mdata, muvec, sdvec, prior)
  
  out <- LabelAss(P0, mdata)
  
  ### Filtering Significant Segments ###
  indSig <- which(out[, 1] != 0)
  
  if (length(indSig) != 0) {
    outSig <- rbind(out[indSig, ])
    P0Sig <- rbind(P0[indSig, ])
    SummaryDataSig <- rbind(SummaryData[indSig, ])
    CNFSig <- 2 * (2 ^ SummaryDataSig[, 4])
    CNSig <- round(CNFSig)
    
    #### FastCall Results in BED format #####
    OutBedSig <- cbind(
      MetaData[SummaryDataSig[, 2], 1],
      MetaData[SummaryDataSig[, 2], 3],
      MetaData[SummaryDataSig[, 3], 4],
      SummaryDataSig[, 4],
      CNFSig,
      CNSig,
      outSig
    )
  }
  if (length(indSig) == 0)  OutBedSig <- c()
  HeaderBed <- c("Chromosome", "Start", "End", "Segment", "CNF", "CN", "Call", "ProbCall")
  OutBedSig <- rbind(HeaderBed, OutBedSig)
  
  FileOutCall <- file.path(args$outputFolder, "Results", ExpTest[zz], paste0("FastCallResults_", ExpTest[zz], ".txt"))
  write.table(OutBedSig, FileOutCall, col.names = F, row.names = F, sep = "\t", quote = F)
  
  #### FastCall Results in VCF format for Regions and Windows #####
  VCFWindowCreate(args$assembly, args$outputFolder, ExpLabelOut, args$targetPath, SummaryData, MetaData, out, DataFilt)
  VCFRegionCreate(args$assembly, args$outputFolder, ExpLabelOut, args$targetPath, SummaryData, MetaData, out)
}
