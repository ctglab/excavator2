#!/usr/bin/env Rscript
library(argparse)
library(yaml)

parser <- ArgumentParser(prog = "EXCAVATORInferenceExome.R", description = "Generates segmentation and calling plots")

parser$add_argument("outputFolder", help = "Path to output folder")
parser$add_argument("samplesFile", help = "Path to samples file")
parser$add_argument("assembly", help = "Reference genome assembly")

args <- parser$parse_args()

#vars.tmp <- commandArgs()
#vars <- vars.tmp[4:length(vars.tmp)]
#split.vars <- unlist(strsplit(vars, ","))

###  Setting input paths for normalized read count and experimental design ###
#DataFolder <- split.vars[1]
#ExperimentalFile <- split.vars[2]

### Load and set experimental design ###
ExpName <- unlist(yaml.load_file(input = args$samplesFile))
LabelName <- names(ExpName)
#ExperimentalTable <- read.table(args$samplesFile, sep = " ", quote = "", header = F)
#LabelName <- as.character(ExperimentalTable[, 1])
#ExpName <- as.character(ExperimentalTable[, 3])

### Create the vector for the experimental design ###
indT <- grep("T", LabelName)
ExpTest <- c(ExpName[indT])

for (zz in 1:length(ExpTest)) {
  FileSeg <- file.path(
    args$outputFolder,
    "Results",
    ExpTest[zz],
    paste0("HSLMResults_", ExpTest[zz], ".txt")
  )
  FileCall <- file.path(
    args$outputFolder,
    "Results",
    ExpTest[zz],
    paste0("FastCallResults_", ExpTest[zz], ".txt")
  )
  PathOut <- file.path(args$outputFolder, "Plots", ExpTest[zz])
  
  ############   ##################
  
  y <- read.table(FileSeg, sep = "\t", quote = "\"", fill = T, header = T)
  
  log2RSeq <- as.numeric(as.character(y[, 5]))
  chrSeq <- as.character(y[, 1])
  SegSeq <- as.numeric(as.character(y[, 6]))
  PositionSeq <- as.numeric(as.character(y[, 2]))
  TargetSeq <- as.character(as.character(y[, 7]))
  
  UniqueChr <- unique(chrSeq)
  
  z <- read.table(FileCall, sep = "\t", quote = "\"", fill = T, header = T)

  chrCall <- as.character(z[, 1])
  StartCall <- as.numeric(as.character(z[, 2]))
  EndCall <- as.numeric(as.character(z[, 3]))
  Call <- as.numeric(as.character(z[, 7]))

  # used for plots (standard and gtrellis)
  sv <- data.frame(
    name = c("2-DEL", "DEL", "AMP", "2-AMP"),
    call = c(-2:-1, 1:2),
    color = c(
      rgb(213, 94, 0, maxColorValue=255), 
      rgb(230, 159, 0, maxColorValue=255), 
      rgb(86, 180, 233, maxColorValue=255), 
      rgb(0, 114, 178, maxColorValue=255)
    ),
    row.names = c("DD","D","A","AA")
  )
  
  for (i in 1:length(UniqueChr)) {
    indSeq <- which(chrSeq == UniqueChr[i])
    PositionSeqC <- PositionSeq[indSeq]
    log2RSeqC <- log2RSeq[indSeq]
    SegSeqC <- SegSeq[indSeq]
    TargetSeqC <- TargetSeq[indSeq]
    
    indIN <- which(TargetSeqC == "IN")
    indOUT <- which(TargetSeqC == "OUT")
    FileName <- paste0("PlotResults_", UniqueChr[i], ".pdf")
    
    FilePlot <- file.path(PathOut, FileName)
    pdf(FilePlot, height = 10, width = 15)
    par(mfrow = c(2, 1))
    if (length(indOUT) != 0) {
      PositionSeqCIN <- PositionSeqC[indIN]
      PositionSeqCOUT <- PositionSeqC[indOUT]
      log2RSeqCIN <- log2RSeqC[indIN]
      log2RSeqCOUT <- log2RSeqC[indOUT]
      SegSeqCIN <- SegSeqC[indIN]
      SegSeqCOUT <- SegSeqC[indOUT]
      
      plot(
        PositionSeqCOUT,
        log2RSeqCOUT,
        ylim = c(-3, 3),
        main = UniqueChr[i],
        pch = 19,
        cex = 0.3,
        xlab = "Position",
        ylab = "log2ratio",
        col = "lightblue"
      )
      points(
        PositionSeqCIN,
        log2RSeqCIN,
        col = "blue",
        pch = 19,
        cex = 0.3
      )
      lines(PositionSeqC, SegSeqC, lwd = 2, col = "red")
      abline(h = 0, lty = 2, lwd = 1, col = "black")
    }
    if (length(indOUT) == 0) {
      plot(
        PositionSeqC,
        log2RSeqC,
        ylim = c(-3, 3),
        main = UniqueChr[i],
        pch = 19,
        cex = 0.3,
        xlab = "Position",
        ylab = "log2ratio",
        col = "blue"
      )
      lines(PositionSeqC, SegSeqC, lwd = 2, col = "red")
      abline(h = 0, lty = 2, lwd = 1, col = "black")
    }
    
    plot(
      PositionSeqC,
      log2RSeqC,
      ylim = c(-3, 3),
      type = "l",
      main = UniqueChr[i],
      lwd = 0.5,
      col = "grey",
      xlab = "Position",
      ylab = "log2ratio"
    )
    
    indCall <- which(chrCall == UniqueChr[i])
    if (length(indCall) != 0) {
      StartCallC <- StartCall[indCall]
      EndCallC <- EndCall[indCall]
      CallC <- Call[indCall]
      for (j in 1:(length(indCall))) {
        if (CallC[j] == 1) {
          rect(
            xleft = StartCallC[j],
            ybottom = 0,
            xright = EndCallC[j],
            ytop = 3,
            density = NA,
            angle = 45,
            col = sv["D",]$color,
            border = sv["D",]$color
          )
        }
        if (CallC[j] == 2) {
          rect(
            xleft = StartCallC[j],
            ybottom = 0,
            xright = EndCallC[j],
            ytop = 3,
            density = NA,
            angle = 45,
            col = sv["DD",]$color,
            border = sv["DD",]$color
          )
        }
        if (CallC[j] == -1) {
          rect(
            xleft = StartCallC[j],
            ybottom = -3,
            xright = EndCallC[j],
            ytop = 0,
            density = NA,
            angle = 45,
            col = sv["A",]$color,
            border = sv["A",]$color
          )
        }
        if (CallC[j] == -2) {
          rect(
            xleft = StartCallC[j],
            ybottom = -3,
            xright = EndCallC[j],
            ytop = 0,
            density = NA,
            angle = 45,
            col = sv["AA",]$color,
            border = sv["AA",]$color
          )
        }
      }
    }
    
    lines(
      PositionSeqC,
      log2RSeqC,
      ylim = c(-3, 3),
      main = UniqueChr[i],
      lwd = 0.5,
      col = "grey"
    )
    abline(h = 0, lwd = 1, col = "black")
    dev.off()
  }
  
  ### gtrellis plot
  
  suppressPackageStartupMessages(library(gtrellis))
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(ComplexHeatmap))
  
  col_fun <- colorRamp2(sv$call, sv$color)
  lgd <- Legend(title = "Class", type = "points", at = sv$call, labels = sv$name, legend_gp = gpar(col = as.character(sv$color)))
  
  pdf(file = file.path(PathOut, "PlotResults.pdf"), width = 13, height = 7)
  
  gtrellis_layout(
    species = tolower(args$assembly),
    n_track = 1,
    ncol = 1,
    track_axis = FALSE,
    xpadding = c(0.05, 0),
    gap = unit(1, "mm"),
    border = FALSE,
    asist_ticks = FALSE,
    add_ideogram_track = TRUE, 
    ideogram_track_height = unit(2, "mm"),
    legend = lgd,
    title = "EXCAVATOR2 Calls"
  )
  
  if (nrow(z) > 0) {
    if (!grepl("^chr", z[1,1])) {
      z[[1]] <- paste0("chr", z[[1]])
    }
  }
  

  add_track(z, panel_fun = function(z) {
    grid.rect(
      x = z[[2]],
      y = unit(0.2, "npc"),
      width = z[[3]] - z[[2]],
      height = unit(0.8, "npc"),
      hjust = 0,
      vjust = 0,
      default.units = "native",
      gp = gpar(fill = col_fun(z[[7]]), col = 'black')
    )
  })
  
  add_track(track = 2, clip = FALSE, panel_fun = function(gr) {
    chr = get_cell_meta_data("name")
    if (chr == "chrY") {
      grid.lines(get_cell_meta_data("xlim"), unit(c(0, 0), "npc"), default.units = "native")
    }
    grid.text(chr, x = 0, y = 0, just = c("left", "bottom"))
  })
  
  dev.off()

}
