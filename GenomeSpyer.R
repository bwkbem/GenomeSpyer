#This software is a rerelease of GenomePlotterGC.R module of the Genome Crawler
#Suite with alterations in the input and ouput format in order to improve
#transition between modules of eMESS.

#No part of this text or software may be used or reproduced in any matter
#whatsoever without written permission, except in the case of brief quotations
#embodied in critical articles of or for reviews.  For information address
#Brian W. Kirk (bwkbem@yahoo.com).
#Software Copyright (C) June 2003  Brian W. Kirk
#Form/Text Copyright (C) June 2003  Brian W. Kirk
#Released for use with R v1.5.1 (C)

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

          
###############################################################################
#                                FUNCTIONS                                    #
###############################################################################


###############################################################################
#                           CONSTANT DECLARATION                              #
###############################################################################

#CONSTANTS AND VARIABLES OF TEMPORARY USE BEGIN WITH "x, y, or z" AND ALL
#USER DEFINED FUNCTIONS BEGIN WITH "f" or "g" IN ORDER TO FACILITATE REMOVABLE
#FROM MEMORY.  WHEN USING VARIABLES WITH THESE RESERVED INTIAL LETTERS,
#NEAREST CLEANUP SHOULD BE IDENTIFIED TO INSURE PROPER PERMANENCY


                 ############################################
                 #######  DEFINE PERMANENT CONSTANTS  #######
                 ############################################

#Define Input Directory and File
ClusterDataIn <- paste(CrawlerDirectory, CrawlerClusterFile, sep="");
PlotDataIn <- paste(CrawlerDirectory, PlotDataFile, sep="");


             #####################################################
             ##########  DEFINE TEMPORARY CONSTANTS  #############
             #####################################################


#Data Positions in CrawlerClusterData output file
xClClusterProbPos <- 3;
xClStartGenePosPos <- 4;
xClEndGenePosPos <- 5;
    

#Data Positions in CrawlerGeneData output file
xPDGenePos <- 1;
xPDExppValuePos <- 2;
xPDFoldpValuePos <- 3;
xPDLogFoldChngPos <- 6;


###############################################################################
#                             VARIABLE DECLARATION                            #
###############################################################################

#CONSTANTS AND VARIABLES OF TEMPORARY USE BEGIN WITH "x", "y", or "z" AND ALL
#USER DEFINED FUNCTIONS BEGIN WITH "f" or "g" IN ORDER TO FACILITATE REMOVABLE
#FROM MEMORY.  WHEN USING VARIABLES WITH THESE RESERVED INTIAL LETTERS,
#NEAREST CLEANUP SHOULD BE IDENTIFIED TO INSURE PROPER PERMANCY


#Input CrawlerClusterData
xCrawlerClusterTable  <- data.frame(read.delim(ClusterDataIn,
                                       skip=CrawlerHeaderLength));
#Input CrawlerGeneData
xPlotDataTable  <- data.frame(read.delim(PlotDataIn,
                                       skip=CrawlerHeaderLength));

                 #####  Create Vectors from Input Tables  #####

##### Data from CrawlerClusterTable
zClStartGenePos <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                      xClStartGenePosPos]));
zClEndGenePos <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                      xClEndGenePosPos]));   
zClClusterProb <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                         xClClusterProbPos]));

###### Data from CrawlerGeneTable
zPDGenes <- as.character(as.vector(xPlotDataTable[, xPDGenePos]));
zPDExppValue <- as.numeric(as.vector(xPlotDataTable[, xPDExppValuePos]));
zPDFoldpValue <- as.numeric(as.vector(xPlotDataTable[,
                                                      xPDFoldpValuePos]));
zPDLogFoldChng <- as.numeric(as.vector(xPlotDataTable[,
                                                        xPDLogFoldChngPos]));

###############################################################################
#                               MAIN BODY                                     #
###############################################################################
zPDIndex <- 1:length(zPDGenes);

PlotStartPos <- zPDIndex[as.logical(match(zPDGenes,
                                             PlotStartGene, nomatch=0))];
PlotEndPos <- zPDIndex[as.logical(match(zPDGenes, PlotEndGene,
                                               nomatch=0))];
    
PlotGenes <- zPDGenes[PlotStartPos:PlotEndPos];


#Set Plot Colors for Genes
if (BGGenes) {
  Color <- rep("grey30", length(zPDGenes));
} else {
  Color <- rep("black", length(zPDGenes));
}


if (ExpGenes) {    
  Color[which(zPDExppValue <= PlotGeneCut)] <-  "blue";   
}


if (FoldGenes) {
  Color[which(zPDFoldpValue <= PlotGeneCut)] <-  "green";
}


                           ##### Plot Data #####

ylo <- trunc(min(zPDLogFoldChng[PlotStartPos:PlotEndPos])) - 1;
if (ylo > -1) {ylo <- -1}

yhigh <- trunc(max(zPDLogFoldChng[PlotStartPos:PlotEndPos])) + 1;
if (yhigh < 1) {yhigh <- 1}

if (OutToScreen) {
  x11(width=PlotWidth, height=PlotHeight);
  par(bg="black", fg="white", col.axis="white", col.lab="white", 
    col.main="white");

  PlotGenesIndex <- 1:length(PlotGenes);
  plot(PlotGenesIndex, zPDLogFoldChng[PlotStartPos:PlotEndPos],
       col=Color[PlotStartPos:PlotEndPos], type='h', bty="n", xlab="Genes",
       xaxt='n', xaxp=c(1,length(PlotGenes),4), ylab="Log Fold Change",
       ylim=c(ylo, yhigh), col.lab="white");
  title(main=PlotTitle, font=1);    
  if (length(PlotGenes) > 15) {
    Interval <- trunc(length(PlotGenes)/5);
    PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
    axis(1, at=c(1, 1:4*Interval, length(PlotGenes)), labels=PlotLabels);

  } else {
    Interval <- length(PlotGenes);
    PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
    axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
  }
  
  if (PlotCluster) {    
    for (i in 1:length(zClClusterProb)) {
      if (zClClusterProb[i] <= PlotClusterCut) {
        rect(zClStartGenePos[i]-PlotStartPos+1,
           -RectHeight,
           zClEndGenePos[i]-PlotStartPos+1,
           RectHeight,
           border="yellow", col="yellow");
      
      }
    }
  }
  abline(h=0, col="white");    
  abline(h=1, col="white", lty=2);
  abline(h=-1, col="white", lty=2);
  if(Interactive) {
    identify(PlotGenesIndex, zPDLogFoldChng[PlotStartPos:PlotEndPos], PlotGenes,
                                         col="white");  
  }
}


if (OutToFile) {

  postscript(file=paste(OutputDirectory, FileTitle, ".eps",
             sep=""), horizontal=TRUE, onefile=TRUE);
  par(bg="black", fg="white", col.axis="white", col.lab="white", 
    col.main="white");

  PlotGenesIndex <- 1:length(PlotGenes);
  plot(PlotGenesIndex, zPDLogFoldChng[PlotStartPos:PlotEndPos],
       col=Color[PlotStartPos:PlotEndPos], type='h', bty="n", xlab="Genes",
       xaxt='n', xaxp=c(1,length(PlotGenes),4), ylab="Log Fold Change",
       ylim=c(ylo, yhigh), col.lab="white");
  title(main=PlotTitle, font=1); 
  if (length(PlotGenes) > 15) {
    Interval <- trunc(length(PlotGenes)/5);
    PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
    axis(1, at=c(1, 1:4*Interval, length(PlotGenes)), labels=PlotLabels);

  } else {
    Interval <- length(PlotGenes);
    PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
    axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
  }

  if (PlotCluster) {    
    for (i in 1:length(zClClusterProb)) {
      if (zClClusterProb[i] <= PlotClusterCut) {
        rect(zClStartGenePos[i]-PlotStartPos+1,
           -RectHeight,
           zClEndGenePos[i]-PlotStartPos+1,
           RectHeight,
           border="yellow", col="yellow");
      
      }
    }
  }
  abline(h=0, col="white");    
  abline(h=1, col="white", lty=2);
  abline(h=-1, col="white", lty=2);
  
  dev.off()
}



rm(list=ls());

