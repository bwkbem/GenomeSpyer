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
CyberTDataIn <- paste(CyberTDirectory, CyberTFile, sep="");


             #####################################################
             ##########  DEFINE TEMPORARY CONSTANTS  #############
             #####################################################


#Data Positions in Crawler output file
xGeClStartIDPos <- 1;
xGeClEndIDPos <- 2;
xGeClStartPosPos <- 3;
xGeClEndPosPos <- 4;
xGeClPermProbPos <- 5;
xGeClLengthPos <- 6;
xGeNamePos <- 7;
xGePositionPos <- 8;
xGetStatPos <- 9; 
xGepValuePos <- 10;
xGeLogFoldChngPos <- 11;
xGeIntensityPos <- 12;
xGeBioSampleNumberPos <- 13; 
xGeArrayNumberPos <- 14;
xKuPos <- 15;
xKsubPos <- 16;
xPgPos <- 17;
xPKsubPos <- 18;
xPKuPos <- 19;
xPKmiPos <- 20;
xPKetPos <- 21; 
xPgKuPos <- 22;
xPgIKuPos <- 23;
xPgEKuPos <- 24;
xPgKsubPos <- 25;
xPgKmiPos <- 26;
xPgKetPos <- 27;

#Data Positions in CyberT .GeCyTDS output file
xGenomeIDPos <- 1;
xGenomePosPos <- 2;
xGenometStatPos <- 3;

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
#Input Genome Data from CyberT
xGenomeDataTable  <- data.frame(read.delim(CyberTDataIn,
                                       skip=CrawlerHeaderLength));

                 #####  Create Vectors from Input Tables  #####

##### Data from CrawlerClusterTable

zClStartID <- as.character(as.vector(xCrawlerClusterTable[,
                                                      xGeClStartIDPos]));
zClEndID <- as.character(as.vector(xCrawlerClusterTable[,
                                                      xGeClEndIDPos]));
zGeClStartPos <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                      xGeClStartPosPos]));
zGeClEndPos <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                      xGeClEndPosPos]));
zGeClPermProb <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                           xGeClPermProbPos]));
zGeClLength <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                         xGeClLengthPos]));
zGeName <- as.character(as.vector(xCrawlerClusterTable[,
                                                      xGeNamePos]));
zGePosition <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                         xGePositionPos]));
zGetStat <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                      xGetStatPos]));
zGepValue <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                       xGepValuePos]));
zGeLogFoldChng <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                      xGeLogFoldChngPos]));
zGeIntensity <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                      xGeIntensityPos]));
zGeBioSampleNumber <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                      xGeBioSampleNumberPos]));
zGeArrayNumber <- as.numeric(as.vector(xCrawlerClusterTable[,
                                                      xGeArrayNumberPos]));
zKu <- as.numeric(as.vector(xCrawlerClusterTable[, xKuPos]));
zKsub <- as.numeric(as.vector(xCrawlerClusterTable[, xKsubPos]));
zPg <- as.numeric(as.vector(xCrawlerClusterTable[, xPgPos]));
zPKsub <- as.numeric(as.vector(xCrawlerClusterTable[, xPKsubPos]));
zPKu <- as.numeric(as.vector(xCrawlerClusterTable[, xPKuPos]));
zPKmi <- as.numeric(as.vector(xCrawlerClusterTable[, xPKmiPos]));
zPKet <- as.numeric(as.vector(xCrawlerClusterTable[, xPKetPos]));
zPgKu <- as.numeric(as.vector(xCrawlerClusterTable[, xPgKuPos]));
zPgIKu <- as.numeric(as.vector(xCrawlerClusterTable[, xPgIKuPos]));
zPgEKu <- as.numeric(as.vector(xCrawlerClusterTable[, xPgEKuPos]));
zPgKsub <- as.numeric(as.vector(xCrawlerClusterTable[, xPgKsubPos]));
zPgKmi <- as.numeric(as.vector(xCrawlerClusterTable[, xPgKmiPos]));
zPgKet <- as.numeric(as.vector(xCrawlerClusterTable[, xPgKetPos]));


###### Data Genome Data from xGenomeDataTable
zGenomeID <- as.character(as.vector(xGenomeDataTable[, xGenomeIDPos]));
zGenomePos <- as.numeric(as.vector(xGenomeDataTable[, xGenomePosPos]));
zGenometStat <- as.numeric(as.vector(xGenomeDataTable[, xGenometStatPos]));


###############################################################################
#                               MAIN BODY                                     #
###############################################################################

##Create filter to select all genes present in the cluster from the cluster file
ClusterFilter1 <- match(zClStartID, ClusterStartGene, nomatch=0);

ClusterFilter2 <- match(zClEndID, ClusterEndGene, nomatch=0);

ClusterFilter <- as.logical(ClusterFilter1 * ClusterFilter2);

ClusterPresentGenes <- zGeName[ClusterFilter];

## From an ordered  genome list of all genes in the genome find all the genes
## in a cluster plus extra genes around for plotting data                  

ClusterStartPos <- zGenomePos[as.logical(match(zGenomeID,
                                            ClusterStartGene, nomatch=0))];
if ((ClusterStartPos-ExtraGenes) > 0) {
 PlotStartPos <- ClusterStartPos-ExtraGenes;
} else {
    PlotStartPos <- 1;
}

ClusterEndPos <- zGenomePos[as.logical(match(zGenomeID, ClusterEndGene,
                                               nomatch=0))];
if ((ClusterEndPos+ExtraGenes) < zGenomePos[length(zGenomePos)]) {
    PlotEndPos <- ClusterEndPos+ExtraGenes;
} else {
    PlotEndPos <- zGenomePos[length(zGenomePos)];
}

#ClusterGenes <- zGenomeID[ClusterStartPos:ClusterEndPos];
PlotGenes <- zGenomeID[PlotStartPos:PlotEndPos];
PlotGenePos <- zGenomePos[PlotStartPos:PlotEndPos];

#### Set Plot Colors for Genes
Color <- rep("grey30", length(PlotGenes));

Color[which(as.logical(match(PlotGenes, ClusterPresentGenes,
                             nomatch=0)))] <- "green";

Color[which(zGenometStat[PlotGenePos] == 0)] <- "red";





############################ Edited to here!!!!!!!!!!!!!!!!!!!!!


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

