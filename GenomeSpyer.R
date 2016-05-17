
####  Run SyerLoad.R first to load files !!!!!!!!!!!!!!!!!!!!!!!

###############################################################################
#                               MAIN BODY                                     #
###############################################################################

OutToFile <- TRUE;

if (OutToScreen) {
    OutToFile <- FALSE;
}

xStringSplit <- unlist(strsplit(ClusterEndGene, "SPy"));
FileTitle <- paste(as.character(ClusterStartGene), "to",
                   as.character(xStringSplit[2]), sep="");
    
##Create filter to select all genes present in the cluster from the cluster file
ClusterFilter1 <- match(zClStartID, ClusterStartGene, nomatch=0);

ClusterFilter2 <- match(zClEndID, ClusterEndGene, nomatch=0);

ClusterFilter <- as.logical(ClusterFilter1 * ClusterFilter2);

if (sum(ClusterFilter) == 0) {
    ClusterFound <- FALSE;
} else {
    ClusterFound <- TRUE;
}

ClusterPresentGenes <- zGeName[ClusterFilter];

ClPermProb <- zGeClPermProb[ClusterFilter];

if (ClusterFound) {
    if (ClPermProb[1] == 0) {
    #Find Total Permutations from .GeClGC file header
        xPermTotScan <- scan(file=ClusterDataIn, what=character(0), sep="\t",
                             nlines=CrawlerHeaderLength);
        xPermTotString <- grep("Total Permutations = ", xPermTotScan,
                               value=TRUE);

        if (length(xPermTotString) != 0) {
            CrawlPermTotFound <- TRUE;
            xStringSplit <- unlist(strsplit(xPermTotString, "= "));
            CrawlPermTot <- as.character(xStringSplit[2]); 
        } else {
            CrawlPermTotFound <- FALSE;
            CrawlPermTot <- "NotFound";
        }
        if (CrawlPermTotFound) {
            TotalPermutations <- as.numeric(CrawlPermTot);
            ClPermProb <- paste("<", 1/TotalPermutations, sep="");
        } else {
            ClPermProb <- "<?";
        }
    } else {
        ClPermProb <- as.character(ClPermProb[1]);
    }
} else {
    ClPermProb <- "Cluster Not Found";
}
## From an ordered  genome list of all genes in the genome find all the genes
## in a cluster plus extra genes around for plotting data                  

ClGenomeStartPos <- zGenomePos[as.logical(match(zGenomeID,
                                            ClusterStartGene, nomatch=0))];

ClGenomeEndPos <- zGenomePos[as.logical(match(zGenomeID, ClusterEndGene,
                                               nomatch=0))];


if ((ClGenomeStartPos-ExtraGenes) > 0) {
 PlotStartPos <- ClGenomeStartPos-ExtraGenes;
} else {
    PlotStartPos <- 1;
}


if ((ClGenomeEndPos+ExtraGenes) < zGenomePos[length(zGenomePos)]) {
    PlotEndPos <- ClGenomeEndPos+ExtraGenes;
} else {
    PlotEndPos <- zGenomePos[length(zGenomePos)];
}

PlotGenes <- zGenomeID[PlotStartPos:PlotEndPos];
PlotGenePos <- zGenomePos[PlotStartPos:PlotEndPos];


#### Set Plot Colors for Genes
Color <- rep("blue", length(PlotGenes));

Color[which(as.logical(match(PlotGenes, ClusterPresentGenes,
                             nomatch=0)))] <- "green";

Color[which(zGenometStat[PlotGenePos] == 0)] <- "red";

                         ##################### 
                         ##### Plot Data #####
                         #####################


######################## Log Fold Change #####################################

PlotLogFoldChng <- 0;
ClusterLogFoldChng <- zGeLogFoldChng[ClusterFilter];


for (i in 1:length(PlotGenes)) {
    if (Color[i] == "green") {
        ClLogFoldFilter <- as.logical(match(ClusterPresentGenes, PlotGenes[i],
                                                                   nomatch=0));
        PlotLogFoldChng <- c(PlotLogFoldChng,
                                        ClusterLogFoldChng[ClLogFoldFilter]);

    } else {
        if (Color[i] == "blue") {
            GeLogFoldFilter <- as.logical(match(zGenomeID, PlotGenes[i],
                                                                   nomatch=0));
            PlotLogFoldChng <- c(PlotLogFoldChng,
                                 zGenomeLogFoldChng[GeLogFoldFilter]);
            } else {
                PlotLogFoldChng <- c(PlotLogFoldChng, 0);
            }
    }
}

PlotLogFoldChng <- PlotLogFoldChng[2:length(PlotLogFoldChng)];

ylo <- trunc(min(PlotLogFoldChng)) - 1;
if (ylo > -1) {ylo <- -1}

yhigh <- trunc(max(PlotLogFoldChng)) + 1;
if (yhigh < 1) {yhigh <- 1}

if (OutToScreen) {


    x11(width=PlotWidth, height=PlotHeight);
    par(bg="black", fg="white", col.axis="white", col.lab="white", 
        col.main="white");
    par(mar=c(5,5,5,5), omi=c(0.5,0.5,0.5,0.5));
    par(mfrow=c(2,4));

  
    PlotGenesIndex <- 1:length(PlotGenes);
    plot(PlotGenesIndex, PlotLogFoldChng,
         col=Color, lwd=4, type='h', bty="n", xlab="Genes", xaxt='n',
         xaxp=c(1,length(PlotGenes),4), ylab="Log Fold Change",
         ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
    title(main=paste("Cluster ", ClusterStartGene, "-", ClusterEndGene,
              "\n PKu(perm) = ", ClPermProb, sep=" "), cex.main=1.5);

    if (length(PlotGenes) > 25) {
        Interval <- trunc(length(PlotGenes)/5);
        PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
        axis(1, at=c(1, 1:4*Interval, length(PlotGenes)), labels=PlotLabels);

    } else {
        Interval <- length(PlotGenes);
        PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
        axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
    }
    
    abline(h=0, col="white");    
    abline(h=1, col="white", lty=2);
    abline(h=-1, col="white", lty=2);

    if (PlotCluster) {    
        rect(ClGenomeStartPos-PlotStartPos+1,
             -RectHeight,
             ClGenomeEndPos-PlotStartPos+1,
             RectHeight,
             border="yellow", col="yellow");
    }

    PlotPos <- 1:length(PlotGenes);
    AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
    AbsentGeneValues <- rep(0, length(AbsentGenesPos));
    points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")
  
    if(Interactive) {
        identify(PlotGenesIndex, PlotLogFoldChng, PlotGenes, col="white");  
    }
}


if (OutToFile) {

    postscript(file=paste(OutputDirectory, FileTitle, ".eps",
                   sep=""), horizontal=TRUE, onefile=FALSE);

    par(bg="black", fg="white", col.axis="white", col.lab="white", 
        col.main="white");

    par(mar=c(5,5,5,5), omi=c(0.5,0.5,0.5,0.5));
    par(mfrow=c(2,4));

    
    PlotGenesIndex <- 1:length(PlotGenes);
    plot(PlotGenesIndex, PlotLogFoldChng,
         col=Color, lwd=2, type='h', bty="n", xlab="Genes", xaxt='n',
         xaxp=c(1,length(PlotGenes),4), ylab="Log Fold Change",
         ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
    title(main=paste("Cluster ", ClusterStartGene, "-", ClusterEndGene,
              "\n PKu(perm) = ", ClPermProb, sep=" "), font=2, cex.main=1.5);


    if (length(PlotGenes) > 25) {
        Interval <- trunc(length(PlotGenes)/5);
        PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
        axis(1, at=c(1, 1:4*Interval, length(PlotGenes)), labels=PlotLabels);

    } else {
        Interval <- length(PlotGenes);
        PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
        axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
    }
  
    abline(h=0, col="white");    
    abline(h=1, col="white", lty=2);
    abline(h=-1, col="white", lty=2);

    if (PlotCluster) {    
        rect(ClGenomeStartPos-PlotStartPos+1,
             -RectHeight,
             ClGenomeEndPos-PlotStartPos+1,
             RectHeight,
             border="yellow", col="yellow");
    }

    PlotPos <- 1:length(PlotGenes);
    AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
    AbsentGeneValues <- rep(0, length(AbsentGenesPos));
    points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

    #dev.off()
}

############# Plot PKu-Pg for a given PKu  ###############

if (ClusterFound) {

    ClusterPKu <-zPKu[ClusterFilter];
    ClusterPg <- zPg[ClusterFilter];
    ClusterPKumPg <- ClusterPKu -ClusterPg;
    
    PlotPKumPg <- 0;

    for (i in 1:length(PlotGenes)) {
        if (Color[i] == "green") {
            TargetFilter <- as.logical(match(ClusterPresentGenes, PlotGenes[i],
                                                                   nomatch=0));
            PlotPKumPg <- c(PlotPKumPg, ClusterPKumPg[TargetFilter]);

        } else {
     
            PlotPKumPg <- c(PlotPKumPg, 0);
        } 
    }
    

    PlotPKumPg <- PlotPKumPg[2:length(PlotPKumPg)];

    ylo <- 0;
    yhigh <- 1;

    if (OutToScreen) {

#        x11(width=PlotWidth, height=PlotHeight);
#        par(bg="black", fg="white", col.axis="white", col.lab="white", 
#            col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPKumPg,
             col=Color, lwd=4, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="PKu - Pg",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-",
                  ClusterEndGene, "\n PKu = ", round(ClusterPKu[1], 5),
                  sep=""), font=2, cex.main=1.5);
        
        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    
        abline(h=1, col="white", lty=1);
        abline(h=0.05, col="white", lty=2);
        abline(h=0.10, col="white", lty=2);
        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
        
        if(Interactive) {
            identify(PlotGenesIndex, PlotPKumPg, PlotGenes, col="white");  
        }
    }


    if (OutToFile) {

#        postscript(file=paste(OutputDirectory, FileTitle, ".eps",
#             sep=""), horizontal=TRUE, onefile=TRUE);
#        par(bg="black", fg="white", col.axis="white", col.lab="white", 
#            col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPKumPg,
             col=Color, lwd=2, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="PKu - Pg",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-", ClusterEndGene,
                  "\n PKu = ", round(ClusterPKu[1], 5), sep=""), font=2,
                  cex.main=1.5);
        
        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    
        abline(h=1, col="white", lty=1);
        abline(h=0.05, col="white", lty=2);
        abline(h=0.10, col="white", lty=2);
        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
 #       dev.off()
    }

} else {

    print(paste("Cluster", ClusterStartGene, "-", ClusterEndGene, "not found.",
                sep=" "));
}


############# Plot PgKu ##############


if (ClusterFound) {

    ClusterPKu <-zPKu[ClusterFilter];
    ClusterPgKu <- zPgKu[ClusterFilter];
    
    PlotPgKu <- 0;

    for (i in 1:length(PlotGenes)) {
        if (Color[i] == "green") {
            TargetFilter <- as.logical(match(ClusterPresentGenes, PlotGenes[i],
                                                                   nomatch=0));
            PlotPgKu <- c(PlotPgKu, ClusterPgKu[TargetFilter]);

        } else {
     
            PlotPgKu <- c(PlotPgKu, 0);
        } 
    }
    

    PlotPgKu <- PlotPgKu[2:length(PlotPgKu)];

    ylo <- 0;
    yhigh <- 1;

    if (OutToScreen) {

#        x11(width=PlotWidth, height=PlotHeight);
#        par(bg="black", fg="white", col.axis="white", col.lab="white", 
#            col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPgKu,
             col=Color, lwd=4, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="P(g | Ku)",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-",
                  ClusterEndGene, "\n PKu = ", round(ClusterPKu[1], 5),
                  sep=""), font=2, cex.main=1.5);
        
        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    
        abline(h=1, col="white", lty=1);
        abline(h=0.95, col="white", lty=2);
        abline(h=0.90, col="white", lty=2);
        abline(h=0.85, col="white", lty=2);
        abline(h=0.80, col="white", lty=2);
        abline(h=0.70, col="white", lty=2);
        abline(h=0.60, col="white", lty=2);
        abline(h=0.50, col="white", lty=2);
        abline(h=0.40, col="white", lty=2);
        abline(h=0.30, col="white", lty=2);
        abline(h=0.20, col="white", lty=2);
        abline(h=0.10, col="white", lty=2);
        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
        
        if(Interactive) {
            identify(PlotGenesIndex, PlotPgKu, PlotGenes, col="white");  
        }
    }


    if (OutToFile) {

#        postscript(file=paste(OutputDirectory, FileTitle, ".eps",
#             sep=""), horizontal=TRUE, onefile=TRUE);
#        par(bg="black", fg="white", col.axis="white", col.lab="white", 
#            col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPgKu,
             col=Color, lwd=2, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="P(g | Ku)",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-", ClusterEndGene,
                  "\n PKu = ", round(ClusterPKu[1], 5), sep=""), font=2,
                  cex.main=1.5);

        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    
        abline(h=1, col="white", lty=1);
        abline(h=0.95, col="white", lty=2);
        abline(h=0.90, col="white", lty=2);
        abline(h=0.85, col="white", lty=2);
        abline(h=0.80, col="white", lty=2);
        abline(h=0.70, col="white", lty=2);
        abline(h=0.60, col="white", lty=2);
        abline(h=0.50, col="white", lty=2);
        abline(h=0.40, col="white", lty=2);
        abline(h=0.30, col="white", lty=2);
        abline(h=0.20, col="white", lty=2);
        abline(h=0.10, col="white", lty=2);
        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
#        dev.off()
    }

} else {

    print(paste("Cluster", ClusterStartGene, "-", ClusterEndGene, "not found.",
                sep=" "));
}

################# Plot PgIKu  ################


if (ClusterFound) {

    ClusterPKu <-zPKu[ClusterFilter];
    ClusterPgIKu <- zPgIKu[ClusterFilter];
    
    PlotPgIKu <- 0;

    for (i in 1:length(PlotGenes)) {
        if (Color[i] == "green") {
            TargetFilter <- as.logical(match(ClusterPresentGenes, PlotGenes[i],
                                                                   nomatch=0));
            PlotPgIKu <- c(PlotPgIKu, ClusterPgIKu[TargetFilter]);

        } else {
     
            PlotPgIKu <- c(PlotPgIKu, 0);
        } 
    }
    

    PlotPgIKu <- PlotPgIKu[2:length(PlotPgIKu)];

    ylo <- 0;
    yhigh <- 1;

    if (OutToScreen) {

#        x11(width=PlotWidth, height=PlotHeight);
#        par(bg="black", fg="white", col.axis="white", col.lab="white", 
#            col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPgIKu,
             col=Color, lwd=4, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="PgIKu",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-",
                  ClusterEndGene, "\n PKu = ", round(ClusterPKu[1], 5),
                  sep=""), font=2, cex.main=1.5);
        
        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    
        abline(h=1, col="white", lty=1);
        abline(h=0.95, col="white", lty=2);
        abline(h=0.90, col="white", lty=2);
        abline(h=0.85, col="white", lty=2);
        abline(h=0.80, col="white", lty=2);
        abline(h=0.70, col="white", lty=2);
        abline(h=0.60, col="white", lty=2);
        abline(h=0.50, col="white", lty=2);
        abline(h=0.40, col="white", lty=2);
        abline(h=0.30, col="white", lty=2);
        abline(h=0.20, col="white", lty=2);
        abline(h=0.10, col="white", lty=2);
        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
        
        if(Interactive) {
            identify(PlotGenesIndex, PlotPgIKu, PlotGenes, col="white");  
        }
    }


    if (OutToFile) {

#        postscript(file=paste(OutputDirectory, FileTitle, ".eps",
#             sep=""), horizontal=TRUE, onefile=TRUE);
#        par(bg="black", fg="white", col.axis="white", col.lab="white", 
#            col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPgIKu,
             col=Color, lwd=2, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="PgIKu",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-", ClusterEndGene,
                  "\n PKu = ", round(ClusterPKu[1], 5), sep=""), font=2,
                  cex.main=1.5);
        
        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    
        abline(h=1, col="white", lty=1);
        abline(h=0.95, col="white", lty=2);
        abline(h=0.90, col="white", lty=2);
        abline(h=0.85, col="white", lty=2);
        abline(h=0.80, col="white", lty=2);
        abline(h=0.70, col="white", lty=2);
        abline(h=0.60, col="white", lty=2);
        abline(h=0.50, col="white", lty=2);
        abline(h=0.40, col="white", lty=2);
        abline(h=0.30, col="white", lty=2);
        abline(h=0.20, col="white", lty=2);
        abline(h=0.10, col="white", lty=2);
        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
#        dev.off()
    }

} else {

    print(paste("Cluster", ClusterStartGene, "-", ClusterEndGene, "not found.",
                sep=" "));
}


################# Plot PgEKu  ################


if (ClusterFound) {

    ClusterPKu <-zPKu[ClusterFilter];
    ClusterPgEKu <- zPgEKu[ClusterFilter];
    
    PlotPgEKu <- 0;

    for (i in 1:length(PlotGenes)) {
        if (Color[i] == "green") {
            TargetFilter <- as.logical(match(ClusterPresentGenes, PlotGenes[i],
                                                                   nomatch=0));
            PlotPgEKu <- c(PlotPgEKu, ClusterPgEKu[TargetFilter]);

        } else {
     
            PlotPgEKu <- c(PlotPgEKu, 0);
        } 
    }
    

    PlotPgEKu <- PlotPgEKu[2:length(PlotPgEKu)];

    ylo <- 0;
    yhigh <- 0.1;

    if (OutToScreen) {

#        x11(width=PlotWidth, height=PlotHeight);
#        par(bg="black", fg="white", col.axis="white", col.lab="white", 
#            col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPgEKu,
             col=Color, lwd=4, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="PgEKu",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-",
                  ClusterEndGene, "\n PKu = ", round(ClusterPKu[1], 5),
                                                 sep=""), font=2, cex.main=1.5);

        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    
        abline(h=0.10, col="white", lty=1);
        abline(h=0.09, col="white", lty=2);
        abline(h=0.08, col="white", lty=2);
        abline(h=0.07, col="white", lty=2);
        abline(h=0.06, col="white", lty=2);
        abline(h=0.05, col="white", lty=2);
        abline(h=0.04, col="white", lty=2);
        abline(h=0.03, col="white", lty=2);
        abline(h=0.02, col="white", lty=2);
        abline(h=0.01, col="white", lty=2);

        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight/10,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
        
        if(Interactive) {
            identify(PlotGenesIndex, PlotPgEKu, PlotGenes, col="white");  
        }
    }


    if (OutToFile) {

#        postscript(file=paste(OutputDirectory, FileTitle, ".eps",
#             sep=""), horizontal=TRUE, onefile=TRUE);
#        par(bg="black", fg="white", col.axis="white", col.lab="white", 
#            col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPgEKu,
             col=Color, lwd=2, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="PgEKu",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-", ClusterEndGene,
                  "\n PKu = ", round(ClusterPKu[1], 5), sep=""), font=2,
                  cex.main=1.5);

        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    
        abline(h=0.10, col="white", lty=1);
        abline(h=0.09, col="white", lty=2);
        abline(h=0.08, col="white", lty=2);
        abline(h=0.07, col="white", lty=2);
        abline(h=0.06, col="white", lty=2);
        abline(h=0.05, col="white", lty=2);
        abline(h=0.04, col="white", lty=2);
        abline(h=0.03, col="white", lty=2);
        abline(h=0.02, col="white", lty=2);
        abline(h=0.01, col="white", lty=2);

        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight/10,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
#        dev.off()
    }

} else {

    print(paste("Cluster", ClusterStartGene, "-", ClusterEndGene, "not found.",
                sep=" "));
}

############# Plot PgKsub ##############


if (ClusterFound) {

    ClusterPKsub <-zPKsub[ClusterFilter];
    ClusterPgKsub <- zPgKsub[ClusterFilter];
    
    PlotPgKsub <- 0;

    for (i in 1:length(PlotGenes)) {
        if (Color[i] == "green") {
            TargetFilter <- as.logical(match(ClusterPresentGenes, PlotGenes[i],
                                                                   nomatch=0));
            PlotPgKsub <- c(PlotPgKsub, ClusterPgKsub[TargetFilter]);

        } else {
     
            PlotPgKsub <- c(PlotPgKsub, 0);
        } 
    }
    

    PlotPgKsub <- PlotPgKsub[2:length(PlotPgKsub)];

    ylo <- 0;
    yhigh <- 1;

    if (OutToScreen) {

 #       x11(width=PlotWidth, height=PlotHeight);
 #       par(bg="black", fg="white", col.axis="white", col.lab="white", 
 #           col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPgKsub,
             col=Color, lwd=4, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="P(g | Ksub)",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-",
                  ClusterEndGene, "\n PKsub = ", round(ClusterPKsub[1], 5),
                                                 sep=""), font=2, cex.main=1.5);

        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    
        abline(h=1, col="white", lty=1);
        abline(h=0.95, col="white", lty=2);
        abline(h=0.90, col="white", lty=2);
        abline(h=0.85, col="white", lty=2);
        abline(h=0.80, col="white", lty=2);
        abline(h=0.70, col="white", lty=2);
        abline(h=0.60, col="white", lty=2);
        abline(h=0.50, col="white", lty=2);
        abline(h=0.40, col="white", lty=2);
        abline(h=0.30, col="white", lty=2);
        abline(h=0.20, col="white", lty=2);
        abline(h=0.10, col="white", lty=2);
        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
        
        if(Interactive) {
            identify(PlotGenesIndex, PlotPgKsub, PlotGenes, col="white");  
        }
    }


    if (OutToFile) {

#        postscript(file=paste(OutputDirectory, FileTitle, ".eps",
#             sep=""), horizontal=TRUE, onefile=TRUE);
#        par(bg="black", fg="white", col.axis="white", col.lab="white", 
#            col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPgKsub,
             col=Color, lwd=2, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="P(g | Ksub)",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-", ClusterEndGene,
                  "\n PKsub = ", round(ClusterPKsub[1], 5), sep=""), font=2,
                  cex.main=1.5);

        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    
        abline(h=1, col="white", lty=1);
        abline(h=0.95, col="white", lty=2);
        abline(h=0.90, col="white", lty=2);
        abline(h=0.85, col="white", lty=2);
        abline(h=0.80, col="white", lty=2);
        abline(h=0.70, col="white", lty=2);
        abline(h=0.60, col="white", lty=2);
        abline(h=0.50, col="white", lty=2);
        abline(h=0.40, col="white", lty=2);
        abline(h=0.30, col="white", lty=2);
        abline(h=0.20, col="white", lty=2);
        abline(h=0.10, col="white", lty=2);
        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
#        dev.off()
    }

} else {

    print(paste("Cluster", ClusterStartGene, "-", ClusterEndGene, "not found.",
                sep=" "));
}



############# Plot PgKmi ##############


if (ClusterFound) {

    ClusterPKmi <-zPKmi[ClusterFilter];
    ClusterPgKmi <- zPgKmi[ClusterFilter];
    
    PlotPgKmi <- 0;

    for (i in 1:length(PlotGenes)) {
        if (Color[i] == "green") {
            TargetFilter <- as.logical(match(ClusterPresentGenes, PlotGenes[i],
                                                                   nomatch=0));
            PlotPgKmi <- c(PlotPgKmi, ClusterPgKmi[TargetFilter]);

        } else {
     
            PlotPgKmi <- c(PlotPgKmi, 0);
        } 
    }
    

    PlotPgKmi <- PlotPgKmi[2:length(PlotPgKmi)];

    ylo <- 0;
    yhigh <- 1;

    if (OutToScreen) {

#        x11(width=PlotWidth, height=PlotHeight);
#        par(bg="black", fg="white", col.axis="white", col.lab="white", 
#            col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPgKmi,
             col=Color, lwd=4, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="P(g | Kmi)",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-",
                  ClusterEndGene, "\n PKmi = ", round(ClusterPKmi[1], 5),
                                                sep=""), font=2, cex.main=1.5);

        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    
        abline(h=1, col="white", lty=1);
        abline(h=0.95, col="white", lty=2);
        abline(h=0.90, col="white", lty=2);
        abline(h=0.85, col="white", lty=2);
        abline(h=0.80, col="white", lty=2);
        abline(h=0.70, col="white", lty=2);
        abline(h=0.60, col="white", lty=2);
        abline(h=0.50, col="white", lty=2);
        abline(h=0.40, col="white", lty=2);
        abline(h=0.30, col="white", lty=2);
        abline(h=0.20, col="white", lty=2);
        abline(h=0.10, col="white", lty=2);
        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
        
        if(Interactive) {
            identify(PlotGenesIndex, PlotPgKmi, PlotGenes, col="white");  
        }
    }


    if (OutToFile) {

#        postscript(file=paste(OutputDirectory, FileTitle, ".eps",
#             sep=""), horizontal=TRUE, onefile=TRUE);
#        par(bg="black", fg="white", col.axis="white", col.lab="white", 
#            col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPgKmi,
             col=Color, lwd=2, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="P(g | Kmi)",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-", ClusterEndGene,
                  "\n PKmi = ", round(ClusterPKmi[1], 5), sep=""), font=2,
                  cex.main=1.5);
        
        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    
        abline(h=1, col="white", lty=1);
        abline(h=0.95, col="white", lty=2);
        abline(h=0.90, col="white", lty=2);
        abline(h=0.85, col="white", lty=2);
        abline(h=0.80, col="white", lty=2);
        abline(h=0.70, col="white", lty=2);
        abline(h=0.60, col="white", lty=2);
        abline(h=0.50, col="white", lty=2);
        abline(h=0.40, col="white", lty=2);
        abline(h=0.30, col="white", lty=2);
        abline(h=0.20, col="white", lty=2);
        abline(h=0.10, col="white", lty=2);
        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
#        dev.off()
    }

} else {

    print(paste("Cluster", ClusterStartGene, "-", ClusterEndGene, "not found.",
                sep=" "));
}

############# Plot PgKet ##############


if (ClusterFound) {

    ClusterPKet <-zPKet[ClusterFilter];
    ClusterPgKet <- zPgKet[ClusterFilter];
    
    PlotPgKet <- 0;

    for (i in 1:length(PlotGenes)) {
        if (Color[i] == "green") {
            TargetFilter <- as.logical(match(ClusterPresentGenes, PlotGenes[i],
                                                                   nomatch=0));
            PlotPgKet <- c(PlotPgKet, ClusterPgKet[TargetFilter]);

        } else {
     
            PlotPgKet <- c(PlotPgKet, 0);
        } 
    }
    

    PlotPgKet <- PlotPgKet[2:length(PlotPgKet)];

    ylo <- 0;
    yhigh <- 0.125;

    if (OutToScreen) {

#        x11(width=PlotWidth, height=PlotHeight);
#        par(bg="black", fg="white", col.axis="white", col.lab="white", 
#            col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPgKet,
             col=Color, lwd=4, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="P(g | Ket)",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-",
                  ClusterEndGene, "\n PKet = ", round(ClusterPKet[1], 5),
                                                sep=""), font=2, cex.main=1.5);

        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    

        abline(h=1, col="white", lty=1);
        abline(h=0.10, col="white", lty=1);
        abline(h=0.09, col="white", lty=2);
        abline(h=0.08, col="white", lty=2);
        abline(h=0.07, col="white", lty=2);
        abline(h=0.06, col="white", lty=2);
        abline(h=0.05, col="white", lty=2);
        abline(h=0.04, col="white", lty=2);
        abline(h=0.03, col="white", lty=2);
        abline(h=0.02, col="white", lty=2);
        abline(h=0.01, col="white", lty=2);

        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight/10,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
        
        if(Interactive) {
            identify(PlotGenesIndex, PlotPgKet, PlotGenes, col="white");  
        }
    }


    if (OutToFile) {

#        postscript(file=paste(OutputDirectory, FileTitle, ".eps",
#             sep=""), horizontal=TRUE, onefile=TRUE);
#        par(bg="black", fg="white", col.axis="white", col.lab="white", 
#            col.main="white");

        PlotGenesIndex <- 1:length(PlotGenes);
        plot(PlotGenesIndex, PlotPgKet,
             col=Color, lwd=2, type='h', bty="n", xlab="Genes", xaxt='n',
             xaxp=c(1,length(PlotGenes),4), ylab="P(g | Ket)",
             ylim=c(ylo, yhigh), col.lab="white", cex.lab=1.5);
        title(main=paste("Cluster ", ClusterStartGene, "-", ClusterEndGene,
                  "\n PKet = ", round(ClusterPKet[1], 5), sep=""), font=2,
                  cex.main=1.5);
        if (length(PlotGenes) > 25) {
            Interval <- trunc(length(PlotGenes)/5);
            PlotLabels <- PlotGenes[c(1,1:4*Interval, length(PlotGenes))];
            axis(1, at=c(1, 1:4*Interval, length(PlotGenes)),
                 labels=PlotLabels);

        } else {
            Interval <- length(PlotGenes);
            PlotLabels <- PlotGenes[c(1:length(PlotGenes))];
            axis(1, at=c(1:length(PlotGenes)), labels=PlotLabels);
        }
  
        abline(h=0, col="white");    
        abline(h=0.10, col="white", lty=1);
        abline(h=0.09, col="white", lty=2);
        abline(h=0.08, col="white", lty=2);
        abline(h=0.07, col="white", lty=2);
        abline(h=0.06, col="white", lty=2);
        abline(h=0.05, col="white", lty=2);
        abline(h=0.04, col="white", lty=2);
        abline(h=0.03, col="white", lty=2);
        abline(h=0.02, col="white", lty=2);
        abline(h=0.01, col="white", lty=2);

        
        if (PlotCluster) {    
            rect(ClGenomeStartPos-PlotStartPos+1,
                 -RectHeight/10,
                 ClGenomeEndPos-PlotStartPos+1,
                 0,
                 border="yellow", col="yellow");
        }

        PlotPos <- 1:length(PlotGenes);
        AbsentGenesPos <- PlotPos[as.logical(match(Color, "red", nomatch=0))];
        AbsentGeneValues <- rep(0, length(AbsentGenesPos));
        points(AbsentGenesPos, AbsentGeneValues, pch=4, lwd=2, cex=2, col="red")

        NonClusterGenesPos <- PlotPos[as.logical(match(Color, "blue",
                                                                nomatch=0))];
        NonClusterGeneValues <- rep(0, length(NonClusterGenesPos));
        points(NonClusterGenesPos, NonClusterGeneValues, pch=4, lwd=2, cex=2,
               col="blue")
  
        dev.off()
    }

} else {

    print(paste("Cluster", ClusterStartGene, "-", ClusterEndGene, "not found.",
                sep=" "));
}

#rm(list=ls());

