PlotClusterCut <- 0.05;   #Must be <= 0.15
PlotGeneCut <- 0.05;      #Must be <= 0.15

PlotStartGene <- "SPy0002";
PlotEndGene <- "SPy2217";

PlotTitle <- "";

OutToScreen <- TRUE;
Interactive <- FALSE;

OutToFile <- FALSE;
FileTitle <- "Spyer";


FoldGenes <- TRUE;      #Green Bars
ExpGenes <- TRUE;       #Blue Bars
BGGenes <- TRUE;        #Grey bars
PlotCluster <- TRUE;    #Yellow Rectangles
RectHeight <- 0.2;

PlotWidth <- 11;
PlotHeight <- 8;

            ################# INPUT / OUTPUT PATHS #################

CrawlerDirectory <- "/home/user/SF370Analysis/Data/GenomeCrawlerSuite/";

CrawlerClusterFile <- "CuratedClusters.BayCl";
PlotDataFile <- "Gene.BayPDM";       
CrawlerHeaderLength <- 9;                #DO NOT EDIT  

OutputDirectory <- "/home/user/SF370Analysis/Data/GenomeCrawlerSuite/";
