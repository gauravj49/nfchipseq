cd /home/rad/users/gaurav/projects/workflows/nfchipseq

# Parameters for the script
species="human"
user="rupert"
projName="mondoAKO"
outdir="/media/rad/HDD1/nfchip"
jobdir="/home/rad/users/gaurav/projects/workflows/nfchipseq"

# Parameters to run the pipeline
projDir="${outdir}/${user}/${projName}"

# Run Diffbind
# https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2020/ChIPSeq/scripts/ChIP_Practical3_DiffBind.html
R

# # In case we do not have these packages
# install.packages("tidyverse")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DiffBind")

# load required packages:
library(DiffBind)
library(tidyverse)
library(rtracklayer)

# 1) Step 1: Reading a peakset
# 1.1) Read a csv file
samples <- read.csv("/media/rad/HDD1/nfchip/rupert/mondoAKO/mondoA_MKO_WT_SampleSheet_DiffBind.csv")
# 1.2) Look at the loaded metadata
names(samples)
# 1.3) The peaksets are read in using the following DiffBind function, that will construct a new DBA object from the sample sheet.
mondoAKO <- dba(sampleSheet=samples)

# 2) Step 2: Take a look at what information gets summarized in the dbObj. How many consensus sites were identified for this dataset? Which sample has a disproportionatley larger number of peaks?
mondoAKO
# 4 Samples, 1416 sites in matrix (9744 total):
#            ID Condition Replicate Caller Intervals
# 1  5B4-MYC_R1       MKO         1 narrow      5821
# 2  5B4-MYC_R2       MKO         2 narrow      2300
# 3 Cas9-MYC_R1  MondoAWT         1 narrow      2705
# 4 Cas9-MYC_R2  MondoAWT         2 narrow      1305

# 3) Step 3: Affinity binding matrix
# 3.1) Take the alignment files and compute count information for each of the peaks/regions in the consensus set. In this step, for each of the consensus regions DiffBind takes the number of aligned reads in the ChIP sample and the input sample, to compute a normalized read count for each sample at every potential binding site. The peaks in the consensus peakset may be re-centered and trimmed based on calculating their summits (point of greatest read overlap) in order to provide more standardized peak intervals.
# - Use the dba.count() function with the following additional parameter:
# - bUseSummarizeOverlaps: to use a more standard counting procedure than the built-in one by default.
mondoAKO <- dba.count(mondoAKO, bUseSummarizeOverlaps=TRUE)
# 4 Samples, 1416 sites in matrix:
#            ID Condition Replicate Caller Intervals
# 1  5B4-MYC_R1       MKO         1 counts      1416
# 2  5B4-MYC_R2       MKO         2 counts      1416
# 3 Cas9-MYC_R1  MondoAWT         1 counts      1416
# 4 Cas9-MYC_R2  MondoAWT         2 counts      1416


# 3.2) With this matrix, the samples can be re-clustered using affinity, rather than occupancy data:
cat("\n- draw a Sample correlation plot using all 1416 consensus sites\n")
SCRplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoA_MKO_WT_sampleCorrelation.jpg", sep='') 
jpeg(filename=SCRplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotHeatmap(mondoAKO)
dev.off()

# 3.3) To see how well the samples cluster with one another, we can draw a PCA plot using all 1621 consensus sites.
cat("\n- draw a PCA plot using all 1416 consensus sites\n")
PCAplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoA_MKO_WT_PCA.jpg", sep='') 
jpeg(filename=PCAplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotPCA(mondoAKO,  attributes=DBA_CONDITION, label=DBA_ID)
dev.off()

# 3.4) Open the file from terminal
eog mondoAKO_PCA.jpg
eog mondoAKO_sampleCorrelation.jpg

# 4) Step 4: Differential binding affinity analysis
# 4.1) 	we have to let DiffBind know how we want to group our samples. 
# 			In our case we will group based on condition. 
# 			This is done using the function, as follows:
mondoAKO <- dba.contrast(mondoAKO, categories=DBA_CONDITION, minMembers=2)
# 1 Contrast:
#   Group1 Members1   Group2 Members2
# 1    MKO        2 MondoAWT        2

# 4.2) The main differential analysis function using DESeq2 is invoked as follows:
mondoAKO$config$th <- 0.2
mondoAKO <- dba.analyze(mondoAKO, method=DBA_ALL_METHODS)

dba.show(mondoAKO , bContrasts=T)
#   Group1 Members1   Group2 Members2 DB.edgeR DB.DESeq2
# 1    MKO        2 MondoAWT        2        5         0

# 4.3) plots
PCAplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoAKO_binding_affinity_heatmap.jpg", sep='') 
jpeg(filename=PCAplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotHeatmap(mondoAKO, ColAttributes = DBA_CONDITION, contrast=1, correlations=TRUE, th=1,method=DBA_EDGER)
dev.off()

# 4.4) Other plots:
# 4.4.1) MA plots
PCAplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoAKO_binding_affinity_MA_plot.jpg", sep='') 
jpeg(filename=PCAplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotMA(mondoAKO, method=DBA_EDGER)
dev.off()
# 4.4.2) Volcano plots
PCAplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoAKO_binding_affinity_volcano_plot.jpg", sep='') 
jpeg(filename=PCAplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotVolcano(mondoAKO, method=DBA_EDGER)
dev.off()
# 4.4.3) PCA plots for contrast=2
PCAplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoAKO_binding_affinity_PCA_MKOvsMondoAWT_plot.jpg", sep='') 
jpeg(filename=PCAplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotPCA(mondoAKO, contrast=1, th=1, method=DBA_EDGER)
dev.off()
# 4.4.4) # Boxplots can give us an idea about the read distribution differences between the classes 
# - in our case the two conditions. The first two boxes show distribution of reads over all 
# differentially bound sites; the middle two show differences on those sites where the affinity 
# increases in Responsive and the two boxes on the right show differences where the affinity increases in Resistant samples.
PCAplotFile     <- paste0("/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoAKO_binding_affinity_box_plot.jpg", sep='') 
jpeg(filename=PCAplotFile, height=900, width=1000, bg="white", res=100); par(mar=c(20,4,1,1));
dba.plotBox(mondoAKO, method=DBA_EDGER)
dev.off()

# 5) Report the differentially bound peak regions, identified by DESeq2. 
# These results files contain the genomic coordinates for all consensus sites and 
# statistics for differential enrichment including fold-change, p-value and FDR.
report <- dba.report(mondoAKO, method=DBA_EDGER)
report
# The value columns are described below:
# 	Chr: Chromosome of binding site
# 	Start: Starting base position of binding site
# 	End: End base position of binding site
# 	Conc: mean read concentration over all the samples (the default calculation uses log2 normalized ChIP read counts with control read counts subtracted)
# 	Conc_group1: Group 1 Concentration
# 	Conc_group2: Group 2 Concentration
# 	Fold: Fold difference – mean fold difference of binding affinity of group 1 over group 2 (Conc1 - Conc2)
# 	p-valule and *FDR statistic indicating significance of difference
# 	Before writing to file we need to convert it to a data frame so that genomic coordinates get written as columns and not GRanges.
report.df <- as.data.frame(report)  
write.table(report.df, "/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoA_MKO_WT_report.csv", sep="\t", quote=F, row.names=F)

# 6) Annotate diffbind results
# Soure: https://gist.github.com/slavailn/a3f94608ac51fc2e6d4cdcaae5d95634
library("ChIPpeakAnno")
library("GenomicRanges")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("biomaRt")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(genomation))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(AnnotationHub))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.ensGene))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

# Read in results table
resFile <- "/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoA_MKO_WT_report.csv"
res <- read.csv(resFile, header=T, sep="\t")

# Convert to GRanges
gr <- makeGRangesFromDataFrame(res, ignore.strand = T, seqnames.field = "seqnames", start.field = "start", end.field = "end")

# Give ranges numeric names in order
names(gr) <- c(1:length(gr))

# Annotate regions
txdb     <- TxDb.Hsapiens.UCSC.hg38.knownGene
andb     <- "org.Hs.eg.db"
peakAnno <- annotatePeak(gr, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb=andb)

annFile <- "/media/rad/HDD1/nfchip/rupert/mondoAKO/analysis/dba/mondoA_MKO_WT_report_annotated.csv"
write.table(peakAnno, annFile, sep = "\t", col.names = T, row.names = F, quote = F)
