library("ChIPpeakAnno")
library("GenomicRanges")
library("org.At.tair.db")
library("TxDb.Athaliana.BioMart.plantsmart28")
library("biomaRt")
# Annotate genomic intervals in bed format using ChIPpeakAnno
# This script was written for Arabidopsis
# Has to be changed to accomodate any other organism available through Ensembl Biomart

setwd("comparisons/")
list.files()

annotateDiffBindResult <- function(resFile, outFile) 
{
  # resFile <- path to resFile
  # name of the out file that contains annotated results
  
  # Read in results table
  res <- read.csv(resFile, header=T)
  
  # Convert to GRanges
  gr <- makeGRangesFromDataFrame(res, ignore.strand = T, seqnames.field = "chr",
                                 start.field = "start", end.field = "end")
  
  # Give ranges numeric names in order
  names(gr) <- c(1:length(gr))
  
  # Create GRanges object with annotations from TxDb database
  annoData <- toGRanges(TxDb.Athaliana.BioMart.plantsmart28, feature="gene")
  
  # Annotate granges with the nearest TSS
  annot <- annotatePeakInBatch(gr, 
                               AnnotationData=annoData, 
                               featureType = "TSS",
                               output="nearestLocation",
                               PeakLocForDistance = "start")
  
  # Load mart
  ensembl <- useMart(biomart = "plants_mart",
                     dataset = "athaliana_eg_gene",
                     host = "plants.ensembl.org")
  
  # Add gene information
  annot <- addGeneIDs(annot, mart = ensembl, feature_id_type = "ensembl_gene_id",
                      IDs2Add = c("entrezgene", "tair_symbol", "description"))
  
  write.table(annot, outFile, sep = "\t", col.names = T, 
              row.names = F, quote = F)
}

annotateDiffBindResult(resFile="DiffBind_result.csv", 
                       outFile="DiffBind_result_annot.txt")
