library(DiffBind)

## 1. Reading in peaksets ##
samples <- dba(sampleSheet="KKCY.csv")
dba.show(samples)

## 2. Occupancy analysis ##
dba.plotHeatmap(samples)

## 3. Counting reads ##
peakset <- dba.peakset(samples, bRetrieve=T, DataType=DBA_DATA_GRANGES)
samples <- dba.count(samples,  peakset, score= DBA_SCORE_TMM_READS_EFFECTIVE, summits=0, bParalle=True)

dba.show(samples)
dba.save(samples, file='counts')

dba.plotPCA(samples)
dev.copy(png,'plotPCA_count.png')
dev.off()


## 4. Differential binding affinity analysis ##
samples$contrasts=NULL
samples <- dba.contrast(samples, samples$masks$WT, samples$masks$YAC, "WT", "YAC")
samples <- dba.analyze(samples, method=c(DBA_EDGER,DBA_DESEQ2))


## 5. plot and report ##
dba.plotHeatmap(samples , contrast=1, method=DBA_DESEQ2, th=0.1, bUsePval=TRUE)
dev.copy(png,'plotHeatmap_DEM_YAC_WT.png')
dev.off()

samples.DB <- dba.report(samples , method=c(DBA_EDGER,DBA_DESEQ2), th=1, bCounts=TRUE, bNormalized=TRUE, bCalled=TRUE, file="MACS_DiffBind_WT_YAC")


## 6. Annotation ##
library(ChIPseeker) 
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peakAnno <- annotatePeak(peakset, TxDb=txdb, level="gene")
write.table(as.data.frame(peakAnno),file="PeakAnnotation.UCSC.mm10.knownGene.txt")
