
# Differntial ChIPseq Peak Analysis

This document assumes [calling of peaks with MACS2](./04-peaks.md) has been completed.

**IF** for some reason it didn't finish, is corrupted or you missed the session, you can copy over a completed copy

```bash
#cp -r /share/biocore/workshops/2020_Epigenetics/ChIPseq/03-Filter /share/workshop/epigenetics_workshop/$USER/chipseq_example/.
#cp -r /share/biocore/workshops/2020_Epigenetics/ATACseq/03-Filter /share/workshop/epigenetics_workshop/$USER/atacseq_example/.
```

### First lets set up our environment and start R

```bash
cd /share/workshop/epigenetics_workshop/msettles/chipseq_example/
mkdir 05-DiffBind
cd 05-DiffBind

module load R
R
```

#### Using a custom library path

we are going to use my library to make sure everyone's env works, but later you can [install your own packages](install_packages).

```r
new_rlib = file.path("/share/workshop/epigenetics_workshop", "msettles","r_lib")

# To use your own
#new_rlib = file.path("/share/workshop/epigenetics_workshop", Sys.getenv("USER") ,"r_lib")

.libPaths("/share/workshop/epigenetics_workshop/msettles/r_lib")
```

And then load your libraries
```r
require(DiffBind)
require(ChIPQC)
library(edgeR)
library(Rsamtools)
library(stringr)
library(ChIPpeakAnno)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
```

### The DiffBind Sample Samplesheet

The DiffBind sample sheet allows you to specify the data and metadata informaton for each samples. The available columns are

* SampleID: Identifier string for sample. Must be unique for each sample.
* Tissue: Identifier string for tissue type
* Factor: Identifier string for factor
* Condition: Identifier string for condition
* Treatment: Identifier string for treatment
* Replicate: Replicate number of sample
* bamReads: file path for bam file containing aligned reads for ChIP sample
* bamControl: file path for bam file containing aligned reads for control sample
* Spikein: file path for bam file containing aligned spike-in reads
* ControlID: Identifier string for control sample
* Peaks: path for file containing peaks for sample. Format determined by PeakCaller field or caller parameter
* PeakCaller: Identifier string for peak caller used. If Peaks is not a bed file, this will determine how the Peaks file is parsed. If missing, will use default peak caller specified in caller parameter. Possible values:
    * “raw”: text file file; peak score is in fourth column
    * “bed”: .bed file; peak score is in fifth column
    * “narrow”: default peak.format: narrowPeaks file
    * “macs”: MACS .xls file
    * “swembl”: SWEMBL .peaks file
    * “bayes”: bayesPeak file
    * “peakset”: peakset written out using pv.writepeakset
    * “fp4”: FindPeaks
* PeakFormat: string indicating format for peak files; see PeakCaller and dba.peakset
* ScoreCol: column in peak files that contains peak scores
* LowerBetter: logical indicating that lower scores signify better peaks
* Counts: file path for externally computed read counts; see dba.peakset (counts parameter)

We are going to use *12* columns even though we don't need all of them. I've got one created all aready, lets load it into R and then visualize

```r
download.file("https://ucdavis-bioinformatics-training.github.io/2020-Epigenetics_Workshop/data_analysis/chip_atac/samplesheets/chipseq_diffbind.tsv","chipseq_diffbind.tsv")
samplesheet = read.table("chipseq_diffbind.tsv",sep="\t", header=T,as.is=T)

samplesheet
```

### We can next perform some quality checks.

1. We'll use ChIPQC to generate a fancy report.

    We are going to use our samplesheet and use all the chromosomes..

    ```r
    experiment = ChIPQC(samplesheet, annotation="mm10", chromosomes=c(paste0("chr",1:19),"chrX","chrY"))

    ChIPQCreport(experiment, reportName="ChIPQC_Epigenetics_Workshop", reportFolder="ChIPQC_allchr_Epigenetics_Workshop")
    ```

    [My ChIPseq QC Report](ChIPseq_QCreport/ChIPQC_allchr_Epigenetics_Workshop)

    * From Encode
        "Library complexity is measured using the Non-Redundant Fraction (NRF) and PCR Bottlenecking Coefficients 1 and 2, or PBC1 and PBC2. Preferred values are as follows: NRF>0.9, PBC1>0.9, and PBC2>10.""

    *Questions*
    1. What are these values in our data? Are they legitamate?
    2. What might we have to do to accurately determine them?

### Differential Peaks using DiffBind and Limma Voom

As you heard already, we tend to prefer Limma Voom over the other techniques out there, due to its model flexibility and speed.


1. Merge and Count
    First thing we need to do if DiffBind to *'merge'* the peaks in our samples and produce a binding affinity matrix (raw counts of reads to merged peaks). There are many algorithms to choose from, the default here should be sufficient:

        1. Counts first establishes summits as the location of maximum overlapping reads.
        1. Recenters the data to the summit ('merged' peak regions)
        1. And then counts the reads that align to the new peak.

    The peaks in the 'merged' peaks regions may be re-centered and trimmed based on the calculated summits (point of greatest read overlap) to provide more standardized peak intervals.

    ```r
    chip_dba = dba(sampleSheet=samplesheet)
    counts = dba.count(chip_dba)
    ```

    *Question/tasks*
    1. What are the different parameters
    1. Look at the object, what are our FRiP scores.
    1. What are the different elements in the resulting list.
    1. How many 'merged peaks to we have'?
    1. Spend a little time getting to know the object.


pdf("DiffPeakPlots_ChIPseq.pdf")
dba.plotPCA(counts,  attributes=DBA_TREATMENT, label=DBA_ID)
plot(counts)
dev.off()

counts$peaks[[1]]$Reads
counts.table <- do.call("cbind",lapply(counts$peaks, function(x)x$Reads))

pdata <- counts$samples
pdata

colnames(counts.table) <- pdata$SampleID

peak.info <- counts$peaks[[1]][,1:3]
rownames(peak.info) <- gsub(" +", "", (apply(peak.info, 1, paste0, collapse = "_")))
rownames(counts.table) <- rownames(peak.info)

counts$names=rownames(peak.info)

write.table(counts.table,file="counts.raw.chipseq.txt",col.names=T,row.names=T,quote=F,sep="\t")

dgelist <- DGEList(counts.table,samples=pdata)
dgelist <- calcNormFactors(dgelist)

### Filtering

dim(dgelist)
cutoff <- 20
drop <- apply(cpm(dgelist), 1, max) < cutoff
dgefilter <- dgelist[!drop,]
dim(dgefilter)

### MDS plots

pdf("chipseq_mds.pdf")
plotMDS(dgefilter)
dev.off()

# Modelling differences

design = model.matrix( ~ 0 + Treatment, dgefilter$samples)
design
colnames(design) <- sub("Treatment","T",colnames(design))
design

contrasts_limma = makeContrasts(T25_weeks-T10_weeks, levels=design)
contrasts_limma

dgevoom <- voom(dgefilter,design=design)

fit <- lmFit(dgevoom,design=design)
fit2 = contrasts.fit(fit, contrasts_limma)
fit2 <- eBayes(fit2)

de.table = topTable(fit2, adjust="BH", sort.by="P", number=Inf)
coords = str_match(rownames(de.table), "(.+?)_(\\d+)_(\\d+)")
colnames(coords) <- c("id", "chr","start","end")
de.table = cbind(coords, de.table)

table(DE=de.table$adj.P.Val<= 0.05, FCpos=de.table$logFC>0)

#FCpos

#DE      FALSE  TRUE
#FALSE 10918 12346
#TRUE      2     3

write.table (de.table, file=paste(gsub(" ", "_", colnames(contrasts_limma)), ".DE_toptable.txt",sep=""), col.names=T, row.names = F, quote = FALSE, sep = "\t" )


annoData <- toGRanges(EnsDb.Mmusculus.v79, feature="gene")
annoData[1:2]

merged_overlaps <- GRanges(peak.info$Chr,IRanges(peak.info$Start, peak.info$End), mcols=data.frame(peakNames=rownames(peak.info)))

overlaps.anno <- annotatePeakInBatch(merged_overlaps,
                                     AnnotationData=annoData,
                                     output="nearestBiDirectionalPromoters",
                                     bindingRegion=c(-2000, 500))


overlaps.anno <- addGeneIDs(overlaps.anno,
                         "org.Mm.eg.db",
                         IDs2Add = "entrez_id")
head(overlaps.anno)

write.csv(as.data.frame(unname(overlaps.anno)),"anno.csv")                                     
