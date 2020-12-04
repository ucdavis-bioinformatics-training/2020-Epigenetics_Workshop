# Calling Peaks using MACS2

This document assumes [calling of peaks with MACS2](./04-peaks.md) has been completed.

**IF** for some reason it didn't finish, is corrupted or you missed the session, you can copy over a completed copy

```bash
#cp -r /share/biocore/workshops/2020_Epigenetics/ChIPseq/03-Filter /share/workshop/epigenetics_workshop/$USER/chipseq_example/.
#cp -r /share/biocore/workshops/2020_Epigenetics/ATACseq/03-Filter /share/workshop/epigenetics_workshop/$USER/atacseq_example/.
```

### First lets set up our environment and start R

```bash
cd /share/workshop/epigenetics_workshop/msettles/atacseq_example/
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
require(ATACseqQC)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
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

We are going to use *12* columns even though we don't need all of them. I've got one created all aready, lets load it into R and then visualize it.

```r
download.file("https://ucdavis-bioinformatics-training.github.io/2020-Epigenetics_Workshop/data_analysis/chip_atac/samplesheets/atacseq_diffbind.tsv","atacseq_diffbind.tsv")
samplesheet = read.table("atacseq_diffbind.tsv",sep="\t", header=T,as.is=T)

samplesheet
```

### WE can next perform some quality checks.

1. First we need to read in the bam files.

    We are going to use our samplesheet to get our bam file paths and then the readGAlignments function to read them in. For the sake of speak we'l only read in "chr1" for our QCs.

    ```r
    bams <- samplesheet$bamReads
    bams.labels <- samplesheet$SampleID

    seqlev <- "chr1" ## subsample data for quick run
    which <- as(seqinfo(Mmusculus)[seqlev], "GRanges")

    gals <- lapply(bams, function(bamfile){
         param <- ScanBamParam( what=c("qname", "flag", "mapq", "isize", "seq", "qual", "mrnm"),  
                                which=which,
                                flag=scanBamFlag(isSecondaryAlignment = FALSE,
                                                  isUnmappedQuery=FALSE,
                                                  isNotPassingQualityControls = FALSE,
                                                  isSupplementaryAlignment = FALSE))
         gal <- readGAlignments(bamfile, param=param)
         header <- scanBamHeader(bamfile, what="text")[[1]]$text
         header <- lapply(header, paste, collapse="\t")
         header <- paste(names(header), unlist(header), sep="\t")
         metadata(gal)$header <- header
         gal
    })
    gals <- GAlignmentsList(gals)
    ```

### Our first QC plots

```r
txs <- transcripts(TxDb.Mmusculus.UCSC.mm9.knownGene)

pt <- PTscore(gals[[1]], txs)
pdf("Promoter-Transcript_body_score.pdf")
plot(pt$log2meanCoverage, pt$PT_score,
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")
dev.off()
```

#### Nucleosome Free Regions (NFR) score

NFR score is a ratio between ATAC cut signal adjacent to TSS and that flanking the corresponding TSS. Each TSS window of 400 bp is first divided into 3 sub-regions: the most upstream 150 bp (n1), the most downstream of 150 bp (n2), and the middle 100 bp (nf). Then the number of fragments with 5’ ends overlapping each region are calculated for each TSS. The NFR score for each TSS is calculated as NFR-score = log2(nf) - log2((n1+n2)/2). A plot can be generated with the NFR scores as Y-axis and the average signals of 400 bp window as X-axis, very like a MA plot for gene expression data.

```r
nfr <- NFRscore(gals[[1]], txs)
pdf("Nucleosome_Free_Regions_score.pdf")
plot(nfr$log2meanCoverage, nfr$NFR_score,
      xlab="log2 mean coverage",
      ylab="Nucleosome Free Regions score",
      main="NFRscore for 200bp flanking TSSs",
      xlim=c(-10, 0), ylim=c(-5, 5))
dev.off()
```

<object data="./ATAC-05-DiffBind/Nucleosome_Free_Regions_score.pdf" type="application/pdf" width="60%"">
    <embed src="http://yoursite.com/the.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="http://yoursite.com/the.pdf">Download PDF</a>.</p>
    </embed>
</object>


```r
tsse <- TSSEscore(gals[[1]], txs)
tsse$TSSEscore
pdf("TSSEscore.pdf")
plot(100*(-9:10-.5), tsse$values, type="b",
     xlab="distance to TSS",
     ylab="aggregate TSS score")
dev.off()
```

### DiffBind and Limma Voom

atac_dba = dba(sampleSheet=samplesheet)

counts = dba.count(atac_dba)

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

### Filtering

dgelist <- DGEList(counts.table,samples=pdata)
dgelist <- calcNormFactors(dgelist)

dim(dgelist)
cutoff <- 20
drop <- apply(cpm(dgelist), 1, max) < cutoff
dgefilter <- dgelist[!drop,]
dim(dgefilter)

pdf("atacseq_mds.pdf")
plotMDS(dgefilter)
dev.off()

design = model.matrix( ~ 0 + Treatment, dgefilter$samples)
design
colnames(design) <- sub("Treatment","",colnames(design))
design

contrasts_limma = makeContrasts(single_treatment-control,double_treatment-control, levels=design)
contrasts_limma

dgevoom <- voom(dgefilter,design=design)

fit <- lmFit(dgevoom,design=design)
fit2 = contrasts.fit(fit, contrasts_limma)
fit2 <- eBayes(fit2)

de.table = topTable(fit2, coef=1, adjust="BH", sort.by="P", number=Inf)
coords = str_match(rownames(de.table), "(.+?)_(\\d+)_(\\d+)")
colnames(coords) <- c("id", "chr","start","end")
de.table = data.frame(coords, de.table)

table(DE=de.table$adj.P.Val<= 0.05, FCpos=de.table$logFC>0)

#FCpos
#DE      FALSE  TRUE
#FALSE 10084  7496
#TRUE   6659  4722

write.table (de.table, file=paste(gsub(" ","_","singleVcontrol"),".DE_toptable.txt",sep=""), row.names = F, quote = FALSE, sep = "\t" )

## Now look at the second contrast


#Now do the Same for the second contrast double_treatment vs control

table(DE=de.table$adj.P.Val<= 0.05, FCpos=de.table$logFC>0)

#FCpos
#DE      FALSE TRUE
#FALSE  8484 6719
#TRUE   8298 5460



# Load libraries

annoData <- toGRanges(EnsDb.Mmusculus.v79, feature="gene")
annoData[1:2]

merged_overlaps <- GRanges(peak.info$Chr,IRanges(peak.info$Start, peak.info$End))
mcols(merged_overlaps)$peakNames=rownames(peak.info)
overlaps.anno <- annotatePeakInBatch(merged_overlaps,
                                     AnnotationData=annoData,
                                     output="nearestBiDirectionalPromoters",
                                     bindingRegion=c(-2000, 500))


overlaps.anno <- addGeneIDs(overlaps.anno,
                         "org.Mm.eg.db",
                         IDs2Add = "symbol")
head(overlaps.anno)

write.csv(as.data.frame(unname(overlaps.anno)),"anno.csv")                                     
de.anno <- as.data.frame(mcols(overlaps.anno))[match(de.table$id, mcols(overlaps.anno)$peakNames),]
de.table = data.frame(de.table, de.anno)

write.table (de.table, file=paste(gsub(" ","_","singleVcontrol"),".anno.DE_toptable.txt",sep=""), row.names = F, quote = FALSE, sep = "\t" )
