#  Differential ATACseq Peak Analysis

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
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
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

    *Questions*
    1. How would you change things to view all chromosomes.


### Our first QC plots

Alot of ATAC-seq QC is dependant on signal around TSS. So we first need to define the location of these sites. We'll use the transcript database (TxDb) from UCSC mm10.

```r
txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
```

#### Promoter/Transcript body (PT) score

PT score is calculated as the coverage of promoter divided by the coverage of its transcript body. PT score will show if the signal is enriched in promoters.


```r
pt <- PTscore(gals[[1]], txs)
pdf("Promoter-Transcript_body_score.pdf")
plot(pt$log2meanCoverage, pt$PT_score,
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")
dev.off()
```

<embed src="./ATAC-05-DiffBind/Promoter-Transcript_body_score.pdf" width="80%" height="80%" frameborder="0" allowfullscreen>

#


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

<embed src="./ATAC-05-DiffBind/Nucleosome_Free_Regions_score.pdf" width="80%" height="80%" frameborder="0" allowfullscreen>


#

### Transcription Start Site (TSS) Enrichment Score

TSS enrichment score is a raio between aggregate distribution of reads centered on TSSs and that flanking the corresponding TSSs. TSS score = the depth of TSS (each 100bp window within 1000 bp each side) / the depth of end flanks (100bp each end). TSSE score = max(mean(TSS score in each window)). TSS enrichment score is calculated according to the definition at https://www.encodeproject.org/data-standards/terms/#enrichment. Transcription start site (TSS) enrichment values are dependent on the reference files used.

```r
tsse <- TSSEscore(gals[[1]], txs)
pdf("TSSEscore.pdf")
plot(100*(-9:10-.5), tsse$values, type="b",
     xlab="distance to TSS",
     ylab="aggregate TSS score")
dev.off()
```

<embed src="./ATAC-05-DiffBind/TSSEscore.pdf" width="80%" height="80%" frameborder="0" allowfullscreen>

#

Encode recommended cutoff values for high quality data are listed below.

*Question*
1. what is the TSSEscore? *hint* its part of the tsse object.

Encode recommends:

|Annotation used | Value| Resulting Data Status |
| :--- | :--- |:--- |
hg19 Refseq TSS annotation|	< 6	|Concerning |
hg19 Refseq TSS annotation| 6 - 10 |	Acceptable |
hg19 Refseq TSS annotation|> 10	| Ideal |
GRCh38 Refseq TSS annotation |	< 5	| Concerning |
GRCh38 Refseq TSS annotation | 5 - 7	| Acceptable |
GRCh38 Refseq TSS annotation |> 7	|Ideal |
mm9 GENCODE TSS annotation |	< 5 |	Concerning |
mm9 GENCODE TSS annotation | 5 - 7 |	Acceptable |
mm9 GENCODE TSS annotation | > 7	| Ideal |
mm10 Refseq TSS annotation |	< 10	| Concerning |
mm10 Refseq TSS annotation | 10 -15	| Acceptable |
mm10 Refseq TSS annotation | > 15	| Ideal |


### Differential Peaks using DiffBind and Limma Voom

As you heard already, we tend to prefer Limma Voom over the other techniques out there, due to its model flexibility and speed.


1. Merge and Count
    First thing we need to do if DiffBind to *'merge'* the peaks in our samples and produce a binding affinity matrix (raw counts of reads to merged peaks). There are many algorithms to choose from, the default here should be sufficient:

        1. Counts first establishes summits as the location of maximum overlapping reads.
        1. Recenters the data to the summit ('merged' peak regions)
        1. And then counts the reads that align to the new peak.

    The peaks in the 'merged' peaks regions may be re-centered and trimmed based on the calculated summits (point of greatest read overlap) to provide more standardized peak intervals.

    ```r
    atac_dba = dba(sampleSheet=samplesheet)
    counts = dba.count(atac_dba)
    ```

    *Question/tasks*
    1. What are the different parameters
    1. Look at the object, what are our FRiP scores.
    1. What are the different elements in the resulting list.
    1. How many 'merged peaks to we have'?
    1. Spend a little time getting to know the object.

    Ideally you want your fraction of reads in called peak regions (FRiP score) to be >0.3, though values greater than 0.2 are acceptable. Note these values are different in that they are using the 'merged' set.

1. We can use DiffBeaks to produce a PCA of samples

    ```r
    pdf("DiffPeakPlots_ATACseq.pdf")
    dba.plotPCA(counts,  attributes=DBA_TREATMENT, label=DBA_ID)
    plot(counts)
    dev.off()
    ```

    <embed src="./ATAC-05-DiffBind/DiffPeakPlots_ATACseq.pdf" width="80%" height="80%" frameborder="0" allowfullscreen>
    #

1. Then we'll extract the reads and build a new counts table for use in other applications (ala Limma).

    ```r
    head(counts$peaks[[1]]$Reads) # what the data look like
    counts.table <- do.call("cbind",lapply(counts$peaks, function(x)x$Reads))

    pdata <- counts$samples
    pdata
    colnames(counts.table) <- pdata$SampleID
    ```

    We'll extract the peak coordinates from the first sample

    ```r
    peak.info <- counts$peaks[[1]][,1:3]
    rownames(peak.info) <- gsub(" +", "", (apply(peak.info, 1, paste0, collapse = "_")))
    rownames(counts.table) <- rownames(peak.info)
    ```

    And now we have an peak abundance table for all samples (counts.table) and experiment metadata table (pdata)

1. Filtering and plotting a MDS

    We need to create our Differential Gene Expression List (DGEList) object that limma uses and Calculate our normalizing factors.

    ```r
    dgelist <- DGEList(counts.table,samples=pdata)
    dgelist <- calcNormFactors(dgelist)

    dim(dgelist)
    ```

    Here we will choose a cutoff of 20 (semi-random choice) to produce a smaller more managable set of peaks. There are strategies using voom and the voom plot to better choose an appropriate cutoff. *BUT* in general we reduce the number of peaks to those that are most likely to be interesting.

    ```r
    cutoff <- 20
    drop <- apply(cpm(dgelist), 1, max) < cutoff
    dgefilter <- dgelist[!drop,]
    dim(dgefilter)

    pdf("atacseq_mds.pdf")
    plotMDS(dgefilter)
    dev.off()
    ```

    <embed src="./ATAC-05-DiffBind/atacseq_mds.pdf" width="80%" height="80%" frameborder="0" allowfullscreen>


1. Finally, lets apply the model and produce our result.

    We use the standard limma modelling after a voom transformation for read count data.

    ```r
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

    write.table (de.table, file=paste(gsub(" ","_","singleVcontrol"),".DE_toptable.txt",sep=""), row.names = F, quote = FALSE, sep = "\t" )
    ```

<!--
#FCpos
#DE      FALSE  TRUE
#FALSE 10084  7496
#TRUE   6659  4722
-->
    *Questions*
    1. How many Differential peaks are there? How many are DE and 'open' in 'single treatment', how many are open in control.
    2. Transfer the DE table to your computer and view it in excel. You can view mine [here](ATAC-05-DiffBind/singleVcontrol.anno.DE_toptable.txt)
    3. Run the test for double_treatment vs control.


<!--
#table(DE=de.table$adj.P.Val<= 0.05, FCpos=de.table$logFC>0)
#FCpos
#DE      FALSE TRUE
#FALSE  8484 6719
#TRUE   8298 5460
-->

### Adding in Annotation using ChIPpeakAnno

We'll use the EnsDb package for Mouse and ChIPpeakAnno. To annotate we merge teh peaks from our data (our peak.info object) with the annoation data and combine with our DE table to produce an annotated to gene DE table.

    ```r
    annoData <- toGRanges(EnsDb.Mmusculus.v79, feature="gene")

    merged_overlaps <- GRanges(peak.info$Chr,IRanges(peak.info$Start, peak.info$End))

    overlaps.anno <- annotatePeakInBatch(merged_overlaps,
                                         AnnotationData=annoData,
                                         output="nearestBiDirectionalPromoters",
                                         bindingRegion=c(-2000, 500))

    # Alot of warnings appear, but none are severe.

    overlaps.anno <- addGeneIDs(overlaps.anno,
                             "org.Mm.eg.db",
                             IDs2Add = "symbol")
    head(overlaps.anno)

    write.csv(as.data.frame(unname(overlaps.anno)),"anno.csv")                                     
    de.anno <- as.data.frame(mcols(overlaps.anno))[match(de.table$id, mcols(overlaps.anno)$peakNames),]
    de.table = data.frame(de.table, de.anno)

    write.table (de.table, file=paste(gsub(" ","_","singleVcontrol"),".anno.DE_toptable.txt",sep=""), row.names = F, quote = FALSE, sep = "\t" )
    ```

*Questions*
1. What does the annotation data look like?
1. What other features can we add (we added Symbol).
1. Copy the new table to your computer and view in excel. My copy is [here](singleVcontrol.anno.DE_toptable.txt)
