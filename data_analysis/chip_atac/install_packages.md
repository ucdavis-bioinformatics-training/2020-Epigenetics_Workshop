# Install the needed R packages


## Setting up a custom library (Optional)

Sometimes its nice (or necessary) to have a local repository for all your packages associated with a project. OR your on a cluster where your home directory isn't very accessible.

```r
.libPaths()
new_rlib = file.path("/share/workshop/epigenetics_workshop", Sys.getenv("USER"),"r_lib")
suppressWarnings(dir.create(new_rlib,  recursive=T))
.libPaths(new_rlib)
.libPaths()
```

## Installing R packages

Then lets install a bunch of packages that can be used ChIP and ATAC Peak analysis.

In the R console run the following commands:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='http://cran.us.r-project.org')

BiocManager::install("ChIPQC")
BiocManager::install("ATACseqQC")

BiocManager::install(c("ChIPpeakAnno", "MotifDb", "GenomicAlignments",
           "BSgenome.Mmusculus.UCSC.mm10", "EpiTxDb.Mm.mm10"))

BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("stringr")

BiocManager::install("Rsamtools")
BiocManager::install("ChIPpeakAnno")
BiocManager::install("EnsDb.Mmusculus.v79")
BiocManager::install("org.Mm.eg.db")
```
