
#lets install to current working directory and then load

.libPaths()
new_rlib = file.path("/share/workshop/epigenetics_workshop", Sys.getenv("USER"),"r_lib")
suppressWarnings(dir.create(new_rlib,  recursive=T))
.libPaths(new_rlib)
.libPaths()

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='http://cran.us.r-project.org')

BiocManager::install("ChIPQC")

require(DiffBind, lib.loc=new_rlib)
require(ChIPQC, lib.loc=new_rlib)

BiocManager::install("ATACseqQC")

require(ATACseqQC, lib.loc=new_rlib)

BiocManager::install(c("ChIPpeakAnno", "MotifDb", "GenomicAlignments",
           "BSgenome.Mmusculus.UCSC.mm10", "EpiTxDb.Mm.mm10"))

BiocManager::install("Rsamtools")
BiocManager::install("ChIPpeakAnno")
BiocManager::install("EnsDb.Mmusculus.v79")
BiocManager::install("org.Mm.eg.db")


BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("stringr")

BiocManager::install("ChIPseeker")
