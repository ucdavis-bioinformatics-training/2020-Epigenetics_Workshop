
#lets install to current working directory and then load

new_rlib = file.path("/share/workshop/epigenetics_workshop", Sys.getenv("USER"),"r_lib")
suppressWarnings(dir.create(new_rlib,  recursive=T))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='http://cran.us.r-project.org', lib=new_rlib)

BiocManager::install("ChIPQC", lib=new_rlib)

require(DiffBind, lib.loc=new_rlib)
require(ChIPQC, lib.loc=new_rlib)

BiocManager::install("ATACseqQC", lib=new_rlib)

require(ATACseqQC, lib.loc=new_rlib)

BiocManager::install(c("ChIPpeakAnno", "MotifDb", "GenomicAlignments",
           "BSgenome.Mmusculus.UCSC.mm10", "EpiTxDb.Mm.mm10"), lib=new_rlib)

BiocManager::install("Rsamtools", lib=new_rlib)
BiocManager::install("ChIPpeakAnno", lib=new_rlib)
BiocManager::install("EnsDb.Mmusculus.v79", lib=new_rlib)
BiocManager::install("org.Mm.eg.db", lib=new_rlib)


BiocManager::install("edgeR", lib=new_rlib)
BiocManager::install("limma", lib=new_rlib)
BiocManager::install("stringr", lib=new_rlib)

BiocManager::install("ChIPseeker", lib=new_rlib)
