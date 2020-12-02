
#lets install to current working directory and then load
suppressWarnings(dir.create("r_lib"))
new_rlib = file.path(getwd(),"r_lib")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='http://cran.us.r-project.org', lib=new_rlib)

BiocManager::install("ChIPQC", lib=new_rlib)
require(DiffBind, lib.loc=new_rlib)
require(ChIPQC, lib.loc=new_rlib)

BiocManager::install("ATACseqQC", lib=new_rlib)
require(ATACseqQC, lib.loc=new_rlib)

BiocManager::install(c("ChIPpeakAnno", "MotifDb", "GenomicAlignments",
           "BSgenome.Mmusculus.UCSC.mm10", "EpiTxDb.Mm.mm10"), lib=new_rlib)
