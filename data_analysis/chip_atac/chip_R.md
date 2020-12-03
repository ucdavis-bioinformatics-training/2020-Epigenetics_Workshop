
cd /share/workshop/epigenetics_workshop/msettles/chipseq_example/
mkdir 05-DiffBind
cd 05-DiffBind

module load R
R


new_rlib = file.path("/share/workshop/epigenetics_workshop", Sys.getenv("USER"),"r_lib")
require(DiffBind, lib.loc=new_rlib)
require(ChIPQC, lib.loc=new_rlib)

#### Load the library

new_rlib = file.path("/share/workshop/epigenetics_workshop", "msettles","r_lib")
require(DiffBind, lib.loc=new_rlib)
require(ChIPQC, lib.loc=new_rlib)


### DiffBind sample samplesheet

download.file("https://ucdavis-bioinformatics-training.github.io/2020-Epigenetics_Workshop/data_analysis/chip_atac/samplesheets/chipseq_diffbind.tsv","chipseq_diffbind.tsv")
samplesheet = read.table("chipseq_diffbind.tsv",sep="\t", header=T,as.is=T)

### Run the ChIPQC Report

experiment = ChIPQC(samplesheet, annotation="mm10")
experiment
ChIPQCreport(experiment)
ChIPQCreport(experiment, reportName="ChIP QC report: Epigenetics Workshop", reportFolder="ChIPseq_QCreport")
