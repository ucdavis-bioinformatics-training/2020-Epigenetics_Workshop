
#cd /share/workshop/epigenetics_workshop/msettles/chipseq_example/
#mkdir 05-DiffBind
#cd 05-DiffBind

#module load R
#R

#### Load the library

new_rlib = file.path("/share/workshop/epigenetics_workshop", "msettles","r_lib")
.libPaths(new_rlib)
require(DiffBind)
require(ChIPQC)


### DiffBind sample samplesheet

download.file("https://ucdavis-bioinformatics-training.github.io/2020-Epigenetics_Workshop/data_analysis/chip_atac/samplesheets/chipseq_diffbind.tsv","chipseq_diffbind.tsv")
samplesheet = read.table("chipseq_diffbind.tsv",sep="\t", header=T,as.is=T)

### Run the ChIPQC Report

experiment = ChIPQC(samplesheet, annotation="mm10", chromosomes=c(paste0("chr",1:19),"chrX","chrY"))

ChIPQCreport(experiment, reportName="ChIPQC_Epigenetics_Workshop", reportFolder="ChIPseq_QCreport/ChIPQC_allchr_Epigenetics_Workshop.html")


#[ChIPseq_QCreport/]

### DiffBind

chip_dba = dba(sampleSheet=samplesheet)
counts = dba.count(chip_dba)

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

library(edgeR)

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
library(stringr, lib.loc=new_rlib)
coords = str_match(rownames(de.table), "(.+?)_(\\d+)_(\\d+)")
colnames(coords) <- c("id", "chr","start","end")
de.table = cbind(coords, de.table)

table(DE=de.table$adj.P.Val<= 0.05, FCpos=de.table$logFC>0)

#FCpos

#DE      FALSE  TRUE
#FALSE 10918 12346
#TRUE      2     3

write.table (de.table, file=paste(gsub(" ", "_", colnames(contrasts_limma)), ".DE_toptable.txt",sep=""), col.names=T, row.names = F, quote = FALSE, sep = "\t" )


library(ChIPpeakAnno)
library(EnsDb.Mmusculus.v79)
annoData <- toGRanges(EnsDb.Mmusculus.v79, feature="gene")
annoData[1:2]

merged_overlaps <- GRanges(peak.info$Chr,IRanges(peak.info$Start, peak.info$End), mcols=data.frame(peakNames=rownames(peak.info)))

overlaps.anno <- annotatePeakInBatch(merged_overlaps,
                                     AnnotationData=annoData,
                                     output="nearestBiDirectionalPromoters",
                                     bindingRegion=c(-2000, 500))


library(org.Mm.eg.db)
overlaps.anno <- addGeneIDs(overlaps.anno,
                         "org.Mm.eg.db",
                         IDs2Add = "entrez_id")
head(overlaps.anno)

write.csv(as.data.frame(unname(overlaps.anno)),"anno.csv")                                     
