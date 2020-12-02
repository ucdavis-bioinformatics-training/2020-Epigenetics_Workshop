library(ggplot2)
library(reshape2)

tss.exprs <- read.table("../05-Deeptools/all.matrix_Genes_expressed.0.tab", skip=3, sep="\t", header=F, stringsAsFactors=F)
tss.noexprs <- read.table("../05-Deeptools/all.matrix_Genes_noexpressed.0.tab", skip=3, sep="\t", header=F, stringsAsFactors=F)

tss.exprs <- apply(tss.exprs, 2, function(x){mean(x, na.rm=TRUE)})
tss.noexprs <- apply(tss.noexprs, 2, function(x){mean(x, na.rm=TRUE)})
s1.exprs <- tss.exprs[1:400]
s2.exprs <- tss.exprs[401:800]
s3.exprs <- tss.exprs[801:1200]

s1.noexprs <- tss.noexprs[1:400]
s2.noexprs <- tss.noexprs[401:800]
s3.noexprs <- tss.noexprs[801:1200]

exprs <- (s1.exprs + s2.exprs + s3.exprs)/3
noexprs <- (s1.noexprs + s2.noexprs + s3.noexprs)/3

x <- 1:400 

data <- data.frame(x=x, exprs=exprs, noexprs=noexprs, stringsAsFactors=F)
data <- melt(data, id="x")

p <- ggplot() +
	geom_line(data=data, aes(x=x, y=value, color=variable)) +
	scale_color_manual(values=c("blue", "red"), labels=c("Expressed", "Not expressed"), name=NULL) +
	theme(legend.position=c(0.5, 0.92)) +
	scale_x_continuous(name="Distance to TSS", breaks=c(1, 100, 200, 300, 400), labels=c("2000", "1000", "0", "1000", "2000"), limits=c(-5, 405)) +
	ylab("Mean methylation percentage") 

pdf("Methylation_TSS.pdf")
p
dev.off()

