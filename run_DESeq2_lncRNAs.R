
module load bioinfo-tools
module load R/4.0.0
module load R_packages/4.0.0
R


library("DESeq2")
library(ComplexHeatmap)
library(ggplot2)
library('PCAtools')
library(ggbiplot, lib.loc = "/proj/uppstore2018034/santhilal/Dakshu_data/R")


rpath <- "/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/";
inpath <- paste0(rpath,"quantification/");
outpath <- paste0(rpath,"DE_analysis/");
plot_path <- paste0(rpath,"DE_analysis/DEplots/");
coldata <- read.table(paste0(rpath,"scripts/sample_table.txt"), sep="\t", row.names=1, header=T)


 
#### data loading (LincRNAs)

annot <- read.table(paste0(rpath,"reference/annotation/gencode.v38.long_noncoding_RNAs_annotation.txt"), sep="\t", header=T, row.names=1)
x <- read.table(paste0(inpath,"gencode.v38_lncRNAs_gcounts.txt"), sep="\t", header=T, row.names=1)

#exprCols <- c(1:14)

exprCols <- c(1,4,8,10,11,14)
#exprCols <- c(1,4,7,8,10,11,14)
coldata <- coldata[c(1,4,8,10,11,14),]
#coldata <- coldata[c(1,4,7,8,10,11,14),]
#tmp <- x[apply(x[, exprCols], 1, function(x) {all(x >= 1)}),c(1,exprCols)]
tmp <- x
cts <- tmp[,exprCols]
cts <- as.matrix(cts)
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))



####### Quality control measures using normalized counts (CPM) - LncRNAs

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ sample_group)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds$condition <- factor(dds$sample_group, levels = levels(coldata$sample_group))
dds <- DESeq(dds)
ddsnorm <- as.data.frame(counts(dds, normalized=T))

pdf(paste0(plot_path,"PCA_plot_lncRNAs.pdf"))

my.pca <- prcomp(ddsnorm, center = TRUE,scale. = TRUE)
ggbiplot(my.pca)


p <- pca(ddsnorm, metadata = coldata, removeVar = 0.1)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p, lab = paste0(p$metadata$OriginalIDs), colby = 'sample_group', hline = 0, vline = 0, legendPosition = 'right')
pairsplot(p)
#eigencorplot(p, metavars = c('sample_type','age_group','age','mutation_group','group','general_group','age_by_median'))
horn <- parallelPCA(ddsnorm)
horn$n
elbow <- findElbowPoint(p$variance)
elbow
screeplot(p,components = getComponents(p, 1:20),vline = c(horn$n, elbow)) + geom_label(aes(x = horn$n + 1, y = 50,label = 'Horn\'s', vjust = -1, size = 8)) + geom_label(aes(x = elbow + 1, y = 50,label = 'Elbow method', vjust = -1, size = 8))
which(cumsum(p$variance) > 80)[1]
biplot(p, lab = paste0(p$metadata$condition), colby = 'sample_group', hline = 0, vline = 0, legendPosition = 'right')

dev.off()

pdf(paste0(plot_path,"complete_correlation_plot_lncRNAs.pdf"))

m1y <- t(apply(ddsnorm, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))
top_ha = HeatmapAnnotation(sampleType=coldata$sample_group, Condition=coldata$condition, Group=coldata$group)
Heatmap(m1y, show_row_names=F, show_column_names=F, top_annotation = top_ha)

top_ha = HeatmapAnnotation(sampleType=coldata$sample_group, Condition=coldata$condition, Group=coldata$group)
row_ha = rowAnnotation(sampleType=coldata$sample_group, Condition=coldata$condition, Group=coldata$group)
Heatmap(cor(ddsnorm), show_row_names=F, show_column_names=F, cluster_columns=T, cluster_rows=T, top_annotation = top_ha, left_annotation = row_ha)

dev.off()


##### multi-group comparisons (lncRNAs)

compGroup <- ~ sample_group
mygroup <- "sample_group"

annot <- read.table(paste0(rpath,"reference/annotation/gencode.v38.long_noncoding_RNAs_annotation.txt"), sep="\t", header=T, row.names=1)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = compGroup )
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds$condition <- factor(coldata[,as.character(mygroup)], levels = levels(factor(coldata[,as.character(mygroup)])) )
dds <- DESeq(dds)
#mysamples <- levels(factor(coldata[,as.character(mygroup)]))


res1 <- results(dds, contrast=c(as.character(mygroup),"CVS","SVS"), alpha=0.05 ) #sample/control for log2FC # HGPS vs young
res1Anot <- merge(as.data.frame(res1),annot,by="row.names")
resSig1 <- subset(res1, abs(res1$log2FoldChange) > 1 & res1$padj < 0.05 )

resSig1 <- resSig1[order(resSig1$padj),]
resSig1Anot <- merge(as.data.frame(resSig1),annot,by="row.names")
colnames(resSig1Anot)[1] <- "Geneid"
			

ddsnorm <- as.data.frame(counts(dds, normalized=T))
ddsnorm[,"Geneid"] <- rownames(ddsnorm)
selectedCols <- subset(rownames(coldata),coldata[,as.character(mygroup)]=="CVS" | coldata[,as.character(mygroup)]=="SVS")
ddsnorm <- ddsnorm[,c("Geneid",as.character(selectedCols))]
			
resSig1AnotExp <- merge(resSig1Anot,ddsnorm,by="Geneid")
resSig1AnotExp <- resSig1AnotExp[order(resSig1AnotExp$padj),]



write.table(resSig1AnotExp, paste0(outpath,"CVS_vs_SVS_lncRNAs_DE.txt"), sep="\t", row.names=F, quote=F)
			
write.table(res1Anot, paste0(outpath,"CVS_vs_SVS_lncRNAs_all_stat.txt"), sep="\t", row.names=F, quote=F)
			
			

pdf(paste0(plot_path,"DE_lncRNAs_heatmap.pdf"))

normDE <-resSig1AnotExp[,12:length(colnames(resSig1AnotExp))]
m1y <- t(apply(normDE, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))
rownames(m1y) <- resSig1AnotExp$GeneSymbol
colx <- coldata[colnames(resSig1AnotExp[,12:length(colnames(resSig1AnotExp))]),]
top_ha = HeatmapAnnotation(sampletype=colx$sample_group, condition=colx$condition)
Heatmap(m1y, show_row_names=T, show_column_names=T, top_annotation = top_ha, column_names_gp = grid::gpar(fontsize = 6), row_names_gp = grid::gpar(fontsize = 6))

dev.off()



pdf(paste0(plot_path,"DE_lncRNAs_volcano_new.pdf"))
plot(res1Anot$log2FoldChange,-log10(res1Anot$padj),col="#ADADAE",bg="#ADADAE", main="CVS vs SVS lncRNAs", xlim=c(-10,10))
points(subset(res1Anot$log2FoldChange,res1Anot$padj < 0.05 & res1Anot$log2FoldChange > 1 ),-log10(subset(res1Anot$padj,res1Anot$padj < 0.05 & res1Anot$log2FoldChange > 1 )),col="#EB4F4F")
points(subset(res1Anot$log2FoldChange,res1Anot$padj < 0.05 & res1Anot$log2FoldChange < -1 ),-log10(subset(res1Anot$padj,res1Anot$padj < 0.05 & res1Anot$log2FoldChange < -1 )),col="#1A21AC")

topFilt <- subset(res1Anot,res1Anot$GeneSymbol %in% c("RP5-841K13.1","EGFLAM-AS1","RP11-43F13.3","AC132217.4","RP11-728F11.4","RP11-234O6.2","TENM3-AS1","AC005754.7","RP11-67L3.5","RP11-67L3.7") )

points(topFilt$log2FoldChange, -log10(topFilt$padj),col="black",bg="#DBA901", pch=21)

abline(h=-log10(0.05),v=c(-1,1), col="#ADADAE", lty = 2)
text(topFilt$log2FoldChange, -log10(topFilt$padj), labels=topFilt$GeneSymbol, pos=4,cex=0.5, col="black")
#points(subset(res1Anot$log2FoldChange,res1Anot$GeneSymbol %in% c("RP5-841K13.1","EGFLAM-AS1","RP11-43F13.3","AC132217.4","RP11-728F11.4","RP11-234O6.2","TENM3-AS1","AC005754.7","RP11-67L3.5","RP11-67L3.7") ),-log10(subset(res1Anot$padj,res1Anot$GeneSymbol %in% c("RP5-841K13.1","EGFLAM-AS1","RP11-43F13.3","AC132217.4","RP11-728F11.4","RP11-234O6.2","TENM3-AS1","AC005754.7","RP11-67L3.5","RP11-67L3.7") )),col="black",bg="#DBA901", pch=21)
#abline(h=-log10(0.05),v=c(-1,1), col="#ADADAE", lty = 2)

#topg <- c("RP5-841K13.1","EGFLAM-AS1","RP11-43F13.3","AC132217.4","RP11-728F11.4","RP11-234O6.2","TENM3-AS1","AC005754.7","RP11-67L3.5","RP11-67L3.7")
#text(subset(res1Anot$log2FoldChange,res1Anot$GeneSymbol %in% c("RP5-841K13.1","EGFLAM-AS1","RP11-43F13.3","AC132217.4","RP11-728F11.4","RP11-234O6.2","TENM3-AS1","AC005754.7","RP11-67L3.5","RP11-67L3.7") ),-log10(subset(res1Anot$padj,res1Anot$GeneSymbol %in% c("RP5-841K13.1","EGFLAM-AS1","RP11-43F13.3","AC132217.4","RP11-728F11.4","RP11-234O6.2","TENM3-AS1","AC005754.7","RP11-67L3.5","RP11-67L3.7") )), labels=topg, pos=4,cex=0.5, col="black")

dev.off()




#### data loading (PCGs)

annot <- read.table(paste0(rpath,"reference/annotation/gencode.v38.basic.annotation.txt"), sep="\t", header=T, row.names=1)
x <- read.table(paste0(inpath,"gencode.v38_all_gcounts.txt"), sep="\t", header=T, row.names=1)



#exprCols <- c(1:14)
#exprCols <- c(1,4,7,8,10,11,14)
exprCols <- c(1,4,8,10,11,14)
#coldata <- coldata[c(1,4,7,8,10,11,14),]
coldata <- coldata[c(1,4,8,10,11,14),]
#tmp <- x[apply(x[, exprCols], 1, function(x) {all(x >= 1)}),c(1,exprCols)]
tmp <- x
cts <- tmp[,exprCols]
cts <- as.matrix(cts)
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

cts <- merge(cts, annot, by=0)
colnames(cts)[1] <- "Geneid"
cts <- subset(cts, cts$Class=="protein_coding")
rownames(cts) <- cts$Geneid
#cts <- cts[,c(2:8)]
cts <- cts[,c(2:7)]
####### Quality control measures using normalized counts (CPM) - LncRNAs

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ sample_group)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds$condition <- factor(dds$sample_group, levels = levels(coldata$sample_group))
dds <- DESeq(dds)
ddsnorm <- as.data.frame(counts(dds, normalized=T))

pdf(paste0(plot_path,"PCA_plot_PCGs.pdf"))
p <- pca(ddsnorm, metadata = coldata, removeVar = 0.1)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p, showLoadings = TRUE, lab = NULL)
pairsplot(p)
#eigencorplot(p, metavars = c('sample_type','age_group','age','mutation_group','group','general_group','age_by_median'))
horn <- parallelPCA(ddsnorm)
horn$n
elbow <- findElbowPoint(p$variance)
elbow
screeplot(p,components = getComponents(p, 1:20),vline = c(horn$n, elbow)) + geom_label(aes(x = horn$n + 1, y = 50,label = 'Horn\'s', vjust = -1, size = 8)) + geom_label(aes(x = elbow + 1, y = 50,label = 'Elbow method', vjust = -1, size = 8))
which(cumsum(p$variance) > 80)[1]
biplot(p, lab = paste0(p$metadata$condition), colby = 'sample_group', hline = 0, vline = 0, legendPosition = 'right')

dev.off()

pdf(paste0(plot_path,"complete_correlation_plot_PCGs.pdf"))

m1y <- t(apply(ddsnorm, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))
#top_ha = HeatmapAnnotation(sampleType=coldata$sample_group, Condition=coldata$condition, Group=coldata$group)
#Heatmap(m1y, show_row_names=F, show_column_names=F, top_annotation = top_ha)

top_ha = HeatmapAnnotation(sampleType=coldata$sample_group, Condition=coldata$condition, Group=coldata$group)
row_ha = rowAnnotation(sampleType=coldata$sample_group, Condition=coldata$condition, Group=coldata$group)
Heatmap(cor(ddsnorm), show_row_names=F, show_column_names=F, cluster_columns=T, cluster_rows=T, top_annotation = top_ha, left_annotation = row_ha)

dev.off()


##### multi-group comparisons (PCGs)

compGroup <- ~ sample_group
mygroup <- "sample_group"

annot <- read.table(paste0(rpath,"reference/annotation/gencode.v38.basic.annotation.txt"), sep="\t", header=T, row.names=1)



dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = compGroup )
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds$condition <- factor(coldata[,as.character(mygroup)], levels = levels(factor(coldata[,as.character(mygroup)])) )
dds <- DESeq(dds)
#mysamples <- levels(factor(coldata[,as.character(mygroup)]))


res1 <- results(dds, contrast=c(as.character(mygroup),"CVS","SVS"), alpha=0.05 ) #sample/control for log2FC # HGPS vs young
res1Anot <- merge(as.data.frame(res1),annot,by="row.names")
resSig1 <- subset(res1, abs(res1$log2FoldChange) > 1 & res1$padj < 0.05 )

resSig1 <- resSig1[order(resSig1$padj),]
resSig1Anot <- merge(as.data.frame(resSig1),annot,by="row.names")
colnames(resSig1Anot)[1] <- "Geneid"
			

ddsnorm <- as.data.frame(counts(dds, normalized=T))
ddsnorm[,"Geneid"] <- rownames(ddsnorm)
selectedCols <- subset(rownames(coldata),coldata[,as.character(mygroup)]=="CVS" | coldata[,as.character(mygroup)]=="SVS")
ddsnorm <- ddsnorm[,c("Geneid",as.character(selectedCols))]
			
resSig1AnotExp <- merge(resSig1Anot,ddsnorm,by="Geneid")
resSig1AnotExp <- resSig1AnotExp[order(resSig1AnotExp$padj),]
			
			
write.table(resSig1AnotExp, paste0(outpath,"CVS_vs_SVS_PCGs_DE.txt"), sep="\t", row.names=F, quote=F)
			
write.table(res1Anot, paste0(outpath,"CVS_vs_SVS_PCGs_all_stat.txt"), sep="\t", row.names=F, quote=F)
		



pdf(paste0(plot_path,"DE_PCGs_heatmap.pdf"))

normDE <-resSig1AnotExp[,12:length(colnames(resSig1AnotExp))]
m1y <- t(apply(normDE, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))
rownames(m1y) <- resSig1AnotExp$GeneSymbol
colx <- coldata[colnames(resSig1AnotExp[,12:length(colnames(resSig1AnotExp))]),]
top_ha = HeatmapAnnotation(sampletype=colx$sample_group, condition=colx$condition)
Heatmap(m1y, show_row_names=T, show_column_names=T, top_annotation = top_ha, column_names_gp = grid::gpar(fontsize = 6), row_names_gp = grid::gpar(fontsize = 1))

dev.off()



pdf(paste0(plot_path,"DE_PCGs_volcano_new.pdf"))
plot(res1Anot$log2FoldChange,-log10(res1Anot$padj),col="#ADADAE",bg="#ADADAE", main="CVS vs SVS PCGs")
points(subset(res1Anot$log2FoldChange,res1Anot$padj < 0.05 & res1Anot$log2FoldChange > 1 ),-log10(subset(res1Anot$padj,res1Anot$padj < 0.05 & res1Anot$log2FoldChange > 1 )),col="#EB4F4F")
points(subset(res1Anot$log2FoldChange,res1Anot$padj < 0.05 & res1Anot$log2FoldChange < -1 ),-log10(subset(res1Anot$padj,res1Anot$padj < 0.05 & res1Anot$log2FoldChange < -1 )),col="#1A21AC")

topFilt <- subset(res1Anot,res1Anot$GeneSymbol %in% c("EGFLAM","CYP2E1","MCOLN2","KCNS3","ENPP2","TBC1D10A","DNAI3","ZC3H7B","OLFML3","BEX1","GJA1","APOBEC3C","RANBP1","CALB2","TMEM184B") )

points(topFilt$log2FoldChange,-log10(topFilt$padj),col="black",bg="#DBA901", pch=21)
abline(h=-log10(0.05),v=c(-1,1), col="#ADADAE", lty = 2)

text(topFilt$log2FoldChange, -log10(topFilt$padj), labels=topFilt$GeneSymbol, pos=4,cex=0.5, col="black")

#points(subset(res1Anot$log2FoldChange,res1Anot$GeneSymbol %in% c("EGFLAM","CYP2E1","MCOLN2","KCNS3","ENPP2","TBC1D10A","DNAI3","ZC3H7B","OLFML3","BEX1","GJA1","APOBEC3C","RANBP1","CALB2","TMEM184B") ),-log10(subset(res1Anot$padj,res1Anot$GeneSymbol %in% c("EGFLAM","CYP2E1","MCOLN2","KCNS3","ENPP2","TBC1D10A","DNAI3","ZC3H7B","OLFML3","BEX1","GJA1","APOBEC3C","RANBP1","CALB2","TMEM184B") )),col="black",bg="#DBA901", pch=21)
#abline(h=-log10(0.05),v=c(-1,1), col="#ADADAE", lty = 2)

#topg <- c("EGFLAM","CYP2E1","MCOLN2","KCNS3","ENPP2","TBC1D10A","DNAI3","ZC3H7B","OLFML3","BEX1","GJA1","APOBEC3C","RANBP1","CALB2","TMEM184B")
#text(subset(res1Anot$log2FoldChange,res1Anot$GeneSymbol %in% c("EGFLAM","CYP2E1","MCOLN2","KCNS3","ENPP2","TBC1D10A","DNAI3","ZC3H7B","OLFML3","BEX1","GJA1","APOBEC3C","RANBP1","CALB2","TMEM184B") ),-log10(subset(res1Anot$padj,res1Anot$GeneSymbol %in% c("EGFLAM","CYP2E1","MCOLN2","KCNS3","ENPP2","TBC1D10A","DNAI3","ZC3H7B","OLFML3","BEX1","GJA1","APOBEC3C","RANBP1","CALB2","TMEM184B") )), labels=topg, pos=4,cex=0.5, col="black")

dev.off()







