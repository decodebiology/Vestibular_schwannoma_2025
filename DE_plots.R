export LC_ALL=C
module load EBModules
module load R/4.1.0-foss-2021a

R

library(circlize)
library(ggplot2)
library(gridExtra)
library(ComplexHeatmap)
library(ggrepel)
library(dplyr)
library(patchwork)



proj_path <- "/grid/beyaz/home/subhash/projects/Vestibular_schwannoma/"
DEpath <- proj_path
plot_path <- paste0(proj_path,"Figures/")
#meta_path <- paste0(proj_path,"scripts/RNA/")
#mdata <- read.table(paste0(meta_path,"sample_table_RNA.txt"), header=T, sep="\t")


PlotDEvolcano <- function(infile, lfc, fdr){

DEg<- read.table(paste0(DEpath,infile), header=T, sep="\t", row.names=1)
#rownames(DEg) <- DEg$GeneSymbol

de <- DEg
de$GeneSymbolx <- make.names(de$GeneSymbol,unique=T)
rownames(de) <- de$GeneSymbolx
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > lfc & de$padj < fdr ] <- "UP"
de$diffexpressed[de$log2FoldChange < -lfc & de$padj < fdr ] <- "DOWN"
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- rownames(de[de$diffexpressed != "NO",])


p <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), fill=diffexpressed, label=delabel)) +
        geom_point(aes(shape=diffexpressed, color=diffexpressed, size=diffexpressed)) + 
        theme(axis.text = element_text(size = 12, color="black"), axis.title = element_text(size = 12),legend.spacing.y = unit(0.3, 'cm'),legend.key.size=unit(0.5, 'cm'),legend.key=element_blank(), legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "#2E2E2E", fill=NA, size=0.5)) +
        geom_text_repel(size=2, fontface="bold.italic", box.padding = 0.2, max.overlaps = 15, segment.size = 0.2, min.segment.length = 0.02, segment.color = '#585858', segment.curvature = 0, segment.ncp = 1, segment.angle = 20, arrow=arrow(angle = 20, length = unit(0.01, "inches"), ends = "last", type = "open")) +
        scale_fill_manual(values=c("UP" = "#BB2649", "DOWN" = "#5F04B4","NO" = "#A4A4A4")) +
		scale_color_manual(values=c("UP" = "#BB2649", "DOWN" = "#5F04B4","NO" = "#A4A4A4")) +
		scale_size_manual(values=c("UP" = 2, "DOWN" = 2,"NO" = 1)) +
        geom_vline(xintercept=c(-lfc, lfc), col="#585858", linetype="longdash", size=0.3) +
        geom_hline(yintercept=-log10(fdr), col="#585858", linetype="longdash", size=0.3) + xlab("log2(fold-change)") + ylab("-log10(FDR)")+
		scale_shape_manual(values=c("UP" = 20, "DOWN" = 20,"NO" = 20))#+ xlim(c(-10,10))
		
		return(p)
}

infile <- "CVS_vs_SVS_PCGs_all_stat.txt"
p1 <- PlotDEvolcano(infile,1,0.05)
x11();p1+ggtitle("DE(CVS/SVS)")+plot_layout(ncol=2,nrow=2)


infile <- "CVS_vs_SVS_lncRNAs_all_stat.txt"
p2 <- PlotDEvolcano(infile,1,0.05)
p2+ggtitle("DE(CVS/SVS)")+plot_layout(ncol=2,nrow=2)


pdf(paste0(plot_path,"CVS_vs_SVS_DE_volcano_highlights.pdf"), width=10)

p1+plot_layout(ncol=3,nrow=2)
p2+plot_layout(ncol=3,nrow=2)

dev.off()

lfc <- 1
fdr <- 0.05
infile1 <- "CVS_vs_SVS_PCGs_DE.txt"
infile2 <- "CVS_vs_SVS_lncRNAs_DE.txt"

mdatax <- data.frame(analysis_names = c("CVS1","CVS4","CVS8","SVS1","SVS2","SVS5"), condition = c("CVS","CVS","CVS","SVS","SVS","SVS"))
rownames(mdatax) <- mdatax$analysis_names

PlotDEheatmap <- function(infile, lfc, fdr, percent){
percent <- percent/100;

oname <- gsub("_DE.txt","",infile)
res <- read.table(paste0(DEpath,infile), header=T, sep="\t", row.names=1)

DEg <- subset(res, abs(res$log2FoldChange)>lfc & res$padj<fdr )
DEg$GeneSymbolx <- DEg$GeneSymbol
DEg <- DEg[order(DEg$padj),]
DEg$GeneSymbol <- make.names(DEg$GeneSymbol, unique=T)
rownames(DEg) <- DEg$GeneSymbol


ghighlight <- DEg$GeneSymbol[1:round(length(DEg$GeneSymbol)*percent)]

matx <- DEg[,c(12:length(colnames(DEg))-1) ]


mat <- t(apply(matx, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))

col_fun = colorRamp2(c(min(mat), 0, max(mat)), c("#08088A", "white", "#B40431")) #c("#0B3861", "white", "#B40431")
row_ha = rowAnnotation(foo = anno_mark(at = c(which(rownames(mat) %in% ghighlight)), labels = subset(rownames(mat),rownames(mat) %in% ghighlight), labels_gp = gpar(fontsize=8,fontface = "italic")), annotation_name_gp= gpar(fontsize = 2) )


ht1<- Heatmap(mat,column_title=paste0(gsub("_DE_allStats.txt","",infile)), cluster_columns = FALSE, show_column_dend = FALSE, heatmap_legend_param = list(legend_direction = "vertical"), show_row_dend = FALSE, row_names_rot = 90, column_names_rot = -90, right_annotation = row_ha, col=col_fun, row_title_gp=gpar(fontsize=1), column_title_gp=gpar(fontsize=8), column_split=mdatax$condition, show_row_names=F, cluster_rows = TRUE, km=2, name=paste0("Expr"), column_names_gp = gpar(fontsize = 8), row_title=NULL)
gb1 = grid.grabExpr(draw(ht1, heatmap_legend_side = "left", annotation_legend_side = "left"))
return(gb1)

}


pdf(paste0(plot_path,"CVS_vs_SVS_lncRNAs_DE_heatmap.pdf"))#, width=10)
p1 <- PlotDEheatmap(infile2,1,0.05,50)
pushViewport(viewport(width = 0.5, height = 1,angle = 90))
grid.draw(p1)

dev.off()
pdf(paste0(plot_path,"CVS_vs_SVS_PCGs_DE_heatmap.pdf"))#, width=10)
p2 <- PlotDEheatmap(infile1,1,0.05,15)
pushViewport(viewport(width = 0.5, height = 1,angle = 90))
grid.draw(p2)

dev.off()




x <- read.table(paste0(DEpath,"HeLa_IER3sh_vs_HeLa_Control_DE_allStats.txt"), header=T, sep="\t")

subset(x, x$gene_name %in% c("IER3","IER3-AS1","ADAM19"))

y <- read.table(paste0(DEpath,"SH5Y5_IER3sh_vs_SH5Y5_Control_DE_allStats.txt"), header=T, sep="\t")

subset(y, y$gene_name %in% c("IER3","IER3-AS1","ADAM19"))

z <- read.table(paste0(DEpath,"SKNB_IER3shvs_SKNB_Control_DE_allStats.txt"), header=T, sep="\t")

subset(z, z$gene_name %in% c("IER3","IER3-AS1","ADAM19"))


subset(x, x$gene_name %in% c("EGR2","ADAM19","COL6A2"))
subset(y, y$gene_name %in% c("EGR2","ADAM19","COL6A2"))
subset(z, z$gene_name %in% c("EGR2","ADAM19","COL6A2"))


####

x <- read.table(paste0(DEpath,"HeLa_IER3sh_vs_HeLa_Control_DE_allStats.txt"), header=T, sep="\t")

y <- read.table(paste0(DEpath,"SH5Y5_IER3sh_vs_SH5Y5_Control_DE_allStats.txt"), header=T, sep="\t")

z <- read.table(paste0(DEpath,"SKNB_IER3shvs_SKNB_Control_DE_allStats.txt"), header=T, sep="\t")


norm <- read.table(paste0(DEpath,"AllSamples_normCounts.txt"), header=T, sep="\t", row.names=1)
cgenes <- read.table(paste0(DEpath,"common_genes.txt"), header=T, sep="\t", row.names=1)


xc <- merge(cgenes,norm,by="gene_name")
xc <- xc[order(xc$gene_name),]



x1 <- merge(cgenes,x,by="gene_name")
y1 <- merge(cgenes,y,by="gene_name")
z1 <- merge(cgenes,z,by="gene_name")

x1c <- x1[,c(1,4,8)]
y1c <- y1[,c(1,4,8)]
z1c <- z1[,c(1,4,8)]

x1c <- x1c[order(x1c$gene_name),]
y1c <- y1c[order(y1c$gene_name),]
z1c <- z1c[order(z1c$gene_name),]

datax <- merge(x1c,y1c,by="gene_name")
colnames(datax)[2:5] <- c("HeLa_LFC","HeLa_padj","SH5Y5_LFC","SH5Y5_padj")
datay <- merge(datax,z1c,by="gene_name")
colnames(datay)[6:7] <- c("SKNB_LFC","SKNB_padj")
datay <- datay[order(datay$gene_name),]

data <- merge(datay,xc,by="gene_name")

rownames(data) <- data$gene_name
xmat <- data[,c(2,4,6)]

Heatmap(xmat,row_km=4)




pv1 <- cor.test(xmat[,"HeLa_LFC"],xmat[,"SH5Y5_LFC"])[[3]]
mcor1 <- cor.test(xmat[,"HeLa_LFC"],xmat[,"SH5Y5_LFC"])[[4]]

pdf(paste0(plot_path,"DE_common_genes_HeLaVsSH5Y5_cor.pdf"))
par(mfrow=c(2,2))
print(plot(xmat[,"HeLa_LFC"],xmat[,"SH5Y5_LFC"], main=paste0("cor: ",mcor1,"\npv< ",pv1),col="#5882FA", pch=16, cex.main=1))
abline(lm(xmat[,"HeLa_LFC"]~xmat[,"SH5Y5_LFC"]), col="#B40431") # regression line (y~x)
print(lines(lowess(xmat[,"HeLa_LFC"],xmat[,"SH5Y5_LFC"]), col="#0B6121")) # lowess line (x,y)
dev.off()


pv2 <- cor.test(xmat[,"HeLa_LFC"],xmat[,"SKNB_LFC"])[[3]]
mcor2 <- cor.test(xmat[,"HeLa_LFC"],xmat[,"SKNB_LFC"])[[4]]

pdf(paste0(plot_path,"DE_common_genes_HeLaVsSKNB_cor.pdf"))
par(mfrow=c(2,2))
print(plot(xmat[,"HeLa_LFC"],xmat[,"SKNB_LFC"], main=paste0("cor: ",mcor2,"\npv< ",pv2),col="#5882FA", pch=16, cex.main=1))
abline(lm(xmat[,"HeLa_LFC"]~xmat[,"SKNB_LFC"]), col="#B40431") # regression line (y~x)
print(lines(lowess(xmat[,"HeLa_LFC"],xmat[,"SKNB_LFC"]), col="#0B6121")) # lowess line (x,y)
dev.off()


pv3 <- cor.test(xmat[,"SH5Y5_LFC"],xmat[,"SKNB_LFC"])[[3]]
mcor3 <- cor.test(xmat[,"SH5Y5_LFC"],xmat[,"SKNB_LFC"])[[4]]

pdf(paste0(plot_path,"DE_common_genes_SH5Y5_LFCVsSKNB_cor.pdf"))
par(mfrow=c(2,2))
print(plot(xmat[,"SH5Y5_LFC"],xmat[,"SKNB_LFC"], main=paste0("cor: ",mcor3,"\npv< ",pv3),col="#5882FA", pch=16, cex.main=1))
abline(lm(xmat[,"SH5Y5_LFC"]~xmat[,"SKNB_LFC"]), col="#B40431") # regression line (y~x)
print(lines(lowess(xmat[,"SH5Y5_LFC"],xmat[,"SKNB_LFC"]), col="#0B6121")) # lowess line (x,y)
dev.off()


dataPCG <- subset(data,data$gene_type=="protein_coding")
dataNCG <- subset(data,!(data$gene_type=="protein_coding"))
xmatPCG <- dataPCG[,c(2,4,6)]
xmatNCG <- dataNCG[,c(2,4,6)]

x11();Heatmap(xmatPCG)
x11();Heatmap(xmatNCG)


datax <- subset(data, (data$HeLa_LFC>0 & data$SH5Y5_LFC<0 & SKNB_LFC<0) | (data$HeLa_LFC<0 & data$SH5Y5_LFC>0 & SKNB_LFC>0))
datay <- subset(data, (data$HeLa_LFC>0 & data$SH5Y5_LFC>0 & SKNB_LFC>0) | (data$HeLa_LFC<0 & data$SH5Y5_LFC<0 & SKNB_LFC<0))

dataPCG <- subset(datax,datax$gene_type=="protein_coding")
dataNCG <- subset(datax,!(datax$gene_type=="protein_coding"))
xmatPCG <- dataPCG[,c(2,4,6)]
xmatNCG <- dataNCG[,c(2,4,6)]
xmat <- datax[,c(2,4,6)]
x11();Heatmap(xmatPCG)
x11();Heatmap(xmatNCG)

pdf(paste0(plot_path,"DE_common_genes_all3lines_heatmapLFC.pdf"), height=6.5, width=2.5)

Heatmap(xmat, row_names_gp = gpar(fontsize = 8), show_column_dend = FALSE, show_row_dend = FALSE)

dev.off()

ymat <- datay[,c(2,4,6)]

pdf(paste0(plot_path,"DE_common_genes_similarExpression_all3lines_heatmapLFC.pdf"), height=5, width=2.5)

Heatmap(ymat, row_names_gp = gpar(fontsize = 8), show_column_dend = FALSE, show_row_dend = FALSE)

dev.off()


