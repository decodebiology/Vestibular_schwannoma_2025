
module load bioinfo-tools
module load R/4.0.0
module load R_packages/4.0.0
R



library(ComplexHeatmap)
library("circlize")

inpath <- "/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/DE_analysis/";
outpath <- "/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/DE_analysis/DEplots/";

x <- read.table(paste0(inpath,"SVS_vs_CVS_lncRNAs_DE.txt"), header=T, sep="\t")
y <- read.table(paste0(inpath,"coexpression/lnc_PC_correlation_stats.txt"), header=T, sep="\t")

xy <- merge(x,y,by="Geneid")

xy <- xy[!xy$GeneSymbol.x == "RP5-841K13.1", ] # outlier


pdf(paste0(outpath,"coexp_scatter.pdf"))
#st <- (1/1837)*100
#end <- (1837/1837)*100
plot(xy$PCGs_coexpressed,-log10(xy$padj),col="#ADADAE",bg="#ADADAE",pch=21, ylab=paste0("DE LncRNAs: Enrichment score (-Log10 FDR)"), xlab="Number of coexpressed DE PCGs")

sigsDown <- subset(xy,xy$GeneSymbol.x=="EGFLAM-AS1" | xy$GeneSymbol.x=="RP11-43F13.3" | xy$GeneSymbol.x=="AC132217.4" | xy$GeneSymbol.x=="RP11-728F11.4")
#sigsDown <- subset(xy,xy$GeneSymbol.x=="RP5-841K13.1" | xy$GeneSymbol.x=="TENM3-AS1" | xy$GeneSymbol.x=="ADIRF-AS1" | xy$GeneSymbol.x=="PCA3")
sigsUp <- subset(xy,xy$GeneSymbol.x=="TENM3-AS1" | xy$GeneSymbol.x=="ADIRF-AS1" | xy$GeneSymbol.x=="PCA3" | xy$GeneSymbol.x=="RP11-108K14.12")
#,col="#D318B4"
points(sigsUp$PCGs_coexpressed,-log10(sigsUp$padj),col="black",bg="#EB4F4F",pch=21)
text(sigsUp$PCGs_coexpressed,-log10(sigsUp$padj), labels=sigsUp$GeneSymbol.x, pos=4,cex=0.8)

#,col="#D318B4"
points(sigsDown$PCGs_coexpressed,-log10(sigsDown$padj),col="black",bg="#1A21AC",pch=21)
text(sigsDown$PCGs_coexpressed,-log10(sigsDown$padj), labels=sigsDown$GeneSymbol.x, pos=4,cex=0.8)

abline(v=20, h=-log10(3e-05), col="#FF00BF", lty = 2, lwd=2)

legend("topright",
       legend=c("CVS Up LncRNAs","CVS Down LncRNAs"),
       pch=c(21, 21), pt.bg=c("#EB4F4F", "#1A21AC"))
	   
text(24,7,">20 coexp PCGs",srt=90,pos=2, col="#FF00BF")
text(85,4.7,"FDR < 3e-05",pos=2, col="#FF00BF")

dev.off()
