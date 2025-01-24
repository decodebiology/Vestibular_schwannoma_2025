
args<-commandArgs(TRUE)

#module load bioinfo-tools
#module load R/3.5.0
#R
#biocLite("ComplexHeatmap", lib="/proj/uppstore2018034/santhilal/R/x86_64-pc-linux-gnu-library")
#biocLite("MethylKit", lib="/proj/uppstore2018034/santhilal/R/x86_64-pc-linux-gnu-library")

#Add new library path
#.libPaths( c( .libPaths(), "/proj/uppstore2018034/santhilal/R/x86_64-pc-linux-gnu-library") )

#check all available library paths 
#.libPaths()

#just keep one library path "/proj/uppstore2018034/santhilal/R/x86_64-pc-linux-gnu-library" and remove remaining

#.libPaths(.libPaths()[3])

#.libPaths()

#library("svMisc")

#expr <- read.table("/proj/uppstore2018034/santhilal/CCM_proj/quantification/ensembl93/R/ensembl_complete_gcounts_new.table.tpm", sep="\t", header=T)


#ENSG00000284642

#ENSG00000161681

#mycorP <- as.numeric(cor.test(t(log2(subset(expr[,c(2:15)]+1,expr$Geneid=="ENSG00000284642"))),t(log2(subset(expr[,c(2:15)]+1,expr$Geneid=="ENSG00000161681"))))[3])
#mycor <- as.numeric(cor.test(t(log2(subset(expr[,c(2:15)]+1,expr$Geneid=="ENSG00000284642"))),t(log2(subset(expr[,c(2:15)]+1,expr$Geneid=="ENSG00000161681"))))[4])


#pdf("corr.pdf")
#plot(t(log2(subset(expr[,c(2:15)]+1,expr$Geneid=="ENSG00000284642"))),t(log2(subset(expr[,c(2:15)]+1,expr$Geneid=="ENSG00000161681"))),xlab=paste0("LncRNA: ",as.character(subset(expr[,"GeneSymbol"],expr$Geneid=="ENSG00000284642"))), ylab=paste0("PCG: ",as.character(subset(expr[,"GeneSymbol"],expr$Geneid=="ENSG00000161681"))), main=paste0("pearson pv<",mycorP,"\n","r=",mycor))


#mycorP <- as.numeric(cor.test(t(log2(subset(expr[,c(2:15)]+1,expr$Geneid=="ENSG00000284642"))),t(log2(subset(expr[,c(2:15)]+1,expr$Geneid=="ENSG00000161681"))),method="spearman")[3])
#mycor <- as.numeric(cor.test(t(log2(subset(expr[,c(2:15)]+1,expr$Geneid=="ENSG00000284642"))),t(log2(subset(expr[,c(2:15)]+1,expr$Geneid=="ENSG00000161681"))),method="spearman")[4])


#plot(t(log2(subset(expr[,c(2:15)]+1,expr$Geneid=="ENSG00000284642"))),t(log2(subset(expr[,c(2:15)]+1,expr$Geneid=="ENSG00000161681"))),xlab=paste0("LncRNA: ",as.character(subset(expr[,"GeneSymbol"],expr$Geneid=="ENSG00000284642"))), ylab=paste0("PCG: ",as.character(subset(expr[,"GeneSymbol"],expr$Geneid=="ENSG00000161681"))), main=paste0("spearman pv<",mycorP,"\n","r=",mycor))


#dev.off()
cat(paste0("Started processing ",args[3],"\n"))

outpath <- "/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/DE_analysis/coexpression/"

lncDEG <- read.table("/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/DE_analysis/SVS_vs_CVS_lncRNAs_DE.txt", sep="\t", header=T)
pcDEG <- read.table("/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/DE_analysis/SVS_vs_CVS_PCGs_DE.txt", sep="\t", header=T)

sink(paste0(outpath,args[3],".txt"))

cat(paste0("lncRNA\tPCG\tcor\tcorP\n"))

sink()

i=0
start <- args[1]
end <- args[2]

for(lnc in lncDEG$Geneid[start:end])
{
	i=i+1
	
	
	#cat(progress(i, progress.bar=TRUE))
	#lnc <- "ENSG00000284642"
	lncsym <- as.character(subset(lncDEG[,"GeneSymbol"],lncDEG$Geneid==lnc))

	for(pc in pcDEG$Geneid)
	{
		pcsym <- as.character(subset(pcDEG[,"GeneSymbol"],pcDEG$Geneid==pc))
		
		
		mycorP <- as.numeric(cor.test(t(log2(subset(lncDEG[,c(12:17)]+1,lncDEG$Geneid==lnc))),t(log2(subset(pcDEG[,c(12:17)]+1,pcDEG$Geneid==pc))),method="spearman")[3])
		mycor <- as.numeric(cor.test(t(log2(subset(lncDEG[,c(12:17)]+1,lncDEG$Geneid==lnc))),t(log2(subset(pcDEG[,c(12:17)]+1,pcDEG$Geneid==pc))),method="spearman")[4])
		if(mycorP<0.05 & mycor> 0.9){
			sink(paste0(outpath,args[3],".txt"), append=TRUE)
		cat(paste0(lnc,"/",lncsym,"\t",pc,"/",pcsym,"\t",mycor,"\t",mycorP,"\n"))
		sink()
		}
	
	}
	
	if(i==length(lncDEG$Geneid[start:end])){cat(paste0("##### Done processing ",args[3],"\n"))}
	
}
cat("\n")

#sink()


#for(i in 0:101)
#{
	
	
#	progress(i)
#	Sys.sleep(0.01)
#	if(i == 101) cat("Done\n")
	
	
#}


