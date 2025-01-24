export LC_ALL=C
module load EBModules
module load R/4.1.0-foss-2021a

R

library(ggplot2)
library(ComplexHeatmap)
library(gridExtra)
library(circlize)
library(data.table)
library(RColorBrewer)
library(data.table)
library(corrplot)
library(stringr)
library(gridExtra)
library(ggpubr)
library(scales)
library(viridis)
library(patchwork)


path <- "/grid/beyaz/home/subhash/projects/Vestibular_schwannoma/functional_analysis/filtered/";
outpath <- "/grid/beyaz/home/subhash/projects/Vestibular_schwannoma/Figures/";

dbs <- c("GO_BP:goa_human","KEGG:hsa","REACTOME:Hs","GO_CC:goa_human","GO_MF:goa_human")#,"Msig_C3_CGP:human"

comparisons <- c("PCGs_coexprd")
#filters <- c("FDR0.05","Pvalue0.05")
#gstatus <- c("up","down","all")

#filter <- "Pvalue0.05"
#comparison <- "WT_Colon_vs_KO_Colon"
#status <- "down"
#db <- "GO_BP:mgi"
num <- 20

	
	for(comparison in comparisons){
			
			for(db in dbs){
				
				cat(paste0("##Comparison:",comparison,"\n","##Database:",db,"\n"))
				
				dbx <- unlist(str_split(db, ":", n = Inf, simplify = FALSE))[1]
				org <- unlist(str_split(db, ":", n = Inf, simplify = FALSE))[2]
				
				x <- read.table(paste0(path,"/",comparison,"_top4_lncRNAs/",comparison,"_top4_lncRNAs.list_",dbx,"_",org,"_functional_classification_filt.txt"), header=T, sep="\t")
				
				colnames(x) <- c("Process.name","DE_genes_involved","TotalGenesInProcess","ProcessCovered","P.value","FDR","Genes")
				
				xn <- subset(x,x$DE_genes_involved>=3)
				
				top <- as.data.frame(head(xn[order(xn$P.value),],n=num))
				rownames(top) <- str_replace(top$Process.name,"~","_")
				mprocess <- top
				
				if(dim(mprocess)[1]>1){
				
				p1 <- ggplot( data=mprocess, aes( x=reorder(stringr::str_wrap(rownames(mprocess),45), -log10(mprocess$P.value) ), y=-log10(mprocess$P.value)) )+
				  geom_bar(stat="identity", position=position_dodge(),color="white", fill="#922329")+
				  geom_text(aes(label=mprocess$DE_genes_involved), hjust=0, color="black", position = position_dodge(0.9), size=5)+
				  labs(title=paste0(comparison," ",status, " genes"," (Top ",num," terms)"), x="", y = "Enrichment Score")+
				  theme(axis.text.x=element_text(angle=90,size=10,hjust=1,color="black"),axis.text.y=element_text(color="black",size=8),plot.title = element_text(hjust = 0.5), panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),panel.background = element_blank())+coord_flip()+ylim( c( 0 ,max(-log10(mprocess$P.value) )+1) )
				 
  				p2 <- ggplot( data=mprocess, aes( x=reorder(stringr::str_wrap(rownames(mprocess),45), -log10(mprocess$P.value) ), y=c(0) ) )+
  				  geom_point(stat="identity", position=position_dodge(), aes(color=-log10(mprocess$P.value), size=mprocess$DE_genes_involved))+coord_flip() + scale_color_gradient(low = "#2AA777", high = "#8904B1",na.value = "#4C0B5F", name=c("Pvalue","Genes") )+ theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_text(color="black",size=8) ,plot.title = element_text(hjust = 0.5), panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),panel.background = element_blank())+xlab("")+  guides(size = guide_legend(name="Genes")) + scale_size(range = c(0, 10))
				  
				  
			
				 
				  pdf(paste0(outpath,"indiv/",comparison,"_top4_lncRNAs",dbx,"_",org,"_","_barplot.pdf"), width=8, height=6)
				  print(p1)
				 
				  dev.off()
				  
				  pdf(paste0(outpath,"indiv/",comparison,"_top4_lncRNAs",dbx,"_",org,"_","_bubbleplot.pdf"), width=8, height=6)
				  print(p2)
				 
				  dev.off()
			  }
								
			}
	
	}
	




#### Combined plots (Up+Down)


path <- "/grid/beyaz/home/subhash/projects/Vestibular_schwannoma/functional_analysis/output/";
outpath <- "/grid/beyaz/home/subhash/projects/Vestibular_schwannoma/Figures/";

dbs <- c("GO_BP:goa_human","KEGG:hsa","REACTOME:Hs","GO_CC:goa_human","GO_MF:goa_human")#,"Msig_C3_CGP:human"

comparisons <- c("CVS")
num <- 20




for(comparison in comparisons){
			
			for(db in dbs){
				
				mprocessx <- 0
				mprocessy <- 0
				dbx <- unlist(str_split(db, ":", n = Inf, simplify = FALSE))[1]
				org <- unlist(str_split(db, ":", n = Inf, simplify = FALSE))[2]
				#MatriGel_vs_ObaGel_DE_PCGs_all.list/MatriGel_vs_ObaGel_DE_PCGs_all.list_KEGG_hsa_functional_classification.tsv
				x <- read.table(paste0(path,"/",comparison,"_up_PCGs.list/",comparison,"_up_PCGs.list_",dbx,"_",org,"_functional_classification.tsv"), header=T, sep="\t")
				y <- read.table(paste0(path,"/",comparison,"_down_PCGs.list/",comparison,"_down_PCGs.list_",dbx,"_",org,"_functional_classification.tsv"), header=T, sep="\t")
				x <- x[,2:7]
				y <- y[,2:7]
				colnames(x) <- c("Process.name","num_of_Genes","gene_group","percentage","Pvalue","FDR")
				colnames(y) <- c("Process.name","num_of_Genes","gene_group","percentage","Pvalue","FDR")
				x$Status <- "Up"
				y$Status <- "Down"
				x$DB <- dbx
				y$DB <- dbx
				x <- x[order(x$Pvalue),]
				y <- y[order(y$Pvalue),]
	
				xn <- subset(x,x$num_of_Genes>=3)
				yn <- subset(y,y$num_of_Genes>=3)
				
				topx <- as.data.frame(head(xn[order(xn$Pvalue),],n=num))
				topy <- as.data.frame(head(yn[order(yn$Pvalue),],n=num))
				
				topx$Process.name <- gsub("_"," ",topx$Process.name)
				
				topy$Process.name <- gsub("_"," ",topy$Process.name)
				
				if(dbx=="KEGG" || dbx=="GO_BP" || dbx == "GO_CC" || dbx == "GO_MF"){
					namesx <- unlist(lapply(strsplit(topx$Process.name, '~', fixed = TRUE), '[', 2)) 
					namesy <- unlist(lapply(strsplit(topy$Process.name, '~', fixed = TRUE), '[', 2))
					
					namesx <- make.names(namesx, unique=T, allow_ = T)
					namesy <- make.names(namesy, unique=T, allow_ = T)
					
					rownames(topx) <- stringr::str_trunc( namesx ,100, "center")
					rownames(topy) <- stringr::str_trunc( namesy ,100, "center")
				
				}else{
					
					topx$Process.name <- make.names(topx$Process.name, unique=T, allow_ = T)
					topy$Process.name <- make.names(topy$Process.name, unique=T, allow_ = T)
					
					rownames(topx) <- stringr::str_trunc( topx$Process.name ,100, "center")
					rownames(topy) <- stringr::str_trunc( topy$Process.name ,100, "center")
				
				}

				mprocessx <- topx
				
				mprocessy <- topy
				
				mprocess <- rbind(mprocessx,mprocessy)
				mprocess$Status <- factor(mprocess$Status,levels=c("Up","Down"))
				x <- as.data.frame(mprocess)
				p1 <- ggplot(x,aes(x=Process.name,y=Status,fill=-log10(Pvalue)),group=Status)+geom_point(stat="identity", position=position_dodge(),aes(fill=-log10(Pvalue),size=num_of_Genes),pch=21)+ggtitle(paste0(comparison))+scale_fill_viridis(direction=-1)+xlab("")+facet_wrap(~DB, scales="free_y")+  guides(size = guide_legend(name="Genes")) + scale_size(range = c(0, 10)) +coord_flip()+ theme(axis.title.x=element_blank(),axis.text.y=element_text(color="black",size=8) ,plot.title = element_text(hjust = 0.5), panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),panel.background = element_blank())
				#+facet_wrap(~Status+DB, scales="free_y")
				pdf(paste0(outpath,comparison,"_",dbx,"_pathways.pdf"), width=8, height=15)
				print(p1+plot_layout(ncol=1,nrow=2))
				dev.off()
				
				scale_neg <- 30
				#combined <- mprocess 
				mprocessx$log_fdr <- -log10(mprocessx$Pvalue)
				mprocessy$log_fdr <- log10(mprocessy$Pvalue)
				combined <- rbind(mprocessx,mprocessy)
				combined$log_fdr <- round(as.numeric(combined$log_fdr),1)
				combined <- as.data.frame(combined)
				
				color <- ifelse(combined$log_fdr < 0, alpha("#0431B4",0.5), alpha("#DF013A",0.5))	

				p1 <- ggplot(combined, aes(x = reorder(rownames(combined), log_fdr), y = log_fdr)) +
			  	geom_bar(stat = "identity",
			           show.legend = FALSE,
			           fill = color,     
			           color = "white") +
			  		 geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
			  	   geom_text(aes(label = rownames(combined), # strtrim(rownames(combined),25), Text with groups,
			                hjust = ifelse(log_fdr < 0, 0, 1),
			                vjust = 0.5), size = 4) +#ifelse(log_fdr < 0, 1, 0)
			  			  xlab("") +
			  			ylab("Log(Pvalue)") +
			  #scale_y_continuous(limits = c(min(combined$log_fdr)-scale_neg,  max(combined$log_fdr)+30) ) +
			  #scale_y_continuous(limits = c(min(combined$log_fdr)-100,  max(combined$log_fdr)+100) ) +
			  	coord_flip() +ggtitle(paste0(comparison))+
			 	 theme_minimal() +
			  	geom_hline(yintercept = c(-1.3,1.3), linetype="dotted", 
			                color = "#585858", size=0.5)+
			  			  theme(axis.line.x = element_blank(),plot.title = element_text(hjust = 0.5, size=12),axis.text.y = element_blank(),axis.title=element_text(size=12),axis.text.x = element_text(colour="black", size=12),
			        	  panel.grid.major.y = element_blank(), panel.border = element_rect(colour = "#2E2E2E", fill=NA, linewidth=0.5),axis.ticks.x=element_line(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) 
				
		    			  pdf(paste0(outpath,comparison,"_",dbx,"_pathways_barplot.pdf"), width=15, height=15)
				
						  #pdf(paste0(outpath,pid,"_",db,"_",species,"_EndoNormaVsCancer_top15.pdf"))
				print( p1+plot_layout(ncol=2,nrow=2) )
				dev.off()
				
				

			}
	
}





