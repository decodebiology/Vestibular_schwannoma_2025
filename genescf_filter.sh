#!/bin/bash

sample=$1;

#CVS_all_PCGs
#CVS_down_PCGs
#CVS_up_PCGs
#PCGs_coexprd_top4_lncRNAs

inpath="/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/DE_analysis/functional_analysis/output";
outpath="/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/DE_analysis/functional_analysis/filtered";

mkdir $outpath/${sample} ;

cat $inpath/${sample}.list/${sample}.list_GO_BP_goa_human_functional_classification.tsv | cut -f1-7 | awk 'BEGIN{FS="\t";OFS="\t"}{if($6<0.05 || $6=="P-value") print $2,$3,$4,$5,$6,$7,$1;}' | sort -g -t$'\t' -k6,6 > $outpath/${sample}/${sample}.list_GO_BP_goa_human_functional_classification_filt.txt ;

cat $inpath/${sample}.list/${sample}.list_GO_CC_goa_human_functional_classification.tsv | cut -f1-7 | awk 'BEGIN{FS="\t";OFS="\t"}{if($6<0.05 || $6=="P-value") print $2,$3,$4,$5,$6,$7,$1;}' | sort -g -t$'\t' -k6,6 > $outpath/${sample}/${sample}.list_GO_CC_goa_human_functional_classification_filt.txt ;

cat $inpath/${sample}.list/${sample}.list_GO_MF_goa_human_functional_classification.tsv | cut -f1-7 | awk 'BEGIN{FS="\t";OFS="\t"}{if($6<0.05 || $6=="P-value") print $2,$3,$4,$5,$6,$7,$1;}' | sort -g -t$'\t' -k6,6 > $outpath/${sample}/${sample}.list_GO_MF_goa_human_functional_classification_filt.txt ;

cat $inpath/${sample}.list/${sample}.list_REACTOME_Hs_functional_classification.tsv | cut -f1-7 | awk 'BEGIN{FS="\t";OFS="\t"}{if($6<0.05 || $6=="P-value") print $2,$3,$4,$5,$6,$7,$1;}' | sort -g -t$'\t' -k6,6 > $outpath/${sample}/${sample}.list_REACTOME_Hs_functional_classification_filt.txt ;


perl /proj/uppstore2018034/santhilal/meena_proj/scripts/mappingIDS.pl $inpath/${sample}.list/${sample}.list_KEGG_hsa_functional_classification.tsv $inpath/${sample}.list/${sample}.list_user_mapped.list | sed 's/Process~name\t\t1/Process~name\tGeneSymbols\tnumGsym/g' | cut -f2 | paste - $inpath/${sample}.list/${sample}.list_KEGG_hsa_functional_classification.tsv > $inpath/${sample}.list/${sample}.list_KEGG_hsa_functional_classification_gsym.tsv ;


cat $inpath/${sample}.list/${sample}.list_KEGG_hsa_functional_classification_gsym.tsv | cut -f1-8 | awk 'BEGIN{FS="\t";OFS="\t"}{if($7<0.05 || $7=="P-value") print $3,$4,$5,$6,$7,$8,$1;}' | sort -g -t$'\t' -k6,6 > $outpath/${sample}/${sample}.list_KEGG_hsa_functional_classification_filt.txt ;


#cat $inpath/${sample}_PCGs_DE.list/${sample}_PCGs_DE.list_MGI_PhenoGenoMap_mouse_functional_classification.tsv | cut -f1-7 | awk 'BEGIN{FS="\t";OFS="\t"}{if($7<0.05 || $7=="Benjamini and Hochberg (FDR)") print $2,$3,$4,$5,$6,$7,$1;}' | sort -g -t$'\t' -k6,6 > $outpath/${sample}/${sample}_PCGs_DE.list_MGI_PhenoGenoMap_mouse_functional_classification_filt.txt ;



