#!/bin/bash -l

outname=$1
sname=$2


module load bioinfo-tools
module load samtools/1.5
module load HISAT2/2.2.1


inpath="/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/fastq"
outpath="/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/aligned"
ht2_index="/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/reference/genome/HISAT2_index/hg38"
spliceSites="/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/genome/HISAT2_index/gencode.v38.basic.annotation.SpliceSites.tx"

#Strand-specificity RF is equivalent to fr-firststrand (dUTP protocol ENCODE)
hisat2 --dta --time --summary-file $outpath/${outname}_alignment_summary.txt --avoid-pseudogene --known-splicesite-infile $spliceSites --phred33 -p 4 --qc-filter -x $ht2_index -1 $inpath/${outname}_P1.fq.gz -2 $inpath/${outname}_P2.fq.gz | samtools view -bS - | samtools sort -n --threads 4 -O BAM -o $outpath/${outname}_sorted.bam - ;
samtools index $outpath/${outname}_sorted.bam $outpath/${outname}_sorted.bam.bai;

wait