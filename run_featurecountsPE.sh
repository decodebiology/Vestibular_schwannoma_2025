#!/bin/bash -l
#SBATCH -A snic2021-22-212
#SBATCH -p core
#SBATCH -n 14
#SBATCH -t 04:00:00
#SBATCH -J quantification
#SBATCH -o quantification.out
#SBATCH -e quantification.error
#SBATCH --mail-user=santhilal.subhash@ki.se

module load bioinfo-tools
module load subread/2.0.0

GTF1="/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/reference/annotation/gencode.v38.long_noncoding_RNAs.gtf";
GTF2="/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/reference/annotation/gencode.v38.basic.annotation.gtf";

path="/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/aligned"
outpath="/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/quantification"

echo "### Quantification Started for lncRNAs `date`";

/home/sant/bin/hESC_project/subread-1.5.3-Linux-x86_64/bin/featureCounts -T 14 -p -s 0 -B -C -a $GTF1 -o $outpath/gencode.v38_lncRNAs_gcounts.out -F GTF -t exon -g gene_id --minOverlap 10 -Q 60 --ignoreDup -J $path/CVS1_sorted.bam $path/CVS2_sorted.bam $path/CVS3_sorted.bam $path/CVS4_sorted.bam $path/CVS5_sorted.bam $path/CVS6_sorted.bam $path/CVS7_sorted.bam $path/CVS8_sorted.bam $path/CVS9_sorted.bam $path/SVS1_sorted.bam $path/SVS2_sorted.bam $path/SVS3_sorted.bam $path/SVS4_sorted.bam $path/SVS5_sorted.bam ;

echo "### Quantification Started for all genes `date`";

/home/sant/bin/hESC_project/subread-1.5.3-Linux-x86_64/bin/featureCounts -T 14 -p -s 0 -B -C -a $GTF2 -o $outpath/gencode.v38_all_gcounts.out -F GTF -t exon -g gene_id --minOverlap 10 -Q 60 --ignoreDup -J $path/CVS1_sorted.bam $path/CVS2_sorted.bam $path/CVS3_sorted.bam $path/CVS4_sorted.bam $path/CVS5_sorted.bam $path/CVS6_sorted.bam $path/CVS7_sorted.bam $path/CVS8_sorted.bam $path/CVS9_sorted.bam $path/SVS1_sorted.bam $path/SVS2_sorted.bam $path/SVS3_sorted.bam $path/SVS4_sorted.bam $path/SVS5_sorted.bam ;


echo "### Quantification done `date`";

wait

