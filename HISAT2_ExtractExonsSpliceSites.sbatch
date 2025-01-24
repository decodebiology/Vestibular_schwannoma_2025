#!/bin/bash -l
#SBATCH -A snic2021-22-212
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 05:00:00
#SBATCH -J hg38ExonsSpliceSites_hisat2
#SBATCH --mail-user=santhilal.subhash@ki.se

module load bioinfo-tools
module load HISAT2/2.2.1


GTF="/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/reference/annotation/gencode.v38.basic.annotation.gtf";

echo "### Started extracting exons `date`";
hisat2_extract_exons.py $GTF > /proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/reference/genome/HISAT2_index/gencode.v38.basic.annotation.Exons.txt ;
echo "### Started extracting splice sites `date`";
hisat2_extract_splice_sites.py $GTF > /proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/reference/genome/HISAT2_index/gencode.v38.basic.annotation.SpliceSites.txt ;

wait
echo "### Completed `date`";
