#!/bin/bash -l
#SBATCH -A snic2021-22-212
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 05:00:00
#SBATCH -J hg38_index_hisat2
#SBATCH -o hg38_index_hisat2.out
#SBATCH -e hg38_index_hisat2.error
#SBATCH --mail-user=santhilal.subhash@gu.se


module load bioinfo-tools
module load HISAT2/2.2.1

# Using splice sites and exons while building genome takes forever (out of memory) - uppmax computation not enough
# --ss /proj/b2010069/santhilal/dbs/hg38/hESC_project/ensembl_genome/HISAT2_index/Homo_sapiens.GRCh38.90.chr.SpliceSites.txt --exon /proj/b2010069/santhilal/dbs/hg38/hESC_project/ensembl_genome/HISAT2_index/Homo_sapiens.GRCh38.90.chr.Exons.txt 
# Alternative approach taken, using --known-splicesite-infile while aligning reads

echo "### Started building index `date`";

hisat2-build -p 8 /proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/reference/genome/GRCh38.primary_assembly.genome.fa /proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/reference/genome/HISAT2_index/hg38 ;

wait
echo "### Completed `date`";