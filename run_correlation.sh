#!/bin/bash
#SBATCH -A snic2021-22-212
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 05:00:00
#SBATCH -J VS_correlation
#SBATCH --mail-user=santhilal.subhash@gu.se


module load bioinfo-tools
module load R/3.5.0

echo "#### Started correlation analysis `date`";

Rscript /proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/scripts/correlation.R 1 65 lnc_PC_correlation &



echo "#### Completed correlation analysis `date`";

wait
