#!/bin/bash -l
#SBATCH -A snic2021-22-212
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 01:30:00
#SBATCH -J SVS3:1_lowrisk
#SBATCH --mail-user=santhilal.subhash@gu.se
echo Aligning HISAT2 SVS3:1_lowrisk `date`;
/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/scripts/run_hisat2PE.sh SVS3 ;

echo Alignment completed SVS3:1_lowrisk `date`;
wait
