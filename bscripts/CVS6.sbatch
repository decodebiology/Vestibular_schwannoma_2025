#!/bin/bash -l
#SBATCH -A snic2021-22-212
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 01:30:00
#SBATCH -J CVS6:2_highrisk
#SBATCH --mail-user=santhilal.subhash@gu.se
echo Aligning HISAT2 CVS6:2_highrisk `date`;
/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/scripts/run_hisat2PE.sh CVS6 ;

echo Alignment completed CVS6:2_highrisk `date`;
wait
