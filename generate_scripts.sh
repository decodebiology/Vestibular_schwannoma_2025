#!/bin/bash

#### Submit multipe jobs
#find . -name '*.sbatch' | xargs -I% sbatch %

### Sbatch parameters

processor=$1
numnode=$2
time=$3
project="snic2021-22-212"



SCRIPT_PATH="/proj/uppstore2018034/santhilal/Vestibular_schwannoma_proj/scripts";


SAMPLES=(`cat $SCRIPT_PATH/sample_table.txt | grep -v "^#"| cut -f1,4 | sed 's/\t/;/g'`); 


for i in ${SAMPLES[@]}

do

	IFS=';' read -ra arr <<< "$i"; 

	echo "#!/bin/bash -l" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch

	echo "#SBATCH -A $project" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch
	echo "#SBATCH -p $processor" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch
	echo "#SBATCH -n $numnode" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch
	echo "#SBATCH -t $time" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch
	echo "#SBATCH -J ${arr[0]}:${arr[1]}" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch
	echo "#SBATCH --mail-user=santhilal.subhash@gu.se" >> bscripts/${arr[0]}.sbatch



	library=${arr[1]};

	echo "echo Aligning HISAT2 ${arr[0]}:${arr[1]} \`date\`;" >> bscripts/${arr[0]}.sbatch;
	#### Aligning (If the fastq files are not cleaned, please include adapter cleaning step before aligning)

			echo "$SCRIPT_PATH/run_hisat2PE.sh ${arr[0]} ;" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch;
			echo "" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch;

	echo "echo Alignment completed ${arr[0]}:${arr[1]} \`date\`;" >> bscripts/${arr[0]}.sbatch;
	echo "wait" >> bscripts/${arr[0]}.sbatch;




done

