
#!/bin/bash

# If your conda environment has a different name than "tupro", change the
# script below accordingly and update the following line:
#   source activate tupro
# to reflect the name of your own environment

if [ -z "$1" ]
then
	echo "Usage: sh `basename $0` CONFIG"
    exit 2
fi
config="$1"

yml=$(grep -A1 "\[YML-OUT\]" $config | tail -n 1)
fastqlocal=$(grep -A1 "\[FASTQ-LOCAL\]" $config | tail -n 1)
temp_out=$(grep -A1 "\[TEMP-OUT\]" $config | tail -n 1)
dropbox=$(grep -A1 "\[OUTPUT-DROPBOX\]" $config | tail -n 1)
cancertype=$(grep -A1 "\[CANCER-TYPE\]" $config | tail -n 1)
user_name=$(grep -A1 "\[USER\]" $config | tail -n 1)

# Get directory of current script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
user="/cluster/home/${user_name}/R/x86_64-pc-linux-gnu-library/3.5"
line="Final process status is success"

# Creating YML files
python ${scriptdir}/../ymlmaker/YMLmaker.py ${yml} ${fastqlocal} ${temp_out} ${dropbox} ${cancertype} ${config}

for file in $(find ${yml} -name \*.yml);do
  jobfiles=${file%.*}.job
	filebase=$(basename $file)
  if [ "$(grep -s "$line" $jobfiles)" == "" ]; then
	    echo "source activate workflow_cwl_tupro_2019; cd ${temp_out}; module load r/3.5.1 libxml2/2.9.4 curl; export R_LIBS_USER=$user; cwltool --preserve-environment R_LIBS_USER ${scriptdir}/../cwltools/pipeline_until_hypoxia.cwl $file" | bsub -n 12 -M 84000 -W 12:00 -R "rusage[mem=7000]" -J ${filebase%.yml} -o ${file%.yml}.job
			echo "Sample ${filebase%.*} is submitted to cluster"
	else
		 echo "Sample ${filebase%.*} is successfully finished runing the job"
	fi
done
