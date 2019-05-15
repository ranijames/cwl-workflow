
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

# Get directory of current script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Creating YML files
python ${scriptdir}/../ymlmaker/YMLmaker.py ${yml} ${fastqlocal} ${temp_out} ${dropbox} ${cancertype} ${config}

# Running the pipeline in cluster
# if not a finished job for each sample_id then submit that job to the cluster
for file in $(find ${yml} -name \*.yml);do
    filebase=$(basename $file)
    echo "source activate tupro; cd ${temp_out}; module load r/3.5.1 libxml2/2.9.4 curl; chmod +R 2770 . ;cwltool --preserve-entire-environment ${scriptdir}/../cwltools/pipeline_part1.cwl $file"| bsub -n 12 -M 84000 -W 12:00 -R "rusage[mem=7000]" -J ${filebase%.yml} -o ${file%.yml}.job
done
