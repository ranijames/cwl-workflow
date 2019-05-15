#!/bin/bash

if [ -z "$1" ]
then
	echo "Usage: sh `basename $0` CONFIG"
    exit 2
fi
config="$1"

# get directory of current script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
outdir=$(grep -A1 "\[TEMP-OUT\]" $config | tail -n 1)

### set paths / cancer type in hypoxia.yml
ct=$(grep -A1 -e "\[CANCER-TYPE\]" $config | tail -n1 | sed -e "s/\//\\\\\//g")
sample_id=$(grep -A1 -e "\[SAMPLES\]" $config | tail -n1 | sed -e "s/\//\\\\\//g")
outdir_sub=$(grep -A1 -e "\[TEMP-OUT\]" $config | tail -n 1 | sed -e "s/\//\\\\\//g")
cat ${scriptdir}/../templates/hypoxia.yml | sed -e "s/OUTPUT_DIR/${outdir_sub}/g" -e "s/CANCER_TYPE/${ct}/g" -e "s/SAMPLES/${sample_id}/g" > ${scriptdir}/../ymlmaker/hypoxia.yml

# run job on cluster
echo "source activate /cluster/home/aalva/software/anaconda/envs/py2/; cd ${outdir}; module load r/3.5.1; 
cwltool --preserve-entire-environment ${scriptdir}/../cwltools/hypoxia.cwl ${scriptdir}/../ymlmaker/hypoxia.yml" |bsub -M 40000 -W 12:00 -R "rusage[mem=40000]" -J hypoxia -o ${outdir}/hypoxia.job
