#!/bin/bash

set -e

if [ "$1" == "-h" -o "$1" == "--help" ]
then
    echo "Usage: $0 [CONFIG]"
    echo
    echo "   If no arguments are provided, the default CONFIG file is used"
    exit 1
fi
if [ -z "$1" ]
then
    config=CONFIG
else
    config=${1}
fi

### set path to expression counting tool
exp_path="$(pwd)/expression/count_expression.py"
[[ ! -f ${exp_path} ]] && echo "script count_expression.py not found" && exit 2
exp_path_sub=$(echo $exp_path | sed -e "s/\//\\\\\//g")
cat templates/count_expression.cwl | sed -e "s/COUNT_EXPRESSION_PATH/${exp_path_sub}/g" > cwltools/count_expression.cwl

### set path to libsize calculation
libsize_path="$(pwd)/expression/compute_libsize.py"
[[ ! -f ${libsize_path} ]] && echo "script compute_libsize.py not found" && exit 2
libsize_path_sub=$(echo $libsize_path | sed -e "s/\//\\\\\//g")
cat templates/libsize_calculator.cwl | sed -e "s/COMPUTE_LIBSIZE_PATH/${libsize_path_sub}/g" > cwltools/libsize_calculator.cwl

### set path to hypoxia R script
pwd_sub=$(pwd | sed -e "s/\//\\\\\//g")
cat templates/project_new_samples_HIF.R | sed -e "s/TUPRO_BASEDIR/${pwd_sub}/g" > hypoxia/R_scripts/project_new_samples_HIF.R
cat templates/tcga_reader_functions.R | sed -e "s/TUPRO_BASEDIR/${pwd_sub}/g" > hypoxia/R_scripts/tcga_reader_functions.R

### set path to hypoxia cwl tool
hypoxia_path="$(pwd)/hypoxia/R_scripts/project_new_samples_HIF.R"
[[ ! -f ${hypoxia_path} ]] && echo "script project_new_samples_HIF.R not found" && exit 2
hypoxia_path_sub=$(echo $hypoxia_path | sed -e "s/\//\\\\\//g")
cat templates/hypoxia.cwl | sed -e "s/HYPOXIA_PATH/${hypoxia_path_sub}/g" > cwltools/hypoxia.cwl

### set path to results mover cwltool
results_path="$(pwd)/ymlmaker/results_to_dropbox.py"
[[ ! -f ${results_path} ]] && echo "script results_to_dropbox.py not found" && exit 2
results_path_sub=$(echo $results_path | sed -e "s/\//\\\\\//g")
cat templates/results_mover.cwl | sed -e "s/RESULTS-MOVE/${results_path_sub}/g" > cwltools/results_mover.cwl

### set path to  md5sum cwl tool
mdsum_path="$(pwd)/hypoxia/mdsum.sh"
[[ ! -f ${mdsum_path} ]] && echo "script mdsum.sh not found" && exit 2
mdsum_path_sub=$(echo $mdsum_path | sed -e "s/\//\\\\\//g")
cat templates/mdsum.cwl | sed -e "s/MDSUM_SHELL/${mdsum_path_sub}/g" > cwltools/mdsum.cwl

### install necessary R packages
#conda_env=$(grep -A1 -e "\[CONDA-ENV\]" $config | tail -n1)
#conda_path=$(dirname $(dirname $(which python) ) )
#source activate $conda_env
#CPPFLAGS="${CPPFLAGS} -I${conda_path}/include" LDFLAGS="${LDFLAGS} -L${conda_path}/lib" Rscript templates/setup_R.R
#source deactivate
