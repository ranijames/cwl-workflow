#!/bin/bash

if [ -z "$1" ]
then
	echo "Usage: sh `basename $0` CONFIG"
    exit 2
fi
config="$1"


dropbox=$(grep -A1 "\[OUTPUT-DROPBOX\]" $config | tail -n 1)
tmpout=$(grep -A1 "\[TEMP-OUT\]" $config | tail -n 1)

# get directory of current script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

sh ${scriptdir}/../hypoxia/mdsum.sh ${config}

# Moving the results and md5sum files to their corresponding folders in dropbox

python ${scriptdir}/../ymlmaker/results_to_dropbox.py ${tmpout} ${dropbox} ${scriptdir}/../${config}
