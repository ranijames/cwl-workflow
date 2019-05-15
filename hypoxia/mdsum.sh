#!/bin/bash

if [ -z "$1" ]
then
	echo "Usage: sh `basename $0` CONFIG"
    exit 2
fi
config="$1"

tmpout=$(grep -A1 "\[TEMP-OUT\]" $config | tail -n 1)
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd ${tmpout}
for file in `ls *`
do
    filebase=$(basename $file)
    md5sum ${filebase} >${filebase}.md5

done