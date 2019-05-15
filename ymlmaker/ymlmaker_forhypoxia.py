import sys
import os
import subprocess
import argparse
from itertools import groupby, chain
import os, shutil, pathlib, fnmatch
from utility_files.utils import Makingdirsandcopyingfastqs

"""

Usage : python YML_hypoxia.py <ymlout> <temp_out> <configFile>
Author: Alva James
Purpose :
    Given the list of sample names and path to the fastq.gz files for those samples, creates YAML input files for the CWL pipeline part1
"""

parser = argparse.ArgumentParser()
parser.add_argument("ymlout", help="output directory for generated yml files")
parser.add_argument("temp_out", help="path to directory containing the fastq.gz files downloaded freashly from openBIS")
parser.add_argument("configFile", help="file containing sample names, one name per line, the paths for dropbox and data_repository, and ymlout")

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(0)

args = parser.parse_args()

## Reading and transfering the samples and participant ID infromation from the CONFIG file to lists
result    = {}
current   = None
with open(args.configFile,'r') as fd: #To close the file automatically
    for line in fd:
        line = line.strip()
        if line.startswith('['):
            current = line.strip('[]')
            result[current] = []
            continue
        if current is None: continue
        if line: result[current].append(line)

#The list for samples and participant ids
samplelistFile    = result['SAMPLES']

## Defining dictionaries for reads within lane1 and lane2 and participant ids
sample_dic_lane1 = []
sample_dict      = dict((sample,[]) for sample in samplelistFile)

#Copying all files from tumpro fastQ folder (downloaded freashly from openBIS) to the data_repository mirrored in leomed cluster
for sam in sample_dict.keys():
    yml =  os.path.join(args.ymlout, sam + '.yml')
    with open(yml, 'w') as ymlFH:
         # Writing sample name
        ymlFH.write("sample_input:\n")
        ymlFH.write(" class: Directory\n")
        ymlFH.write(" path: "  + args.temp_out + "\n")
        # Star input for BAM name
        ymlFH.write("outputdir: " +  args.temp_out + "\n")
        ymlFH.write("cancertype: SKCM\n")
        ymlFH.write("sample_id: " + sam   + "\n")
