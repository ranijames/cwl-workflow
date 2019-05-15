import sys
import os
import subprocess
import argparse
from itertools import groupby, chain
import os, shutil, pathlib, fnmatch
#from utils import Makingdirsandcopyingfastqs


"""

Usage : python YML_hypoxia.py <ymlout> <data_repository> <FASTQlocal> <temp_out> <dropbox> <cancertype> <configFile>
Author: Alva James
Purpose :
    Given the list of sample names and path to the fastq.gz files for those samples, creates YAML input files for the CWL complete workflow for Hypoxia score generation
"""

parser = argparse.ArgumentParser()

parser.add_argument("FASTQlocal", help="path to directory containing the fastq.gz files downloaded freshly from openBIS")
parser.add_argument("dropbox", help="path to directory containing the fastq.gz files downloaded freshly from openBIS")
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

## The list for samples and participant ids
samplelistFile    = result['SAMPLES']
participant_ID    = result['PARTI_ID']

## Defining dictionaries for reads within lane1 and lane2 and participant ids
sample_dic_lane1 = []
patid_dict       = []
sample_dic_lane1 = dict((sample,[]) for sample in samplelistFile)
sample_dic_lane2 = dict((sample,[]) for sample in samplelistFile)
patid_dict       = dict((pid,[]) for pid in participant_ID)

## Defining the local fastq directory and the data_repository
fastqlocal       = args.FASTQlocal
raw              = args.dropbox

## Calling the function from previusly defined class
for sam in sample_dic_lane1.keys():
    for pids in patid_dict:
        PID = pids.split('_')[1]
        if sam in pids:
            destination = raw +  "/bkRNA/" + sam + "/raw/"
            if not os.path.isdir(destination):
                pathlib.Path(destination).mkdir(parents=True, exist_ok=True)
            for files in fnmatch.filter(os.listdir(fastqlocal), pat="*.tsv"):
                if sam in files:
                    #print(os.path.join(destination,files))
                    shutil.copy(os.path.join(fastqlocal,files), os.path.join(destination,files))
