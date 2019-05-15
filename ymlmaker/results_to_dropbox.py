import os
import os, shutil, fnmatch
import sys
import os
import subprocess
import glob
import fnmatch
import re
import argparse
from utils import Moving_results

"""
Usage : python results_to_dropbox.py <tmpoutdir> <dropbox> <configFile>
Author: Alva James
Purpose :
    Given the list of sample names and path to the fastq.gz files for those samples, creates YAML input files for the CWL pipeline part1
"""

parser = argparse.ArgumentParser()

parser.add_argument("tmpoutdir", help="where all the results are currently saved in")
parser.add_argument("dropbox", help="path where all the output files should be saved into")
parser.add_argument("configFile", help="file containing sample names, one name per line, the paths for dropbox and data_repository, and ymlout")

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(0)

args = parser.parse_args()

## Reading and transfering the samples and participant ID infromation from the CONFIG file to lists
result    = {}
current   = None

#Open the CONFIG file here
with open(args.configFile,'r') as fd:
    for line in fd:
        line = line.strip()
        if line.startswith('['):
            current = line.strip('[]')
            result[current] = []
            continue
        if current is None: continue
        if line: result[current].append(line)

samplelistFile     = result['SAMPLES']
samples_dictionary = dict((sample,[]) for sample in samplelistFile)
dropbox_mover      = Moving_results()
tempoutdir         = args.tmpoutdir
dropbox            = args.dropbox
dropbox_mover.move_results_to_dropbox(samples_dictionary,tempoutdir,dropbox)
