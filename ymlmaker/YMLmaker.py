import sys
import subprocess
import argparse
from itertools import groupby, chain
import os, shutil, fnmatch
from utils import Makingdirsandcopyingfastqs

"""

Usage : python YML_hypoxia.py <ymlout> <FASTQlocal> <temp_out> <dropbox> <cancertype> <configFile>
Author: Alva James
Purpose :
    Given the list of sample names and path to the fastq.gz files for those samples, creates YAML input files for the CWL complete workflow for Hypoxia score generation
"""

parser = argparse.ArgumentParser()
parser.add_argument("ymlout", help="output directory for generated yml files")
parser.add_argument("FASTQlocal", help="path to directory containing the fastq.gz files downloaded freshly from openBIS")
parser.add_argument("temp_out", help="path to directory containing the fastq.gz files downloaded freshly from openBIS")
parser.add_argument("dropbox", help="path to directory containing the fastq.gz files downloaded freshly from openBIS")
parser.add_argument("cancertype", help="path to directory containing the fastq.gz files downloaded freshly from openBIS")
parser.add_argument("configFile", help="file containing sample names, one name per line, the paths for dropbox and data_repository, and ymlout")

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(0)
args = parser.parse_args()

## Reading and transfering the samples and the labkey prefix from the CONFIG file to lists
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

# The list for samples and participant ids
samplelistFile    = result['SAMPLES']
cancer_type       = result['CANCER-TYPE']
labkey_prefix     = result['LABKEY_PREFIX']
dropbox           = args.dropbox
fastqlocal        = args.FASTQlocal

# Defining dictionaries for reads within lane1 and lane2 and labkey prefix
labkey_prefix_dict    = []
sample_dic_lane1      = []
sample_dic_lane2      = []

# Populating keys to the dictionaries from the input sample_id and labkey prefix
labkey_prefix_dict    = dict((pid,[]) for pid in labkey_prefix)
sample_dic_lane1      = dict((sample,[]) for sample in samplelistFile)
sample_dic_lane2      = dict((sample,[]) for sample in samplelistFile)

# Populating the values to the labkey_prefix_dict dictionary from the sample_id, and sorting the dictionary
for ids in labkey_prefix_dict:
    for names in sample_dic_lane1:
        if names in ids:
            labkey_prefix_dict[ids].append(names)
for prefix in labkey_prefix_dict:
    labkey_prefix_dict[prefix].sort()

# Calling the function for filemanagement from previously defined class
filemanagementinstance     = Makingdirsandcopyingfastqs()
filemanagementinstance.fastqmover_reads_for_lanes(sample_dic_lane1, sample_dic_lane2, dropbox, fastqlocal)

# Start writing the YML files
# The labkey prefix are taken from the above variable and used for in yml file for each sample(s)
for sam in sample_dic_lane1.keys():
    yml =  args.ymlout + "/" + sam + '.yml'
    for lkeyids in labkey_prefix_dict:
        labkey_prefix = lkeyids.split('_',1)[1]
        sample_id = lkeyids.split('_',1)[0]
        if sam == sample_id:

# defining paths where raw fastq files are located for each samples
            for root,sub,file in os.walk(dropbox):
                for folders in sub:
                    if folders.startswith("raw"):
                        fastqdir = dropbox + "/bkRNA/" + sam + "/" + "raw/"
            with open(yml, 'w') as ymlFH:
                ymlFH.write("reads1: [\n")
                ln1=len(sample_dic_lane1[sam])
                ct1=0
                for R1 in sorted(sample_dic_lane1[sam]):
                    ct1+=1
                    if ct1 < ln1:
                        ymlFH.write(" {class: File, path: "+ fastqdir + R1 + "},\n")
                    elif ct1 == ln1 :
                        ymlFH.write(" {class: File, path: "+ fastqdir  + R1 + "}\n")
                        ymlFH.write("]\n")

                # Inputs for star tool
                ymlFH.write("reads2: [\n")
                ln2=len(sample_dic_lane2[sam])
                ct2=0
                for R2 in sorted(sample_dic_lane2[sam]):
                    ct2+=1
                    if ct2 < ln2:
                        ymlFH.write(" {class: File, path: "+ fastqdir  + R2 + "},\n")
                    elif ct2 == ln2 :
                        ymlFH.write(" {class: File, path: "+ fastqdir  + R2 + "}\n")
                        ymlFH.write("]\n")

                # Writing sample name
                ymlFH.write("sample: " + sam + labkey_prefix + "__" + "\n")
                ymlFH.write("outSAMattrRGline: " + "ID::" + sam + labkey_prefix + "__" + "\n")

                # Adding the genome Directory
                ymlFH.write("genomeDir: \n")
                ymlFH.write(" class: Directory\n")
                ymlFH.write(" location: /cluster/work/grlab/share/databases/genomes/H_sapiens/STARgenomes/Hsap_GRCh37_plus_phiX174_NC_001422.1_overhang80\n")
                ymlFH.write("genome:\n")
                ymlFH.write(" class: File\n")
                ymlFH.write(" path: /cluster/work/tumorp/share/bulkRNA_references/genome.hs37d5_plus_Enterobacteria_phage_phiX174.fa\n")

               # Input for count expresion tool
                ymlFH.write("annotation:\n")
                ymlFH.write(" class: File\n")
                ymlFH.write(" path: /cluster/work/grlab/projects/TCGA/PanCancer/annotation/gencode.v19.annotation.hs37d5_chr.gtf\n")
                ymlFH.write("filterM: -m\n")
                ymlFH.write("filterB: -B\n")
                ymlFH.write("filterV: -v\n")
                ymlFH.write("bam:\n")
                ymlFH.write(" class: File\n")
                ymlFH.write(" path: "   + sam + labkey_prefix + "__Aligned.sortedByCoord.out.bam\n")
                ymlFH.write("exp_out: "   + sam + labkey_prefix + "__counts.tsv\n")

                # Inputs for libsize cwl tool
                ymlFH.write("samples:\n")
                ymlFH.write(" class: File\n")
                ymlFH.write(" path: "  + sam + labkey_prefix +  "__counts.tsv\n")
                ymlFH.write("outputfile: "  + sam + labkey_prefix + "__libsize.tsv\n")

                # Input for Hypoxia score
                ymlFH.write("count_expression:\n")
                ymlFH.write(" class: File\n")
                ymlFH.write(" path: "  +  args.temp_out + "/" + sam + labkey_prefix +  "__counts.tsv\n")
                ymlFH.write("libsize:\n")
                ymlFH.write(" class: File\n")
                ymlFH.write(" path: "  +  args.temp_out  + "/" + sam + labkey_prefix + "__libsize.tsv\n")
                ymlFH.write("result_dir:\n")
                ymlFH.write(" class: Directory\n")
                ymlFH.write(" path: " +  args.temp_out + "\n")
                ymlFH.write("cancertype: " + args.cancertype + "\n")
                ymlFH.write("sample_id: " + sam   + "\n")
                ymlFH.write("labkey_prefix: "  + labkey_prefix   + "\n")

                # Input for the results_mover and mdsum_generator
                ymlFH.write("config:\n")
                ymlFH.write(" class: File\n")
                ymlFH.write(" path: "  + os.path.abspath("CONFIG") + "\n")
                ymlFH.write("tempout:\n")
                ymlFH.write(" class: Directory\n")
                ymlFH.write(" path: "  + args.temp_out + "\n")
                ymlFH.write("dropbox:\n")
                ymlFH.write(" class: Directory\n")
                ymlFH.write(" path: "  + args.dropbox + "\n")
                ymlFH.write("hypoxia_newsample:\n")
                ymlFH.write(" class: File\n")
                ymlFH.write(" path: "  + sam + labkey_prefix + "__hypoxia_scores.tsv\n")
                ymlFH.write("hypoxia_bg:\n")
                ymlFH.write(" class: File\n")
                ymlFH.write(" path: "  + sam + labkey_prefix + "__tcga_gtex_hypoxia_scores.tsv\n")
