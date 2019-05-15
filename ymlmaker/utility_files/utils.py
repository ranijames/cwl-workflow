import os
import os, shutil, pathlib, fnmatch

__author__ = 'Alva James'

class Makingdirsandcopyingfastqs():
    def __init__(self, sample_dic_lane1=[], sample_dic_lane2=[], patid_dict=[], raw=None, fastqlocal=None, outdir=None):
        self.raw = raw
        self.fastqlocal = fastqlocal
        self.sample_dic_lane1= sample_dic_lane1
        self.sample_dic_lane2= sample_dic_lane2
        self.patid_dict = patid_dict

    def fastqmover_reads_for_lanes(self, sample_dic_lane1, sample_dic_lane2, patid_dict, raw, fastqlocal):
        for sam in sample_dic_lane1.keys():
            for pids in patid_dict:
                PID = pids.split('_')[1]
                if sam in pids:
                    destination =  raw + "/data_repository/" + "study/" +  PID  + "/" + sam + "/bkRNA/"+ "raw/"
                    if not os.path.isdir(destination):
                        pathlib.Path(destination).mkdir(parents=True)
                    for files in fnmatch.filter(os.listdir(fastqlocal), pat="*.gz"):
                        if sam in files:
                            shutil.copy(os.path.join(fastqlocal,files), os.path.join(destination,files))
                    for root, sub, file in os.walk(destination):
                        for name in file:
                            if name.endswith(".fastq.gz"):
                                R1_list =  ([f for f in os.listdir(root) if fnmatch.fnmatch(f, "*_R1_001_MM_1.fastq.gz")])
                                R2_list =  ([f for f in os.listdir(root) if fnmatch.fnmatch(f, "*_R2_001_MM_1.fastq.gz")])
                    for foreads1 in sample_dic_lane1:
                        for reads in R1_list:
                            if foreads1 in reads:
                                sample_dic_lane1[foreads1].append(reads)
                    for foreads2 in sample_dic_lane2:
                        for read2 in R2_list:
                            if foreads2 in read2:
                                sample_dic_lane2[foreads2].append(read2)
                            for x in sample_dic_lane1:
                                    sample_dic_lane1[x].sort()
                            for x in sample_dic_lane2:
                                    sample_dic_lane2[x].sort()
        return Makingdirsandcopyingfastqs(sample_dic_lane1,sample_dic_lane2)

class Moving_results():
    def __init__(self, samples_dictionary=[],tempoutdir=None, dropbox=None):
        self.samples_dictionary    = samples_dictionary
        self.dropbox    = dropbox
        self.tempoutdir = tempoutdir

    def move_results_to_dropbox(self,samples_dictionary, tempoutdir, dropbox):
        for samples in samples_dictionary.keys():
            for root,sub,file in os.walk(tempoutdir):
                for results in file:
                    if samples in results:
                        samples_dictionary[samples].append(results)
                        for sampleids in samples_dictionary.keys():
                             if sampleids in results:
                                derivedpath = dropbox +  "/bkRNA/"  + "/" + sampleids + "/derived/"
                                if not os.path.isdir(derivedpath):
                                    pathlib.Path(derivedpath).mkdir(parents=True)
                                for derivedfiles in fnmatch.filter(os.listdir(tempoutdir), pat="*"):
                                    if sampleids in derivedfiles:
                                        dropbox_results=shutil.move(os.path.join(tempoutdir,derivedfiles), os.path.join(derivedpath,derivedfiles))
        return Moving_results(dropbox_results)
