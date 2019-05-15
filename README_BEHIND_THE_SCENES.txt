
   #### NEW ADDITIONS TO THE PREVIOUS PIPELINE #####

AS OF APRIL 12-2019/ 
__author__ = 'Alva James'

The complete workflow generates all output files from alignment, quantification of expression, lib size calculation, hypoxia score steps defined in the workflow to a temporary directory. From the `md5sum_generator` step defined in the workflow generates `.md5` file for all result files/output files in the temporary directory. Finally, the `move_results`, in the complete workflow moves all result files to their corresponding directories/subfolders in Dropbox.


(1) Added a python class with 2 functions:

 Currently, the `ymlmaker.py`, searches for all the fastq files in the local directory (the place where the fastqs are downloaded directly from the openBIS (`/cluster/work/grlab/projects/tumor_profiler/data/fastq`, which is termed as FASTQ-LOCAL in CONFIG file). After that, the script copies all fastq files for the given input samples from the local directory to the data_repository created in leomed. There the function generates folder and subfolders for all given samples based on their basename (sample_ID), technology (bkRNA) and rootname (bkRNA/hospitalID).  Finally, for all steps in the main workflow the `ymlmaker.py` generates parameters and file paths.

All these are accomplished by using the following functions. These functions are defined in `ymlmaker/utility_files/utils.py` as a global class.

    (i) Function1:  For moving/copying all `fastq` files from local folder to the mirrored data-repository in leomed
                 Firstly, the function searches for `fastq` files for each given sample_id(s) in the local fastq folder.
                 Secondly, the part of copies all the reads to their respective directories and subfolders in data_repository. The function makes all directory and subfolders by its own if they are not existing
                 Thirdly, for each sample, the part of the function will search for reads (READ1 and READ2) and make a list with the reads from both lanes, if 2 lane exists, otherwise for single lane
                 Finally, that the list is populated to the predefined dictionaries which is having its key: the sample_id, and the lists would be populated as their items
                 
      The `YMLmaker.py` uses this information to generate reads parameter for the STAR tool.
                    
    (ii) Function2: This function looks at the global temporary output directory where all output files are collected from all steps
                    It collects those output files and copy them to the respective and subfolder within each directory created for each sample.
                    The function generates both directories (with sample_IDs) and subfolders (/derived). 

(2) Generating the corresponding directories and subfolders and the `md5sum` files
      Then, moving the result files and `md5sum` files to their corresponding folders or dictionaries in dropbox
      The mover pythons script: `ymlmaker/results_to_dropbox.py`

(3) Developed cwltools for md5sum and result mover script. The results mover is a Python script that does the data management based on sample metadata from all previous steps. Here, mainly it moves the files to dropbox based on sample metadata, from the temp folder

(4) Result mover script takes the function from the utility class and moves all files from temp folder to the DROPBOX

(5) Created additional two cwl tools

  1. moving the results
  2. generating the md5sum for each results

(6) Added additional inputs for hypoxia.cwl tool, mdsum.cwl, results_mover.cwl, All these are added in order to wire a relationship the last steps of the workflow so that they would wait until the all previous steps finishes. This is added mainly because CWL has only explicit catches, and no catch all mechanisms from the previous steps. We need to catch all result files. Here the best logic was to wire a connection from all the previous steps to result mover and mdsum maker tools. 


     #### Other New additions ######

(1) Changed the filenames for output from hypoxia generator R script

(2) Added a java expression to `star.cwl` to capture the log files, which gives the library size or  count

# April -15-2019

(3) Added a if loop conditions to write the sample_ids to the column sample_id in hypoxia dataframe, for each given input. For example, whether or not single sample input or multiple samples, the sample_IDS have to be capture correctly.

In addition to that, there was version conflict for cwltool2018 and 2019, in terms of generating log files from STAR aligner, so I have removed the secondary files section from star.cwl's outputs and moved it to the glob session.

