# TuPro Workflow (CWL tools)

This README describes how to compute hypoxia scores for one or many samples, provided in fastq format. The requirements and necessary steps are explained in the following.

### Prerequisites   ###

You need the following tools available in the PATH of your environment:

 * **cwltool (v1.0.20181012180214 or higher)**
 
    * Please refer to the [CWL documentation](https://github.com/common-workflow-language/cwltool) if you are interested in further details.

 * **R-version 3.5 or higher**
    * Required packages will be installed automatically using an `install` script

 * **Python 2.7**
    * The following packages are needed: pip nodejs pysam h5py scipy openssl

We recommend using a virtual environment offered by package managers such as Anaconda. The following steps are an example to setup a working environment (named `workflow_cwl_tupro_2019`) with Anaconda:

```
conda create -n workflow_cwl_tupro_2019 -c default python=2.7 pip nodejs pysam h5py scipy
source activate  workflow_cwl_tupro_2019
conda install -c default openssl=1.0
pip install cwltool==1.0.20181012180214
pip install pathlib

```

After installing all the requirements, please clone this git repository to your local machine:

```
git clone git@github.com:ratschlab/projects2018_tumorprofiler.git
```
Then make sure to check out the branch `Tupro_CWL_workflow_2019`:

```
 git checkout Tupro_CWL_workflow_2019
```

For the following, we will assume, that you will be using the `workflow_cwl_tupro_2019` environment. If you are using a different environment, please edit the section `[CONDA-ENV]` in the file `CONFIG` accordingly.

As a last preparatory step, please run the installation script, to set internal paths according to your local setup (this needs to be done only once, but every time you update the code from the repository).


```
./install.sh

```

## CWL workflow ##

The pipelines in this repository are run using CWL, a workflow management system. The workflow pipeline consists of two components:

   i) A `.cwl` file, main workflow that contains the logical steps of the workflow and the names of all necessary parameters
  ii) A `.yml` file that provides the actual values for each parameter

The `.yml` are generated using the `YMLmaker.py` script. In addition to generating yml files for all input samples, the script
also copies all fastq files from the local directory to the corresponding locations in `dropbox` directory.


The pipeline we provide consists of four main parts. 

1) Copying raw fastqs files to dropbox/raw folder
2) Alignment and quantification
3) Hypoxia score computation
4) Move all the derived result files to dropbox/derived folder

Before we get started, please customize the `CONFIG` file that contains all the necessary input parameters that you need to provide:

* `[YML-OUT]`        : The directory where you want your yml files written to (CWL-operational folder)
* `[OUTPUT-DROPBOX]` : The directory where the final results would be written to
* `[SAMPLES]`        : The IDs of the samples you are analyzing (they should be a part of the fastq file names); provide one ID per line
* `[LABKEY_PREFIX]`     : The labkey prefix for each patient as the suffix to sample_id defined in `[SAMPLES]`,
                       separated with   `_`; if in doubt please refer to the current `CONFIG` file for examples.
                       In order to  find the Labkey prefixes please refer the [labkey](https://tp-labkey.ethz.ch/labkey/Tumor%20Profiler%20-%20Melanoma/login-login.view?returnUrl=%2Flabkey%2FTumor%2520Profiler%2520-%2520Melanoma%2Fproject-begin.view%3F).
* `[CANCER-TYPE]`   : The cancer type of the samples you want to analyze
* `[CONDA-ENV]`     : The name of the conda environment you are using, if it differs from `workflow_cwl_tupro_2019` (optional)
* `[TEMP-OUT]`      : The results are saved here temporarily  (CWL-operational folder)
* `[USER]`          : Plese provide your Leomed user name here


Currently, we agreed on stable paths for `[OUTPUT-DROPBOX]`, `[YML-OUT]` and `[TEMP-OUT]` are defined previously for Tumour profiler project, therefore we request not to change them for Tupro project unless you are testing.

**Please Note for testing**: If you are interested in testing the setup on test fastq files please refer the [test CONFIG](https://github.com/ratschlab/projects2018_tumorprofiler/blob/Tupro_CWL_workflow_2019/config_test) for all paths and test datasets.  


From here, we describe how to get this workflow run on your samples. The workflow is currently in two parts, as described below.

## Part 1.  Alignment to hypoxia score  ##

The steps from alignment to computing hypoxia score is now embedded into one single complete workflow. The complete workflow is described in [pipeline_until_hypoxia.cwl](https://github.com/ratschlab/projects2018_tumorprofiler/blob/Tupro_CWL_workflow_2019/cwltools/pipeline_until_hypoxia.cwl), which calls the following steps or cwltools:

`* Alignment (star.cwl)`

`* Expression quantification (count_expression.cwl)`

`* Library size calculation (libsize_calculation.cwl)`

`* Hypoxia score calculation (hypoxia.cwl)`

To run the workflow (we assume your in `Tupro_CWL_workflow_2019` environment), you can call the run script directly from the top-level directory of the repository:

```
bash hypoxia/Compute_hypoxia_score.sh CONFIG
```

You can either adapt the CONFIG file or adapt the above call to use a custom config file (that follows the same format).
The workflow will take approximately 4-5 hours, depending upon the sample size and to generate the following set of output files for each of the input sample(s):

* `.bam file`
* `.Log.final.out`
* `.Log.out`
* `.SJ.out.tab`
* `.count file`
* `.libsize file`
* `*.hypoxia_scores.tsv`


## Part 2. Data management ##

The second part of the workflow is datamangement, where the result files (defined above) would be moved to their respective folders and subfodlers based on the sample metadata (sample_id). The secound part generates `md5` files for each result samples and move the results to the dropbox. Inorder to run the part 2 please call the script from the top-level directory of the repository. Please make sure you are still in `Tupro_CWL_workflow_2019` environment while runing the below script.

``` 
bash hypoxia/Results_mover.sh CONFIG 

```



