cwlVersion: v1.0
class: CommandLineTool
doc: "Script to calculate the hypoxia score"

requirements:
 InlineJavascriptRequirement: {}

baseCommand: [Rscript, /cluster/home/aalva/Projects/Kjo_proct/projects2018_tumorprofiler/hypoxia/R_scripts/project_new_samples_HIF.R]

inputs:
 result_dir:
  type: Directory
  inputBinding:
   position: 1
 cancertype:
  type: string?
  inputBinding:
   position: 2
 sample_id:
  type: string?
  inputBinding:
   position: 3
 count_expression:
  type: File
  inputBinding:
   position: 4
 libsize:
  type: File
  inputBinding:
   position: 5
 labkey_prefix:
  type: string
  inputBinding:
   position: 6


outputs:
 Hypoxiaresult:
  type: Directory
  outputBinding:
   glob: $(runtime.outdir)
