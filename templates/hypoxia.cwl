cwlVersion: v1.0
class: CommandLineTool

requirements:
 InlineJavascriptRequirement: {}

baseCommand: [Rscript, HYPOXIA_PATH]

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
