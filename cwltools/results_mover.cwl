cwlVersion: v1.0
class: CommandLineTool


baseCommand: [python,  /cluster/home/aalva/Projects/Kjo_proct/projects2018_tumorprofiler/ymlmaker/results_to_dropbox.py]

inputs:
 tempout:
  type: Directory
  inputBinding:
   position: 1
 dropbox:
  type: Directory
  inputBinding:
   position: 2
 config:
  type: File
  inputBinding:
   position: 3

outputs:
 results_move:
  type: Directory
  outputBinding:
   glob: $(runtime.outdir)
