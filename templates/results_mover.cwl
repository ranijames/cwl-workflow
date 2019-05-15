cwlVersion: v1.0
class: CommandLineTool


baseCommand: [python,  RESULTS-MOVE]

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
