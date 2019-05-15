cwlVersion: v1.0
class: CommandLineTool

doc: "The tool for calculating the libsize for a given set of input samples"

baseCommand: [python, COMPUTE_LIBSIZE_PATH]

inputs:
 #sample: string
 samples:
  type: File
  inputBinding:
   position: 1
 outputfile:
  type: string
  inputBinding:
   position: 2

outputs:
 libsize:
  type: File
  outputBinding:
   glob: $(inputs.outputfile)
