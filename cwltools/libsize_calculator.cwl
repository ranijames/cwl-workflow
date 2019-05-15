cwlVersion: v1.0
class: CommandLineTool

doc: "The tool for calculating the lib-size for a given set of input samples"

baseCommand: [python, /cluster/home/aalva/Projects/Kjo_proct/projects2018_tumorprofiler/expression/compute_libsize.py]

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
