cwlVersion: v1.0
class: CommandLineTool



baseCommand: [sh, /cluster/home/aalva/Projects/Kjo_proct/projects2018_tumorprofiler/hypoxia/mdsum.sh]


inputs:
 config:
  type: File
  inputBinding:
   position: 1

outputs:
 outFile:
  type: File[]
  outputBinding:
   glob: '*'
