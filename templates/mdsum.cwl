cwlVersion: v1.0
class: CommandLineTool



baseCommand: [sh, MDSUM_SHELL]


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
