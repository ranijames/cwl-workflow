#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
doc: "Script to calculate the expression counts"

#### The tool count expression is dependent on python version 2.7
baseCommand: [python, COUNT_EXPRESSION_PATH]

inputs:
 filterM:
  type: string
  inputBinding:
   position: 1
 filterB:
  type: string
  inputBinding:
   position: 2
 filterV:
  type: string
  inputBinding:
   position: 3
 sample:
  type: string
 annotation:
  type: File
  inputBinding:
   prefix: -a
 bam:
  type: File
  inputBinding:
   prefix: -A
 exp_out:
  type: string
  inputBinding:
   prefix: -o


outputs:
 expression_out:
  type: File
  outputBinding:
   glob: $(inputs.exp_out)
