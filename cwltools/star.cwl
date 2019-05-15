#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [/cluster/work/grlab/share/modules/packages/star/2.5.3a/bin/STAR]
doc: "STAR: Alignment"

requirements:
 ResourceRequirement:
  coresMax: $(inputs.runThreadN)
  ramMax: 80000
 InlineJavascriptRequirement: {}
 #ShellCommandRequirement: {}
 #InitialWorkDirRequirement: {}


inputs:
 genomeDir:
  type: Directory
  inputBinding:
   position: 1
   prefix: --genomeDir
 reads1:
  type: File[]
  inputBinding:
   position: 2
   itemSeparator: ','
   prefix: --readFilesIn
 reads2:
  type: File[]
  inputBinding:
   position: 3
   itemSeparator: ','
 runThreadN:
  type: int?
  inputBinding:
   position: 24
   prefix: --runThreadN
 outFilterMultimapScoreRange:
  type: int?
  inputBinding:
   position: 5
   prefix: --outFilterMultimapScoreRange
 outFilterMultimapNmax:
  type: int?
  inputBinding:
   position: 6
   prefix: --outFilterMultimapNmax
 outFilterMismatchNmax:
  type: int?
  inputBinding:
   position: 7
   prefix: --outFilterMismatchNmax
 alignIntronMax:
  type: int?
  inputBinding:
   position: 8
   prefix: --alignIntronMax
 alignMatesGapMax:
  type: int?
  inputBinding:
   position: 9
   prefix: --alignMatesGapMax
 sjdbScore:
  type: int?
  inputBinding:
   position: 10
   prefix: --sjdbScore
 alignSJDBoverhangMin:
  type: int?
  inputBinding:
   position: 11
   prefix: --alignSJDBoverhangMin
 genomeLoad:
  type: string?
  inputBinding:
   position: 12
   prefix: --genomeLoad
 limitBAMsortRAM:
  type: long?
  inputBinding:
   position: 13
   prefix: --limitBAMsortRAM
 readFilesCommand:
  type: string[]
  inputBinding:
   position: 14
   prefix: --readFilesCommand
 outFilterMatchNminOverLread:
  type: float?
  inputBinding:
   position: 15
   prefix: --outFilterMatchNminOverLread
 outFilterScoreMinOverLread:
  type: float?
  inputBinding:
   position: 16
   prefix: --outFilterScoreMinOverLread
 sjdbOverhang:
  type: int?
  inputBinding:
   position: 17
   prefix: --sjdbOverhang
 #outSAMstrandField:
  #type: string?
  #inputBinding:
   #position: 18
   #prefix: --outSAMstrandField
 outSAMattributes:
  type: string[]
  inputBinding:
   position: 18
   prefix: --outSAMattributes
   shellQuote: false
 #sjdbGTFfile:
  #type: File?
  #inputBinding:
   #position: 19
   #prefix: --sjdbGTFfile
 #limitSjdbInsertNsj:
  #type: int?
  #inputBinding:
   #position: 20
   #prefix: --limitSjdbInsertNsj
 outSAMunmapped:
  type: string?
  inputBinding:
   position: 19
   prefix: --outSAMunmapped
 outSAMtype:
  type: string[]
  inputBinding:
   position: 20
   prefix: --outSAMtype
  doc: |
    strings: type of BAM output
    1st word:
    BAM  ... output BAM without sorting
    2nd, 3rd:
    Unsorted           ... standard unsorted
    SortedByCoordinate ... sorted by coordinate. This option will allocate extra memory for sorting which can be specified by --limitBAMsortRAM.
 outSAMheaderHD:
  type: string[]
  inputBinding:
   position: 21
   prefix: --outSAMheaderHD
 outSAMattrRGline:
  type: string?
  inputBinding:
   position: 22
   prefix: --outSAMattrRGline
 #twopassMode:
  #type: string?
  #inputBinding:
   #position: 25
   #prefix: --twopassMode
 #outSAMmultNmax:
  #type: int?
  #inputBinding:
   #position: 26
   #prefix: --outSAMmultNmax
 outFileNamePrefix:
  type: string
  inputBinding:
   position: 23
   prefix: --outFileNamePrefix
 OutsamMode:
  type: string
  inputBinding:
   position: 24
   prefix: --outSAMmode

outputs:
 star_bam:
  type: File
  outputBinding:
   glob:
     ${
       var p = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
     if (inputs.outSAMtype.indexOf("SortedByCoordinate") > -1) {
        return p+"Aligned.sortedByCoord.out.bam";
     }
     }
  secondaryFiles: |
    ${
         var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
         return [
           {"path": p+"Log.final.out", "class":"File"},
           {"path": p+"SJ.out.tab", "class":"File"},
           {"path": p+"Log.out", "class":"File"}
         ];
     }
