#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

doc: "Workflow Components: alignment -> count expression -> lib-size -> hypoxia score"

requirements:
 - class: ScatterFeatureRequirement
 - class: SubworkflowFeatureRequirement
 - class: InlineJavascriptRequirement

inputs:
 reads1:
  type: File[]?
 reads2:
  type: File[]?
 sample:
  type: string
 genome:
  type: File
 genomeDir:
  type: Directory
 outSAMattrRGline: string
 annotation:
  type: File
 bam:
  type: File
 exp_out:
  type: string
 samples:
  type: File
 outputfile:
  type: string
 filterM:
  type: string
 filterB:
  type: string
 filterV:
  type: string
 result_dir:
  type: Directory
 cancertype:
  type: string
 sample_id:
  type: string
 count_expression:
  type: File
 libsize:
  type: File
 labkey_prefix:
  type: string


outputs:
 alignment_out:
  type: File
  outputSource: star/star_bam
 expression_out:
  type: File
  outputSource: expressioncount/expression_out
 libsize:
  type: File
  outputSource: libsize_calculator/libsize
 Hypoxiaresult:
  type: Directory
  outputSource: hypoxia_calculator/Hypoxiaresult


steps:
  star:
    run: star.cwl
    in:
     genomeDir: genomeDir
     reads1: reads1
     reads2: reads2
     outFileNamePrefix: sample
     runThreadN:
      default: 24
     outFilterMultimapScoreRange:
      default: 1
     outFilterMultimapNmax:
      default: 20
     outFilterMismatchNmax:
      default: 10
     alignIntronMax:
      default: 500000
     alignMatesGapMax:
      default: 1000000
     sjdbScore:
      default: 2
     alignSJDBoverhangMin:
      default: 1
     genomeLoad:
      default: NoSharedMemory
     limitBAMsortRAM:
      default: 70000000000
     readFilesCommand:
      default: [zcat]
     outFilterMatchNminOverLread:
      default: 0.33
     outFilterScoreMinOverLread:
      default: 0.33
     sjdbOverhang:
      default: 80
     outSAMstrandField:
      default: intronModif
     outSAMattributes:
      default: [NH, HI, NM, MD, AS, XS]
     #limitSjdbInsertNsj:
      #default: 2000000
     outSAMunmapped:
      default: Within
     outSAMtype:
      default: [BAM, Unsorted, SortedByCoordinate]
     outSAMheaderHD:
      default: ["@HD", "VN:1.4"]
     OutsamMode:
      default: Full
    out: [star_bam]
  expressioncount:
   run: count_expression.cwl
   in:
    sample: sample
    filterM: filterM
    filterB: filterB
    filterV: filterV
    annotation: annotation
    bam: star/star_bam
    exp_out: exp_out
   out: [expression_out]
  libsize_calculator:
   run: libsize_calculator.cwl
   in:
    samples: expressioncount/expression_out
    outputfile: outputfile
   out: [libsize]
  hypoxia_calculator:
   run: hypoxia.cwl
   in:
    count_expression: expressioncount/expression_out
    libsize: libsize_calculator/libsize
    result_dir: result_dir
    cancertype: cancertype
    sample_id: sample_id
    labkey_prefix: labkey_prefix
   out: [Hypoxiaresult]
