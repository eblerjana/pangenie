cwlVersion: v1.2
class: CommandLineTool
label: "PanGenie"

baseCommand: PanGenie

requirements:
  - class: DockerRequirement
    dockerPull: images.sb.biodatacatalyst.nhlbi.nih.gov/andrew.blair/pangenie:from_sif_v4.1.0_samtools
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement

inputs:

  fastq:
    type: File
    inputBinding:
      prefix: -i
      separate: true

  reference:
    type: File
    inputBinding:
      prefix: -r
      separate: true

  vcf:
    type: File
    inputBinding:
      prefix: -v
      separate: true

  sample_id:
    type: string
    inputBinding:
      prefix: -s
      separate: true

  threads_kmer:
    type: int
    inputBinding:
      prefix: -j
      separate: true

  threads_genotyping:
    type: int
    inputBinding:
      prefix: -t
      separate: true

  output_prefix:
    type: string
    inputBinding:
      prefix: -o
      separate: true

outputs:
  output_files:
    type: File[]
    outputBinding:
      glob: "$(inputs.output_prefix)*"
