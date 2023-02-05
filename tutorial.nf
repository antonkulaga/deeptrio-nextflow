#!/usr/bin/env nextflow

params.reads = "data/*[_-rR]{1,2}.fastq"
params.outdir = 'results'


def samples_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

process FASTQC {
  

  tag "FASTQC on $sample_id"
  publishDir params.outdir, mode:'copy'
  container 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'

  input:
    tuple val(sample_id), path(reads)

  output:
      path "fastqc_${sample_id}_logs" 


  script:
  """
  mkdir fastqc_${sample_id}_logs
  fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
  """
}


workflow {
  
  FASTQC(samples_ch)
  //return samples_ch.view()
}