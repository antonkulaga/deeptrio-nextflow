#!/usr/bin/env nextflow

params.genome = "GRCh38.p13"
params.genome_folder = "data/reference"

int cores = Runtime.getRuntime().availableProcessors();

params.threads = cores -1 //let's keep one thread unused for other stuff
params.deepvariant_version = "deeptrio-1.4.0"
params.model = "WGS" //[WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]**
params.input_dir = "test"
params.reads_child = "HG001.chr20.10_10p1mb_sorted.bam"
params.reads_parent1 = "NA12891.chr20.10_10p1mb_sorted.bam"
params.reads_parent2 = "NA12892.chr20.10_10p1mb_sorted.bam"
params.reference_genome = "hs37d5.chr20.fa.gz"
params.output_dir = "./results"
params.family_name = "test_family"
params.child_name = "test_child"
params.parent1_name = "test_parent1"
params.parent2_name = "test_parent2"
//params.dry_run = false

process download_genome {
    
    publishDir params.genome_folder

    container "quay.io/biocontainers/genomepy:0.14.0--pyh7cba7a3_2"
    input: 
        val genome
        path genome_folder
    
    script:
    """
    genomepy install --genomes_dir ${genome_folder} $genome
    """
    output:
        path "$genome"
        path "$genome/${params.genome}.fa"
}

process deeptrio {

  container "google/deepvariant:$params.deepvariant_version"

  publishDir params.output_dir

  tag "deeptrio on $params.family_name"

  input:
    path input_dir
    path reference_genome

  output:
    path "${results}"

  script:
  """
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
  --model_type=${params.model} \
  --ref=${reference_genome} \
  --reads_child=${input_dir}/${params.reads_child} \
  --reads_parent1=${input_dir}/${params.reads_parent1} \
  --reads_parent2=${input_dir}/${params.reads_parent2} \
  --sample_name_child=${params.child_name} \
  --sample_name_parent1=${params.parent1_name} \
  --sample_name_parent2=${params.parent2_name} \
  --output_vcf_child ${results}/${params.child_name}.vcf \
  --output_vcf_parent1 ${results}/${params.parent1_name}.vcf \
  --output_vcf_parent2 ${results}/${params.parent2_name}.vcf \
  --output_gvcf_child ${results}/${params.child_name}.g.vcf \
  --output_gvcf_parent1 ${results}/${params.parent1_name}.g.vcf \
  --output_gvcf_parent2 ${results}/${params.parent2_name}.g.vcf \
  --output_gvcf_merged ${results}/${params.family_name}.g.vcf \
  --num_shards=${params.cores} \
  --logging_dir ${results}/logs \
  --runtime_report \
  --vcf_stats_report \
  --intermediate_results_dir ${results}/intermediate
  """
  // --postprocess_variants_extra_args flag_name=flag_value
}

workflow {
  main:
    download_genome(Channel.value(params.genome), Channel.fromPath(params.genome_folder))
    input_dir_ch = Channel.fromPath(params.input_dir)
    deeptrio(input_dir_ch, download_genome.out[1])

  emit:
    deeptrio.out
}


/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nWorkflow succeeded, you can check the results!" : "The workflow failed, please check for errors!" )
}