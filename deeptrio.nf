#!/usr/bin/env nextflow

params.genome = "hs37d5.chr20" //"GRCh38.p13" 
params.reference_folder = "data/reference"
params.genome_fasta_name = "hs37d5.chr20.fa.gz"

int cores = Runtime.getRuntime().availableProcessors();
params.threads = cores -1 //let's keep one thread unused for other stuff
params.deepvariant_version = "deeptrio-1.4.0"
params.model = "WGS" //[WGS,WES,PACBIO,HYBRID_PACBIO_ILLUMINA]**
params.input_dir = "test"
params.reads_child = "HG001.chr20.10_10p1mb_sorted.bam"
params.reads_parent1 = "NA12891.chr20.10_10p1mb_sorted.bam"
params.reads_parent2 = "NA12892.chr20.10_10p1mb_sorted.bam"

params.output_dir = "./results"
params.family_name = "test_family"
params.child_name = "test_child"
params.parent1_name = "test_parent1"
params.parent2_name = "test_parent2"



process download_genome {

    container "quay.io/biocontainers/genomepy:0.14.0--pyh7cba7a3_2"
    input: 
        val genome
        path genome_folder
    
    script:
    """
    genomepy install --genomes_dir ${genome_folder} $genome
    """
    output:
        path "${genome_folder}/${genome}"
}

process deeptrio {

  debug true

  container "google/deepvariant:$params.deepvariant_version"

  publishDir params.output_dir

  tag "deeptrio on ${params.family_name}"

  input:
    path input_dir
    path reference_dir
    val genome_fasta_name
    

  output:
    path params.family_name


  script:
  """
  echo "starting deeptrio variant calling"
  echo "for the family ${params.family_name} "
  echo "for the reference ${reference_dir}"
  echo "with the input dir ${input_dir}"
  echo "with the input genome fasta name ${genome_fasta_name}"
  mkdir ${params.family_name}
  mkdir intermediate
  echo "the process folders are:"
  echo \$(ls)
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
  --model_type=${params.model} \
  --ref=${reference_dir}/${genome_fasta_name} \
  --reads_child=${input_dir}/${params.reads_child} \
  --reads_parent1=${input_dir}/${params.reads_parent1} \
  --reads_parent2=${input_dir}/${params.reads_parent2} \
  --sample_name_child=${params.child_name} \
  --sample_name_parent1=${params.parent1_name} \
  --sample_name_parent2=${params.parent2_name} \
  --output_vcf_child ./${params.family_name}/${params.child_name}.vcf \
  --output_vcf_parent1 ./${params.family_name}/${params.parent1_name}.vcf \
  --output_vcf_parent2 ./${params.family_name}/${params.parent2_name}.vcf \
  --output_gvcf_child ./${params.family_name}/${params.child_name}.g.vcf \
  --output_gvcf_parent1 ./${params.family_name}/${params.parent1_name}.g.vcf \
  --output_gvcf_parent2 ./${params.family_name}/${params.parent2_name}.g.vcf \
  --output_gvcf_merged ./${params.family_name}/${params.family_name}.g.vcf \
  --num_shards=${params.threads} \
  --logging_dir logs \
  --intermediate_results_dir intermediate \
  --runtime_report \
  --vcf_stats_report
  """
  // --postprocess_variants_extra_args flag_name=flag_value
}

workflow {
  main:
    input_dir_ch = Channel.fromPath(params.input_dir)
    def genome_folder = params.reference_folder + "/" + params.genome 
    def expected_fasta = file(genome_folder + "/" + params.genome_fasta_name)
    def reference_provided = (params.genome_fasta_name != null) && expected_fasta.isFile()

    if(reference_provided){
      genome_dir_ch = Channel.fromPath(genome_folder)
      genome_fasta_name_ch = Channel.value(params.genome_fasta_name)
      deeptrio(input_dir_ch, genome_dir_ch, genome_fasta_name_ch)
    } else {
      println("reference ${expected_fasta} does not seem to exist!")
      println("trying to download the genome with genomepy")
      genome_ch = Channel.value(params.genome)
      reference_folder_ch = Channel.fromPath(params.reference_folder)

      download_genome(genome_ch,  reference_folder_ch)
      genome_path = download_genome.out
      println("DOWNLOADED the genome, starting the deeptrio process!")
      genome_dir_ch = Channel.fromPath(genome_folder)
      genome_fasta_name_ch = Channel.value(params.genome + ".fa") //it is rully ugly, but should work with Ensembl
      deeptrio(input_dir_ch, genome_dir_ch, genome_fasta_name_ch)
      
    }

    /*
    if(params.reference_genome_fasta && file(params.reference_genome_fasta))
    {
      ref_genome = Channel.fromPath("${params.input_dir}/${params.genome}.fasta")
      ref_genome_fai = Channel.fromPath("${params.input_dir}/${params.genome}.fai")
      deeptrio(input_dir_ch, ref_genome, ref_genome_fai)
    }
    else {
      download_genome(Channel.value(params.genome), Channel.fromPath(params.genome_folder))
      ref_genome = Channel.fromPath("${download_genome.out}/${params.genome}.fasta")
      ref_genome_fai = Channel.fromPath("${download_genome.out}/${params.genome}.fai")
      deeptrio(input_dir_ch, ref_genome, ref_genome_fai)
    }
    */

  //emit:
  //  deeptrio.out
}


/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nWorkflow succeeded, you can check the results!" : "The workflow failed, please check for errors!" )
}