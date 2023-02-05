
//WORK IN PROGRESS, need to make  ped (pedgree files) and make sure that read groups correspond with samples
process whatshap{
    container: "quay.io/biocontainers/whatshap:1.7--py310h30d9df9_0"

    input:
        path reference_genome
        path input_vcf
        path bam

    script:
    """
    whatshap phase --ped ${pedigree} -o phased.vcf --reference=${reference_genome} ${input_vcf} ${bam}
    """

    output:
        path "phased.vcf"
}

workflow {

}