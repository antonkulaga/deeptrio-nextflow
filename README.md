Just deep-trio-nextflow pipeline
=================================
This pipeline should process DeepTrio samples. You have to provide reference fasta, genomes (as sorted indexed bam) of the child and its two parents.
Names of the child and parents are also used to make files, so try to use underscore instead of space and latin characters.

To run the pipeline you need nextflow and singularity (see https://sylabs.io/singularity/) installed and run the following command:
```
nextflow run deeptrio.nf <parameters>
```
To download the latest human genome you can use prepare.nf workflow.

You can read more about nextflow at https://www.nextflow.io/ and watch introductory tutorial at https://www.youtube.com/watch?v=wbtMbJTo1xo&list=PLPZ8WHdZGxmUVZRUfua8CsjuhjZ96t62R (note many examples there are for DSL 1, while here I am using DSL 2)

You can also read general deeptrio docs at
https://github.com/google/deepvariant/blob/r1.4/docs/deeptrio-details.md 


KNOWN ISSUES
------------

In deepvariant parameters they expicetly ask for bam and fasta files but not for bai, fai files. At the same time they check if those files exist in the current folder. I cd-ed the folder to input so deepvariant can get those files, but it created other troubles - it started writing results to input. Probably it makes sence to symlink the files to the folder instead or use another workaround.

Phasing
-------

Looks like deeptrio does not produce phased vcf-s. Maybe something is wrong with parameters or maybe additional tool is needed for phasing.
I put whatshap.nf for phasing with whatsap tool (not tested yet), maybe something else is needed.