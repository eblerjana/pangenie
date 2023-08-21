# Running PanGenie given a reference panel VCF #

This pipeline prepares a given input VCF to be used as input for PanGenie, and then runs PanGenie to compute genotypes for the input variants.

## Required input data

* a reference genome in FASTA format
* the input VCF (fully-phased, multi-sample VCF representing a reference panel)
* short sequencing reads in FASTA/FASTQ format

The pipeline can be used to run PanGenie on multiple samples. Paths to reads for each sample can be specified in the config file.

## How to run

* insert paths to the input data (reference genome, VCF, reads, PanGenie executables) into the config file (`` contig.yaml ``)
* run pipeline: `` snakemake -j <n_cores> --use-conda `` 

## Output

The final genotypes will be written to:

* `` <outdir>/genotypes/<sample>-genotypes.vcf ``
