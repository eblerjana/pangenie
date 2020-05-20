# PanGenie

A genotyper for various types of genetic variants (such as SNPs, indels and structural variants). Genotypes are computed based on read k-mer counts and a panel of known haplotypes.

## Requirements
* gcc 4.7+
* cmake
* jellyfish

## Installation
`` git clone https://jana_ebler@bitbucket.org/jana_ebler/pangenie.git``  
`` cd pangenie``  
``mkdir build; cd build; cmake .. ; make``

## Installing into a conda environment
`` git clone https://jana_ebler@bitbucket.org/jana_ebler/pangenie.git``  
`` cd pangenie``  
`` conda env create -f environment.yml``  
`` conda activate pangenie``   
``mkdir build; cd build; cmake .. ; make``

## Usage
Genotyping can be run as shown below. Required arguments are short read sequencing reads in FASTQ-format, the reference genome in FASTA-format and a fully-phased, multisample VCF-file that contains a panel of known haplotypes.
Per default, the program will run genotyping (Forward-Backward algorithm) and phasing (Viterbi algorithm) and produces a VCF-file with genotypes for all variants in the input VCF (named `` <outname>_genotyping.vcf``) and a separate
VCF-file with haplotypes (named `` <outname>_phasing.vcf``). Flags `` -g `` and `` -p `` can be used to run only genotyping or phasing, respectively.


```bat

usage: PanGenie [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf>

options:
	-c	count all read kmers instead of only those located in graph.
	-d	do not add reference as additional path.
	-g	only run genotyping (Forward backward algorithm)
	-i VAL	sequencing reads in FASTA/FASTQ format or Jellyfish database in jf format (required).
	-j VAL	number of threads to use for kmer-counting (default: 1).
	-k VAL	kmer size (default: 31).
	-m VAL	regularization constant for copynumber probabilities (default: 0.001).
	-n VAL	effective population size (default: 0.00001).
	-o VAL	prefix of the output files (default: result).
	-p	only run phasing (Viterbi algorithm)
	-r VAL	reference genome in FASTA format (required).
	-s VAL	name of the sample (will be used in the output VCFs) (default: sample).
	-t VAL	number of threads to use for core algorithm. Largest number of threads possible is the number of chromosomes given in the VCF (default: 1).
	-u	output genotype ./. for variants not covered by any unique kmers.
	-v VAL	variants in VCF format (required).
```

Note that for the current version of this algorithm, the number of haplotypes (which is twice the number of samples) given in the input VCF should not be larger than 30. For larger panels, first subsample a smaller panel and use that as input for genotyping.


