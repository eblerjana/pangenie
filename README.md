# PanGenie

A genotyper for various types of genetic variants (such as SNPs, indels and structural variants). Genotypes are computed based on read k-mer counts and a panel of known haplotypes. A description of the method can be found here: https://www.biorxiv.org/content/10.1101/2020.11.11.378133v1.

## Requirements
* gcc 4.9+
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

Note: we have observed that sometimes it is necessary to set ``export PKG_CONFIG_PATH="<path-to-miniconda>/miniconda3/envs/pangenie/lib/pkgconfig"``   before building the software.

## Runtime and memory usage

Runtime and memory usage depend on the number of variants genotyped and the number of haplotypes present in the graph.

With the data described here: https://www.biorxiv.org/content/10.1101/2020.11.11.378133v1, PanGenie ran in 1 hour and 25 minutes walltime using 22 cores (16 CPU hours) and used 68 GB RAM.
The largest dataset that we have tested contained around 16M variants, 64 haplotypes and around 30x read coverage. Using 24 cores, PanGenie run in 1 hour and 46 minutes (24 CPU hours) and used 120 GB of RAM.

## Usage
Genotyping can be run as shown below. Required arguments are short read sequencing reads in FASTQ-format, the reference genome in FASTA-format and a fully-phased, multisample VCF-file that contains a panel of known haplotypes.
Per default, the program will run genotyping (Forward-Backward algorithm) and phasing (Viterbi algorithm) and produces a VCF-file with genotypes for all variants in the input VCF (named `` <outname>_genotyping.vcf``) and a separate
VCF-file with haplotypes (named `` <outname>_phasing.vcf``). Flags `` -g `` and `` -p `` can be used to run only genotyping or phasing, respectively.


```bat


program: PanGenie - genotyping and phasing based on kmer-counting and known haplotype sequences.
author: Jana Ebler

usage: PanGenie [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf>

options:
	-a VAL	sample subsets of paths of this size. (default: 0).
	-c	count all read kmers instead of only those located in graph.
	-d	do not add reference as additional path.
	-e VAL	size of hash used by jellyfish. (default: 3000000000).
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


## Pan-Genome graph

PanGenie expects a directed and acyclic pangenome graph as input (-v option) which is represented as fully-phased, multisample VCF file. Variant records represent variant bubbles and each haplotype defines one path through this graph. We generate such graphs from haplotype-resolved assemblies using this pipeline: https://bitbucket.org/jana_ebler/vcf-merging. However, any fully-phased, multisample VCF file can be used as input. Important is, that variant records with overlapping coordinates are merged into a single, multiallelic VCF record, since otherwise PanGenie will filter them out. 


## Demo

The typical use case is to run PanGenie on a whole genome dataset. The following example is just a little demo illustrating how to run PanGenie. 

We run PanGenie given a pangenome graph (VCF file,``test-variants.vcf``), sequencing reads (FASTA/FASTQ file, ``test-reads.fa``) and a reference sequence (FASTA file, ``test-reference.fa``) provided in the ``demo/`` folder. After installation, PanGenie's genotyping algorithm can be run using the following command:

`` ./build/src/PanGenie -i test-reads.fq -r test-reference.fasta -v test-variants.vcf -g -o test -e 100000 ``

The result will be a VCF file named `` test_genotyping.vcf `` containing the same variants as the input VCF with additional genotype predictions, genotype likelihoods and genotype qualities.

To run only phasing, use the -p option instead to obtain phasing in a file that will be named `` test_phasing.vcf ``.

`` ./build/src/PanGenie -i test-reads.fq -r test-reference.fasta -v test-variants.vcf -p -o test -e 100000 ``

If `` -g `` and `` -p `` are omitted, PanGenie will run both genotying and phasing (which will take a few seconds for this example). Parameter `` -e `` sets the hash size used by Jellyfish for k-mer counting. When running PanGenie on a whole genome dataset, this parameter can be omitted (so that PanGenie uses the default value).

Per default, PanGenie uses a single thread. The number of threads used for k-mer counting and genotyping/phasing can be set via parameters ``-j`` and ``-t``, respectively. 
