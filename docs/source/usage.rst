======
Usage
======

There are two ways of running PanGenie. The first way (recommended) is to first run a preprocessing step with ``PanGenie-index`` and then ``PanGenie`` with option ``-f``::

 PanGenie-index -v <variants.vcf> -r <reference.fa> -t <number of threads> -o <outfile-prefix>
 PanGenie -f <outfile-prefix>` -i <reads.fa/fq>  -s <sample-name> -j <nr threads kmer-counting> -t <nr threads genotyping>


The second way is to skip the preprocessing step and just run ``PanGenie`` with options ``-v`` and ``-r``, which will produce in the same genotyping results, but uses more memory::


 PanGenie -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf> -s <sample-name> -j <nr threads kmer-counting> -t <nr threads genotyping>


Both ways will produce the same end results, but running PanGenie in two separate steps is especially useful in cases where one wants to genotype the same set of variants across multiple samples. In such a case, ``PanGenie-index`` allows to do all preprocessing of the variant data only once instead of doing it over and over again for each sample. So when genotyping multiple samples, one needs to run ``PanGenie-index`` only a single time, and then runs ``PanGenie`` separately on each sample re-using the precomputed data. This reduces memory usage and runtime.

Below, details on these commands are provided.


-----------------------------------------------
Preprocessing step (optional, but recommended)
-----------------------------------------------

During preprocessing, steps **unrelated to the genotyped sample(s)** are performed, like processing the input variants and determining unique k-mers in the graph. In a setting in which the same set of input variants are genotyped across multiple samples, the advantage is that this preprocessing step **needs to be run only once**. The preprocessing step can be run using the command ``PanGenie-index``::


 PanGenie-index -v <variants.vcf> -r <reference.fa> -t <number of threads> -o <outfile-prefix>


The full list of options is previded below::


 program: PanGenie - genotyping based on kmer-counting and known haplotype sequences.
 author: Jana Ebler

 version: v4.2.1
 usage:
 PanGenie-index [options] -r <reference.fa> -v <variants.vcf> -o <index-prefix>

 options:
        -e VAL  size of hash used by jellyfish (default: 3000000000).
        -k VAL  kmer size (default: 31).
        -o VAL  prefix of the output files. NOTE: the given path must not include non-existent folders.
        -r VAL  reference genome in FASTA format. NOTE: INPUT FASTA FILE MUST NOT BE COMPRESSED.
        -t VAL  number of threads to use for kmer-counting (default: 1).
        -v VAL  variants in VCF format. NOTE: INPUT VCF FILE MUST NOT BE COMPRESSED.



The pre-proccessing step will result in a set of files (listed below) that can be used by ``PanGenie`` in order to genotype a specific sample:

* ``<outfile-prefix>_<chromosome>_Graph.cereal`` (one for each chromosome) serialization of Graph object
* ``<outfile-prefix>_<chromosome>_kmers.tsv.gz`` (one for each chromosome) containing unique k-mers
* ``<outfile-prefix>_UniqueKmersMap.cereal`` serialization of UniqueKmersMap object
* ``<outfile-prefix>_path_segments.fasta`` containing all reference and allele sequences of the graph

You don't need to understand what any of these files represent. They mainly contain information important to the subsequent genotyping step and ``PanGenie`` automatically processes them while running. So the only important thing is to not delete them prior to running ``PanGenie``.


---------------------
Genotyping step
---------------------

After preprocessing is completed, the genotyping step can be run in order to genotype a specific sample. If multiple samples shall be genotyped, this step needs to be run on each of these samples separately (while the preprocessing needs to be done only once). Based on the sequencing reads of a sample and the pre-computed files, genotyping is run using the command ``PanGenie`` with option ``-f``::


 PanGenie -f <outfile-prefix> -i <reads.fa/fq> -s <sample-name> -j <nr threads kmer-counting> -t <nr threads genotyping>


The full list of options is provided below::


 program: PanGenie - genotyping based on kmer-counting and known haplotype sequences.
 author: Jana Ebler

 version: v4.2.1
 usage:
 PanGenie [options] -f <index-prefix> -i <reads.fa/fq> -o <outfile-prefix>
 PanGenie [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf> -o <outfile-prefix>

 options:
        -a VAL  sample subsets of paths of this size (default: 0).
        -b VAL  effective population size for sampling step. (default: 0.01).
        -c      count all read kmers instead of only those located in graph
        -d      write sampled panel to additional output VCF.
        -e VAL  size of hash used by jellyfish (default: 3000000000).
        -f VAL  Filename prefix of files computed by PanGenie-index (i.e. option -o used with PanGenie-index).
        -g      run genotyping (Forward backward algorithm, default behaviour)
        -i VAL  sequencing reads in FASTA/FASTQ format or Jellyfish database in jf format. NOTE: INPUT FASTA/Q FILE MUST NOT BE COMPRESSED.
        -j VAL  number of threads to use for kmer-counting (default: 1).
        -k VAL  kmer size (default: 31).
        -o VAL  prefix of the output files. NOTE: the given path must not include non-existent folders (default: result).
        -p      run phasing (Viterbi algorithm). Experimental feature
        -r VAL  reference genome in FASTA format. NOTE: INPUT FASTA FILE MUST NOT BE COMPRESSED.
        -s VAL  name of the sample (will be used in the output VCFs) (default: sample).
        -t VAL  number of threads to use for core algorithm. Largest number of threads possible is the number of chromosomes given in the VCF (default: 1).
        -u      output genotype ./. for variants not covered by any unique kmers
        -v VAL  variants in VCF format. NOTE: INPUT VCF FILE MUST NOT BE COMPRESSED.
        -x VAL  to which size the input panel shall be reduced. (default: 15).
        -y VAL  Penality used for already selected alleles in sampling step. (default: 5).




The result will be a VCF file containing genotypes of the sample for the variants provided in the input VCF. Per default, the name of the output VCF is ``result_genotyping.vcf``. You can specify the prefix of the output file using option ``-o <prefix>``, i.e. the output file will be named as ``<prefix>_genotyping.vcf``.
The full list of options is provided below.

If you want to genotype the same set of variants across more than one sample, run the command above separately on each sample. The preprocessing step only needs to be run once (as long as the VCF does not change).

-----------------------------
Optimize compute resources
-----------------------------

The genotyping command itself can also be run in two steps, separating the genotyping step from the step that writes the final VCF. Unlike genotyping, writing the output VCF is always done using only a single core (regardless of ``-j`` and ``-t`` parameters). Therefore, running the two steps separately can be useful to optimize resource usage. These are the commands to use::


 PanGenie -f <outfile-prefix> -i <reads.fa/fq> -s <sample-name> -j <nr threads kmer-counting> -t <nr threads genotyping> -w -o <result-prefix>
 PanGenie-vcf -f <outfile-prefix> -z <result-prefix>_genotyping.cereal -s <sample-name> -o <result-prefix>


Note that the only difference for the genotyping command is the additional flag ``-w``. This will make PanGenie produce a ``<result-prefix>_genotyping.cereal`` file (as before, output prefix can be set using option ``-o <result-prefix>``) instead of an output VCF. The second command then converts this file into a VCF.

----------------------------------------
Running PanGenie with a single command
----------------------------------------

We also provide the option of running  ``PanGenie`` without running the preprocessing step first. This can be done by running it with parameters ``-v`` and ``-r`` instead of ``-f``. This will automatically do the preprocessing steps. In contrast to ``PanGenie-index``, it does not write as many files to disk during preprocessing to save time, but needs more RAM (similar to previous release v2.1.1). Running PanGenie like this might be useful in cases where one wants to genotype a single sample only, or to save some disk space.

As mentioned before, especially when genotyping more than one sample, it is beneficial to run both steps separately, since the preprocessing needs to be run only once for all samples, while the genotyping step needs to be run separately on each sample. Running PanGenie with a single command works as follows::

 PanGenie -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf> -s <sample-name> -j <nr threads kmer-counting> -t <nr threads genotyping> ``
