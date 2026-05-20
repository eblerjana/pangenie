=====
Demo
=====

The typical use case is to run PanGenie on a whole genome dataset. The following example is just a little demo illustrating how to run PanGenie.

We run PanGenie given a pangenome graph (VCF file,``test-variants.vcf``), sequencing reads (FASTA/FASTQ file, ``test-reads.fa``) and a reference sequence (FASTA file, ``test-reference.fa``) provided in the ``demo/`` folder of the PanGenie repository which can be found `here <https://github.com/eblerjana/pangenie/tree/master/demo>`_. After installation, PanGenie's genotyping algorithm can be run using the following commands (which will take a few seconds for this example)::

 PanGenie-index -r test-reference.fa -v test-variants.vcf -o preprocessing -e 100000
 PanGenie -f preprocessing -i test-reads.fa -o test -e 100000


The result will be a VCF file named ``test_genotyping.vcf`` containing the same variants as the input VCF with additional genotype predictions, genotype likelihoods and genotype qualities.

Parameter ``-e`` sets the hash size used by Jellyfish for k-mer counting. When running PanGenie on a whole genome dataset, this parameter can be omitted (so that PanGenie uses the default value).

Per default, PanGenie uses a single thread. The number of threads used for k-mer counting and genotyping/phasing can be set via parameters ``-j`` and ``-t``, respectively.
