PanGenie
===================================

A short-read genotyper for various types of genetic variants (such as SNPs, indels and structural variants) represented in a pangenome graph. Genotypes are computed based on read k-mer counts and a panel of known, fully assembled haplotypes. A description of the method can be found in our `publication <https://doi.org/10.1038/s41588-022-01043-w>`_.

Citation
---------

| J. Ebler, P. Ebert, W. E. Clarke, T. Rausch, P. A. Audano, T. Houwaart, Y. Mao, J. Korbel, E. E. Eichler, M. C. Zody, A. T. Dilthey, and T. Marschall.
| *Pangenome-based genome inference*
| Nature Genetics, 54(4):518–525, 2022
| doi: `<https://doi.org/10.1038/s41588-022-01043-w>`_


Remarks
-------

* PanGenie is designed for whole genome genotyping, i.e. using the full set of variants as input rather than restricting to certain regions and/or variant types. In case you want to genotype a certain genomic region only (which is not the ideal use case for PanGenie), make sure the provided reference genome as well as the provided reads only contain data for these respective regions.


* PanGenie is a re-genotyping tool, **not a variant caller**. It re-genotypes the variants provided in the input pangenome, meaning the output VCF will contain the same variant records as the input, just with genotypes added for the sample PanGenie is run on. It cannot detect new variant alleles not represented in the pangenome graph.


Limitations
------------

* The runtime of PanGenie gets slow as the number of haplotype paths increases. Due to technical reasons, the current implementation of PanGenie cannot handle more than 65,534 input haplotypes (32,767 diploid samples).
* PanGenie can only genotype diploid genomes. It cannot be used for polyploid samples.



Issue tracker
--------------

Please use your `Issue Tracker <https://github.com/eblerjana/pangenie/issues>`_  for bug reports and feature requests.


Contents
--------

.. toctree::

   installation
   input
   usage
   resources
   data
   nested
