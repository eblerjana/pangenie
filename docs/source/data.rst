====================
Data and genotypes
====================


We have already produced input reference panels for several datasets from high-quality, haplotype-resolved assemblies that can be used as input to PanGenie. These files were used to produce genotyping results for the HGSVC and HPRC projects. Genotypes for 3,202 samples from the 1000 Genomes Project produced based on these VCFs are also linked below.

**Note**: results produced by different versions of PanGenie are not directly comparable, since newer versions of PanGenie produce more accurate genotyping results.

----------------------
PanGenie v1.0.0
----------------------

.. csv-table::
  :header: "Dataset", "PanGenie input VCF", "Callset VCF", "1000G Genotypes (n=3,202)"
  :width: 100%
  :widths: auto
  :align: center

  "HGSVC-GRCh38 (freeze3, 64 haplotypes)",  `bubble-VCF <https://zenodo.org/record/7763717/files/pav-panel-freeze3.vcf.gz?download=1>`_ ,  `callset-VCF <https://zenodo.org/record/7763717/files/pav-calls-freeze3.vcf.gz?download=1>`_, `1000G-VCF <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/PanGenie_results/pangenie_merged_bi_all.vcf.gz>`_  (PanGenie v1.0.0)
  "HGSVC-GRCh38 (freeze4, 64 haplotypes)",  `bubble-VCF <https://zenodo.org/record/7763717/files/pav-panel-freeze4.vcf.gz?download=1>`_ , `callset-VCF <https://zenodo.org/record/7763717/files/pav-calls-freeze4.vcf.gz?download=1>`_,  `1000G-VCF <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/PanGenie_results/20201217_pangenie_merged_bi_all.vcf.gz>`_ (PanGenie v1.0.0)
  "HPRC-GRCh38 (88 haplotypes)", `bubble-VCF <https://zenodo.org/record/6797328/files/cactus_filtered_ids.vcf.gz?download=1>`_,  `callset-VCF <https://zenodo.org/record/6797328/files/cactus_filtered_ids_biallelic.vcf.gz?download=1>`_ , `1000G-VCF <https://zenodo.org/record/6797328/files/all-samples_bi_all.vcf.gz?download=1>`_  (PanGenie v1.0.0)



----------------
PanGenie v2.1.1
----------------

.. csv-table::
  :header: "Dataset", "PanGenie input VCF", "Callset VCF", "1000G Genotypes (n=3,202)"
  :width: 100%
  :widths: auto
  :align: center

  "HPRC-CHM13 (88 haplotypes)",  `bubble-VCF <https://zenodo.org/record/7839719/files/chm13_cactus_filtered_ids.vcf.gz?download=1>`_ , `callset-VCF <https://zenodo.org/record/7839719/files/chm13_cactus_filtered_ids_biallelic.vcf.gz?download=1>`_,  `1000G-VCF <https://zenodo.org/record/7839719/files/chm13_all-samples_bi_all.vcf.gz?download=1>`_ (PanGenie v2.1.1)

---------------------
PanGenie v3.1.0
---------------------

.. csv-table::
  :header: "Dataset", "PanGenie input VCF", "Callset VCF", "1000G Genotypes (n=3,202)"
  :width: 100%
  :widths: auto
  :align: center

   "HGSVC3+HPRC-CHM13 (214 haplotypes)", `bubble-VCF <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Genotyping_1kGP/PanGenie-genotypes/1.0/MC_hgsvc3-hprc_chm13_filtered_bubbles.vcf.gz>`_, `callset-VCF <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Genotyping_1kGP/PanGenie-genotypes/1.0/MC_hgsvc3-hprc_chm13_filtered_decomposed.vcf.gz>`_,   `1000G-VCF <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Genotyping_1kGP/PanGenie-genotypes/1.0/pangenie_chm13_all_decomposed_lenient.vcf.gz>`_ (PanGenie v3.1.0)

-----------------
PanGenie v4.2.1
-----------------

.. csv-table::
  :header: "Dataset", "PanGenie input VCF", "Callset VCF", "1000G Genotypes (n=3,202)"
  :width: 100%
  :widths: auto
  :align: center

  "HPRC2-CHM13 (462 haplotypes)", `bubble-VCF <https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2026_03_30_pangenie/mc_filtered_ids.vcf.gz>`_, `callset-VCF <https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2026_03_30_pangenie/mc_filtered_ids_biallelic.vcf.gz>`_ , `1000G-VCF <https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2026_03_30_pangenie/pangenie_all-samples_filtered.vcf.gz>`_ (PanGenie v4.2.1)



In all cases, the graph-VCFs provided in the second column were given as input to PanGenie. The callset-VCFs (third column) were used to convert the genotyped VCFs into a biallelic, callset representation using the following command::

 cat <pangenie-output> | python3 convert-to-biallelic.py <callset-VCF>  > callset-genotypes.vcf

The script ``convert-to-biallelic.py`` can be found `here <https://github.com/eblerjana/pangenie/blob/master/pipelines/run-from-callset/scripts/convert-to-biallelic.py>`_.
