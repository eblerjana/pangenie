# Filtering cohort genotypes based on a support vector regression approach

This pipeline can be used to filter PanGenie genotypes computed for a larger cohort of samples. It is similar to the filtering approaches used to filter the genotypes computed for the 3,202 1000 Genomes samples in the HGSVC and HPRC projects. Besides the cohort genotypes, this pipeline requires genotypes for the graph samples (the panel samples that were present in the input VCF to PanGenie) as input.


## Required input data

* a biallelic, multi-sample VCF file with genotypes of all cohort samples and genotypes of the graph samples
* a PED file specifying trio relationships. Unrelated samples must be listed as well. The file must be in the following format:

```bat


#FamilyID	SampleID	FatherID	MotherID
BB46	HG02449	0	0
BB46	HG02450	0	0
BB46	HG02451	HG02449	HG02450

```

## How to run

* insert paths to the input data into the config file (`` contig.yaml ``)
* run pipeline: `` snakemake -j <n_cores> --use-conda `` 

## Output

The final filtered genotypes will be written to:

* `` <outdir>/filtered-vcfs/filtered_lenient.vcf.gz ``
