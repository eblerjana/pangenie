# PanGenie Vignettes

**Author:** Andrew Blair

This directory contains vignettes and tutorials for running PanGenie on paired WGS and RNA-seq cohorts.

## GEUVADIS GRCh38 Analysis

A step-by-step guide for genotyping GEUVADIS samples using PanGenie with the GRCh38 reference genome.

### Prerequisites

#### Download Reference Genome

```bash
wget -O GRCh38_full_analysis_set_plus_decoy_hla.fa \
  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

#### Get Pangenome VCF

```bash
# Download HPRC pangenome reference
wget hprc-v2.0-mc-grch38.pgin.vcf.gz
```

#### Build PanGenie Index

```bash
singularity exec pangenie.sif PanGenie-index \
  -v hprc-v2.0-mc-grch38.pgin.vcf.gz \
  -r GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -o pangenie_index/index
```

### Running the Pipeline

#### Prepare Sample Manifest

Sample manifest `sample_data.csv` with GEUVADIS samples:
```
sample_id,cram_url
NA12778,ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239484/NA12778.final.cram,b03ae320c1a3b13c750d23d28a8dbc13
```


---

## Contributing

Additional vignettes and improvements are welcome.
