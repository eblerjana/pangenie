# PanGenie Vignettes

**Author:** Andrew Blair

This directory contains vignettes and tutorials for running PanGenie on paired WGS and RNA-seq cohorts.

## GEUVADIS 

### CHM13 Analysis

### GRCh38 Analysis

Documentation and examples for genotyping GEUVADIS samples using PanGenie with the GRCh38 reference genome.

Reference index
```
rsync -avz \
  rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
  .
```

IGSR reference
```
wget -O GRCh38_full_analysis_set_plus_decoy_hla.fa \
  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

VCF
```
hprc-v2.0-mc-grch38.pgin.vcf.gz
```

### Contents


## Contributing

Additional vignettes and improvements are welcome.
