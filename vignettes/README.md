# PanGenie Vignettes

**Author:** Andrew Blair

This directory contains vignettes and tutorials for running PanGenie on paired WGS and RNA-seq cohorts.

## GEUVADIS GRCh38 Analysis

A step-by-step guide for genotyping GEUVADIS samples using PanGenie with the GRCh38 reference genome.

### Prerequisites

#### Download Reference Genome

```bash
rsync -avz \
  rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
  .
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
  -r hg38.fa.gz \
  -o pangenie_index/index
```

**Note:** You'll also need the IGSR reference for CRAM to FASTQ conversion:
```bash
wget -O GRCh38_full_analysis_set_plus_decoy_hla.fa \
  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```



### Running the Workflow

#### Prepare Sample Manifest

Sample manifest `sample_data.csv` with GEUVADIS samples:
```
sample_id,cram_url
NA12778,ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239484/NA12778.final.cram,b03ae320c1a3b13c750d23d28a8dbc13
```

```
aria2c \
  -x 8 -s 8 -k 1M \
  --file-allocation=none \
  --continue=true \
  --summary-interval=30 \
  --max-tries=20 \
  --retry-wait=30 \
  --timeout=60 \
  --connect-timeout=30 \
  -o "${SAMPLE_ID}.cram" \
  "${CRAM_URL}"

samtools fastq \
  -@ "${THREADS}" \
  --reference "${REFERENCE}" \
  -1 "${SAMPLE_ID}_R1.fastq.gz" \
  -2 "${SAMPLE_ID}_R2.fastq.gz" \
  -s /dev/null \
  -0 /dev/null \
  -n \
  "${SAMPLE_ID}.cram"

singularity exec \
  -B "${SCRATCH}:/mnt/work" \
  -B "${INDEX_DIR}:/mnt/pangenie_index" \
  "${SIF}" \
  PanGenie \
    -f /mnt/pangenie_index/index \
    -i /mnt/work/"${MERGED_FASTQ}" \
    -s "${SAMPLE_ID}" \
    -j "${THREADS}" \
    -t "${THREADS}" \
    -o /mnt/work/"${SAMPLE_ID}"
```

---

## Contributing

Additional vignettes and improvements are welcome.
