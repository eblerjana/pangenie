# CLAUDE.md

**Author:** Andrew Blair

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PanGenie is a pangenome-based genotyping tool for diploid genomes that uses short-read sequencing data and k-mer counting to genotype variants (SNPs, indels, and structural variants) represented in a pangenome graph. It uses Hidden Markov Models (HMM) with known haplotype sequences to compute genotype likelihoods.

Version: v4.2.1

## Building the Project

### Using Singularity (Recommended)
```bash
cd container/
sudo singularity build pangenie.sif pangenie.def
```

Run commands through the container:
```bash
singularity exec pangenie.sif PanGenie-index <PARAMETERS>
singularity exec pangenie.sif PanGenie <PARAMETERS>
```

### Using Conda
```bash
conda env create -f environment.yml
conda activate pangenie
mkdir build && cd build
cmake ..
make
```

### Manual Build (requires jellyfish-2.0)
```bash
mkdir build && cd build
cmake ..
make
```

### Running Tests
```bash
cd build
make check
# or
ctest
```

Build system uses CMake with C++20 standard. Main dependencies: Jellyfish 2.0, Cereal (serialization), ZLIB.

## Core Executables

- **PanGenie-index**: Preprocessing step - creates indices from reference FASTA + VCF (required for multi-sample genotyping)
- **PanGenie**: Main genotyping tool - computes genotypes for a sample using k-mer counts
- **PanGenie-vcf**: Converts serialized genotyping results to VCF format
- **PanGenie-sampling**: Performs haplotype panel downsampling
- **Analyze-UK**: Analyzes unique k-mer distributions

## Architecture Overview

### Data Flow: Two-Stage Workflow

**Stage 1: Indexing (PanGenie-index)** - Run once per VCF panel
```
Reference FASTA + VCF → GraphBuilder → Graph objects (per chromosome)
                      ↓
                  UniqueKmerComputer → K-mer indices
                      ↓
              Serialize to disk:
              - *_Graph.cereal (per chromosome)
              - *_kmers.tsv.gz (per chromosome)
              - *_UniqueKmersMap.cereal
              - *_path_segments.fasta
```

**Stage 2: Genotyping (PanGenie)** - Run once per sample
```
Reads (FASTA/FASTQ) + Pre-computed indices
          ↓
    JellyfishCounter → K-mer counts from reads
          ↓
    Histogram → Coverage estimation
          ↓
    HMM (Forward-Backward) → Genotype likelihoods
          ↓
    Output VCF with genotypes
```

### Key Components

**Core Data Structures** (src/)
- `variant.hpp/cpp`: Variant representation with alleles, paths, flanking sequences
- `graph.hpp/cpp`: Graph of variant bubbles organized by chromosome; merges nearby variants
- `dnasequence.hpp/cpp`: Space-efficient DNA sequence (2 bases per byte)
- `genotypingresult.hpp/cpp`: Stores genotype likelihoods and haplotype info

**K-mer Management**
- `uniquekmers.hpp`: Abstract base class for k-mer storage
- `biallelicuniquekmers.hpp/cpp`: Optimized for biallelic variants (uses KmerPath16)
- `multiallelicuniquekmers.hpp/cpp`: Handles multiallelic variants (uses full KmerPath)
- `uniquekmercomputer.hpp/cpp`: Computes unique k-mers distinguishing variant alleles
- `jellyfishcounter.hpp/cpp`: Wraps Jellyfish library for parallel k-mer counting

**HMM & Genotyping**
- `hmm.hpp/cpp`: Core HMM implementation
  - Forward-Backward algorithm for genotyping (computes posterior probabilities)
  - Viterbi algorithm for phasing (finds most likely path)
- `emissionprobabilitycomputer.hpp/cpp`: P(observed k-mers | genotype)
- `transitionprobabilitycomputer.hpp/cpp`: P(genotype_i → genotype_j) using recombination rate
- `copynumber.hpp/cpp`: Binomial/beta-binomial models for k-mer count distributions

**Command Interface**
- `commands.hpp/cpp`: High-level workflows (`run_index_command`, `run_genotype_command`, etc.)
- `commandlineparser.hpp/cpp`: Argument parsing and validation
- `pangenie-genotype.cpp`, `pangenie-index.cpp`, etc.: Executable entry points

### Architectural Patterns

**Polymorphic K-mer Storage**: Uses abstract `UniqueKmers` base class with specialized implementations for biallelic vs multiallelic variants, selected at runtime based on variant properties.

**Serialization**: Heavy use of Cereal library for object persistence. Pre-computed indices are serialized as `.cereal` files and can be reused across samples.

**Memory Optimization**:
- `DnaSequence`: 2-bit encoding (4 bases per byte)
- `KmerPath16`: 16-bit optimization for common biallelic cases
- Lazy deletion of variant data after processing

**Threading**: `threadpool.hpp` provides multi-threaded k-mer counting and per-chromosome HMM computation.

## Important Input Requirements

### VCF Requirements
PanGenie requires specially formatted input VCFs:
1. **Multi-sample** with phased haplotypes
2. **Fully phased** in a single block (start to end)
3. **Non-overlapping variants** - overlapping alleles must be merged into multi-allelic records
4. **Sequence-resolved** - no symbolic alleles like `<DEL>`

For Minigraph-Cactus VCFs, use vcfbub filtering:
```bash
vcfbub -l 0 -r 100000 --input input.vcf > pangenie-ready.vcf
```

### Reads
- Must be short reads (k-mer based approach)
- Single FASTA or FASTQ file (uncompressed)
- Can also accept pre-computed Jellyfish database (.jf format)

### Reference
- FASTA format (uncompressed)
- Must match VCF coordinate system

## Common Workflows

### Single Sample Genotyping (Recommended)
```bash
# Preprocessing (once per VCF panel)
PanGenie-index -v variants.vcf -r reference.fa -t 8 -o index_prefix

# Genotyping (once per sample)
PanGenie -f index_prefix -i reads.fq -s sample_name -j 8 -t 8 -o output_prefix
# Output: output_prefix_genotyping.vcf
```

### Multiple Sample Genotyping
```bash
# Run PanGenie-index once
PanGenie-index -v variants.vcf -r reference.fa -t 8 -o index_prefix

# Run PanGenie separately for each sample
for sample in sample1 sample2 sample3; do
    PanGenie -f index_prefix -i ${sample}.fq -s ${sample} -j 8 -t 8 -o ${sample}_output
done
```

### Single-Command Genotyping (Higher Memory)
```bash
# Skip preprocessing step (uses more RAM, suitable for single samples)
PanGenie -i reads.fq -r reference.fa -v variants.vcf -s sample_name -j 8 -t 8 -o output_prefix
```

### Optimized Resource Usage (Separate VCF Writing)
```bash
# Genotyping only (multi-threaded)
PanGenie -f index_prefix -i reads.fq -s sample_name -j 8 -t 8 -w -o output_prefix
# Output: output_prefix_genotyping.cereal

# VCF writing (single-threaded, can run with fewer resources)
PanGenie-vcf -f index_prefix -z output_prefix_genotyping.cereal -s sample_name -o output_prefix
# Output: output_prefix_genotyping.vcf
```

## Code Organization

### Source Structure (src/ - 93 files)
```
src/
├── pangenie-*.cpp           # Executable entry points (5 files)
├── commands.{hpp,cpp}       # High-level command workflows
├── commandlineparser.{hpp,cpp}
├── variant.{hpp,cpp}        # Variant representation
├── graph.{hpp,cpp}          # Variant graph per chromosome
├── variantreader.{hpp,cpp}  # VCF parsing
├── graphbuilder.{hpp,cpp}   # Graph construction
├── hmm.{hpp,cpp}            # HMM algorithms (Forward-Backward, Viterbi)
├── *uniquekmers.{hpp,cpp}   # K-mer storage (biallelic/multiallelic)
├── uniquekmercomputer.{hpp,cpp}  # K-mer computation
├── jellyfish*.{hpp,cpp}     # K-mer counting interface
├── *probability*.{hpp,cpp}  # Probability computations
├── copynumber.{hpp,cpp}     # Statistical distributions
├── haplotypesampler.{hpp,cpp}    # Panel downsampling
└── utilities (dnasequence, timer, threadpool, etc.)
```

### Test Structure (tests/ - 31 files)
Uses Catch2 framework. Major test files:
- `HMMTest.cpp` (35KB) - Comprehensive HMM algorithm tests
- `GraphBuilderTest.cpp` (18KB) - Graph construction tests
- `VariantReaderTest.cpp`, `UniqueKmersTest.cpp`, etc.

Run with: `make check` or `ctest` in build directory

### Helper Scripts (scripts/)
Python utilities for analysis and validation:
- `genotype-concordance*.py` - Compare genotyping results
- `plot-*.py` - Visualization tools
- `fit-mixture-mode.py` - Distribution fitting

## Key Parameters

**K-mer size** (`-k`): Default 31. Affects specificity vs sensitivity.

**Threading**:
- `-j`: Threads for k-mer counting (Jellyfish)
- `-t`: Threads for genotyping (per-chromosome parallelism)

**Hash size** (`-e`): Jellyfish hash size. Default 3000000000 (3GB). Reduce for small datasets (demo uses 100000).

**Panel reduction** (`-x`): Maximum panel size to keep (default 15). Reduces computational cost.

## Important Notes

**Input Format**: PanGenie cannot handle compressed input files. All FASTA/FASTQ and VCF inputs must be uncompressed.

**Variant Representation**: PanGenie represents variants as bubbles in a pangenome graph. It is a re-genotyping tool, not a variant caller - it only genotypes variants present in the input VCF.

**Chromosome Names**: Reference FASTA chromosome names must exactly match those in the VCF.

**Limitations**:
- Maximum 65534 input haplotypes (32767 diploid samples)
- Only diploid genomes supported
- Runtime increases with number of haplotypes in panel
- Designed for whole-genome genotyping (not targeted regions)

**Memory Management**: Using PanGenie-index preprocessing reduces memory usage. Single-command mode (with `-v` and `-r`) uses more RAM but saves disk space.

## HMM Algorithm Details

**Emission Probabilities**: Given a genotype, compute P(observed k-mer counts | genotype) using copy number models (binomial or beta-binomial distributions based on estimated coverage).

**Transition Probabilities**: Use recombination rate to model P(switching between haplotype pairs). Constant recombination rate across genome (not position-specific).

**Forward-Backward Algorithm**: Computes posterior genotype probabilities at each variant position. Produces genotype likelihoods (GL field in VCF).

**Viterbi Algorithm** (optional with `-p`): Finds most likely sequence of genotypes across all positions. Used for phasing and panel downsampling.

**Coverage Estimation**: Uses k-mer count histogram to identify heterozygous peak and estimate sequencing coverage.

## Working with HPRC Data

The `hprc-README.md` file contains details about HPRC v2.0 Minigraph-Cactus graphs:
- Pre-processed VCFs available for PanGenie (`.pgin.vcf.gz` files)
- Biallelic decomposed VCFs (`.pgbi.vcf.gz` files) for downstream analysis
- Multiple reference builds: CHM13, GRCh38, GRCh37
- Benchmark graphs with held-out samples (HG002, HG005, NA19240)

PanGenie-ready HPRC VCFs were processed using:
https://github.com/eblerjana/genotyping-pipelines/tree/main/prepare-vcf-MC

## Vignettes

The `vignettes/` directory contains practical examples and tutorials for running PanGenie on specific datasets.

### GEUVADIS GRCh38 Analysis

A complete workflow for genotyping GEUVADIS samples using PanGenie with the GRCh38 reference genome. This vignette demonstrates:

- Processing 1000 Genomes Project GEUVADIS RNA-seq samples
- Converting CRAM files to FASTQ format using the correct GRCh38 reference (`GRCh38_full_analysis_set_plus_decoy_hla.fa`)
- Running PanGenie on short-read data aligned to GRCh38
- Integration with WDL-based workflow management using Toil

**Key considerations for GEUVADIS data:**
- CRAM files require the specific 1000 Genomes GRCh38 reference with decoy sequences and HLA alleles
- Reference mismatch will cause CRAM decompression failures with MD5 checksum errors
- FASTQ extraction must use `samtools bam2fq` with the matching reference genome

See `vignettes/` for detailed documentation and example workflows.

## Serialization Format

Uses Cereal library for C++ object serialization. Pre-computed index files:
- `.cereal` files: Binary serialized objects (Graph, UniqueKmersMap, GenotypingResult)
- `_kmers.tsv.gz`: Compressed k-mer tables per chromosome
- `_path_segments.fasta`: Reference and allele sequences

These files enable fast loading and reuse across samples without recomputing variant processing.
