# Create PanGenie-ready VCF from haplotype resolved assemblies #

This pipeline calls variants from haplotype-resolved assemblies of human samples and produces a fully-phased VCF file representing a pangenome graph.

Note: this pipeline was previously available under: https://bitbucket.org/jana_ebler/vcf-merging/src/master/pangenome-graph-from-assemblies/

## Required input data

* a reference genome in FASTA format
* a FASTA/FASTQ file with assembly contigs for each haplotype of each sample

## How to run

The provided config file (`` contig.json ``) looks as follows:

``` bat
{
	"reference": 
		{
			"filename": "/path/to/reference.fasta"
		},
	"assemblies":
		{
			"sample1" : ["/path/to/hap-1.fasta", "/path/to/hap-2.fasta"]
		},
	"trios":
		{
			"child_sample" : ["father_sample", "mother_sample"]
		},
	"scripts": "scripts",
	"outdir" : "results/"
}

```

It needs to be modified in the following way. 

* provide the path to your reference FASTA file in field `` "filename" ``.
* provide the paths to the haplotype-resolved assemblies for each sample in field `` "assemblies" ``. Give the name of a sample and paths to the assemblies of haplotype 1 and haplotype 2. Multiple samples can be listed here.

* if trio data is available, add them to the "trios" section in the config file using the child sample as key. An example is shown below. Otherwise just leave the "trios" section empty.

	```bat
 
   	"trios":
    	    {   
    	        "HG00733" : ["HG00731", "HG00732"]    
    	    },
	```

Once the config file is ready, you can run the pipeline as follows:
	
* run pipeline: `` snakemake -j <n_cores> --use-conda`` 

## Output

The pipeline first detects variants from the haplotype-resolved assemblies and then constructs a pangenome
representation of the respective callset by inserting the variant calls into the reference genome and merging overlapping alleles into
bubbles. The pangenome graph is represented in terms of a fully-phased, multi-sample VCF file containing one record per bubble.

* a multisample VCF containing all detected variant alleles will be written to: `` <outdir>/multisample-vcfs/callset-filtered.vcf ``
* a multisample VCF containing the pangenome graph representation will be written to: `` <outdir>/multisample-vcfs/graph-filtered.vcf ``

## Notes

* the pipeline assumes data of a diploid sample
