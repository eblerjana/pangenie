import gzip

# at least this fraction of haplotypes needs to cover a variant allele with a non '.' allele
min_frac = 0.8


# correct genotypes on sex chromosomes (human)
rule correct_sex_chromosomes:
	input:
		vcf = lambda wildcards: CALLSETS[wildcards.caller]['vcf'],
		sample_info = lambda wildcards: CALLSETS[wildcards.caller]['sample_info']
	output:
		"{results}/vcf/{caller}/{caller}_corrected-sex-chromosomes.vcf.gz"
	log:
		"{results}/vcf/{caller}/{caller}_corrected-sex-chromosomes.log"
	conda:
		"../envs/whatshap.yml"
	resources:
		mem_mb=20000,
		walltime="01:59:00"
	shell:
		"""
		python3 workflow/scripts/correct-sex-chromosomes.py {input.vcf} {input.sample_info} 2> {log} | bcftools +fill-tags -Oz -o {output} -- -t AN,AC,AF
		"""


# remove sites that:
# - are covered by missing alleles in more than min_an haplotypes
# - contain Ns in sequence
rule filter_vcf:
	input:
		lambda wildcards: "{results}/vcf/{caller}/{caller}.vcf.gz" if not CALLSETS[wildcards.caller]['sample_info'] else "{results}/vcf/{caller}/{caller}_corrected-sex-chromosomes.vcf.gz"
	output:
		temp("{results}/vcf/{caller}/{caller}_filtered.vcf")
	conda:
		"../envs/genotyping.yml"
	log:
		"{results}/vcf/{caller}/{caller}_filtered.log"
	resources:
		mem_mb=20000,
		walltime="01:59:00"
	params:
		exclude = ','.join(CALLSETS[wildcards.caller]['exclude'])
	shell:
		"bcftools view --samples ^{params.exclude} --force-samples  {input} | bcftools view --min-ac 1 | python3 workflow/scripts/filter-vcf.py {min_frac} 2> {log} 1> {output}"


# remove alternative alleles that are not covered by any haplotype
rule trim_alt_alleles:
	input:
		"{results}/vcf/{caller}/{caller}_filtered.vcf"
	output:
		temp("{results}/vcf/{caller}/{caller}_filtered_trim.vcf")
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb=20000,
		walltime="00:30:00"
	shell:
		"bcftools view --trim-alt-alleles {input} > {output}"


# annotate the VCF and create an equivalent biallelic version of it
rule annotate_vcf:
	input:
		vcf="{results}/vcf/{caller}/{caller}_filtered_trim.vcf",
		gfa=lambda wildcards: CALLSETS[wildcards.caller]['gfa']
	output:
		multi="{results}/vcf/{caller}/{caller}_filtered_ids.vcf",
		multi_tmp=temp("{results}/vcf/{caller}/{caller}_filtered_ids-tmp.vcf"),
		biallelic="{results}/vcf/{caller}/{caller}_filtered_ids_biallelic.vcf",
		bi_tmp=temp("{results}/vcf/{caller}/{caller}_filtered_ids-tmp_biallelic.vcf")
	log:
		"{results}/vcf/{caller}/{caller}_filtered_ids.log"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_mb=200000,
		walltime="08:59:00"
	params:
		outname = "{results}/vcf/{caller}/{caller}_filtered_ids-tmp"
	shell:
		"""
		python3 workflow/scripts/annotate_vcf.py -vcf {input.vcf} -gfa {input.gfa} -o {params.outname} &> {log}
		cat {output.multi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' > {output.multi}
		cat {output.bi_tmp} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n\"}}' > {output.biallelic}
		"""

# compress annotated
rule compress_annotated:
	input:
		"{results}/vcf/{caller}/{caller}_filtered_{which}.vcf"
	output:
		"{results}/vcf/{caller}/{caller}_filtered_{which}.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bgzip -c {input} > {output}
		tabix -p vcf {output}
		"""
