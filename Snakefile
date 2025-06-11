import os
import subprocess
from snakemake.utils import min_version
import pandas as pd
import gzip
from cyvcf2 import VCF

###################
## Helper functions
###################

def count_variants(vcf_file_path):
    """
    Count the number of variants in a VCF file using cyvcf2.
    """
    vcf_reader = VCF(vcf_file_path)
    variant_count = 0

    for variant in vcf_reader:
        variant_count += 1

    vcf_reader.close()
    return variant_count


#########################
## Pipeline configuration
##########################
configfile: "config/config.yaml"

wildcard_constraints:
  dataset="[Aa-Zz0-9]+"

# directories
WORKING_DIR = config["temp_dir"]
RES_DIR = config["result_dir"]
if not os.path.exists(RES_DIR):
    os.makedirs(RES_DIR)
if not os.path.exists(WORKING_DIR):
    os.makedirs(WORKING_DIR)


# get list of samples
samples_df = pd.read_csv(config["samples"], sep = "\t", index_col = 0)
SAMPLES = samples_df.index.values.tolist()

# get fastq file
def get_vcf_file(wildcards):
    vcf_file = samples_df.loc[wildcards.sample,"vcf"]
    return vcf_file

####################
## Desired outputs
####################
FILTERED_VCF = expand(RES_DIR + "filtered/{sample}.vcf.gz", sample = SAMPLES)
SNP_COUNTS= expand(RES_DIR + "counts/{sample}.n_snps.txt", sample = SAMPLES)
GENOTYPES = expand(RES_DIR + "genotypes/{sample}.genotypes.txt", sample = SAMPLES)

rule all:
    input:
        FILTERED_VCF, 
        SNP_COUNTS,
        GENOTYPES
    message:
        "All done!"
    shell:
        "rm -r {WORKING_DIR}/"


########################
## Original VCF metrics
#######################

rule count_original_snps:
    input:
        vcf = get_vcf_file
    output:
        n_snps = RES_DIR + "counts/{sample}.n_snps.txt"
    message:
        "Counting initial number SNPs in {wildcards.sample} VCF file"
    threads: 1
    run:
        number_of_variants = count_variants(input.vcf)
        print(f"The number of variants in {input.vcf} is: {number_of_variants}")
        with open(output[0], "w") as f:
            f.write(f"step0:\t{number_of_variants}\n")

###########################
## Keep only biallelic SNPs 
###########################

rule filter_to_keep_biallelic_snps:
    input:
        vcf = get_vcf_file
    output:
        WORKING_DIR + "{sample}.biallelic.vcf.gz"
    message:
        "Keeping only biallelic SNPs in {wildcards.sample} VCF file"
    threads: 20
    shell:
        "bcftools view --max-alleles 2 "
        "--with-header -Oz --output {output} "
        "--threads {threads} "
        "{input} "

###########################################
## SNP filters (applied to all individuals)
## Filters on SNP quality, depth, etc. 
###########################################

rule first_filters_on_snp_sites:
    input:
        WORKING_DIR + "{sample}.biallelic.vcf.gz"
    output:
        WORKING_DIR + "{sample}.biallelic.qc1.vcf.gz"
    message:
        "Step1: filtering {wildcards.sample} VCF file on SNP depth and SNP quality"
    params:
        min_snp_quality = config["bcftools"]["snp"]["min_snp_quality"],
        min_snp_depth = config["bcftools"]["snp"]["min_snp_depth"],
        mean_depth_per_site = config["bcftools"]["snp"]["mean_depth_per_snp_site"]
    threads: 20
    shell:
        "bcftools view --include 'QUAL >= {params.min_snp_quality} && (FORMAT/DP) >= {params.min_snp_depth} && MEAN(FMT/DP) >= {params.mean_depth_per_snp_site}' "
        "--threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}"

rule filter_on_fraction_missing_per_snp:
    input:
        vcf = WORKING_DIR + "{sample}.biallelic.qc1.vcf.gz"
    output:
        vcf = WORKING_DIR + "{sample}.biallelic.qc2.vcf.gz",
        miss = WORKING_DIR + "{sample}.missing_snp.txt",
        snps_to_keep = WORKING_DIR + "{sample}.snp_sites_to_keep.txt"
    message:
        "Step2: filtering {wildcards.sample} biallelic VCF file on percentage of missing data for SNP sites"
    params:
        max_missing_fraction_per_snp = config["bcftools"]["snp"]["max_missing_fraction_per_snp_site"]
    threads: 20
    shell:
        "vcftools --gzvcf {input.vcf} --missing-site --out {output.miss}; "
        "awk '{{if ($5 < {params.max_missing_fraction_per_snp}) print $1, $2}}' {output.miss} > {output.snps_to_keep}; "
        "bcftools view -R {output.snps_to_keep} "
        "{input.vcf} "
        "-Oz "
        "-o {output.vcf}"

#############################################################
## Individuals filters on selected list of genotypes
## Filter on MAF, F_MISSING etc
#############################################################

if config["individuals"] == "all":
    pass  # No need to filter individuals, keep all
# check if file exists
elif config["individuals"].endswith('.tsv'):
    rule select_individuals:
        input:
            WORKING_DIR + "{sample}.biallelic.qc2.vcf.gz"
        output:
            WORKING_DIR + "{sample}.biallelic.qc2.selected.vcf.gz"
        message:
            "Selecting individuals from {wildcards.sample} biallelic VCF file"
        params:
            individuals = config["individuals"] # List of individuals to keep ("all" to keep all)
        threads: 20
        shell:
            "bcftools view -S {params.individuals} --threads {threads} "
            "{input} "
            "-Oz "
            "-o {output}"
else: 
    raise ValueError("The 'individuals' parameter in the config file must be either 'all' or a path to a TSV file with individuals to keep.")

rule filter_on_maf: 
    input:
        WORKING_DIR + "{sample}.biallelic.qc2.selected.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.biallelic.qc2.selected.maf.vcf.gz"
    message:
        "Filtering {wildcards.sample} biallelic VCF file on MAF"
    params:
        min_maf = config["bcftools"]["individuals"]["min_maf"]
    threads: 20
    shell:
        "bcftools view -i 'MAF > {params.min_maf}' --threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}"

rule filter_on_fraction_missing:
    input:
         WORKING_DIR + "filtered/{sample}.qc.biallelic.maf.vcf.gz"
    output:
        RES_DIR + "filtered/{sample}.vcf.gz"
    message:
        "Filtering {wildcards.sample} biallelic VCF file on percentage of missing genotype calls"
    params:
        missing = config["bcftools"]["max_fraction"]
    threads: 20
    shell:
        "bcftools view -i 'F_MISSING < {params.missing}' --threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}"

rule extract_genotypes:
    input:
        RES_DIR + "filtered/{sample}.vcf.gz"
    output:
        RES_DIR + "genotypes/{sample}.genotypes.txt"
    message:
        "Extracting genotypes from {wildcards.sample} VCF file"
    threads: 20
    shell:
        "bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' {input} > {output} "

# Extract allele frequencies from the filtered VCF file
rule extract_allele_frequencies:
    input:
        RES_DIR + "filtered/{sample}.vcf.gz"
    output:
        RES_DIR + "allele_frequencies/{sample}.allele_frequencies.txt"
    message:
        "Extracting allele frequencies from {wildcards.sample} VCF file"
    threads: 20
    shell:
        "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' {input} > {output}"

#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' vfaba_hedin_peamust_SNP_GATK.vcf.gz > allele_frequencies.txt

# Filter to keep sample sequenced to a certain depth
"""     shell:
        "bcftools view --include 'MEAN(FMT/DP) >= {params.min_sample_depth}' "
        "--threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}" """