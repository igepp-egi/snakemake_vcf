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

def calculate_heterozygosity_rate(genotype_file_path):
    """
    Calculate the heterozygosity rate from a genotype file.
    """
    with open(genotype_file_path, 'r') as f:
        lines = f.readlines()

    total_snps = len(lines)
    heterozygous_snps = sum(1 for line in lines if line.strip().split('\t')[2] == '0/1')

    heterozygosity_rate = heterozygous_snps / total_snps if total_snps > 0 else 0
    return heterozygosity_rate

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
        print(f"The number of variants is: {number_of_variants}")
        with open(output[0], "w") as f:
            f.write(f"step0:\t{number_of_variants}\n")

rule calculate_heterozygosity_rate_in_original_vcf:
    input:
        vcf = get_vcf_file
    output:
        freq_het = RES_DIR + "metrics/raw/{sample}.freq_het.txt"
    message:
        "Counting initial number SNPs in {wildcards.sample} VCF file"
    threads: 1
    run:
        number_of_variants = count_variants(input.vcf)
        print(f"The number of variants is: {number_of_variants}")
        with open(output[0], "w") as f:
            f.write(f"step0:\t{number_of_variants}\n") 

####################################################
## Basic QC filters: min quality, min mean depth etc
####################################################

rule first_qc_filters:
    input:
        vcf = get_vcf_file
    output:
        WORKING_DIR + "{sample}.qc.vcf.gz"
    message:
        "Filtering {wildcards.sample} VCF file on quality"
    params:
        min_site_quality = config["bcftools"]["min_site_quality"],
        min_sample_depth = config["bcftools"]["min_sample_depth"],
        min_snp_depth = config["bcftools"]["min_snp_depth"],
    threads: 20
    shell:
        "bcftools view --include 'QUAL >= {params.min_site_quality} && (FORMAT/DP) >= {params.min_snp_depth} && MEAN(FMT/DP) >= {params.min_sample_depth}' "
        "--threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}"

#############################################################
## Keep only biallelic SNPs, filter on MAF, F_MISSING etc
#############################################################

rule filter_to_keep_biallelic_snps:
    input:
        WORKING_DIR + "{sample}.qc.vcf.gz"
    output:
        WORKING_DIR + "{sample}.qc.biallelic.vcf.gz"
    message:
        "Keeping only biallelic SNPs in {wildcards.sample} VCF file"
    threads: 20
    shell:
        "bcftools view --max-alleles 2 "
        "--with-header -Oz --output {output} "
        "--threads {threads} "
        "{input} "

rule filter_on_maf: 
    input:
        WORKING_DIR + "{sample}.qc.biallelic.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.qc.biallelic.maf.vcf.gz"
    message:
        "Filtering {wildcards.sample} biallelic VCF file on MAF"
    params:
        maf = config["bcftools"]["maf"]
    threads: 20
    shell:
        "bcftools view -i 'MAF > {params.maf}' --threads {threads} "
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
