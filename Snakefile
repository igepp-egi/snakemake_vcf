import os
import subprocess
from snakemake.utils import min_version
import pandas as pd
import gzip
from cyvcf2 import VCF

# import Python functions from utils directory
from utils.count_variants_and_individuals import count_variants_by_type, count_individuals
from utils.parse_genotype_and_indiv import create_genotype_tsv
from utils.calculate_heterozygosity_rates import calculate_snp_heterozygosity_rates, calculate_individuals_heterozygosity_rates

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

# get VCF file
def get_vcf_file(wildcards):
    vcf_file = samples_df.loc[wildcards.sample,"vcf"]
    return vcf_file

# Small functions that reads a chromosome list file and returns a list of chromosomes separated by commas

def get_chromosomes():
    with open(config["chromosomes"], "r") as f:
        chromosomes = f.read().strip().split("\n")
    return ",".join(chromosomes)

####################
## Desired outputs
####################

if config["impute_genotypes"] == "yes":
    FILTERED_VCF = expand(RES_DIR + "filtered/{sample}_filtered_imputed_maf.vcf.gz", sample = SAMPLES) 
elif config["impute_genotypes"] == "no":
    FILTERED_VCF = expand(RES_DIR + "filtered/{sample}_filtered_not_imputed.vcf.gz", sample = SAMPLES) 
else: 
    raise ValueError("The 'impute_genotypes' parameter in the config file must be either 'yes' or 'no'.")

ALL_COUNTS =  expand(RES_DIR + "counts/counts_merged.{type_of_counts}.csv", type_of_counts=["snp","ind"])
GENOTYPES = expand(RES_DIR + "genotypes/{sample}.genotypes.tsv", sample = SAMPLES)
TABLES = expand(RES_DIR + "tables/table01_snp_density/{sample}.snp_density.{type}.tsv", sample = SAMPLES, type=["raw","filtered"])

if config["keep_temp_dir"] == True:
    rule all:
        input:
            GENOTYPES,
            FILTERED_VCF, 
            ALL_COUNTS,
            TABLES
        message: "All done! Temporary directory will be kept."
        shell:
            "cp config/config.yaml {RES_DIR}/; "
            "cp config/samples.tsv {RES_DIR}/; "
else:       
    rule all:
        input:
            FILTERED_VCF, 
            ALL_COUNTS,
            GENOTYPES,
            TABLES
        message: "All done!"
        shell:
            "rm -r {WORKING_DIR}/;"
            "cp config/config.yaml {RES_DIR}/; "
            "cp config/samples.tsv {RES_DIR}/; "




#####################################
## Step 1: keep selected individuals 
#####################################

rule step1_select_individuals:
    input:
        vcf = get_vcf_file
    output:
        WORKING_DIR + "filtered/{sample}.selected.vcf.gz"
    message:
        "Step 1: selecting individuals from {wildcards.sample} raw VCF file"
    params:
        individuals = config["individuals"]
    threads: 
        config["threads"]
    shell:
        "bcftools view -S {params.individuals} "
        "--threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}"

###################################
## Step 2: select chromosomes
###################################

rule step2_select_chromosomes:
    input:
        vcf =  WORKING_DIR + "filtered/{sample}.selected.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.selected.chr.vcf.gz"
    message:
        "Step 2: Selecting chromosomes in {wildcards.sample} VCF file"
    params:
        chromosomes = get_chromosomes() 
    threads: 
        config["threads"]
    shell:
        # index file
        "bcftools index {input};"
        # select chromosomes
        "bcftools view -r {params.chromosomes} "
        "--threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}"

###################################
## Step 3: keep only biallelic SNPs 
####################################

rule step3_keep_biallelic_snps:
    input:
        vcf =  WORKING_DIR + "filtered/{sample}.selected.chr.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.vcf.gz"
    message:
        "Step 3: Keeping only biallelic SNPs in {wildcards.sample} VCF file"
    threads: config["threads"]
    shell:
        "bcftools view --max-alleles 2 "
        "--with-header -Oz --output {output} "
        "--threads {threads} "
        "{input} "

###########################################################
## Step 4: SNP quality filters (applied to all individuals)
## Filters on SNP quality, depth, etc. 
###########################################################

rule step4_filter_snp_sites:
    input:
        WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.vcf.gz"
    message:
        "Step 4: filtering {wildcards.sample} VCF file on SNP depth and SNP quality"
    params:
        min_snp_quality = config["bcftools"]["snp"]["min_snp_quality"],
        min_snp_depth = config["bcftools"]["snp"]["min_snp_depth"],
        mean_depth_per_snp_site = config["bcftools"]["snp"]["mean_depth_per_snp_site"],
        max_snp_depth = config["bcftools"]["snp"]["max_snp_depth"]
    threads: config["threads"]
    shell:
        "bcftools view --include 'QUAL >= {params.min_snp_quality} && (FORMAT/DP) >= {params.min_snp_depth} && MEAN(FMT/DP) >= {params.mean_depth_per_snp_site} && (FORMAT/DP) <= {params.max_snp_depth}' "
        "--threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}"

#######################################################################################
## Step 5: filter on SNP missing rate (exclude if fraction of missing SNPs > threshold)
#######################################################################################

rule step5_filter_on_fraction_missing_per_snp:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.vcf.gz"
    output:
        vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.vcf.gz"
    message:
        "Step 5: filtering {wildcards.sample} biallelic VCF file to exclude SNP with a missing rate higher or equal to {params.max_missing_fraction_per_snp}"
    params:
        max_missing_fraction_per_snp = config["bcftools"]["snp"]["max_missing_fraction_per_snp"]
    threads: 
        config["threads"]
    shell:
        "bcftools view --exclude 'F_MISSING >= {params.max_missing_fraction_per_snp}' "
        "--threads {threads} "
        "{input.vcf} "
        "-Oz "
        "-o {output.vcf}"

###################################################################
## Step 6: filter on Minor Allele Frequency (MAF) before imputation
###################################################################

rule step6_filter_on_maf_before_imputation: 
    input:
        WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.maf.vcf.gz"
    message:
        "Step 6: filtering {wildcards.sample} biallelic VCF file on MAF higher than {params.min_maf}."
    params:
        min_maf = config["bcftools"]["individuals"]["min_maf_before_imputation"]
    threads: config["threads"]
    shell:
        "bcftools filter --include 'MAF > {params.min_maf}' --threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}"

############################################################################
## Step 7: filter on missing genotype calls using bcftools 'F_PASS(GT="mis")'
## Keep only genotypes with a maximum percentage of missing genotype calls
############################################################################

rule step7_filter_fraction_missing_per_genotype:
    input:
        WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.maf.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.maf.miss.vcf.gz"
    message:
        "Step 7: exclude genotypes from {wildcards.sample} VCF file if percentage of missing genotype calls is higher than {params.missing}"
    params:
        missing = config["bcftools"]["individuals"]["max_missing_fraction_per_genotype"]
    threads: config["threads"]
    shell:
        "bcftools view --exclude 'F_PASS(GT=\"mis\") > {params.missing}' "
        "{input} "
        "-Oz "
        "-o {output}"

######################################################################
## Step 8: filter on heterozygosity rate (keep only if het < threshold)
## 8.1: make a table of individuals + their genotype
## 8.2: calculate their heterozygosity rates
## 8.3: filter on heterozygoty excess
#####################################################################

rule step8_extract_individuals_and_genotypes:
    input: 
        vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.maf.miss.vcf.gz"
    output: 
        individuals = WORKING_DIR + "genotypes/{sample}.list_individuals.csv",
        genotypes = WORKING_DIR + "genotypes/{sample}.genotypes.GT.FORMAT" 
    message:
        "Step 8: extracting genotypes and individuals for {wildcards.sample} from {input.vcf} VCF file (step 8.1)"
    params: 
        out_prefix = WORKING_DIR + "genotypes/{sample}.genotypes"
    threads: 1
    shell:
        "bcftools query --list-samples {input.vcf} > {output.individuals};"
        "vcftools --gzvcf {input.vcf} --extract-FORMAT-info GT --out {params.out_prefix} "
            
rule step8_add_individuals_to_genotypes_tsv:
    input:
        individuals = WORKING_DIR + "genotypes/{sample}.list_individuals.csv",
        genotypes = WORKING_DIR + "genotypes/{sample}.genotypes.GT.FORMAT"
    output:
        WORKING_DIR + "genotypes/{sample}.genotypes.tsv"
    message:
        "Step 8: add individuals to the {wildcards.sample} genotype tsv file (step 8.1)"
    params:
        output_file_path = WORKING_DIR + "genotypes/{sample}.genotypes.tsv"
    threads: 1
    run:
        create_genotype_tsv(
            genotypes_tsv=input.genotypes, 
            individuals_csv=input.individuals, 
            output_tsv=params.output_file_path)

rule step8_calculate_individuals_heterozygosity_rates:
    input:
        genotypes = WORKING_DIR + "genotypes/{sample}.genotypes.tsv"
    output:
        ind_het = WORKING_DIR + "genotypes/{sample}.individuals.heterozygosity_rates.tsv"
    message:
        "Step 8: calculating heterozygosity rates for {wildcards.sample} genotypes (step 8.2)"
    threads: 1
    params:
        out_individuals_het_rates = WORKING_DIR + "genotypes/{sample}.individuals.heterozygosity_rates.tsv"
    run:
        calculate_individuals_heterozygosity_rates(
            input_genotypes_tsv=input.genotypes,
            out_individuals_het_rates_tsv=params.out_individuals_het_rates)

rule step8_extract_list_of_individuals_below_het_threshold:
    input:
        ind_het = WORKING_DIR + "genotypes/{sample}.individuals.heterozygosity_rates.tsv"
    output:
        ind_het_below = WORKING_DIR + "genotypes/{sample}.individuals_below_heterozygosity_threshold.list"
    message:
        "Step 8: extracting individuals that are below heterozygosity threshold from {wildcards.sample} genotypes (step 8.3)"
    params:
        max_het = config["bcftools"]["individuals"]["max_heterozygosity"],
        out_path = WORKING_DIR + "genotypes/{sample}.individuals_below_heterozygosity_threshold.list"
    threads: 1
    run:
        het = pd.read_csv(input.ind_het, sep="\t", index_col=None)
        max_het_rate_allowed = params.max_het
        filtered_het = het[het['heterozygosity_rate'] < float(params.max_het)]
        individuals_to_keep = filtered_het.individual_id
        individuals_to_keep.to_csv(params.out_path, index=False, header=False)

rule step8_filter_on_heterozygosity_rate:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.maf.miss.vcf.gz",
        ind_het_below = WORKING_DIR + "genotypes/{sample}.individuals_below_heterozygosity_threshold.list"
    output:
        vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.maf.miss.het.vcf.gz"
    message:
        "Step 8: filtering {wildcards.sample} VCF file to keep individuals below {params.het_rate_threshold} heterozygosity rate (step 8.3)"
    params:
        het_excess = WORKING_DIR + "genotypes/{sample}.individuals_below_heterozygosity_threshold.list",
        het_rate_threshold = config["bcftools"]["individuals"]["max_heterozygosity"]
    threads: config["threads"]
    shell:
        "bcftools view -S {input.ind_het_below} --threads {threads} "
        "{input.vcf} "
        "-Oz "
        "-o {output.vcf}"

###########################################
## Impute missing genotypes using Beagle v5
###########################################

if config["impute_genotypes"] == "yes":
    rule step9_impute_missing_genotypes:
        input:
            vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.maf.miss.het.vcf.gz"
        output:
            vcf = WORKING_DIR + "filtered/{sample}_filtered_imputed.vcf.gz"
        message:
            "Step 9: imputing missing genotypes in {wildcards.sample} VCF file using Beagle v5"
        params:
            beagle_memory = config["beagle"]["memory"],
            beagle_impute = config["beagle"]["impute"],
            output_prefix = WORKING_DIR + "filtered/{sample}_filtered_imputed"
        threads: config["threads"]
        shell:
            "beagle gt={input.vcf} "
            "out={params.output_prefix} "
            "nthreads={threads}"
elif config["impute_genotypes"] == "no":
    rule step9_no_imputation:
        input:
            vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.maf.miss.het.vcf.gz"
        output:
            vcf = RES_DIR + "filtered/{sample}_filtered_not_imputed.vcf.gz"
        message:
            "Step 9: not imputing missing genotypes in {wildcards.sample} VCF file; copying file from step 6 in {RES_DIR}"
        threads: config["threads"]
        shell:
            "cp {input.vcf} {output.vcf}"
else: 
    raise ValueError("The 'impute_genotypes' parameter in the config file must be either 'yes' or 'no'.")

##################################
## Filter on MAF after imputation 
##################################

if config["impute_genotypes"] == "yes":
    rule step10_filter_on_maf_after_imputation: 
        input:
            vcf = WORKING_DIR + "filtered/{sample}_filtered_imputed.vcf.gz"
        output:
            vcf = RES_DIR + "filtered/{sample}_filtered_imputed_maf.vcf.gz"
        message:
            "Step 10: filtering {wildcards.sample} imputed VCF file on MAF higher than {params.min_maf}."
        params:
            min_maf = config["bcftools"]["individuals"]["min_maf_after_imputation"]
        threads: 
            config["threads"]
        shell:
            "bcftools filter --exclude 'MAF < {params.min_maf}' --threads {threads} "
            "{input.vcf} "
            "-Oz "
            "-o {output.vcf}"
elif config["impute_genotypes"] == "no":
    pass
else: 
    raise ValueError("The 'impute_genotypes' parameter in the config file must be either 'yes' or 'no'.")

#####################
## SNP counts metrics
#####################

include: "subworkflows/count_snps_at_each_step.smk"

############################
## Individuals counts metrics
############################

include: "subworkflows/count_individuals_at_each_step.smk"

##################################
# All counts: SNPs and individuals
##################################

include: "subworkflows/count_snps_and_individuals_all_steps.smk"

#########################################################
# Extract final genotypes in tsv format from filtered VCF        
##########################################################

include: "subworkflows/extract_final_genotypes.smk"

#########################################################
# Tables for plots
# Table 1: SNP counts per chromosome and per bin (raw and final filtered VCF)
# Table 2: Fst per chromosome and per bin (raw and final filtered VCF)        
##########################################################

rule calculate_snp_density:
    input:
        raw_vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.maf.miss.het.vcf.gz",
        filtered_vcf = RES_DIR + "filtered/{sample}_filtered_imputed_maf.vcf.gz" if config["impute_genotypes"] == "yes" else RES_DIR + "filtered/{sample}_filtered_not_imputed.vcf.gz"
    output:
        raw_density = RES_DIR + "tables/table01_snp_density/{sample}.snp_density.raw.tsv",
        filtered_density = RES_DIR + "tables/table01_snp_density/{sample}.snp_density.filtered.tsv"
    message:
        "Calculating SNP density in {wildcards.sample} VCF raw and filtered files"
    params:
        bin_size = config["snp_density_bin_size"],
        raw_out_prefix = RES_DIR + "tables/table01_snp_density/{sample}.snp_density.raw",
        filtered_out_prefix = RES_DIR + "tables/table01_snp_density/{sample}.snp_density.filtered"
    threads: 1        
    shell:
        "mkdir -p {RES_DIR}/tables/table01_snp_density/; "
        "vcftools --gzvcf {input.raw_vcf} --SNPdensity {params.bin_size} --out {params.raw_out_prefix};"
        "vcftools --gzvcf {input.filtered_vcf} --SNPdensity {params.bin_size} --out {params.filtered_out_prefix};"
        "mv {params.raw_out_prefix}.snpden {output.raw_density};"
        "mv {params.filtered_out_prefix}.snpden {output.filtered_density};"