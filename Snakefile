import os
import subprocess
from snakemake.utils import min_version
import pandas as pd
import gzip
from cyvcf2 import VCF

# import Python functions from utils directory
from utils.count_variants import count_variants_by_type
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

# get fastq file
def get_vcf_file(wildcards):
    vcf_file = samples_df.loc[wildcards.sample,"vcf"]
    return vcf_file

####################
## Desired outputs
####################

if config["impute_genotypes"] == "yes":
    FILTERED_VCF = expand(RES_DIR + "filtered/{sample}_filtered_imputed.vcf.gz", sample = SAMPLES) 
elif config["impute_genotypes"] == "no":
    FILTERED_VCF = expand(RES_DIR + "filtered/{sample}_filtered_no_imputed.vcf.gz", sample = SAMPLES) 
else: 
    raise ValueError("The 'impute_genotypes' parameter in the config file must be either 'yes' or 'no'.")

ALL_COUNTS =  RES_DIR + "counts/counts_merged.csv"

if config["keep_temp_dir"] == True:
    rule all:
        input:
            FILTERED_VCF, 
            ALL_COUNTS
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
            BED
        message: "All done!"
        shell:
            "rm -r {WORKING_DIR}/;"
            "cp config/config.yaml {RES_DIR}/; "
            "cp config/samples.tsv {RES_DIR}/; "



#############################################################
## Individuals filters on selected list of genotypes
#############################################################

if config["individuals"] == "all":
    pass  # No need to filter individuals, keep all
# check if file exists
elif config["individuals"].endswith('.tsv'):
    rule step0_select_individuals:
        input:
            vcf = get_vcf_file
        output:
            WORKING_DIR + "filtered/{sample}.selected.vcf.gz"
        message:
            "Selecting individuals from {wildcards.sample} raw VCF file"
        params:
            individuals = config["individuals"] # List of individuals to keep ("all" to keep all)
        threads: config["threads"]
        shell:
            "bcftools view -S {params.individuals} --threads {threads} "
            "{input} "
            "-Oz "
            "-o {output}"
else: 
    raise ValueError("The 'individuals' parameter in the config file must be either 'all' or a path to a TSV file with individuals to keep.")

###########################
## Keep only biallelic SNPs 
###########################

rule step1_keep_biallelic_snps:
    input:
        vcf =  WORKING_DIR + "filtered/{sample}.selected.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.selected.biallelic.vcf.gz"
    message:
        "Step1: Keeping only biallelic SNPs in {wildcards.sample} VCF file"
    threads: config["threads"]
    shell:
        "bcftools view --max-alleles 2 "
        "--with-header -Oz --output {output} "
        "--threads {threads} "
        "{input} "

###########################################
## SNP filters (applied to all individuals)
## Filters on SNP quality, depth, etc. 
###########################################

rule step2_filter_snp_sites:
    input:
        WORKING_DIR + "filtered/{sample}.selected.biallelic.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.vcf.gz"
    message:
        "Step2: filtering {wildcards.sample} VCF file on SNP depth and SNP quality"
    params:
        min_snp_quality = config["bcftools"]["snp"]["min_snp_quality"],
        min_snp_depth = config["bcftools"]["snp"]["min_snp_depth"],
        mean_depth_per_snp_site = config["bcftools"]["snp"]["mean_depth_per_snp_site"]
    threads: config["threads"]
    shell:
        "bcftools view --include 'QUAL >= {params.min_snp_quality} && (FORMAT/DP) >= {params.min_snp_depth} && MEAN(FMT/DP) >= {params.mean_depth_per_snp_site}' "
        "--threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}"

rule step3_filter_on_fraction_missing_per_snp:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.vcf.gz"
    output:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.vcf.gz"
    message:
        "Step3: filtering {wildcards.sample}  biallelic VCF file to keep SNP with a minimum call rate percentage of {params.max_missing_fraction_per_snp}"
    params:
        vcf_file_prefix = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2",
        vcf_file_complete_name = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.vcf",
        max_missing_fraction_per_snp = config["vcftools"]["snp"]["max_missing_fraction_per_snp_site"] # more a calling rate with 0 (allow completely missing) to 1 (no missing)
    threads: config["threads"]
    shell:
        "vcftools --max-missing {params.max_missing_fraction_per_snp} --gzvcf {input.vcf} --recode --recode-INFO-all --out {params.vcf_file_prefix} ;"
        # rename to remove the recode extension and produce a .gz file
        "mv {params.vcf_file_prefix}.recode.vcf {params.vcf_file_prefix}.vcf ;"
        # compress using gzip
        "gzip {params.vcf_file_complete_name} ;"

############################
## Filter on MAF, F_MISSING 
############################

rule step4_filter_on_maf_before_imputation: 
    input:
        WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.vcf.gz"
    message:
        "Step4: filtering {wildcards.sample} biallelic VCF file on MAF higher than {params.min_maf}; before imputation"
    params:
        min_maf = config["bcftools"]["individuals"]["min_maf_before_imputation"]
    threads: config["threads"]
    shell:
        "bcftools filter --include 'MAF > {params.min_maf}' --threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}"

rule step5_filter_fraction_missing_per_genotype:
    input:
        WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.miss.vcf.gz"
    message:
        "Step5: exclude genotypes from {wildcards.sample} VCF file if percentage of missing genotype calls is higher than {params.missing}"
    params:
        missing = config["bcftools"]["individuals"]["max_missing_fraction_per_genotype"]
    threads: config["threads"]
    shell:
        "bcftools view --exclude 'F_PASS(GT=\"mis\") > {params.missing}' "
        "{input} "
        "-Oz "
        "-o {output}"

###########################################
## Extract genotypes and allele frequencies
## Filtering on HWE and heterozygote excess
###########################################

rule extract_individuals_and_genotypes:
    input: 
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.miss.vcf.gz"
    output: 
        individuals = WORKING_DIR + "genotypes/{sample}.list_individuals.csv",
        genotypes = WORKING_DIR + "genotypes/{sample}.genotypes.GT.FORMAT" 
    message:
        "Extracting genotypes and individuals for {wildcards.sample} from {input.vcf} VCF file"
    params: 
        out_prefix = WORKING_DIR + "genotypes/{sample}.genotypes"
    threads: 1
    shell:
        "bcftools query --list-samples {input.vcf} > {output.individuals};"
        "vcftools --gzvcf {input.vcf} --extract-FORMAT-info GT --out {params.out_prefix} "
            
rule add_individuals_to_genotypes_tsv:
    input:
        individuals = WORKING_DIR + "genotypes/{sample}.list_individuals.csv",
        genotypes = WORKING_DIR + "genotypes/{sample}.genotypes.GT.FORMAT"
    output:
        WORKING_DIR + "genotypes/{sample}.genotypes.tsv"
    message:
        "Create the {wildcards.sample} genotype tsv file"
    params:
        output_file_path = WORKING_DIR + "genotypes/{sample}.genotypes.tsv"
    threads: 1
    run:
        create_genotype_tsv(
            genotypes_tsv=input.genotypes, 
            individuals_csv=input.individuals, 
            output_tsv=params.output_file_path)

rule calculate_individuals_heterozygosity_rates:
    input:
        genotypes = WORKING_DIR + "genotypes/{sample}.genotypes.tsv"
    output:
        ind_het = WORKING_DIR + "genotypes/{sample}.individuals.heterozygosity_rates.tsv"
    message:
        "Calculating heterozygosity rates for {wildcards.sample} genotypes"
    threads: 1
    params:
        out_individuals_het_rates = WORKING_DIR + "genotypes/{sample}.individuals.heterozygosity_rates.tsv"
    run:
        calculate_individuals_heterozygosity_rates(
            input_genotypes_tsv=input.genotypes,
            out_individuals_het_rates_tsv=params.out_individuals_het_rates)

rule extract_list_of_individuals_below_het_threshold:
    input:
        ind_het = WORKING_DIR + "genotypes/{sample}.individuals.heterozygosity_rates.tsv"
    output:
        ind_het_below = WORKING_DIR + "genotypes/{sample}.individuals_below_heterozygosity_threshold.list"
    message:
        "Extracting individuals that are below heterozygosity threshold from {wildcards.sample} genotypes"
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

rule step6_filter_on_heterozygosity_rate:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.miss.vcf.gz",
        ind_het_below = WORKING_DIR + "genotypes/{sample}.individuals_below_heterozygosity_threshold.list"
    output:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.miss.het.vcf.gz"
    message:
        "Step6: filtering {wildcards.sample} VCF file to keep individuals below {params.het_rate_threshold} heterozygosity rate"
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
    rule step7_impute_missing_genotypes:
        input:
            vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.miss.het.vcf.gz"
        output:
            vcf = RES_DIR + "filtered/{sample}_filtered_imputed.vcf.gz"
        message:
            "Step7: imputing missing genotypes in {wildcards.sample} VCF file using Beagle v5"
        params:
            beagle_memory = config["beagle"]["memory"],
            beagle_impute = config["beagle"]["impute"]
        threads: config["threads"]
        shell:
            "beagle gt={input.vcf} "
            "out={output.vcf} "
            "nthreads={threads}"
elif config["impute_genotypes"] == "no":
    rule step7_no_imputation:
        input:
            vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.miss.het.vcf.gz"
        output:
            vcf = RES_DIR + "filtered/{sample}_filtered_no_imputed.vcf.gz"
        message:
            "Step7: not imputing missing genotypes in {wildcards.sample} VCF file; copying file from step 6 in {RES_DIR}"
        threads: config["threads"]
        shell:
            "cp {input.vcf} {output.vcf}"
else: 
    raise ValueError("The 'impute_genotypes' parameter in the config file must be either 'yes' or 'no'.")

#####################
## SNP counts metrics
#####################

rule count_at_step0_raw_snps:
    input:
        vcf = get_vcf_file
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step0.csv"
    message:
        "Counting initial number SNPs in {wildcards.sample} VCF file (all types, SNPs and indels)"
    threads: 1
    params: 
        step_name = "step0: raw file"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_after_step1_biallelic_snps: 
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step1.csv" 
    message:
        "Counting number of biallelic SNPs in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "step1: biallelic SNPs"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=threads, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_after_step2_snp_filters:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step2.csv"
    message:
        "Counting number of SNPs after first filters in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "step2: filters on SNPs"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_after_step3_frac_missing_per_snp: 
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step3.csv"
    message:
        "Counting number of SNPs after filtering on fraction missing per SNP in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "step3: filter on fraction missing per SNP"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_after_step4_maf_filter: 
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step4.csv"
    message:
        "Counting number of SNPs after MAF filter in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "step4: MAF filter"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_after_step5_frac_missing_per_genotype: 
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.miss.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step5.csv"
    message:
        "Counting number of SNPs after filtering on fraction missing per genotype in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "step5: filter on fraction missing per genotype"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_after_step6_het_excess_filter: 
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.miss.het.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step6.csv"
    message:
        "Counting number of SNPs after filtering on heterozygosity excess in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "step6: filter on heterozygosity excess"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_after_step7_imputation: 
    input:
        vcf = RES_DIR + "filtered/{sample}_filtered_imputed.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step7.csv"
    message:
        "Counting number of SNPs after imputation in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "step7: after imputation"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule merge_all_step_counts: 
    input:
        expand(WORKING_DIR + "counts/{sample}.step{step}.csv",
        sample=SAMPLES, 
        step=[0, 1, 2, 3, 4, 5, 6])
    output:
        RES_DIR + "counts/counts_merged.csv"
    message:
        "Merging all counts from different steps into a summary file"
    params: 
        out_path = RES_DIR + "counts/counts_merged.csv"
    threads: 1
    run:
        counts_df = []
        for f in input:
            df = pd.read_csv(f, index_col=0).head() 
            counts_df.append(df)
        counts_df = pd.concat(counts_df, axis=0)
        counts_df.to_csv(path_or_buf=params.out_path, index=True)

