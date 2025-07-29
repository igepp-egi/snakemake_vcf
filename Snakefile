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
FILTERED_VCF = expand(RES_DIR + "filtered/{sample}.vcf.gz", sample = SAMPLES)
ALL_COUNTS =  RES_DIR + "counts/counts_merged.csv"
GENOTYPES = expand(RES_DIR + "genotypes/{sample}.genotypes.tsv", sample = SAMPLES)
HET_RATES = expand(RES_DIR + "genotypes/{sample}.{type}.heterozygosity_rates.tsv",
                  sample = SAMPLES, 
                  type = ["individuals"])

#BED = expand(WORKING_DIR + "plink/{sample}_{status}.lmiss", sample = SAMPLES, status = ["raw", "filtered"])

if config["keep_temp_dir"] == True:
    rule all:
        input:
            FILTERED_VCF, 
            ALL_COUNTS,
            GENOTYPES, 
            HET_RATES
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
        threads: 20
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
        vcf = WORKING_DIR + "filtered/{sample}.selected.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.selected.biallelic.vcf.gz"
    message:
        "Step1: Keeping only biallelic SNPs in {wildcards.sample} VCF file"
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
    threads: 20
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
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc2.vcf.gz"
    message:
        "Step3: filtering {wildcards.sample}  biallelic VCF file to keep SNP with less than {params.max_missing_fraction_per_snp} fraction missing"
    params:
        max_missing_fraction_per_snp = config["bcftools"]["snp"]["max_missing_fraction_per_snp_site"],
        miss_prefix = WORKING_DIR + "filtered/{sample}.missing_snp",
        miss_file = WORKING_DIR + "filtered/{sample}.missing_snp.lmiss",
        snps_to_keep = WORKING_DIR + "filtered/{sample}.snp_sites_to_keep.txt"
    threads: 20
    shell:
        "vcftools --gzvcf {input.vcf} --missing-site --out {params.miss_prefix}; "
        "awk -v OFS='\t' '{{if ($5 < {params.max_missing_fraction_per_snp}) print $1, $2}}' {params.miss_file} > {params.snps_to_keep}; "
        # filter the VCF file
        "bcftools index {input.vcf}; "
        "bcftools view -R {params.snps_to_keep} "
        "{input.vcf} "
        "-Oz "
        "-o {output.vcf}"

############################
## Filter on MAF, F_MISSING 
############################

rule step4_filter_on_maf_before_imputation: 
    input:
        WORKING_DIR + "filtered/{sample}.selected.biallelic.qc2.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.selected.biallelic.qc2.maf.vcf.gz"
    message:
        "Step4: filtering {wildcards.sample} biallelic VCF file on MAF higher than {params.min_maf}; before imputation"
    params:
        min_maf = config["bcftools"]["individuals"]["min_maf_before_imputation"]
    threads: 20
    shell:
        "bcftools view -i 'MAF > {params.min_maf}' --threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}"

rule step5_filter_fraction_missing_per_genotype:
    input:
        WORKING_DIR + "filtered/{sample}.selected.biallelic.qc2.maf.vcf.gz"
    output:
        WORKING_DIR + "filtered/{sample}.selected.biallelic.qc2.maf.miss.vcf.gz"
    message:
        "Step5: filtering {wildcards.sample} biallelic VCF file on percentage of missing genotype calls"
    params:
        missing = config["bcftools"]["individuals"]["max_missing_fraction_per_genotype"]
    threads: 20
    shell:
        "bcftools view -i 'F_MISSING < {params.missing}' --threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}"

###########################################
## Extract genotypes and allele frequencies
## Filtering on HWE and heterozygote excess
###########################################

rule extract_individuals_and_genotypes:
    input: 
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc2.maf.miss.vcf.gz"
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

rule extract_list_of_individuals_heterozygosity_excess:
    input:
        ind_het = WORKING_DIR + "genotypes/{sample}.individuals.heterozygosity_rates.tsv"
    output:
        ind_het_excess = WORKING_DIR + "genotypes/{sample}.individuals.heterozygosity_excess.list"
    message:
        "Extracting individuals with heterozygosity excess from {wildcards.sample} genotypes"
    params:
        max_het = config["bcftools"]["individuals"]["max_heterozygosity"]
    threads: 20
    shell:
        "awk -v OFS='\\t' '{{if ($3 < {params.max_het}) print $1}}' {input.ind_het} > {output.ind_het_excess}"

rule step6_filter_on_heterozygosity_excess:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc2.maf.miss.vcf.gz",
        ind_het_excess = WORKING_DIR + "genotypes/{sample}.individuals.heterozygosity_excess.list"
    output:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc2.maf.miss.het.vcf.gz"
    message:
        "Step6: filtering {wildcards.sample} VCF file on heterozygosity excess"
    params:
        het_excess = WORKING_DIR + "genotypes/{sample}.individuals.heterozygosity_excess.list"
    threads: 20
    shell:
        "bcftools view -S {input.ind_het_excess} --threads {threads} "
        "{input.vcf} "
        "-Oz "
        "-o {output.vcf}"

###########################################
## Impute missing genotypes using Beagle v5
###########################################

rule step7_impute_missing_genotypes:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc2.maf.miss.het.vcf.gz"
    output:
        vcf = RES_DIR + "filtered/{sample}.filtered.imputed.vcf.gz"
    message:
        "Step7: imputing missing genotypes in {wildcards.sample} VCF file using Beagle v5"
    params:
        beagle_memory = config["beagle"]["memory"]
        beagle_impute = config["beagle"]["impute"]
    threads: 20
    shell:
        "beagle gt={input.vcf} "
        "out={output.vcf} "
        "jar={params.beagle_jar} "
        "{params.memory} "
        "--threads {threads}"

#####################
## SNP counts metrics
#####################

rule count_original_snps_by_types:
    input:
        vcf = get_vcf_file
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step0.csv"
    message:
        "Counting initial number SNPs in {wildcards.sample} VCF file (all types, SNPs and indels)"
    threads: 10
    params: 
        step_name = "step0: raw file"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_biallelic_snps: 
    input:
        vcf = WORKING_DIR + "filtered/{sample}.biallelic.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step1.csv"
    message:
        "Counting number of biallelic SNPs in {wildcards.sample} VCF file"
    threads: 10
    params: 
        step_name = "step1: biallelic SNPs"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=threads, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_after_first_filters:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.biallelic.qc1.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step2.csv"
    message:
        "Counting number of SNPs after first filters in {wildcards.sample} VCF file"
    threads: 10
    params: 
        step_name = "step2: filters on SNPs"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_after_filter_fraction_missing_per_snp: 
    input:
        vcf = WORKING_DIR + "filtered/{sample}.biallelic.qc2.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step3.csv"
    message:
        "Counting number of SNPs after filtering on fraction missing per SNP in {wildcards.sample} VCF file"
    threads: 10
    params: 
        step_name = "step3: filter on fraction missing per SNP"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_after_maf_filter: 
    input:
        vcf = WORKING_DIR + "filtered/{sample}.biallelic.qc2.selected.maf.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step4.csv"
    message:
        "Counting number of SNPs after MAF filter in {wildcards.sample} VCF file"
    threads: 10
    params: 
        step_name = "step4: MAF filter"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_after_fraction_missing_per_genotype: 
    input:
        vcf = RES_DIR + "filtered/{sample}.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step5.csv"
    message:
        "Counting number of SNPs after filtering on fraction missing per genotype in {wildcards.sample} VCF file"
    threads: 10
    params: 
        step_name = "step5: filter on fraction missing per genotype"
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
        step=[0, 1, 2, 3, 4, 5])
    output:
        RES_DIR + "counts/counts_merged.csv"
    message:
        "Merging all counts from different steps into a summary file"
    params: 
        out_path = RES_DIR + "counts/counts_merged.csv"
    run:
        counts_df = []
        for f in input:
            df = pd.read_csv(f, index_col=0).head() 
            counts_df.append(df)
        counts_df = pd.concat(counts_df, axis=0)
        counts_df.to_csv(path_or_buf=params.out_path, index=True)

###############################
# QC before and after filtering
###############################

rule convert_vcf_files_to_plink_format:
    input:
        raw_vcf = get_vcf_file,
        filtered_vcf = RES_DIR + "filtered/{sample}.vcf.gz"
    output:
        plink_raw = WORKING_DIR + "plink/{sample}_raw.bed",
        plink_filtered = WORKING_DIR + "plink/{sample}_filtered.bed"
    message:
        "Preparing PLINK files from {wildcards.sample} VCF files"
    threads: 20
    params: 
        plink_prefix_raw = WORKING_DIR + "plink/{sample}_raw",
        plink_prefix_filtered = WORKING_DIR + "plink/{sample}_filtered"
    shell:
        "mkdir -p {WORKING_DIR}/plink/; "
        "plink --vcf {input.raw_vcf}      --make-bed --out {params.plink_prefix_raw} --allow-extra-chr; "
        "plink --vcf {input.filtered_vcf} --make-bed --out {params.plink_prefix_filtered} --allow-extra-chr;"

rule plink_compute_snp_and_genotype_missing_rate: 
    input:
        plink = WORKING_DIR + "plink/{sample}_{status}.bed"
    output:
        snp_missing = WORKING_DIR + "plink/{sample}_{status}.lmiss",
        geno_missing = WORKING_DIR + "plink/{sample}_{status}.imiss"
    message:
        "Computing SNP and genotype missing rates for {wildcards.sample} {wildcards.status} PLINK files"
    params: 
        plink_prefix = WORKING_DIR + "plink/{sample}_{status}"
    threads: 20
    shell:
        "plink --bfile {params.plink_prefix} --missing --out {params.plink_prefix} --allow-extra-chr"

rule parse_and_plot_snp_and_genotype_missing_rates: 
    input:
        snp_missing = WORKING_DIR + "plink/{sample}_{status}.lmiss",
        geno_missing = WORKING_DIR + "plink/{sample}_{status}.imiss"
    output:
       snp_plot = RES_DIR + "plots/{sample}_{status}.snp_missing.png",
       geno_plot = RES_DIR + "plots/{sample}_{status}.genotype_missing.png"
    message:
        "Parsing SNP and genotype missing rates for {wildcards.sample} {wildcards.status} PLINK files"
    threads: 20
    params: 
        outdir = RES_DIR + "plots/"
    shell:
        "awk -v OFS='\t' '{{print $1, $2, $3, $4, $5}}' {input.snp_missing} > {output.snp_missing_parsed}; "
        "awk -v OFS='\t' '{{print $1, $2, $3, $4, $5}}' {input.geno_missing} > {output.geno_missing_parsed}; "
        "python utils/plot_snp_missing_rates.py --outdir {params.outdir} --input {input.snp_missing} {output.snp_plot}; "
        "python utils/plot_genotype_missing_rates.py --outdir {params.outdir} --input {input.geno_missing} {output.geno_plot}; "


#####################
# Extra draft section
#####################

#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' vfaba_hedin_peamust_SNP_GATK.vcf.gz > allele_frequencies.txt

# Filter to keep sample sequenced to a certain depth
"""     shell:
        "bcftools view --include 'MEAN(FMT/DP) >= {params.min_sample_depth}' "
        "--threads {threads} "
        "{input} "
        "-Oz "
        "-o {output}" """

# Convert VCF to PLINK format
# plink --make-bed --vcf vfaba_hedin_peamust_SNP_GATK.vcf.gz --out vfaba --allow-extra-chr

# Extract allele frequencies from the VCF file
""" CHR	Chromosome code
SNP	Variant identifier
A1	Allele 1 (usually minor)
A2	Allele 2 (usually major)
C(HOM A1)	A1 homozygote count
C(HET)	Heterozygote count
C(HOM A2)	A2 homozygote count
C(HAP A1)	Haploid A1 count (includes male X chromosome)
C(HAP A2)	Haploid A2 count
C(MISSING)	Missing genotype count """
# plink --freqx --out vfaba_freq --bfile vfaba --allow-extra-chr