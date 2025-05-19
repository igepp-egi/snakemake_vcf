import os
import subprocess
from snakemake.utils import min_version
import pandas as pd
from datetime import datetime 

from helpers import check_samples_tsv_file
from helpers import create_df_of_seq_length_distributions

###############################
# OS and related configurations
###############################

##### set minimum snakemake version #####
min_version("5.4.3")

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

#########################
## Pipeline configuration
##########################
configfile: "config/config.yaml"

wildcard_constraints:
  dataset="[Aa-Zz0-9]+"

# directories
WORKING_DIR = config["temp_dir"]
RES_DIR = config["result_dir"]
CURRENT_TIME = datetime.now().strftime("%Y-%m-%d_%H-%M")  # time when pipeline started (will be used to rename result directory)


# Samples: verify 
# get list of samples
# The first functions verify that the provided sample file 
samples_df = check_samples_tsv_file(sample_tsv_file = "config/samples.tsv")
SAMPLES = samples_df.index.values.tolist()

# get fastq file
def get_fastq_file(wildcards):
    fastq_file = samples_df.loc[wildcards.sample,"fastq"]
    return fastq_file

# ShortStack parameters
SHORTSTACK_PARAMS = " ".join(config["shortstack"].values())

######################
# Local rule execution
######################

# Some rules are simple jobs which should not be submitted as jobs
# The following rules will be run locally = on the head node on a HPC cluster.

localrules: all
localrules: fastp
localrules: multiqc_report
localrules: read_length_distribution
localrules: shortstack

####################
## Desired outputs
####################
QC = RES_DIR + "qc/multiqc_report.html"

SEQ_DISTRI = RES_DIR + "seq_length_distribution.tsv"

SHORTSTACK = expand(RES_DIR + "shortstack/{sample}/Results.txt",sample = SAMPLES)

rule all:
    input:
        QC,
	SEQ_DISTRI,
	SHORTSTACK
    message:
        "All done! Removing intermediate files in {WORKING_DIR}. Adding date and current time to {RES_DIR} folder name"
    params:
        new_result_dir_name =  CURRENT_TIME + "_" + RES_DIR
    shell:
        "rm -rf {WORKING_DIR};" # removes unwanted intermediate files
        "mv {RES_DIR} {params.new_result_dir_name};"


######################
## Shortstack analysis
######################

rule shortstack:
    input:
        reads =  RES_DIR + "trimmed/{sample}.trimmed.fastq"
    output:
        RES_DIR + "shortstack/{sample}/Results.txt"
    message:"Shortstack analysis of {wildcards.sample} using {params.genome} reference"
    params:
        resdir = RES_DIR + "shortstack/{sample}/",
        genome = lambda wildcards: samples_df.loc[wildcards.sample,"genome"]
    threads: 20
    shell:
        "ShortStack "
        "--outdir {wildcards.sample} "
        "--threads {threads} "
        "{SHORTSTACK_PARAMS} "
        "--readfile {input.reads} "
	"--dn_mirna "
        "--genome {params.genome};"
        "cp -r {wildcards.sample}/* {params.resdir};"
        "rm -r {wildcards.sample};"

###########################################################
## Get read length distribution (before and after trimming)
##########################################################

rule read_length_distribution:
    input: 
        expand(RES_DIR + "trimmed/{sample}.trimmed.fastq", sample = SAMPLES)
    output:
        RES_DIR + "seq_length_distribution.tsv"
    message: 
        "Computing sequence length distribution for all samples"
    params:
        path_to_fastq_files = RES_DIR + "trimmed/"
    run:
        create_df_of_seq_length_distributions(path_to_fastq_files =  params.path_to_fastq_files, outfile = output[0])

###########################################
## Trim reads for all samples and QC report
###########################################
rule multiqc_report:
    input:
        expand(WORKING_DIR + "trimmed/{sample}_fastp.json", sample = SAMPLES)
    output:
        RES_DIR + "qc/multiqc_report.html"
    message:
        "Compiling QC reports from fastp"
    params:
        input_directory = WORKING_DIR + "trimmed/",
        output_directory = RES_DIR + "qc/"
    shell:
        "multiqc "
        "--force " # force directory to be created
        "--outdir {params.output_directory} "
        "{params.input_directory}"

rule fastp:
    input:
        get_fastq_file
    output:
        fastq = RES_DIR + "trimmed/{sample}.trimmed.fastq",
        json = WORKING_DIR + "trimmed/{sample}_fastp.json",
        html = WORKING_DIR + "trimmed/{sample}_fastp.html"
    message:
        "trimming {wildcards.sample} reads on quality and adapter presence"
    threads: 10
    resources:
        mem_mb = 1000
    params:
        adapters_fasta =                config["fasta_adapters"],
        qualified_quality_phred =       config["fastp"]["qualified_quality_phred"],
        average_quality =               config["fastp"]["average_quality"]
    shell:
        "fastp -i {input} "
        "--stdout "
        "--json {output.json} "
        "--html {output.html} "
        "--qualified_quality_phred {params.qualified_quality_phred} "
        "--average_qual {params.average_quality} "
        "--adapter_fasta {params.adapters_fasta} > {output.fastq}"