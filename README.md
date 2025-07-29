# Snakemake VCF

A small pipeline to filter a VCF file and produce QC plots before/after filterings

## Usage

1. Create the conda environment to install all required dependencies.   

`conda create --name snakevcf --file environment.yaml`

2. Activate the environment

`conda activate snakevcf`  

3. On the SURM cluster, run: 

`screen -S vcf` # Make sure the session does not log out. 

`srun --time=24:00:00 --cpus-per-task=10 --pty bash -i` (connect to a compute node + allocate resources). 

`snakemake -j 10` # start the pipeline.   

# Softwares/dependencies 

## bcftools 

See https://samtools.github.io/bcftools/bcftools.html

## vcftools 

See https://manpages.ubuntu.com/manpages/xenial/man1/vcftools.1.html

