# Snakemake VCF

A small pipeline to filter a VCF file and produce QC plots before/after filterings


## TODOs
1. Add imputation with Beagle v5
2. Add a filtering step on the MAF after imputation

## Configuration

All configuration files are in the `config/` folder.

1. **VCF files:** the vcf files to be processed should be indicated in `samples.tsv`.  
2. **Individuals:** the `individuals.tsv` file in the `config/` folder contains the list of individuals to keep.  
3. **Chromosomes:** the `chromosomes.tsv` file in the `config/` folder contains the list of chromosomes to keep with a one-column format: CHROM.    
4. **Parameters for Snakemake, vcftools and bcftools** can be found in [`config.yaml`](config/config.yaml). Details are provided in the file itself. 

### Chromosomes 

To extract the list of chromosomes and contigs from your VCF file, use: `bcftools query -f '%CHROM\n' yourfile.vcf.gz | sort -u`.

### Genotype imputation with Beagle v5

It is possible to specify two additional steps by specifying "yes" in the `impute_genotype` field in the `config.yaml` file (only "yes" or "no" is allowed).   
When `impute_genotype: "yes"` is specified, the VCF file will undergo two additional steps:   
- Genotype imputation with Beagle (related parameters are available in the config file in the "beagle" section).  
- An additional filter on Minor Allele Frequency (e.g. "MAF > 0.05") is performed.
The pipeline performs 10 steps in total when imputing genotypes and 8 steps when no imputation is performed. 


## Usage

1. Create the conda environment to install all required dependencies: `conda create --name snakevcf --file environment.yaml`
2. Activate the environment: `conda activate snakevcf`  
3. On the SURM cluster, run: 
- `screen -S vcf` # To make sure the pipeline does not crash because of a lost connection. 
- `srun --time=24:00:00 --cpus-per-task=10 --pty bash -i` (connect to a compute node + allocate resources). 
- `snakemake -j 10` # start the pipeline.   

## Softwares/dependencies 

1. bcftools: see https://samtools.github.io/bcftools/bcftools.html
2. vcftools: see https://manpages.ubuntu.com/manpages/xenial/man1/vcftools.1.html
3. beagle 5: see https://faculty.washington.edu/browning/beagle/beagle.html

## Pipeline structure

The pipeline steps can be vizualised using: `snakemake --rulegraph | dot -Tpdf > dag.pdf`

[DAG graph](./dag.pdf)
