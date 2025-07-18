import pandas as pd

def calculate_heterozygosity_rates(input_genotypes_tsv, out_snp_het_rates_tsv, out_individuals_het_rates_tsv):
    """
    Calculate the heterozygosity rates for each SNP in the genotype TSV file.
    
    Parameters:
    - genotypes_tsv: Path to the genotype TSV file.
    - snp_het_rates_tsv: Path to save the SNP heterozygosity rates TSV file.
    - individuals_het_rates_tsv: Path to save the individuals heterozygosity rates TSV file.

    Returns:
    - A .tsv file containing the heterozygosity rates for each SNP.
    - A .tsv file containing the heterozygosity rates for each individual.

    Example input for genotype TSV file format (5 rows x 5 columns):
        CHROM	POS	    IVIP100034	IVIP100081	IVIP100120	IVIP100121	IVIP100123	
        chr1S	165485	1/1 0/0	0/0	0/0	0/0	
        chr1S	165530	0/0	1/1	0/0	1/1	0/0	
        chr1S	165533	0/0	0/0	0/0	0/0	0/0	
        chr1S	167677	0/0	0/0	0/0	0/0	0/0	
    
    Example output for SNP heterozygosity TSV file format:
        snp_id	heterozygosity_rate
        chr1S_165485	0.0
        chr1S_165530	0.5
        chr1S_165533	0.0
        chr1S_167677	0.0

    Example output for individuals heterozygosity TSV file format:
        individuals	heterozygosity_rate
        IVIP100034	0.0
        IVIP100081	0.5
        IVIP100120	0.0
        IVIP100121	0.5
        IVIP100123	0.0
    
    """
    # Read the genotype TSV file
    genotypes_df = pd.read_csv(input_genotypes_tsv, sep="\t", index_col=0)

    # If cell content is equal to 0/0 or 0|0 # then replace it with 0
    # if cell content is equal to 1/1 or 1|1 # then replace it with 2
    # if cell content is equal to 0/1 or 1/0 or 1|0 or 0|1 # then replace it with 1
    genotypes_df = genotypes_df.replace(
        {"0/0": 0,
         "0|0": 0, 
         "1/1": 2, 
         "1|1": 2, 
         "0/1": 1, 
         "1/0": 1, 
         "1|0": 1, 
         "0|1": 1})

    # For each SNP, sum the number of individuals that are homozygous (0 or 2) and heterozygous (1); calculate the heterozygosity rate
    genotypes_df['homozygous'] = genotypes_df.apply(lambda x: (x == 0).sum() + (x == 2).sum(), axis=1)
    genotypes_df['heterozygous'] = genotypes_df.apply(lambda x: (x == 1).sum(), axis=1)
    genotypes_df['heterozygosity_rate'] = genotypes_df['heterozygous'] / (genotypes_df['homozygous'] + genotypes_df['heterozygous'])    
    snp_het_rates = genotypes_df.drop(columns=genotypes_df.columns[:-3])  # remove the indivuals columns
    snp_het_rates.to_csv(out_snp_het_rates_tsv, sep="\t")
    
    # calculate the heterozygosity rate for each individual
    individuals_het_rates = genotypes_df.drop(columns=['homozygous', 'heterozygous', 'heterozygosity_rate'])

    # per column, calculate the heterozygosity rate
    # axis = 0 means apply function to each column
    individuals_het_rates["homozygous"] = individuals_het_rates.apply(lambda x: (x == 0).sum() + (x == 2).sum(), axis=0)
    individuals_het_rates["heterozygous"] = individuals_het_rates.apply(lambda x: (x == 1).sum(), axis=0)
    # calculate the heterozygosity rate
    individuals_het_rates["heterozygosity_rate"] = individuals_het_rates['heterozygous'] / (individuals_het_rates['homozygous'] + individuals_het_rates['heterozygous'])

    individuals_het_rates.to_csv(out_individuals_het_rates_tsv, sep="\t")