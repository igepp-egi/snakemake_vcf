import pandas as pd

def create_genotype_tsv(genotypes_tsv, individuals_csv, output_tsv):
    """
    Reads a genotype TSV file and an individuals CSV file, processes the data,
    and returns a DataFrame with the genotype information as a tab-sep file (.tsv).

    Parameters:
    - genotypes_tsv: Path to the genotype TSV file.
    - individuals_csv: Path to the individuals CSV file.

    Returns:
    - A .tsv DataFrame containing the processed genotype data as a .tsv file with individuals in column names and genotypes as rows.

    Example genotype TSV file format (5 rows x 5 columns):
        CHROM	POS	    IVIP100034	IVIP100081	IVIP100120	IVIP100121	IVIP100123	
        chr1S	165485	1/1 0/0	0/0	0/0	0/0	
        chr1S	165530	0/0	1/1	0/0	1/1	0/0	
        chr1S	165533	0/0	0/0	0/0	0/0	0/0	
        chr1S	167677	0/0	0/0	0/0	0/0	0/0	
    """
    # Read the genotypes TSV file
    genotypes_df = pd.read_csv(genotypes_tsv, sep="\t")
    
    # Concatenate the first two columns (CHROM and POS) to create a snp_id column and set as index
    genotypes_df['snp_id'] = genotypes_df.iloc[:, 0].astype(str) + "_" + genotypes_df.iloc[:, 1].astype(str)
    genotypes_df = genotypes_df.drop(["CHROM", "POS"], axis=1)
    genotypes_df.set_index('snp_id', inplace=True)

    # Read the individuals CSV file and set the column names to individuals
    individuals_df = pd.read_csv(individuals_csv, names=["individuals"], header=None)
    genotypes_df.columns = individuals_df.individuals.values

    # Write the processed genotype DataFrame to a .tsv file
    genotypes_df.to_csv(output_tsv, sep="\t")

 