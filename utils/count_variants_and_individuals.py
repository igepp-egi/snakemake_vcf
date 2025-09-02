from cyvcf2 import VCF
import pandas as pd

def count_variants_by_type(vcf_file_path, n_threads=1, step_name="step01"):
    """
    Count the number of variants for each type in a VCF file using cyvcf2.
    Args:
        vcf_file_path (str): Path to the VCF file.
        n_threads (int): Number of threads to use for reading the VCF file.
        step_name (str): Name of the step used in the column "step".
    Returns:
        Pandas DataFrame: A DataFrame with the counts of SNPs and InDels.
    """
    vcf_reader = VCF(vcf_file_path, threads=n_threads)
    snp_variant_count = 0
    indel_variant_count = 0
    sv_variant_count = 0

    for variant in vcf_reader:
        # Check if the variant is a SNP
        if len(variant.REF) == 1 and len(variant.ALT[0]) == 1:
            snp_variant_count += 1
        # Check if the variant is an indel
        elif len(variant.REF) != len(variant.ALT[0]):
            indel_variant_count += 1
        else:   
            pass # SVs are not counted in this function

    vcf_reader.close()

    # Create a DataFrame to hold the counts
    variant_counts = {
        "step": [step_name],
        "snp_variant_count": [snp_variant_count],
        "indel_variant_count": [indel_variant_count],
        "sv_variant_count": [sv_variant_count]
    }
    variant_counts_df = pd.DataFrame.from_dict(variant_counts, orient='columns')
    return variant_counts_df

# Count the number of individuals in a VCF file using cyvcf2
def count_individuals(vcf_file, n_threads=1, step_name="step01"):
    """
    Count the number of individuals (samples) in a VCF file using cyvcf2
    
    Args:
        vcf_file (str): Path to the VCF file
    
    Returns:
        int: Number of individuals/samples in the VCF file
    """
    vcf_reader = VCF(vcf_file, threads=n_threads)
    # create a pandas df to hope the results
    res_dict = {'step': step_name, 'n_individuals': len(vcf_reader.samples)}
    df = pd.DataFrame.from_records([res_dict])
    return df