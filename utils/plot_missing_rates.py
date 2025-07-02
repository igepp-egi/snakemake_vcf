import pandas as pd
import seaborn as sns   
import matplotlib.pyplot as plt

def plot_snp_missing_rates(outdir, raw_vcf, filtered_vcf, output_plot):
    """
    Plot the fraction of missing SNPs from a TSV file.
    
    Args:
        outdir (str): Output directory for the plot.
        raw_vcf (str): Path to the raw VCF file.
        filtered_vcf (str): Path to the filtered VCF file.
        output_plot (str): Path to save the output plot.

    Returns:
        None (a plot is saved to the specified output path)
    """
    # Read the VCF files
    raw = pd.read_csv(raw_vcf, sep="\t")
    filtered = pd.read_csv(filtered_vcf, sep="\t")
    
    # Plot the distribution of missing SNPs side by side raw and filtered
    fig, axes = plt.subplots(1, 2, figsize=(18, 10))

    # first plot for raw unfiltered SNPs
    sns.histplot(raw['F_MISS'], bins=50, kde=True, ax=axes[0], color='blue')
    axes[0].set_title('Raw SNPs')
    axes[0].set_xlabel('F_MISS')
    axes[0].set_xlim(0, 1)

    # second plot for filtered SNPs
    sns.histplot(filtered['F_MISS'], bins=50, kde=True, ax=axes[1], color='orange')
    axes[1].set_title('Filtered SNPs')
    axes[1].set_xlabel('F_MISS')
    axes[1].set_xlim(0, 1)

    plt.tight_layout() 
    fig.savefig(output_plot, dpi=300)
    plt.close(fig)