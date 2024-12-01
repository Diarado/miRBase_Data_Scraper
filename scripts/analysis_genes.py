import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid")
summary = pd.read_csv('gh38_genes.csv')

seq_lengths = summary['length']

mean_seq_len = seq_lengths.mean()
median_seq_len = seq_lengths.median()

# Create a figure with appropriate size
plt.figure(figsize=(10, 7))

# Plot histogram with logarithmic x-axis
plt.hist(seq_lengths, bins=100, color='skyblue', edgecolor='black', alpha=0.7, log=True)

# Title and labels
plt.title('Histogram of Human Gene Lengths (hg38)', fontsize=16)
plt.xlabel('Gene Length (Base Pairs)', fontsize=14)
plt.ylabel('Number of Genes (Log Scale)', fontsize=14)

# Add mean and median lines
plt.axvline(mean_seq_len, color='red', linestyle='dashed', linewidth=2, label=f'Mean: {mean_seq_len:,.0f} bp')
plt.axvline(median_seq_len, color='green', linestyle='dashed', linewidth=2, label=f'Median: {median_seq_len:,.0f} bp')

# Add legend
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('hg38_gene_length_histogram.png', dpi=300)
plt.show()

