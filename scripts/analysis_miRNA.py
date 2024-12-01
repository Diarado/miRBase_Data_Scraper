import pandas as pd
import matplotlib.pyplot as plt

# Load data
summary = pd.read_csv('miRNA_human.csv')
# summary = pd.read_csv('miRNA_Ecoli.csv')

print("First few rows of the data:")
print(summary.head())

# Check the data types of the columns
print("\nData types before conversion:")
print(summary.dtypes)

# Convert 'End' and 'Start' columns to numeric types
# `errors='coerce'` will convert non-convertible values to NaN
summary['End'] = pd.to_numeric(summary['End'], errors='coerce')
summary['Start'] = pd.to_numeric(summary['Start'], errors='coerce')

# Check for any NaN values introduced by conversion
nan_end = summary['End'].isna().sum()
nan_start = summary['Start'].isna().sum()

if nan_end > 0 or nan_start > 0:
    print(f"\nWarning: Found {nan_end} NaN values in 'End' and {nan_start} NaN values in 'Start'. These rows will be dropped.")
    # Drop rows with NaN in 'End' or 'Start'
    summary = summary.dropna(subset=['End', 'Start'])

# Verify row count after dropping NaNs
print(f"\nRow count after dropping NaNs: {len(summary)}")
# assert len(summary) <= 1917, "Row count exceeds expected after dropping NaNs."

# Proceed with calculations
ai_cnt = summary['Summary'].notna().sum()
non_ai_cnt = len(summary) - ai_cnt

# Calculate sequence lengths
seq_lengths = summary['End'] - summary['Start']
mean_seq_len = seq_lengths.mean()
median_seq_len = seq_lengths.median()

# Verify row count matches expected (optional, adjust as needed)
# assert len(summary) <= 1917, "Row count mismatch after processing."

# Print statistics
print(f"\nAI count (rows with summaries): {ai_cnt}")    
print(f"Non-AI count (rows without summaries): {non_ai_cnt}")
print(f"Mean sequence length: {mean_seq_len:.2f}")
print(f"Median sequence length: {median_seq_len:.2f}")

# Plot histogram of sequence lengths
plt.figure(figsize=(8, 6))
plt.hist(seq_lengths, bins=30, color='skyblue', edgecolor='black')
plt.title('Histogram of miRNA Sequence Lengths (human)')
plt.xlabel('Sequence Length')
plt.ylabel('Frequency')

# Optional: Add lines for mean and median
plt.axvline(mean_seq_len, color='red', linestyle='dashed', linewidth=1, label=f'Mean: {mean_seq_len:.2f}')
plt.axvline(median_seq_len, color='green', linestyle='dashed', linewidth=1, label=f'Median: {median_seq_len:.2f}')
plt.legend()

# Show the plot
plt.tight_layout()
plt.savefig('human_miRNA_distribution.png', dpi=300)
plt.show()

