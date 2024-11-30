import pandas as pd

# Load data
summary = pd.read_csv('start_end_summary_data.csv')

# Count summaries that are non-empty (using notna() since empty strings appear as '')
ai_cnt = summary['Summary'].notna().sum()
non_ai_cnt = len(summary) - ai_cnt

# Calculate sequence lengths
seq_lengths = summary['End'] - summary['Start']
mean_seq_len = seq_lengths.mean()
median_seq_len = seq_lengths.median()

# Verify row count
assert len(summary) == 1917

print(f"AI count (rows with summaries): {ai_cnt}")    
print(f"Non-AI count (rows without summaries): {non_ai_cnt}")
print(f"Mean sequence length: {mean_seq_len}")
print(f"Median sequence length: {median_seq_len}")