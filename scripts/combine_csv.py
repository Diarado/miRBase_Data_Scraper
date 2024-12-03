import pandas as pd

# Define organisms and column names
organisms = ['human', 'mouse', 'rat', 'fly', 'worm', 'Arabidopsis']
columns = ['Organism', 'Name', 'Start', 'End', 'Sequence', 'Seq_len', 
           'Summary', 'Summary_len', 'Chr', 'Strand']

# Initialize empty list to store dataframes
all_dfs = []

# Process each organism's CSV file
for organism in organisms:
    csv_file = f'results/miRNA_{organism}_with_sequence.csv'
    try:
        # Read CSV and add organism column
        df = pd.read_csv(csv_file)
        df.insert(0, 'Organism', organism)
        all_dfs.append(df)
    except FileNotFoundError:
        print(f"Warning: {csv_file} not found")

# Combine all dataframes
if all_dfs:
    combined_df = pd.concat(all_dfs, ignore_index=True)
    # Ensure columns are in correct order
    combined_df = combined_df[columns]
    # Save combined data
    combined_df.to_csv('miRNA_combined.csv', index=False)
    print("Successfully created miRNA_combined.csv")
else:
    print("No data files were found to combine")
    
    