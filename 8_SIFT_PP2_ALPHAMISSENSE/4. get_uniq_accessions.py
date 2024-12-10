import pandas as pd

# Open results
sp = pd.read_csv('./Final_results/Evaluated_vars_from_strict_pairs.txt')
sp = sp[['Entry']]
si = pd.read_csv('./Final_results/Evaluated_vars_from_similar_initial_AA_pairs.txt')
si = si[['Entry']]
sf = pd.read_csv('./Final_results/Evaluated_vars_from_similar_final_AA_pairs.txt')
sf = sf[['Entry']]
sif = pd.read_csv('./Final_results/Evaluated_vars_from_similar_i+f_AA_pairs.txt')
sif = sif[['Entry']]

# Step 1: Concatenate all DataFrames
combined_df = pd.concat([sp, si, sf, sif])

# Step 2: Get unique entries in the 'Entry' column
unique_entries_df = combined_df.drop_duplicates(subset='Entry').reset_index(drop=True)

# Display the resulting DataFrame
print(unique_entries_df)

# Saving the merged DataFrame to a CSV file
unique_entries_df.to_csv('./Alphamissense/Unique_accession_from_pairs.txt', index=False)