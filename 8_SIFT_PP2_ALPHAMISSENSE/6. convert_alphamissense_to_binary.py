import pandas as pd
import numpy as np

# Load the datasets
sp = pd.read_csv('./Final_results/Evaluated_vars_from_strict_pairs(alphamissense).txt')
si = pd.read_csv('./Final_results/Evaluated_vars_from_similar_initial_AA_pairs(alphamissense).txt')
sf = pd.read_csv('./Final_results/Evaluated_vars_from_similar_final_AA_pairs(alphamissense).txt')
sif = pd.read_csv('./Final_results/Evaluated_vars_from_similar_i+f_AA_pairs(alphamissense).txt')

def add_verify_sift_column(df):
    # Define conditions
    conditions = [
        (df['ALPHAMISSENSE_prediction'] == 'pathogenic') & (df['db'] == 'gnomAD'),
        (df['ALPHAMISSENSE_prediction'] == 'pathogenic') & (df['db'] == 'ClinVar'),
        (df['ALPHAMISSENSE_prediction'] == 'benign') & (df['db'] == 'ClinVar'),
        (df['ALPHAMISSENSE_prediction'] == 'benign') & (df['db'] == 'gnomAD'),
        (df['ALPHAMISSENSE_prediction'] == ' '),
        (df['ALPHAMISSENSE_prediction'] == 'ambiguous') & (df['db'] == 'ClinVar'),
        (df['ALPHAMISSENSE_prediction'] == 'ambiguous') & (df['db'] == 'gnomAD'),
    ]

    # Define corresponding values
    values = [0, 1, 0, 1, np.nan, 0, 0]

    # Apply conditions to create the 'Verify_ALPHAMISSENSE' column
    df['Verify_ALPHAMISSENSE'] = np.select(conditions, values, default=np.nan)
    return df

# List of DataFrames and their file paths
datasets = {
    './Final_results/Evaluated_vars_from_strict_pairs(alphamissense).txt': sp,
    './Final_results/Evaluated_vars_from_similar_initial_AA_pairs(alphamissense).txt': si,
    './Final_results/Evaluated_vars_from_similar_final_AA_pairs(alphamissense).txt': sf,
    './Final_results/Evaluated_vars_from_similar_i+f_AA_pairs(alphamissense).txt': sif
}

# Apply the function and save each DataFrame to its respective file
for file_path, df in datasets.items():
    updated_df = add_verify_sift_column(df)
    updated_df.to_csv(file_path, index=False)
    print(f"Processed and saved DataFrame to {file_path}.")
