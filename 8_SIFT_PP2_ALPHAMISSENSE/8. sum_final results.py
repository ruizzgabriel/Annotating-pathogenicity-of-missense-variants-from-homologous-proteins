import pandas as pd

# Open results
sp = pd.read_csv('./Final_results/Evaluated_vars_from_strict_pairs(alphamissense).txt')
si = pd.read_csv('./Final_results/Evaluated_vars_from_similar_initial_AA_pairs(alphamissense).txt')
sf = pd.read_csv('./Final_results/Evaluated_vars_from_similar_final_AA_pairs(alphamissense).txt')
sif = pd.read_csv('./Final_results/Evaluated_vars_from_similar_i+f_AA_pairs(alphamissense).txt')

# Function to calculate sums and total rows for each DataFrame
def calculate_sums_and_totals(df):
    return {
        'Sum of Verify_SIFT': df['Verify_SIFT'].sum(skipna=True),
        '% Success SIFT': (df['Verify_SIFT'].sum(skipna=True) / len(df)) * 100,
        'Sum of Verify_POLYPHEN': df['Verify_POLYPHEN'].sum(skipna=True),
        '% Success POLYPHEN': (df['Verify_POLYPHEN'].sum(skipna=True) / len(df)) * 100,
        'Sum of Verify_ALPHAMISSENSE': df['Verify_ALPHAMISSENSE'].sum(skipna=True),
        '% Success ALPHAMISSENSE': (df['Verify_ALPHAMISSENSE'].sum(skipna=True) / len(df)) * 100,
        'Total': len(df),  # Total rows in the DataFrame
    }

# Calculate sums and totals for each DataFrame
results = {
    'Strict Pairs': calculate_sums_and_totals(sp),
    'Similar Initial AA Pairs': calculate_sums_and_totals(si),
    'Similar Final AA Pairs': calculate_sums_and_totals(sf),
    'Similar i+f AA Pairs': calculate_sums_and_totals(sif)
}

# Create a summary DataFrame
summary_df = pd.DataFrame(results).T  # Transpose to make the DataFrame more readable

# Round numeric values to 2 decimal places
summary_df = summary_df.round(2)

# Display the summary DataFrame
print(summary_df)

# Save the summary DataFrame as a CSV file
summary_df.to_csv('./Final_results/summary_results.txt', index=True)