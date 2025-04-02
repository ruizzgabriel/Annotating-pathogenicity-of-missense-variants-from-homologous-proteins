"""
This script compares Pfam domains between pathogenic and non-pathogenic genes.
It handles multiple Pfam codes per gene (separated by semicolons) and creates
a detailed comparison output including gene lists and counts.
"""

import pandas as pd

# File paths
PATHO_FILE = './Domains/patho_genes.txt'
NON_PATHO_FILE = './Domains/non_patho.txt'
OUTPUT_FILE = './Domains/pfam_comparison.xlsx'

def load_and_preview_data(file_path, sep='\t'):
    """Load a tab-separated file and display its first few rows."""
    df = pd.read_csv(file_path, sep=sep)
    print(f"Preview of {file_path}:")
    print(df.head())
    return df

def get_unique_pfam_codes(df, column='Pfam'):
    """Extract unique Pfam codes from a column, handling semicolon-separated values."""
    # Split semicolon-separated values and get unique codes
    pfam_series = df[column].str.split(';').explode().str.strip()
    return set(pfam_series.dropna())

def count_pfam_occurrences(df, pfam, column='Pfam'):
    """Count occurrences of a specific Pfam code in a dataframe."""
    # Count rows where pfam appears in the semicolon-separated list
    return df[column].fillna('').str.contains(pfam, regex=False).sum()

def get_genes_with_pfam(df, pfam, gene_col, pfam_col='Pfam'):
    """Get list of genes containing a specific Pfam code."""
    mask = df[pfam_col].fillna('').str.contains(pfam, regex=False)
    genes = df.loc[mask, gene_col].dropna()
    return ';'.join(genes)

def main():
    """Main function to compare Pfam domains between pathogenic and non-pathogenic genes."""
    # Load data
    patho = load_and_preview_data(PATHO_FILE)
    non_patho = load_and_preview_data(NON_PATHO_FILE)

    # Get unique Pfam codes from both datasets
    patho_pfam = get_unique_pfam_codes(patho)
    non_patho_pfam = get_unique_pfam_codes(non_patho)

    # Find common Pfam codes
    common_pfam = patho_pfam.intersection(non_patho_pfam)

    # Create results dictionary
    pfam_counts = {}
    for pfam in common_pfam:
        pfam_counts[pfam] = {
            'patho_count': count_pfam_occurrences(patho, pfam),
            'non_patho_count': count_pfam_occurrences(non_patho, pfam),
            'patho_genes': get_genes_with_pfam(patho, pfam, 'Gene Names'),
            'non_patho_genes': get_genes_with_pfam(non_patho, pfam, 'Gene Names')
        }

    # Create and format result dataframe
    result_df = pd.DataFrame.from_dict(pfam_counts, orient='index')
    result_df.index.name = 'Pfam_code'
    result_df = result_df.reset_index()
    result_df = result_df.sort_values('Pfam_code')
    result_df = result_df[['Pfam_code', 'patho_count', 'non_patho_count', 
                          'patho_genes', 'non_patho_genes']]

    # Save to Excel
    result_df.to_excel(OUTPUT_FILE, index=False)

    # Print summary
    print(f"\nExcel file '{OUTPUT_FILE}' has been created with {len(result_df)} common Pfam codes.")
    print("Columns in the output:")
    print("- Pfam_code: The PFAM identifier")
    print("- patho_count: Number of occurrences in pathogenic proteins")
    print("- non_patho_count: Number of occurrences in non-pathogenic proteins")
    print("- patho_genes: Semicolon-separated list of pathogenic genes containing this Pfam")
    print("- non_patho_genes: Semicolon-separated list of non-pathogenic genes containing this Pfam")

if __name__ == "__main__":
    main()