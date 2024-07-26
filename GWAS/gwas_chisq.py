import pandas as pd
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import fdrcorrection

# Load the dataframe
df = pd.read_csv('~/Downloads/filtered_gt_counts.csv')

# Replace '["AT","GT"]' with '["A","G"]' and '["CG","TG"]' with '["C","T"]' in 'alleles_y' column
df['alleles_y'] = df['alleles_y'].replace({'["AT","GT"]': '["A","G"]', '["CG","TG"]': '["C","T"]'})

# Filter the dataframe to only those rows where alleles_x equals alleles_y
filtered_df = df[df['alleles_x'] == df['alleles_y']]

# Columns to be used for the chi-square test
columns_ced = ['wild_type_count_ced', 'heterozygous_count_ced', 'homozygous_count_ced']
columns_control = ['wild_type_count_control', 'heterozygous_count_control', 'homozygous_count_control']

#potential update for minor allele flip
filtered_df['flip wild_type_count_ced']=filtered_df['homozygous_count_ced']
filtered_df['flip homozygous_count_ced']=filtered_df['wild_type_count_ced']
#allele freq

filtered_df['ced_maf']=(2*filtered_df['homozygous_count_ced']+filtered_df['heterozygous_count_ced'])/(2*filtered_df[columns_ced].sum(axis=1))
filtered_df['control_maf']=(2*filtered_df['homozygous_count_control']+filtered_df['heterozygous_count_control'])/(2*filtered_df[columns_control].sum(axis=1))

#now flip
filtered_df.loc[filtered_df['AF_alt']>0.5,'homozygous_count_ced']=filtered_df.loc[filtered_df['AF_alt']>0.5,'flip wild_type_count_ced']
filtered_df.loc[filtered_df['AF_alt']>0.5,'wild_type_count_ced']=filtered_df.loc[filtered_df['AF_alt']>0.5,'flip homozygous_count_ced']

#delta calc
filtered_df.loc[filtered_df['AF_alt']>0.5,'ced_maf']=1-filtered_df.loc[filtered_df['AF_alt']>0.5,'ced_maf']
filtered_df.loc[filtered_df['AF_alt']>0.5,'control_maf']=1-filtered_df.loc[filtered_df['AF_alt']>0.5,'control_maf']
filtered_df['delta']=filtered_df['ced_maf']-filtered_df['control_maf']
filtered_df.loc[filtered_df['AF_alt']>0.5,'AF_alt']=1-filtered_df.loc[filtered_df['AF_alt']>0.5,'AF_alt']
# Function to perform chi-square test and return chi2 statistic and p-value, handling failure
def chi_square_test_with_stat(row):
    try:
        table = [
            [row[columns_ced[0]], row[columns_ced[1]], row[columns_ced[2]]],
            [row[columns_control[0]], row[columns_control[1]], row[columns_control[2]]]
        ]
        chi2, p, dof, ex = chi2_contingency(table)
        return chi2, p
    except:
        return None, None

# Apply the chi-square test to each row and get both chi2 statistic and p-value
filtered_df[['chi_square_statistic', 'chi_square_p_value']] = filtered_df.apply(
    lambda row: pd.Series(chi_square_test_with_stat(row)), axis=1
)

# Perform FDR correction
p_values = filtered_df['chi_square_p_value'].dropna()  # Remove None (NaN) values
_, p_values_corrected = fdrcorrection(p_values)

# Assign the corrected p-values back to the dataframe
filtered_df.loc[filtered_df['chi_square_p_value'].notna(), 'chi_square_fdr_corrected'] = p_values_corrected
filtered_df['Variant ID']=filtered_df['locus'].str.replace(':','-')+'-'+filtered_df['ref']+'-'+filtered_df['alt']
filtered_df.loc[filtered_df['AF_alt']>0.5,'Variant ID']=filtered_df.loc[filtered_df['AF_alt']>0.5,'locus'].str.replace(':','-')+'-'+filtered_df.loc[filtered_df['AF_alt']>0.5,'alt']+'-'+filtered_df.loc[filtered_df['AF_alt']>0.5,'ref']
filtered_df=filtered_df[~pd.isna(filtered_df['DISEASE/TRAIT'])]

filtered_df.to_csv('Filtered_Chi-Square_Test_with_Statistics_flip.csv',index=False)
