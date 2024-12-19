import pandas as pd
import numpy as np
df=pd.read_csv('ancestry_count.csv')

#sort by OR of eur against others
df['alt_ct']=df[['afr', 'amr'  ]].sum(axis=1)
df['eur other_OR']=df['eur']*(df.loc[44,'alt_ct']-df['alt_ct'])/(df['alt_ct']*(df.loc[44,'eur']-df['eur']))
df=df[((df['pvalue']<=0.05)|(df['amr pvalue']<=0.05))&(df['Unnamed: 0']!='genetic_female')]
print(df.sort_values(['amr OR','afr OR'])['amr OR'])
condition_order=df.sort_values(['amr OR','afr OR'])['Unnamed: 0']#.drop(44)
# Load the cleaned dataframe
df_cleaned = df[['Unnamed: 0', 'afr OR', 'amr OR',  'afr pvalue', 'amr pvalue']].copy()
df_cleaned.columns = ['condition', 'afr OR', 'amr OR', 'afr pvalue', 'amr pvalue']
df_cleaned['eur pvalue'] = 1  # Set EUR p-value to 1
df_cleaned['eur OR'] = 1  # Set EUR OR to 1
# Melt the dataframe to have ancestry in one column
df_melted_no_female = df_cleaned.melt(id_vars=['condition', 'afr pvalue', 'amr pvalue', 'eur pvalue'], 
                                      var_name='ancestry', value_name='OR')

# Map the correct p-values based on ancestry
pvalue_map = {
    'afr OR': 'afr pvalue',
    'amr OR': 'amr pvalue',
    'eur OR': 'eur pvalue'
}

# Apply the mapping to add p-values
df_melted_no_female['pvalue'] = df_melted_no_female.apply(lambda row: row[pvalue_map[row['ancestry']]], axis=1)

# Convert p-values to -log10(p-values)
df_melted_no_female['log_pvalue'] = -np.log10(df_melted_no_female['pvalue'])
df_melted_no_female['ancestry'].replace(['eur OR','amr OR','afr OR'  ],['eur (reference)','amr','afr'  ],inplace=True)
# Custom order for x-axis
custom_order = ['eur (reference)','amr','afr'  ]
# Use .loc to set the custom order for 'Impact label'
df_melted_no_female = df_melted_no_female.set_index('ancestry').loc[custom_order].reset_index()
df_melted_no_female = df_melted_no_female.set_index('condition').loc[condition_order].reset_index()
#remove reference group from plot
df_melted_no_female = df_melted_no_female[df_melted_no_female['ancestry']!='eur (reference)']

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
def shiftedColorMap(cmap, start=0.5, midpoint=1, stop=2.5, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

shifted_cmap = shiftedColorMap(matplotlib.cm.seismic,start=0, midpoint=1, stop=1.84, name='shifted')
# Set the font type to 42 (TrueType)
plt.rcParams['pdf.fonttype'] = 42

# Re-filter the data including 'eur' percentages and conditions where p-values < 0.05
df_filtered_including_eur = df_melted_no_female[
    ~df_melted_no_female['condition'].str.contains('total', case=False, na=False)
]
df_filtered_including_eur.loc[df_filtered_including_eur['log_pvalue']<-np.log10(0.05),'OR']=1
# Recreate the plot including the 'eur' percentages and save as PDF
plt.figure(figsize=(4.5, 6.5))  # Maintain the same figure size for PDF

# Plot the filtered data including 'eur'
plt.scatter(
    x=df_filtered_including_eur['ancestry'],
    y=df_filtered_including_eur['condition'],
    s=df_filtered_including_eur['log_pvalue'] * 70+10,  # Reduce bubble size to 90% of the current size
    c=df_filtered_including_eur['OR'],           # Color based on -log10(pvalue)
    cmap='shifted',                                         # Use 'Reds' colormap for white to red scale
    alpha=1.0,                                           # Make circles fully opaque
    edgecolor='black',                                   # Black border for better contrast
    linewidth=0.5,                                         # Thicker black borders for circles
    marker='o'                                           # Use 'o' marker for filled circles
)

# Add a color bar to represent the -log10(p-value)
cbar = plt.colorbar()
cbar.set_label('Odds ratio', rotation=270, labelpad=15)

# Add axis labels and title
plt.xlabel('Ancestry')
plt.ylabel('Condition')
plt.title('Clinical Conditons differing by ancestry')

# Rotate x-axis labels for readability
plt.xticks(rotation=45)

# Add the corrected legend for circle sizes
size_labels=[0.05,0.01,0.001]
sizes = -np.log10(size_labels)  # Example percentages for the legend
size_legend = [Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=np.sqrt(size * 70), 
                      markeredgecolor='black', markeredgewidth=1) for size in sizes]
plt.legend(size_legend, [f'{size}' for size in size_labels], scatterpoints=1, frameon=True, labelspacing=0.5, 
           title='pvalue', bbox_to_anchor=(1.5, 0.5), loc='center left', borderpad=1, handletextpad=1)

# Save the plot to a PDF file with fonttype 42
pdf_file_path = 'bubble_plot_ancestry_condition_ORs_noref_1104.pdf'
plt.tight_layout()
plt.savefig(pdf_file_path,bbox_inches='tight', format='pdf')
#plt.show()
# Provide the download link
pdf_file_path

