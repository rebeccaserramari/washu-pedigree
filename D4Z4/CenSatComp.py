# load packages
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import argparse 
import numpy as np

# parse inputs 
parser = argparse.ArgumentParser()
parser.add_argument('--bed', '-i', type=str, action='store', help='bed file goes here')
parser.add_argument('--outName', '-o', type=str, action='store', help='output file name prefix')
parser.add_argument('--leaf_order', '-l', type=str, action='order of chromosome ids to plot')
parser.add_argument('--pLAM', '-p', type=str, action='bed file with pLAM locations')

args = parser.parse_args()

bed = args.bed
outName = args.outName
leaf_list = args.leaf_order
pLAM = args.pLAM


# Load the BED file with the annotations
df = pd.read_csv(bed, sep="\t", header=None)

# load the bed file with the plam 
pLAMdf = pd.read_csv(pLAM, sep="\t", header=None)
pLAMdf.columns = ['chrom', 'start', 'end']

# Extract relevant columns
df.columns = ["chrom", "start", "end", "name", "score", "strand", "start2", "end2", "color"]

# calculate the length
df['length'] = df['end']-df['start']

# assign colors based on the size of the repeat Annotation
log_norm = (np.log1p(df["length"]) - np.log1p(df["length"]).min()) / \
           (np.log1p(df["length"]).max() - np.log1p(df["length"]).min())

cmap = plt.cm.viridis 
df['color'] = log_norm.map(lambda x: cmap(x))

# Reverse the order of chromosomes to plot in correct order 
chrom_labels = df['chrom'].unique()[::-1]  



#use an existing list to sort the order 
leaves = pd.read_csv(leaf_list, sep="\t", header=None)
leaves.columns = ['Samplehaplotype']

plot_order = leaves['Samplehaplotype'].tolist()[::-1]

# create a dictionary for the y positions 
y_pos_dict = {chrom: i for i, chrom in enumerate(plot_order)}  


# Normalize the start and end positions for each chromosome so they start at 0
normalizer = df.groupby('chrom')['start'].min().to_dict()
print(normalizer)


for dataframe in [df, pLAMdf]:
    dataframe[['normalized_start', 'normalized_end']] = dataframe.apply(
    lambda row: [
        row['start'] - normalizer[row['chrom']],
        row['end'] - normalizer[row['chrom']],
    ], axis=1, result_type='expand')


# Determine the maximum normalized end for x-axis limits
max_normalized_end = df['normalized_end'].max()


# Set x-axis limits with 500,000 base pairs padding
x_min = 0
x_max = max_normalized_end


# Plotting the annotations
fig, ax = plt.subplots(figsize=(20, len(df['chrom'].unique()) * 0.5))

for chrom in chrom_labels:
    short_chrom = "#".join(chrom.split('#')[:2])
    if short_chrom in y_pos_dict:
        chrom_data = df[df['chrom'] == chrom]
        plam_data = pLAMdf[pLAMdf['chrom'] == chrom]
        y_pos = y_pos_dict[short_chrom]  # Get the correct y position from the dictionary
        for _, row in chrom_data.iterrows():
            start = row['normalized_start']
            end = row['normalized_end']
            color = row['color']
            ax.barh(y_pos, end - start, left=start, color=color, height=0.5, edgecolor='black')
        for _, row in plam_data.iterrows():
            start = row['normalized_start']
            end = row['normalized_end']
            ax.barh(y_pos-0.08, end - start, left=start, color='blue', alpha=0.3, height=1, edgecolor='blue')
        
        
    else:
        continue

# Set the y-ticks with chromosome names and increase font size
ax.set_yticks(range(len(plot_order)))
ax.set_yticklabels(plot_order, fontsize=14)

# Adjust the y-axis limits to reduce space at the top and bottom
ax.set_ylim(-0.5, len(plot_order) - 0.5)

# Set plot title with larger font size
ax.set_title(outName, fontsize=18)


# Remove x-ticks and x-tick labels
ax.xaxis.set_ticks([])
ax.xaxis.set_ticklabels([])

# Save the plot without the key
plt.savefig(outName+".jpg", dpi=100, bbox_inches='tight')