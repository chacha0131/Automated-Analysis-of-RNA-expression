import os
import numpy as np
import pandas as pd

# Set up experimental variables
path = "/Users/emilycha/Desktop/imgs_bagot/redo/PVT/out/quant/"
files = os.listdir(path)

# Create an empty list to store the data for each file
data = []

# Loop through each file and append the data to the list
for file in files:
    print(file)
    df = pd.read_csv(path + file + '/cells_quant.csv')
    df['file'] = file  # Add a new column for the file name
    data.append(df)

# Concatenate the data for all files into a single DataFrame
all_data = pd.concat(data)

# Group the data by file and type, and count the cells in each group
counts = all_data.groupby(['file']).apply(lambda x: pd.Series({
    'sdk1_puncta': x['sdk1'].sum(),  # Add a new column for the sum of sdk1
    **x['type'].value_counts(normalize=False).to_dict()
})).unstack(fill_value=0)

# Reset the index to remove the multi-level column index
counts = counts.reset_index()

# Save the counts to a CSV file
counts.to_csv(path + 'cell_counts_all_files.csv', index=False)
