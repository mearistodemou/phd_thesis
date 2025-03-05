############################################
# Brain Atlas plots for Regional MSSD
############################################

# load packages
import os
import nilearn
import nilearn.plotting
import pandas as pd
import numpy as np
from nilearn import input_data, image, datasets, plotting, masking
import matplotlib.pyplot as plt
import nibabel as nib

######################################################################
# Cortical Atlas
######################################################################

# load cortical brain atlas
# Load Harvard-Oxford atlas
# dataset = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr25-2mm') # subcortical regions
dataset = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm') # cortical regions
atlas_filename = dataset.maps
labels = dataset.labels

# load dataframe
# effect size
mssd_per_region = pd.read_csv(r'C:\Users\micha\OneDrive\Documents\CatVar\Visualization\MPH_delta_cort.csv') # MPH
mssd_per_region = pd.read_csv(r'C:\Users\micha\OneDrive\Documents\CatVar\Visualization\SUL_delta_cort.csv') # SUL
# pvalue
mssd_per_region = pd.read_csv(r'C:\Users\micha\OneDrive\Documents\CatVar\Visualization\MPH_pval_cort.csv') # MPH
mssd_per_region = pd.read_csv(r'C:\Users\micha\OneDrive\Documents\CatVar\Visualization\SUL_pval_cort.csv') # SUL

# Convert the DataFrame into a Series
mssd_per_region = mssd_per_region.set_index('Region')["Value"]
# Ensure the data type is float32 to match the original Series
mssd_per_region = mssd_per_region.astype('float32')
# Check the structure of the converted Series
mssd_per_region.info()

# Remove 'Background' from the list of labels
labels = [label for label in labels if label != 'Background']

# Display the updated labels of regions
print("Updated labels:", labels)

# Now you can calculate the number of regions excluding the background
num_regions = len(labels)  # Number of regions in the atlas
mssd_values = mssd_per_region  # Replace this with actual MSSD values

# Step 2: Create a new image based on the MSSD values
atlas_data = image.load_img(atlas_filename).get_fdata()
mssd_values_for_voxels = np.zeros(atlas_data.shape)

# Check the length of mssd_values to avoid indexing errors
if len(mssd_values) != num_regions:
    raise ValueError(f"Number of MSSD values ({len(mssd_values)}) does not match number of regions ({num_regions}).")

# Map MSSD values to the corresponding voxels in the atlas
for idx, region_name in enumerate(labels):
    # Find the label in the atlas that corresponds to the current region
    region_label = idx + 1  # Adjust if your region labels start from 1
    # Check if region_label exists in atlas_data before indexing
    if region_label in np.unique(atlas_data):
        mssd_values_for_voxels[atlas_data == region_label] = mssd_values[idx]
# space

# Step 3: Create a Nifti image from the MSSD values
mssd_img = image.new_img_like(atlas_filename, mssd_values_for_voxels)

# Create a figure with specified size and DPI
plt.figure(figsize=(10, 10), dpi=300)  # You can adjust figsize as needed

# Step 4: Plot the atlas with the MSSD heatmaplabels
plotting.plot_stat_map(
    mssd_img,
    title='Significance of MSSD change under SUL in cortical regions',
    display_mode='ortho',  # Use 'ortho' for 3D views
    draw_cross=False,
    cmap='brown_blue',  # Change to a suitable colormap
    threshold=0,  # Adjust as needed
    colorbar=True
)

# Show the plot
plt.savefig('C:/Users/micha/OneDrive/Documents/CatVar/Visualization/HO_pval_cort_SUL.png', bbox_inches='tight', dpi=300)
# Show the plot
plotting.show()


######################################################################
# Subcortical Atlas
######################################################################

# load cortical brain atlas
# Load Harvard-Oxford atlas
dataset = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr25-2mm') # subcortical regions
# dataset = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm') # cortical regions
atlas_filename = dataset.maps
labels = dataset.labels

# load dataframe
# effect size
mssd_per_region = pd.read_csv(r'C:\Users\micha\OneDrive\Documents\CatVar\Visualization\MPH_delta_subcort.csv') # MPH
mssd_per_region = pd.read_csv(r'C:\Users\micha\OneDrive\Documents\CatVar\Visualization\SUL_delta_subcort.csv') # SUL
# pvalue
mssd_per_region = pd.read_csv(r'C:\Users\micha\OneDrive\Documents\CatVar\Visualization\MPH_pval_subcort.csv') # MPH
mssd_per_region = pd.read_csv(r'C:\Users\micha\OneDrive\Documents\CatVar\Visualization\SUL_pval_subcort.csv') # SUL

# Convert the DataFrame into a Series
mssd_per_region = mssd_per_region.set_index('Region')["Value"]
# Ensure the data type is float32 to match the original Series
mssd_per_region = mssd_per_region.astype('float32')
# Check the structure of the converted Series
mssd_per_region.info()

# Remove 'Background' from the list of labels
labels = [label for label in labels if label != 'Background']

# Display the updated labels of regions
print("Updated labels:", labels)

# Now you can calculate the number of regions excluding the background
num_regions = len(labels)  # Number of regions in the atlas
mssd_values = mssd_per_region  # Replace this with actual MSSD values

# Step 2: Create a new image based on the MSSD values
atlas_data = image.load_img(atlas_filename).get_fdata()
mssd_values_for_voxels = np.zeros(atlas_data.shape)

# Check the length of mssd_values to avoid indexing errors
if len(mssd_values) != num_regions:
    raise ValueError(f"Number of MSSD values ({len(mssd_values)}) does not match number of regions ({num_regions}).")

# Map MSSD values to the corresponding voxels in the atlas
for idx, region_name in enumerate(labels):
    # Find the label in the atlas that corresponds to the current region
    region_label = idx + 1  # Adjust if your region labels start from 1
    # Check if region_label exists in atlas_data before indexing
    if region_label in np.unique(atlas_data):
        mssd_values_for_voxels[atlas_data == region_label] = mssd_values[idx]
# space

# Step 3: Create a Nifti image from the MSSD values
mssd_img = image.new_img_like(atlas_filename, mssd_values_for_voxels)

# Create a figure with specified size and DPI
plt.figure(figsize=(10, 10), dpi=300)  # You can adjust figsize as needed

# Step 4: Plot the atlas with the MSSD heatmaplabels
plotting.plot_stat_map(
    mssd_img,
    title='Significance of MSSD change under SUL in subcortical regions',
    display_mode='ortho',  # Use 'ortho' for 3D views
    draw_cross=False,
    cmap='brown_blue',  # Change to a suitable colormap
    threshold=0,  # Adjust as needed
    colorbar=True
)

# Show the plot
plt.savefig('C:/Users/micha/OneDrive/Documents/CatVar/Visualization/HO_pval_subcort_SUL.png', bbox_inches='tight', dpi=300)
# Show the plot
plotting.show()