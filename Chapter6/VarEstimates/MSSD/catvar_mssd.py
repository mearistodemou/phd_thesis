# Calculate MSSD per voxel and region
# Author: Michael Aristodemou

############ Appendix #############
#0. import packages & functions

#1. Voxel-wise analyses
#1.1. create dataframe for drug1
#1.2. create dataframe for drug2
#1.3. create dataframe for drug3

#2. Regional analyses (cortical)
#2.1. create dataframe for drug1
#2.2. create dataframe for drug2
#2.3. create dataframe for drug3

#3. Regional analyses (subcortical)
#3.1. create dataframe for drug1
#3.2. create dataframe for drug2
#3.3. create dataframe for drug3 
#####################################

#################################
#0.1 import packages
#################################

import os
import nilearn
import nilearn.plotting
import pandas as pd
import numpy as np
from nilearn import input_data, image, datasets, plotting, masking
import matplotlib.pyplot as plt
import nibabel as nib
# from sklearn.decomposition import PCA

######################################################
#0.2 Functions
######################################################

# Define the function to compute MSSD per voxel
def mean_successive_squared_differences(time_series):
    differences = np.diff(time_series, axis=0)
    squared_differences = np.square(differences)
    mean_squared_differences = np.mean(squared_differences, axis=0)
    return mean_squared_differences

###########################################################################
#1. Estimate voxel-wise MSSD
###########################################################################

########################################################
#1.1 Create dataframe for drug1
#########################################################

# List of participant IDs (example for 148 participants)
participant_ids = [f'sub-{i:03d}' for i in range(1, 101)]

# Base directory path
base_dir = "C:/Users/micha/OneDrive/Documents/CatVar/Data/NII"

# Create an empty dataframe to store the results
all_mssd_df = pd.DataFrame()

# Loop over each participant
for participant in participant_ids:
    # Construct the file path
    file_name = f"{participant}_ses-drug1_task-MID_run-4_bold.nii"
    fmri_file = os.path.join(base_dir, file_name, file_name, file_name)
    
    # Check if the file exists before loading
    if os.path.exists(fmri_file):
        # Load fMRI data
        fmri_img = image.load_img(fmri_file)
        
        # Extract timeseries from each voxel
        masker = input_data.NiftiMasker(standardize=True, detrend=True, smoothing_fwhm=6)
        time_series = masker.fit_transform(fmri_img)

        # Estimate MSSD
        mssd = mean_successive_squared_differences(time_series)
        
        # Create a dataframe to store the MSSD for the current participant
        mssd_df = pd.DataFrame(mssd, columns=["MSSD"])
        mssd_df["Participant"] = participant
        
        # Append the current participant's MSSD dataframe to the main dataframe
        all_mssd_df = pd.concat([all_mssd_df, mssd_df], ignore_index=True)
    else:
        print(f"File for participant {participant} does not exist at {fmri_file}")

# Save the combined MSSD DataFrame to a CSV file
all_mssd_df.to_csv("C:/Users/micha/OneDrive/Documents/CatVar/Data/MSSD/mssd_ses-drug1.csv", index_label="Voxel")

########################################################
#1.2. Create dataframe for drug2
#########################################################

# List of participant IDs (example for 148 participants)
participant_ids = [f'sub-{i:03d}' for i in range(1, 101)]

# Base directory path
base_dir = "C:/Users/micha/OneDrive/Documents/CatVar/Data/NII"

# Create an empty dataframe to store the results
all_mssd_df = pd.DataFrame()

# Loop over each participant
for participant in participant_ids:
    # Construct the file path
    file_name = f"{participant}_ses-drug2_task-MID_run-4_bold.nii"
    fmri_file = os.path.join(base_dir, file_name, file_name, file_name)
    
    # Check if the file exists before loading
    if os.path.exists(fmri_file):
        # Load fMRI data
        fmri_img = image.load_img(fmri_file)
        
        # Extract timeseries from each voxel
        masker = input_data.NiftiMasker(standardize=True, detrend=True, smoothing_fwhm=6)
        time_series = masker.fit_transform(fmri_img)
        
        # Estimate MSSD
        mssd = mean_successive_squared_differences(time_series)
        
        # Create a dataframe to store the MSSD for the current participant
        mssd_df = pd.DataFrame(mssd, columns=["MSSD"])
        mssd_df["Participant"] = participant
        
        # Append the current participant's MSSD dataframe to the main dataframe
        all_mssd_df = pd.concat([all_mssd_df, mssd_df], ignore_index=True)
    else:
        print(f"File for participant {participant} does not exist at {fmri_file}")

# Save the combined MSSD DataFrame to a CSV file
all_mssd_df.to_csv("C:/Users/micha/OneDrive/Documents/CatVar/Data/MSSD/mssd_ses-drug2.csv", index_label="Voxel")

########################################################
#1.3. Create dataframe for drug3
#########################################################

# List of participant IDs (example for 148 participants)
participant_ids = [f'sub-{i:03d}' for i in range(1, 101)]

# Base directory path
base_dir = "C:/Users/micha/OneDrive/Documents/CatVar/Data/NII"

# Create an empty dataframe to store the results
all_mssd_df = pd.DataFrame()

# Loop over each participant
for participant in participant_ids:
    # Construct the file path
    file_name = f"{participant}_ses-drug3_task-MID_run-4_bold.nii"
    fmri_file = os.path.join(base_dir, file_name, file_name, file_name)
    
    # Check if the file exists before loading
    if os.path.exists(fmri_file):
        # Load fMRI data
        fmri_img = image.load_img(fmri_file)
        
        # Extract timeseries from each voxel
        masker = input_data.NiftiMasker(standardize=True, detrend=True, smoothing_fwhm=6)
        time_series = masker.fit_transform(fmri_img)
        
        # Estimate MSSD
        mssd = mean_successive_squared_differences(time_series)
        
        # Create a dataframe to store the MSSD for the current participant
        mssd_df = pd.DataFrame(mssd, columns=["MSSD"])
        mssd_df["Participant"] = participant
        
        # Append the current participant's MSSD dataframe to the main dataframe
        all_mssd_df = pd.concat([all_mssd_df, mssd_df], ignore_index=True)
    else:
        print(f"File for participant {participant} does not exist at {fmri_file}")

# Save the combined MSSD DataFrame to a CSV file
all_mssd_df.to_csv("C:/Users/micha/OneDrive/Documents/CatVar/Data/MSSD/mssd_ses-drug3.csv", index_label="Voxel")

#####################################################################################
#2. Estimate regional MSSD (cortical)
#####################################################################################

# Load Harvard-Oxford atlas
# dataset = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr25-2mm') # subcortical regions
dataset = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm') # cortical regions
atlas_filename = dataset.maps
labels = dataset.labels

# Create a figure with specified size and DPI
plt.figure(figsize=(10, 10), dpi=300)  # You can adjust figsize as needed
# Plot atlas (48 regions including background)
plotting.plot_roi(atlas_filename, title='Harvard-Oxford Atlas Cortical Regions', draw_cross=False)
# Show the plot
plt.savefig('C:/Users/micha/OneDrive/Documents/CatVar/Visualization/HO_atlas.png', bbox_inches='tight', dpi=300)
plotting.show()

# List of participant IDs (example for 148 participants)
participant_ids = [f'sub-{i:03d}' for i in range(1, 101)]

# Base directory path
base_dir = "C:/Users/micha/OneDrive/Documents/CatVar/Data/NII"


######################################################################
#2.1. Estimate regional MSSD for session 1
######################################################################

# Create an empty dataframe to store the results
all_mssd_df = pd.DataFrame()

# Loop over each participant
for participant in participant_ids:
    # Construct the file path
    file_name = f"{participant}_ses-drug1_task-MID_run-4_bold.nii"
    fmri_file = os.path.join(base_dir, file_name, file_name, file_name)
    
    # Check if the file exists before loading
    if os.path.exists(fmri_file):
        # Load fMRI data
        fmri_img = image.load_img(fmri_file)
        
        # Extract timeseries from each voxel
        masker = input_data.NiftiLabelsMasker(labels_img=atlas_filename, labels= labels, standardize=True, detrend=True, smoothing_fwhm=6)
        time_series = masker.fit_transform(fmri_img)

        # Estimate MSSD
        mssd = mean_successive_squared_differences(time_series)
        
        # Create a dataframe to store the MSSD for the current participant
        mssd_df = pd.DataFrame(mssd, columns=["MSSD"])
        mssd_df["Participant"] = participant
        mssd_df["Region"] = masker.region_names_
        
        # Append the current participant's MSSD dataframe to the main dataframe
        all_mssd_df = pd.concat([all_mssd_df, mssd_df], ignore_index=True)
    else:
        print(f"File for participant {participant} does not exist at {fmri_file}")

# Save the combined MSSD DataFrame to a CSV file
all_mssd_df.to_csv("C:/Users/micha/OneDrive/Documents/CatVar/Data/MSSD/mssd_ses-drug1_regions.csv", index_label="Voxel")









######################################################################
#2.2. Estimate regional MSSD for session 2
######################################################################

# Create an empty dataframe to store the results
all_mssd_df = pd.DataFrame()

# Loop over each participant
for participant in participant_ids:
    # Construct the file path
    file_name = f"{participant}_ses-drug2_task-MID_run-4_bold.nii"
    fmri_file = os.path.join(base_dir, file_name, file_name, file_name)
    
    # Check if the file exists before loading
    if os.path.exists(fmri_file):
        # Load fMRI data
        fmri_img = image.load_img(fmri_file)
        
        # Extract timeseries from each voxel
        masker = input_data.NiftiLabelsMasker(labels_img=atlas_filename, labels= labels, standardize=True, detrend=True, smoothing_fwhm=6)
        time_series = masker.fit_transform(fmri_img)

        # Estimate MSSD
        mssd = mean_successive_squared_differences(time_series)
        
        # Create a dataframe to store the MSSD for the current participant
        mssd_df = pd.DataFrame(mssd, columns=["MSSD"])
        mssd_df["Participant"] = participant
        mssd_df["Region"] = masker.region_names_
        
        # Append the current participant's MSSD dataframe to the main dataframe
        all_mssd_df = pd.concat([all_mssd_df, mssd_df], ignore_index=True)
    else:
        print(f"File for participant {participant} does not exist at {fmri_file}")

# Save the combined MSSD DataFrame to a CSV file
all_mssd_df.to_csv("C:/Users/micha/OneDrive/Documents/CatVar/Data/MSSD/mssd_ses-drug2_regions.csv", index_label="Voxel")

######################################################################
#2.3. Estimate regional MSSD for session 3
######################################################################

# Create an empty dataframe to store the results
all_mssd_df = pd.DataFrame()

# Loop over each participant
for participant in participant_ids:
    # Construct the file path
    file_name = f"{participant}_ses-drug3_task-MID_run-4_bold.nii"
    fmri_file = os.path.join(base_dir, file_name, file_name, file_name)
    
    # Check if the file exists before loading
    if os.path.exists(fmri_file):
        # Load fMRI data
        fmri_img = image.load_img(fmri_file)
        
        # Extract timeseries from each voxel
        masker = input_data.NiftiLabelsMasker(labels_img=atlas_filename, labels= labels, standardize=True, detrend=True, smoothing_fwhm=6)
        time_series = masker.fit_transform(fmri_img)

        # Estimate MSSD
        mssd = mean_successive_squared_differences(time_series)
        
        # Create a dataframe to store the MSSD for the current participant
        mssd_df = pd.DataFrame(mssd, columns=["MSSD"])
        mssd_df["Participant"] = participant
        mssd_df["Region"] = masker.region_names_
        
        # Append the current participant's MSSD dataframe to the main dataframe
        all_mssd_df = pd.concat([all_mssd_df, mssd_df], ignore_index=True)
    else:
        print(f"File for participant {participant} does not exist at {fmri_file}")

# Save the combined MSSD DataFrame to a CSV file
all_mssd_df.to_csv("C:/Users/micha/OneDrive/Documents/CatVar/Data/MSSD/mssd_ses-drug3_regions.csv", index_label="Voxel")

#####################################################################################
#3. Estimate regional MSSD (subcortical)
#####################################################################################

# Load Harvard-Oxford atlas
dataset = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr25-2mm') # subcortical regions
# dataset = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm') # cortical regions
atlas_filename = dataset.maps
labels = dataset.labels

# Plot atlas (48 regions inc background)
plotting.plot_roi(atlas_filename)
plotting.show()

# List of participant IDs (example for 148 participants)
participant_ids = [f'sub-{i:03d}' for i in range(1, 101)]

# Base directory path
base_dir = "C:/Users/micha/OneDrive/Documents/CatVar/Data/NII"


######################################################################
#3.1. Estimate regional MSSD for session 1
######################################################################

# Create an empty dataframe to store the results
all_mssd_df = pd.DataFrame()

# Loop over each participant
for participant in participant_ids:
    # Construct the file path
    file_name = f"{participant}_ses-drug1_task-MID_run-4_bold.nii"
    fmri_file = os.path.join(base_dir, file_name, file_name, file_name)
    
    # Check if the file exists before loading
    if os.path.exists(fmri_file):
        # Load fMRI data
        fmri_img = image.load_img(fmri_file)
        
        # Extract timeseries from each voxel
        masker = input_data.NiftiLabelsMasker(labels_img=atlas_filename, labels= labels, standardize=True, detrend=True, smoothing_fwhm=6)
        time_series = masker.fit_transform(fmri_img)

        # Estimate MSSD
        mssd = mean_successive_squared_differences(time_series)
        
        # Create a dataframe to store the MSSD for the current participant
        mssd_df = pd.DataFrame(mssd, columns=["MSSD"])
        mssd_df["Participant"] = participant
        mssd_df["Region"] = masker.region_names_

        # Remove any duplicates in mssd_df (within the same participant)
        mssd_df.drop_duplicates(subset=["Region"], inplace=True)
        
        # Append the current participant's MSSD dataframe to the main dataframe
        all_mssd_df = pd.concat([all_mssd_df, mssd_df], ignore_index=True)
    else:
        print(f"File for participant {participant} does not exist at {fmri_file}")

# Save the combined MSSD DataFrame to a CSV file
all_mssd_df.to_csv("C:/Users/micha/OneDrive/Documents/CatVar/Data/MSSD/mssd_ses-drug1_subreg.csv", index_label="Voxel")


######################################################################
#3.2. Estimate regional MSSD for session 2
######################################################################

# Create an empty dataframe to store the results
all_mssd_df = pd.DataFrame()

# Loop over each participant
for participant in participant_ids:
    # Construct the file path
    file_name = f"{participant}_ses-drug2_task-MID_run-4_bold.nii"
    fmri_file = os.path.join(base_dir, file_name, file_name, file_name)
    
    # Check if the file exists before loading
    if os.path.exists(fmri_file):
        # Load fMRI data
        fmri_img = image.load_img(fmri_file)
        
        # Extract timeseries from each voxel
        masker = input_data.NiftiLabelsMasker(labels_img=atlas_filename, labels= labels, standardize=True, detrend=True, smoothing_fwhm=6)
        time_series = masker.fit_transform(fmri_img)

        # Estimate MSSD
        mssd = mean_successive_squared_differences(time_series)
        
        # Create a dataframe to store the MSSD for the current participant
        mssd_df = pd.DataFrame(mssd, columns=["MSSD"])
        mssd_df["Participant"] = participant
        mssd_df["Region"] = masker.region_names_

        # Remove any duplicates in mssd_df (within the same participant)
        mssd_df.drop_duplicates(subset=["Region"], inplace=True)
        
        # Append the current participant's MSSD dataframe to the main dataframe
        all_mssd_df = pd.concat([all_mssd_df, mssd_df], ignore_index=True)
    else:
        print(f"File for participant {participant} does not exist at {fmri_file}")

# Save the combined MSSD DataFrame to a CSV file
all_mssd_df.to_csv("C:/Users/micha/OneDrive/Documents/CatVar/Data/MSSD/mssd_ses-drug2_subreg.csv", index_label="Voxel")

######################################################################
#3.3. Estimate regional MSSD for session 3
######################################################################

# Create an empty dataframe to store the results
all_mssd_df = pd.DataFrame()

# Loop over each participant
for participant in participant_ids:
    # Construct the file path
    file_name = f"{participant}_ses-drug3_task-MID_run-4_bold.nii"
    fmri_file = os.path.join(base_dir, file_name, file_name, file_name)
    
    # Check if the file exists before loading
    if os.path.exists(fmri_file):
        # Load fMRI data
        fmri_img = image.load_img(fmri_file)
        
        # Extract timeseries from each voxel
        masker = input_data.NiftiLabelsMasker(labels_img=atlas_filename, labels= labels, standardize=True, detrend=True, smoothing_fwhm=6)
        time_series = masker.fit_transform(fmri_img)

        # Estimate MSSD
        mssd = mean_successive_squared_differences(time_series)
        
        # Create a dataframe to store the MSSD for the current participant
        mssd_df = pd.DataFrame(mssd, columns=["MSSD"])
        mssd_df["Participant"] = participant
        mssd_df["Region"] = masker.region_names_
        
        # Remove any duplicates in mssd_df (within the same participant)
        mssd_df.drop_duplicates(subset=["Region"], inplace=True)
        
        # Append the current participant's MSSD dataframe to the main dataframe
        all_mssd_df = pd.concat([all_mssd_df, mssd_df], ignore_index=True)
    else:
        print(f"File for participant {participant} does not exist at {fmri_file}")

# Save the combined MSSD DataFrame to a CSV file
all_mssd_df.to_csv("C:/Users/micha/OneDrive/Documents/CatVar/Data/MSSD/mssd_ses-drug3_subreg.csv", index_label="Voxel")

############################################################
# plotting mssd using atlas
############################################################

import matplotlib.pyplot as plt
from nilearn import plotting
from nilearn import image

######################################################################
# Step 4: Map MSSD to the Brain Atlas and Visualize
######################################################################

# Now we need to average the MSSD per region across participants.
mssd_per_region = all_mssd_df.groupby("Region")["MSSD"].mean()

from nilearn import image, plotting
import numpy as np
import matplotlib.pyplot as plt

# Step 1: Map the MSSD values to the brain regions in the atlas
# Load the atlas as a Nifti image
atlas_img = image.load_img(atlas_filename)

# Step 2: Extract the atlas data (3D array with integer labels for each region)
atlas_data = atlas_img.get_fdata()

# Step 3: Create an empty array to hold MSSD values for each voxel
# Initialize the array with the same shape as the atlas
mssd_values_for_voxels = np.zeros_like(atlas_data)

# Map MSSD values to brain regions in the atlas
for idx, region_name in enumerate(masker.region_names_):
    if region_name in mssd_per_region.index:
        region_label = masker.labels_[idx]
        mssd_value = mssd_per_region[region_name]

        print(f"Processing region: {region_name} with label: {region_label} and MSSD value: {mssd_value}")

        # Assign the MSSD value to all voxels that correspond to the current region label
        # Ensure that region_label matches the values in atlas_data
        mssd_values_for_voxels[atlas_data == region_label] = mssd_value

        # Debug output to check how many voxels were assigned
        num_assigned = np.sum(atlas_data == region_label)
        print(f"Number of voxels assigned for {region_name}: {num_assigned}")

# Check if any values were assigned
print("Values in mssd_values_for_voxels after assignment:")
print(np.unique(mssd_values_for_voxels))  # Should show non-zero values if the assignment was successful

# Create a new Nifti image with MSSD values
mssd_img = image.new_img_like(atlas_img, mssd_values_for_voxels)

# Step 6: Plot the resulting MSSD image using nilearn's plot_stat_map
plotting.plot_stat_map(
    mssd_img,
    title="Mean MSSD per Region",
    display_mode="ortho",  # You can use 'x', 'y', 'z', or 'ortho' for different views
    draw_cross=False,
    cmap="coolwarm",  # You can choose a colormap that suits the data range
    threshold=0.01,  # Set a threshold to remove background or irrelevant values
    colorbar=True
)

# Show the plot
plt.show()