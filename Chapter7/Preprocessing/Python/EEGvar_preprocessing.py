##########################################################################################################################
##############################################           Appendix           ##############################################
# 1. Cloning data using datalad
# 2. Installing data from participants we want to preprocess
# 3. Visualizing data with BIDS
# 4.  Preprocessing data with MNE BIDS
##########################################################################################################################

##############################################################
# 1. Cloning data using datalad
##############################################################
# Data cloning and extraction for PEERS experiments
# Change working directory to save clone data
# import datalad

import os
os.chdir(r'C:\Users\micha\OneDrive\Documents\EEGvar')
os.getcwd()

# Clone dataset to PEERS folder
# datalad install https://github.com/OpenNeuroDatasets/ds004395.git

# Import specific data files for preprocessing
# ^Z #exit python to enter command line
# mkdir df_eegvar #create folder for local data
# cd ds004395
# cd C:\Users\micha\OneDrive\Documents\EEGvar\bids_chirstoph

################################################################
# 2. Download data from participants we want to preprocess
################################################################

# Get
# datalad get -d^ 'sub-003\ses-1\eeg\sub-003_ses-1_task-multiple_run-01_eeg.eeg'
# Drop
# datalad drop -d^ 'sub-002\ses-1\eeg\sub-002_ses-1_task-multiple_run-01_eeg.eeg'

##################################################################
# 3. Visualizing data with BIDS
###################################################################

# BIDS tutorial
import os
import os.path as op
import openneuro
import mne_bids
# pip install #openneuro-py #run this in powershell


from mne.datasets import sample
from mne_bids import (
    BIDSPath,
    read_raw_bids,
    print_dir_tree,
    make_report,
    find_matching_paths,
    get_entity_vals,
)

# Get to bids_root
dataset = "bids_chirstoph"
subject = "sub-151"
sub = "151"

# Set root directory for mne
bids_root = op.join(op.dirname(r'C:\Users\micha\OneDrive\Documents\EEGvar\l'), dataset)
# print(bids_root)

# openneuro.download(dataset=dataset, target_dir=bids_root, include=[f"sub-{subject}"])


# Explore the contents of the dataset
# print_dir_tree(bids_root, max_depth= 4)

# Prepare dataset for reading
session = "1"
datatype = "eeg"
run = "01"
bids_path = BIDSPath(root=bids_root, session=session, run=run, datatype=datatype)
# print(bids_path.match(ignore_json=True))

# print(bids_path)

# Create path with all the info required to actually go view EEG data
task = "multiple"
suffix = "eeg"
subject = sub

bids_path = bids_path.update(subject=subject, task=task, suffix=suffix)
# print(bids_path)

# Read BIDS data
raw = read_raw_bids(bids_path=bids_path, verbose=False)
print(raw.info["subject_info"])
print(raw.info["line_freq"])
print(raw.info["sfreq"])
print(raw.annotations)

# Plot raw data
# raw.plot()


################################################################################
# 4. Preprocessing data by José C. García Alanis
################################################################################


"""
=============
Preprocessing
==============

Extracts relevant data and removes artefacts.

Authors: José C. García Alanis <alanis.jcg@gmail.com>
Michael E. Aristodemou <michael.aristodemou@radboudumc.nl>

License: BSD (3-clause)
"""
# %%
# imports
import sys
import os

import warnings

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from mne import events_from_annotations, concatenate_raws, Report
from mne.filter import notch_filter
from mne.preprocessing import ICA
from mne.utils import logger

from mne_bids import BIDSPath, read_raw_bids

from pyprep.prep_pipeline import NoisyChannels
from mne_icalabel import label_components

# %%
# create bids path for import
FNAME = BIDSPath(root=bids_root,
                 subject=sub,
                 task='multiple',
                 session="1",
                 run="01",
                 datatype='eeg',
                 extension = '.vhdr')

# if not os.path.exists(FNAME):
#     warnings.warn(MISSING_FPATH_BIDS_MSG.format(FNAME))
#     sys.exit()

# %%
# get the data
raw = read_raw_bids(FNAME)
raw.load_data()

# Set reference to average
raw = raw.set_eeg_reference('average')

#Drop Cz
# raw = raw.drop_channels(ch_names = ['Cz'], on_missing='raise')

# get sampling rate
# sfreq = raw.info['sfreq']
# raw.resample(300, npad = 'auto')  # resample to 300Hz, do not do this
sfreq = raw.info['sfreq'] # get new samling frequency


# remove line noise
# jobs = 1
line_freq = [50.0,100.0]
raw._data = notch_filter(raw.get_data(),
                         Fs=sfreq,
                         freqs=line_freq,
                         method="spectrum_fit",
                         mt_bandwidth=2,
                         p_value=0.05,
                         filter_length="10s",
                         n_jobs=4)

# raw.compute_psd(fmax = 150.0).plot(average = True)

# fig, ax = plt.subplots(1,1)
# raw.compute_psd(n_jobs=jobs).plot(axes=ax)
# ax.set_title('Subejct %s' % subject)

# %%
# extract relevant parts of the recording

################################################################
# Subset to Sternberg task
################################################################

# eeg events
custom_mapping = {"Stimulus/S 16": 16, "Stimulus/S 17": 17}
(cue_evs, event_dict) = events_from_annotations(
    raw, event_id=custom_mapping
)

events = events_from_annotations(raw)
cue_evs

#Define start and end of Sternberg task
tmin = (cue_evs[0,0])/sfreq  # Start of Sternberg task
tmax = (cue_evs[1,0])/sfreq  # End of Sternberg task

#Crop Sternberg task based on start and end
raw_b = raw.copy().crop(tmin=tmin, tmax=tmax)

# Plot the raw datasets (fill raw_b[] with block number)
# raw_b.plot()

#############################################################
# Apply ICA to copy of data
#############################################################
# %%
# bad channel detection and interpolation

# make a copy of the data in question
raw_copy = raw_b.copy()
raw_copy.load_data()

# apply a 100Hz low-pass filter to data
raw_copy = raw_copy.filter(l_freq=None, h_freq=100.0,
                           picks=['eeg', 'eog'],
                           filter_length='auto',
                           l_trans_bandwidth='auto',
                           h_trans_bandwidth='auto',
                           method='fir',
                           phase='zero',
                           fir_window='hamming',
                           fir_design='firwin',
                           n_jobs=4)

# find bad channels !remove bad channels from copy data
noisy_dectector = NoisyChannels(raw_copy, random_state=42, do_detrend=True)
noisy_dectector.find_all_bads(ransac=False)

# create summary for PyPrep output
bad_channels = {'bads_by_deviation:': noisy_dectector.bad_by_deviation,
                'bads_by_hf_noise:': noisy_dectector.bad_by_hf_noise,
                'bads_by_correlation:': noisy_dectector.bad_by_correlation,
                'bads_by_SNR:': noisy_dectector.bad_by_SNR}

# interpolate the identified bad channels
raw_b.info['bads'] = noisy_dectector.get_bads()
raw_b.load_data()
raw_b.interpolate_bads(mode='accurate')

# %%
# prepare ICA

# set eeg reference
# raw_bl = raw_bl.set_eeg_reference('average')

# set ICA parameters
method = 'infomax'
reject = dict(eeg=250e-6)
ica = ICA(n_components=0.95,
          method=method,
          fit_params=dict(extended=True))

# make copy of raw with 1Hz high-pass filter
# raw_4_ica = raw_b.copy().filter(l_freq=1.0, h_freq=100.0, n_jobs=1)

raw_4_ica = raw_b.copy().filter(l_freq=1.0, h_freq=100.0, n_jobs=4)

# run ICA
ica.fit(raw_4_ica,
        reject=reject,
        reject_by_annotation=True)

# %%
# find bad components using ICA label
ic_labels = label_components(raw_4_ica, ica, method="iclabel")

labels = ic_labels["labels"]
exclude_idx = [idx for idx, label in
               enumerate(labels) if label not in ["brain"]]

logger.info(f"Excluding these ICA components: {exclude_idx}")

# exclude the identified components and reconstruct eeg signal
ica.exclude = exclude_idx
ica.apply(raw_b)

# clean up
del raw_4_ica

# %%
# apply filter to data, 40Hz
raw_concatenated = raw_b.filter(l_freq=0.01, h_freq=40.0,
                       picks=['eeg', 'eog'],
                       filter_length='auto',
                       l_trans_bandwidth='auto',
                       h_trans_bandwidth='auto',
                       method='fir',
                       phase='zero',
                       fir_window='hamming',
                       fir_design='firwin',
                       n_jobs=4)


############################################################
# Create report and save data
############################################################

# %%
# create path for preprocessed data
FPATH_DERIVATIVES = r'C:\Users\micha\OneDrive\Documents\EEGvar' # set path
FPATH_PREPROCESSED = os.path.join(
    FPATH_DERIVATIVES,
    'preprocessing',
    'sub-%s' % sub,
    'eeg',
    'sub-%s_task-%s_preprocessed-raw.fif' % (sub, 'dpx')) # changed subject

# check if directory exists
if not Path(FPATH_PREPROCESSED).exists():
    Path(FPATH_PREPROCESSED).parent.mkdir(parents=True, exist_ok=True)

# save file
# overwrite = True
# if overwrite:
#    logger.info("`overwrite` is set to ``True`` ")

raw_concatenated.save(FPATH_PREPROCESSED)

# , overwrite=overwrite)
report = True
# %%
if report:
    # make path
    FPATH_REPORT = os.path.join(
        FPATH_DERIVATIVES,
        'preprocessing',
        'sub-%s' % sub,
        'report') # changed from report
    # create path on the fly
    if not Path(FPATH_REPORT).exists():
        Path(FPATH_REPORT).mkdir(parents=True, exist_ok=True)

# create data report
bidsdata_report = Report(title='Subject %s' % sub)
bidsdata_report.add_raw(raw=raw_concatenated, title='Raw data',
                            butterfly=False,
                            replace=True,
                            psd=True)

# add bad channels
bads_html = """
    <p>Bad channels identified by PyPrep:</p>
    <p>%s</p> 
    """ % '<br> '.join([key + ' ' + str(val)
                        for key, val in bad_channels.items()])
bidsdata_report.add_html(title='Bad channels',
                             tags='bads',
                             html=bads_html,
                             replace=True)

# add ica
fig = ica.plot_components(show=False, picks=np.arange(ica.n_components_))
plt.close('all')

bidsdata_report.add_figure(
        fig=fig,
        tags='ica',
        title='ICA cleaning',
        caption='Bad components identified by ICA Label: %s' % ', '.join(
            str(ix) + ': ' + labels[ix] for ix in exclude_idx),
        image_format='PNG',
        replace=True
    )

for rep_ext in ['hdf5', 'html']:
        FPATH_REPORT_O = os.path.join(
            FPATH_REPORT,
            'sub-%s_task-%s_preprocessing_report.%s'
            % (sub, 'dpx', rep_ext))

bidsdata_report.save(FPATH_REPORT_O,
                                 # overwrite=overwrite,
                                 open_browser=False)

###########################################################################################################################
# PREPROCESSING END HERE
###########################################################################################################################