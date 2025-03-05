###################################################################
###############     Extract relevant epochs     ###################
###################################################################

"""
Extracts relevant data and removes artefacts.

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)

"""
# %%
# imports
import sys
import os

import warnings

from pathlib import Path

import numpy as np
import pandas as pd

from mne import events_from_annotations, Epochs

from mne.io import read_raw_fif
from mne.utils import logger

#from config import (
#    FPATH_DERIVATIVES,
#    MISSING_FPATH_BIDS_MSG,
#    SUBJECT_IDS,
#    BAD_SUBJECTS_SES_01,
#    BAD_SUBJECTS_SES_02
#)

# from utils import parse_overwrite

############################################
# Set default dictionary
############################################

# %%
# default settings (use subject 1, don't overwrite output files)
subject = 52
session = 1
task = 'multiple'
overwrite = False
report = False
jobs = 4

#############################################
# Load preprocessed Sternberg data 
#############################################
subject = 151
sub = "151"
task = "multiple"

FPATH_DERIVATIVES = r'C:\Users\micha\OneDrive\Documents\EEGvar'
# create path for preprocessed data
FPATH_PREPROCESSED = os.path.join(FPATH_DERIVATIVES,
                                  'Preprocessing',
                                  'sub-%s' %sub)

if overwrite:
    logger.info("`overwrite` is set to ``True`` ")

# %%
#  create path for import
FNAME = os.path.join(
    FPATH_PREPROCESSED,
    'eeg',
    'sub-%s_task-%s_preprocessed-raw.fif' % (sub, 'dpx')
)

# %%
# get the data
raw = read_raw_fif(FNAME)
raw.load_data()

# raw.plot()

# only keep eeg channels
# raw.pick(['eeg'])
# sampling rate
sfreq = raw.info['sfreq']

# events = events_from_annotations(raw)

############################################################
# Code to only keep memory sets
############################################################

custom_mapping = {"Stimulus/S 61": 61, "Stimulus/S 62": 62}
(cue_evs, event_dict) = events_from_annotations(
    raw, event_id=custom_mapping
)

# cue_evs

# Create 
stims = cue_evs[(cue_evs[:,2] == 61) | (cue_evs[:,2] == 62)]

from mne import Epochs

cue_type = list()
for i in range(stims.shape[0]):
    if stims[i, 2] == 62:
        cue_type.append("lure")
    else:
        cue_type.append('target')

import pandas as pd
# construct epochs metadata

metadata = {'subject': subject,
            'cue_type': cue_type,
            }
metadata = pd.DataFrame(metadata)

cue_epochs = Epochs(raw, stims, baseline=None, tmin=0, tmax=1.0, metadata=metadata)

# cue_epochs.plot()

############################################################
# Save data to preprocessing
############################################################

# create path for preprocessed data
FPATH_EPOCHS = os.path.join(FPATH_DERIVATIVES,
                            'epochs',
                            'sub-%s' % sub)

FPATH = os.path.join(FPATH_EPOCHS,
                            'sub-%s_task-%s_%scue-epo.fif' % (sub, task, '001'))

# check if directory exists
if not Path(FPATH).exists():
    Path(FPATH).parent.mkdir(parents=True, exist_ok=True)

# save file
# if overwrite:
    logger.info("`overwrite` is set to ``True`` ")

#save to path
cue_epochs.save(FPATH)

###############################################################################################################
# EPOCH END HERE
###############################################################################################################
