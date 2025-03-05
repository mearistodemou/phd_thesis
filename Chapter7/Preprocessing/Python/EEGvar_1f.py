########################################################################
###############     Estiamation of 1/fa noise        ###################
# 1. I am missing the arange command
# 2.
# 3.
# 4. 
#########################################################################

"""
===================================
Compute neural variability measures
===================================

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

from scipy.signal import periodogram

from mne import read_epochs
from mne.utils import logger

# from utils import parse_overwrite # does not work
from signal_variability import parallel_analysis

# from config import (
#    FPATH_DERIVATIVES,
#    MISSING_FPATH_BIDS_MSG,
#    SUBJECT_IDS,
#    BAD_SUBJECTS_SES_01,
#    BAD_SUBJECTS_SES_02
# )



######################################################################################
# Start here
######################################################################################

# %%
#  create path for import
subject = 151
sub = "151"
task = "multiple"

FPATH_DERIVATIVES = r'C:\Users\micha\OneDrive\Documents\EEGvar'
# create path for preprocessed data
FPATH_PREPROCESSED = os.path.join(FPATH_DERIVATIVES,
                                  'epochs',
                                  'sub-%s' %sub)

overwrite = True
if overwrite:
    logger.info("`overwrite` is set to ``True`` ")

# %%
#  create path for import
FNAME = os.path.join(
    FPATH_PREPROCESSED,
    'sub-%s_task-%s_%scue-epo.fif' % (sub, task, '001')
)

# %%
# get the data
epochs = read_epochs(FNAME, preload=True)
# epochs.plot()

# information needed for analyses
sfreq = epochs.info['sfreq']
n_epochs = len(epochs)
n_channels = len(epochs.ch_names)
metadata = epochs.metadata

# %%
psd_adds = np.array(
    ['subject', 'epoch', 'cue_type']
)

# %%
# run initial frequency analysis

# get EEG signal (this used to be epochs)
signal = epochs.get_data()

# compute power spectral density (PSD)
freqs, psd = periodogram(signal,
                         detrend='constant',
                         fs=sfreq,
                         nfft=sfreq * 2,
                         window='hamming')

# normalize psd (epoch- and sensor-wise)
for epo in range(psd.shape[0]):
    for ch in range(psd[epo, :, :].shape[0]):
        psd[epo, ch, :] = psd[epo, ch, :] / psd[epo, ch, :].sum(keepdims=True)

# %%
# make psd results dataframe
e, c, f = psd.shape
psd_st = np.column_stack(
    (
        np.repeat(metadata.subject, c),
        np.repeat(epochs.selection, c),
        np.tile(np.arange(0, 34), e), # I changed this
        psd.reshape(e*c, -1)
    )
)
psd_results = pd.DataFrame(
    psd_st,
    columns=np.concatenate(
        (psd_adds, np.array(['f_' + str(f) for f in freqs])))
)

# only keep relevant frequencies
psd_results.drop(
    psd_results.columns[np.arange(209, 1004)], # I made this change, because of an out of range error
    axis=1,
    inplace=True
)

# make column types for pandas dataframe
types_subj = {
    "subject": int,
    "epoch": int,
    "cue_type": int
}
types_fp = [np.float64 for i in np.arange(0, 201)]
types_fq = {'f_' + str(key): val
            for key, val
            in zip(freqs[freqs <= 100], types_fp)}
types_psd_results = types_subj | types_fq

# set column types
psd_results = psd_results.astype(types_psd_results)

# %%
# export PSD results to .tsv
FPATH_PSD = os.path.join(
    FPATH_DERIVATIVES,
    'power_spectral_density',
    'sub-%s_task-%s_psd.tsv' % (sub, task)
)

if not Path(FPATH_PSD).parent.exists():
    Path(FPATH_PSD).parent.mkdir(parents=True, exist_ok=True)

if os.path.exists(FPATH_PSD) and not overwrite:
    raise RuntimeError(
        f"'{FPATH_PSD}' already exists; consider setting 'overwrite' to True"
    )

psd_results.to_csv(FPATH_PSD, index=False, sep='\t', float_format='%.5f')

# tidy up
del psd_results, psd_st, e, c, f, epo, ch

# %%
# compute signal variability measures

# in the frequency domain
frequency_results, measures_fq = parallel_analysis(
    inst=psd, freqs=freqs, jobs=4)
# in the amplitude domain
amplitude_results, measures_amp = parallel_analysis(
    inst=signal, freqs=None, jobs=jobs)

# %%
# save frequency domain results

# make frequency results dataframe
meas, e, c = frequency_results.shape

fq = np.column_stack(
    (
        np.repeat(metadata.subject, c),
        np.repeat(epochs.selection, c),
        np.tile(np.arange(0, 34), e),
        frequency_results[0].reshape(e*c),
        frequency_results[1].reshape(e*c),
        frequency_results[2].reshape(e*c),
    )
)
fq_res = pd.DataFrame(
    fq,
    columns=np.concatenate(
        (psd_adds, measures_fq)
    )
)

# export frequency results results to .tsv
FPATH_FQ_VAR = os.path.join(
    FPATH_DERIVATIVES,
    'signal_variability',
    'sub-%s_task-%s_freq_var_single_trial.tsv'
    % (sub, task)
)

if not Path(FPATH_FQ_VAR).parent.exists():
    Path(FPATH_FQ_VAR).parent.mkdir(parents=True, exist_ok=True)

fq_res = fq_res.astype(
    {"subject": int, "epoch": int,
     "1f_offset": np.float64,
     "1f_exponent": np.float64,
     "spectral_entropy": np.float64,
     }
)

if os.path.exists(FPATH_FQ_VAR) and not overwrite:
    raise RuntimeError(
        f"'{FPATH_FQ_VAR}' already exists; consider setting 'overwrite' to True"
    )

fq_res.to_csv(FPATH_FQ_VAR, index=False, sep='\t', float_format='%.4f')

# tidy up
del fq_res, fq, meas, e, c

####################################################################################################################
# 1/f END HERE
####################################################################################################################
