'''
NWB parameters common across all animals and recording sessions
Parameters that may change:
  dataset
'''

import pathlib
import os
import re

global rawDataFolder, derivedDataFolder, derivedDataFolderNWB

projectName = 'Brainwide Infraslow Activity Dynamics'
experimenter = 'Martynas Dervinis'
institution = 'University of Leicester'
publications = []
lab = 'Michael Okun lab'
dataset = 'neuropixels'
videoFrameRate = 25.0 # Hz

# Find the repository root folder
rootFolderName = 'infraslow-dynamics'
rootFolderCandidate = pathlib.Path(__file__).parent.resolve()
endInd = re.search(rootFolderName, str(rootFolderCandidate)).end()
rootFolder = str(rootFolderCandidate)[:endInd]

# Assign data folders
if 'neuronexus' in dataset:
  rawDataFolder = os.path.join(rootFolder, '03_data', '001_uol_neuronexus_exp_raw_derived'); # Raw and certain derived Neuronexus data for all animals
  derivedDataFolder = os.path.join(rootFolder, '04_data_analysis', '001_uol_neuronexus_exp', 'runNSG'); # Neuroscience Gateway analysis scripts, functions, and downloaded processed Neuronexus data for all animals
  derivedDataFolderNWB = os.path.join(rootFolder, '03_data', '101_uol_neuronexus_exp_derived_data', 'NWB'); # Derived Neuronexus data for all animals placed in NWB format
elif 'neuropixels' in dataset:
  rawDataFolder = os.path.join(rootFolder, '03_data', '002_uol_neuropixels_exp_raw_derived'); # Raw and certain derived Neuropixels data for all animals
  derivedDataFolder = os.path.join(rootFolder, '04_data_analysis', '002_uol_neuropixels_exp', 'runNSG'); # Neuroscience Gateway analysis scripts, functions, and downloaded processed Neuropixels data for all animals
  derivedDataFolderNWB = os.path.join(rootFolder, '03_data', '102_uol_neuropixels_exp_derived_data', 'NWB'); # Derived Neuropixels data for all animals placed in NWB format