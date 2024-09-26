'''
Session-specific NWB parameters
Parameters that are different for different animals are:
  sessionID
  sessionDescription
  sessionNotes
  areas
  endCh
  probeInserted
  electrodeName
  nChannelsPerShank
  electrodeCoordinates
  electrodeImplantationType
'''

import numpy as np
import datetime
import os
import h5py
import re
from localFunctions import electrodeLocations

# General info
sessionID = ['20191106163201']
sessionDescription = ['anaesthesia']
sessionNotes = ['4200 um in D; isoflurane at 1% throughout the recording.']
endCh = [[None] * 2 for i in range(len(sessionID))]
endCh[0][0] =  np.array([]) # Corresponding probe end channels starting from the tip of the probe. Corresponding and previous end channels are used to work out probe channels that reside in the corresponding brain area.
endCh[0][1] = np.array([30, 98, 138, 220, 306, 384])

# Initialise variables
sessionStartTime = [[None] for i in range(len(sessionID))]
probeInserted = [[None] * 2 for i in range(len(sessionID))]
electrodeName = [[None] * 2 for i in range(len(sessionID))]
electrodeDescription = [[None] * 2 for i in range(len(sessionID))]
electrodeManufacturer = [[None] * 2 for i in range(len(sessionID))]
electrodeFolder = [[None] * 2 for i in range(len(sessionID))]
nShanks = [[None] * 2 for i in range(len(sessionID))]
nChannelsPerShank = [[None] * 2 for i in range(len(sessionID))]
nCh = [[None] * 2 for i in range(len(sessionID))]
electrodeLocation = [[None] * 2 for i in range(len(sessionID))]
electrodeCoordinates = [[None] * 2 for i in range(len(sessionID))]
electrodeLabel = [[None] * 2 for i in range(len(sessionID))]
electrodeImplantationType = [[None] * 2 for i in range(len(sessionID))]

# Probes info
for iSess in range(len(sessionID)):
  sessionStartTime[iSess] = datetime.datetime(int(sessionID[iSess][0:4]), int(sessionID[iSess][4:6]), int(sessionID[iSess][6:8]), \
    int(sessionID[iSess][8:10]), int(sessionID[iSess][10:12]), int(sessionID[iSess][12:14]))
  
  # Probe #1 info
  ref = 0 # probe reference
  probeInserted[iSess][ref] = False # If the probe used at all
  if probeInserted[iSess][ref] and endCh[iSess][ref].any():
    electrodeName[iSess][ref] = 'Neuropixels 1.0'
    electrodeDescription[iSess][ref] = 'Single shank high density probe in position ' + str(ref+1)
    electrodeManufacturer[iSess][ref] = 'imec'
    electrodeFolder[iSess][ref] = os.path.join(animalRawDataFolder, sessionID[iSess] + '1')
    if not os.path.isdir(electrodeFolder[iSess][ref]): # in case a probe is missing
      electrodeFolder[iSess][ref] = os.path.normpath(electrodeFolder[iSess][ref][0:-1])
    if re.search('neuropixels', electrodeName[iSess][ref], re.IGNORECASE):
      electrodeMap = os.path.join(electrodeFolder[iSess][ref], 'forPRB_Neuropixels.mat');
    else:
      electrodeMap = os.path.join(electrodeFolder[iSess][ref], 'forPRB_' + electrodeName[iSess][ref] + '.mat')
    nShanks[iSess][ref] = 1
    nChannelsPerShank[iSess][ref] = 384
    nCh[iSess][ref] = nChannelsPerShank[iSess][ref]*nShanks[iSess][ref] # total number of probe channels
    areas = ['VB','LP','LGN','CA3','S1'] # brain areas that this probe spans
    electrodeLocation[iSess][ref] = electrodeLocations(areas, endCh[iSess][ref], nCh[iSess][ref]) # Brain area assigned to each recording channel.
    electrodeCoordinates[iSess][ref] = [-1.8, -2.5, 0] # Electrode insertion location on the cortical surface in Paxinos coords: AP (posterior negative), ML (left negative), DV (recording site position starting with the tip of the probe.
    electrodeCoordinates[iSess][ref] = [electrodeCoordinates[iSess][ref]] * nCh[iSess][ref] # Coordinates of each probe recording channel (the probe rotation angle is not taken into account). Y coordinates are relative to the tip of the probe.
    f = h5py.File(electrodeMap,'r')
    ycoords = np.array(f.get('ycoords'))/1000
    electrodeLabel[iSess][ref] = 'probe' + str(ref+1)
    electrodeImplantationType[iSess][ref] = 'acute'
  else: # The case when the probe #1 is missing
    probeInserted[iSess][ref] = False
    electrodeName[iSess][ref] = []
    electrodeDescription[iSess][ref] = []
    electrodeManufacturer[iSess][ref] = []
    electrodeFolder[iSess][ref] = []
    nShanks[iSess][ref] = []
    nChannelsPerShank[iSess][ref] = []
    nCh[iSess][ref] = []
    electrodeLocation[iSess][ref] = []
    electrodeCoordinates[iSess][ref] = []
    electrodeLabel[iSess][ref] = []
    electrodeImplantationType[iSess][ref] = []

  # Probe #2 info
  ref = 1 # probe reference
  probeInserted[iSess][ref] = True # If the probe used at all
  if probeInserted[iSess][ref] and endCh[iSess][ref].any():
    electrodeName[iSess][ref] = 'Neuropixels 1.0'
    electrodeDescription[iSess][ref] = 'Single shank high density probe in position ' + str(ref+1)
    electrodeManufacturer[iSess][ref] = 'imec'
    electrodeFolder[iSess][ref] = os.path.join(animalRawDataFolder, sessionID[iSess] + '26')
    if not os.path.isdir(electrodeFolder[iSess][ref]) or len(os.listdir(electrodeFolder[iSess][ref])) <= 0: # in case the raw data folder name is shorter
      electrodeFolder[iSess][ref] = os.path.normpath(electrodeFolder[iSess][ref][0:-1])
    if not os.path.isdir(electrodeFolder[iSess][ref]): # in case a probe is missing
      electrodeFolder[iSess][ref] = os.path.normpath(electrodeFolder[iSess][ref][0:-1])
    if re.search('neuropixels', electrodeName[iSess][ref], re.IGNORECASE):
      electrodeMap = os.path.join(electrodeFolder[iSess][ref], 'forPRB_Neuropixels.mat');
    else:
      electrodeMap = os.path.join(electrodeFolder[iSess][ref], 'forPRB_' + electrodeName[iSess][ref] + '.mat')
    nShanks[iSess][ref] = 1
    nChannelsPerShank[iSess][ref] = 384
    nCh[iSess][ref] = nChannelsPerShank[iSess][ref]*nShanks[iSess][ref]
    areas = ['VB','Po','LP','DG','CA1','RSC']
    electrodeLocation[iSess][ref] = electrodeLocations(areas, endCh[iSess][ref], nCh[iSess][ref])
    electrodeCoordinates[iSess][ref] = [-1.8, -0.5, 0]
    electrodeCoordinates[iSess][ref] = [electrodeCoordinates[iSess][ref]] * nCh[iSess][ref]
    f = h5py.File(electrodeMap,'r')
    ycoords = np.array(f.get('ycoords'))/1000
    electrodeLabel[iSess][ref] = 'probe' + str(ref+1)
    electrodeImplantationType[iSess][ref] = 'acute'
  else: # The case when the probe #2 is missing
    probeInserted[iSess][ref] = False
    electrodeName[iSess][ref] = []
    electrodeDescription[iSess][ref] = []
    electrodeManufacturer[iSess][ref] = []
    electrodeFolder[iSess][ref] = []
    nShanks[iSess][ref] = []
    nChannelsPerShank[iSess][ref] = []
    nCh[iSess][ref] = []
    electrodeLocation[iSess][ref] = []
    electrodeCoordinates[iSess][ref] = []
    electrodeLabel[iSess][ref] = []
    electrodeImplantationType[iSess][ref] = []