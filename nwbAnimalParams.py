'''
NWB parameters for an animal
Parameters that are likely to change:
  animalID
  dob
  firstProcedureDate
  sex
  weight
  description
'''

import datetime
import os

animalID = 'M191106_MD'
dob = '20190905'
dob = datetime.datetime(int(dob[0:4]), int(dob[4:6]), int(dob[6:8])) # Convert to datetime format
firstProcedureDate = '20191106'
firstProcedureDate = datetime.datetime(int(firstProcedureDate[0:4]), int(firstProcedureDate[4:6]), int(firstProcedureDate[6:8]))
ageInDays = (firstProcedureDate - dob).days
age = 'P' + str(ageInDays) + 'D' # Convert to ISO8601 format: https://en.wikipedia.org/wiki/ISO_8601#Durations
strain = 'C57BL/6J'
sex = 'M'
species = 'Mus musculus'
weight = []
description = '013' # Animal testing order.
animalRawDataFolder = os.path.join(rawDataFolder, animalID) # Raw and certain derived data for the animal. Probe map and unit waveforms files are stored here.
animalDerivedDataFile = os.path.join(derivedDataFolder, 'runNSG_' + animalID, animalID + '.mat') # Spiking, behavioural, and certain analysis data produced after running analyses routines on the Neuroscience Gateway Portal are stored here in a single mat file.
animalDerivedDataFolderNWB = os.path.join(derivedDataFolderNWB, animalID) # The output folder: The location where derived data after converting to the NWB format is stored.