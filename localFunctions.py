'''
Local functions for convert2nwb
'''

import numpy as np
import h5py
from scipy import sparse
from pynwb.core import VectorData


def electrodeLocations(areas, endCh, nCh):
  '''
  electrodeLocationVec = electrodeLocations(areas, endCh, nCh)
  Function assigns brain area to a channel
  '''

  electrodeLocationVec = []
  for iCh in range(nCh):
    areaInd = np.nonzero(endCh >= iCh)[0].min()
    electrodeLocationVec.append(areas[areaInd])

  return electrodeLocationVec


def createElectrodeTable(nwb, input):
  '''
  (tbl, columnLabels, columnDescription) = createElectrodeTable(nwb, input)

  Function creates an electrode table with the following columns:
    channel_id: A unnique probe channel ID formed by combining session ID,
                probe reference number, and channel number relative to the
                tip of the probe.
    channel_local_index: Channel index relative to the tip of the probe.
                        Channel indices are only unique within a probe.
    x: Channel AP brain surface coordinate (probe inisertion location mm).
    y: Channel ML brain surface coordinate (probe inisertion location mm).
    z: Channel location relative to the tip of the probe in mm.
    imp: Channel impedance.
    location: Channel brain area location.
    filtering: Type of LFP filtering applied.
    group: Channel electrode group (e.g., shank 1). NWB documentation on
          ElectrodeGroup datatype is provided in the following links:
          https://nwb-schema.readthedocs.io/en/latest/format.html#electrodegroup
          https://nwb-schema.readthedocs.io/en/latest/format.html#sec-electrodegroup-src
    channel_label
    probe_label.
  The rows of the table correspond to individual recording channels.

    Input: nwb - a nwb object.
          input - a dictionary with the following keys:
            iElectrode: electrode reference number.
            electrodeDescription: a cell array (n probes) with probe
                                  desciptions.
            electrodeManufacturer: a cell array of electrode manufacturers.
            nShanks: a cell array of number of shanks.
            nChannelsPerShank: a cell array of electrode number of
                                recording channels per shank.
            electrodeLocation: a cell array (n channels) of channel brain
                                area locations.
            electrodeCoordinates: a cell array (n probes) with recording
                                  channel coordinate arrays (n channels by
                                  3).
            sessionID: a string with the session ID.
            electrodeLabel: a cell array (n probes) with probe labels.

    Output: tbl - a numpy array object with rows and columns as described
                  above.
            columnLabels - a list of column labels of the output table.
            columnDescription - a list of descriptions of data in the columns.
  '''

  # Parse input
  iEl = input["iElectrode"]
  nSh = input["nShanks"]
  nCh = input["nChannelsPerShank"]

  # Create a table with given column labels
  columnLabels = ['channel_id', 'channel_local_index', 'x', 'y', 'z', 'imp', 'location', 'filtering', 'group', 'channel_label', 'probe_label']
  columnDescription = [
    'A unnique probe channel ID formed by combining session ID, probe reference number, and channel number relative to the tip of the probe.',
    'Channel index relative to the tip of the probe. Channel indices are only unique within a probe.',
    'Channel AP brain surface coordinate (probe inisertion location mm).',
    'Channel ML brain surface coordinate (probe inisertion location mm).',
    'Channel location relative to the tip of the probe in mm.',
    'Channel impedance.',
    'Channel brain area location.',
    'Type of LFP filtering applied.',
    'Channel electrode group (e.g., shank 1).',
    'Channel_label.',
    'Probe_label.'
  ]
  tbl = np.empty((0, len(columnLabels)))

  # Register the probe device
  device = nwb.create_device(
    name = 'probe' + str(iEl+1),
    description = input["electrodeDescription"][iEl],
    manufacturer = input["electrodeManufacturer"][iEl])

  for iShank in range(nSh[iEl]):
    
    # Register a shank electrode group (only one because this is a single shank probe)
    electrode_group = nwb.create_electrode_group(
      name = 'probe' + str(iEl+1),
      description = 'electrode group for probe' + str(iEl+1),
      location = input["electrodeLocation"][iEl][-1],
      device = device,
      position = input["electrodeCoordinates"][iEl][-1])
    
    # Populate the electrode table
    for iCh in range(nCh[iEl]):
      if iCh < 10-1:
        channelID = str(input["sessionID"] + str(iEl+1) + '00' + str(iCh+1))
      elif iCh < 99-1:
        channelID = str(input["sessionID"] + str(iEl+1) + '0' + str(iCh+1))
      else:
        channelID = str(input["sessionID"] + str(iEl+1) + str(iCh+1))
      channel_label = 'probe' + str(iEl+1) + 'shank' + str(iShank+1) + 'elec' + str(iCh+1)
      tbl = np.append(tbl, np.matrix([
        channelID, iCh+1, input["electrodeCoordinates"][iEl][iCh][0], input["electrodeCoordinates"][iEl][iCh][1], input["electrodeCoordinates"][iEl][iCh][2],
        np.nan, input["electrodeLocation"][iEl][iCh], 'unknown', electrode_group, channel_label, input["electrodeLabel"][iEl]]), axis=0)
    
  return np.array(tbl), columnLabels, columnDescription


def array2list(tbl, columnLabels, columnDescription):
# Function definition

  columnList = []
  if tbl.any():
    tblColumns=tbl.transpose().tolist()
    for iCol in range(len(columnLabels)):
      columnList.append(VectorData(
        name=columnLabels[iCol],
        data=tblColumns[iCol],
        description=columnDescription[iCol]
      ))
  
  return columnList


def concatenateMat(mat1, mat2, method='vertical'):
  '''
  concatenatedMat = concatenateMat(mat1, mat2, method)

  Function concatenates two 2-D numpy arrays by appending the
  smaller matrix with trailing zeros or NaNs.
  Input: mat1.
         mat2.
         method - concatenate either vertically ('vertical') or
                  horizontally ('horizontal'). This input is optional and
                  the default method is vertical. In case when trailing
                  NaNs are needed rather than zeros, corresponding methods
                  are 'verticalnan' and 'horizontalnan', respectively.
  Output: concatenatedMat - a concatenated numpy array.
  '''

  if mat1.any() and mat2.any():
    if method == 'vertical':
      diff = mat1.shape[1] - mat2.shape[1]
      if diff > 0:
        trailingMat = np.zeros([mat2.shape[0], abs(diff)])
        mat2 = np.append(mat2, trailingMat, axis=1)
      elif diff < 0:
        trailingMat = np.zeros([mat1.shape[0], abs(diff)])
        mat1 = np.append(mat1, trailingMat, axis=1)
      concatenatedMat = np.append(mat1, mat2, axis=0)
    elif method == 'horizontal':
      diff = mat1.shape[0] - mat2.shape[0]
      if diff > 0:
        trailingMat = np.zeros([abs(diff), mat2.shape[1]])
        mat2 = np.append(mat2, trailingMat, axis=0)
      elif diff < 0:
        trailingMat = np.zeros([abs(diff), mat1.shape[1]])
        mat1 = np.append(mat1, trailingMat, axis=0)
      concatenatedMat = np.append(mat1, mat2, axis=1)
    elif method == 'verticalnan':
      diff = mat1.shape[1] - mat2.shape[1]
      if diff > 0:
        trailingMat = np.nan([mat2.shape[0], abs(diff)])
        mat2 = np.append(mat2, trailingMat, axis=1)
      elif diff < 0:
        trailingMat = np.nan([mat1.shape[0], abs(diff)])
        mat1 = np.append(mat1, trailingMat, axis=1)
      concatenatedMat = np.append(mat1, mat2, axis=0)
    elif method == 'horizontalnan':
      diff = mat1.shape[0] - mat2.shape[0]
      if diff > 0:
        trailingMat = np.nan([abs(diff), mat2.shape[1]])
        mat2 = np.append(mat2, trailingMat, axis=0)
      elif diff < 0:
        trailingMat = np.nan([abs(diff), mat1.shape[1]])
        mat1 = np.append(mat1, trailingMat, axis=0)
      concatenatedMat = np.append(mat1, mat2, axis=1)    
  elif mat1.any() and not mat2.any():
    concatenatedMat = mat1
  elif not mat1.any() and mat2.any():
    concatenatedMat = mat2
  else:
    concatenatedMat = np.array()
      
  return concatenatedMat


def getSpikes(animalDerivedDataFile, animalID, sessionID, electrodeTbl):
  '''
  (spikes, metadataTbl, derivedData, columnLabels, columnDescription) = getSpikes(animalDerivedDataFile, animalID, sessionID, electrodeTbl)

  Function loads Neuronexus spiking data from a MAT file with a custom data
  structure. Input:
    animalDerivedDataFile - a string with derived data file name or already
                            loaded data.
    animalID - an animal ID string.
    sessionID - a session of interest ID string.
    electrodeTbl - a numpy array with electrode information generated by
                   the function createElectrodeTable.
  Output: spikes - a 1-by-n list of numpy arrays (n units) with unit spike
                   times in seconds.
          metadataTbl - a numpy array with rows corresponding to
                        individual clusters (units) and columns to:
            cluster_id: a unique cluster ID formed by combining session
                        ID, probe reference number, and unit cluster ID.
            local_cluster_id: a unit cluster ID. This is only unique
                              within the probe.
            type - activity type: single unit (unit) or multi-unit (mua).
            peak_channel_index: recording channel with the highest unit peak
                          index relative to the tip of the probe.
            peak_channel_id: a corresponding unnique probe channel ID formed by
                        combining session ID, probe reference number, and
                        channel number relative to the tip of the probe.
            local_peak_channel_id: a corresponding channel index relative to the
                              tip of the probe. Channel indices are only
                              unique within a probe.
            rel_horz_position: relative horizontal position in um.
            rel_vert_position: probe tip-relative vertical position in um.
            isi_violations: interspike interval violations, a cluster
                            quality measure.
            isolation_distance: cluster isolation distance, a cluster
                            quality measure.
            area: unit brain area location.
            probe_id: probe label.
            electrode_group: channel electrode group (e.g., shank 1). NWB
                            documentation on ElectrodeGroup datatype is
                            provided in the following links:
                            https://nwb-schema.readthedocs.io/en/latest/format.html#electrodegroup
                            https://nwb-schema.readthedocs.io/en/latest/format.html#sec-electrodegroup-src
          derivedData - animal data loaded from the MAT derived data file.
          columnLabels - a list of column labels of the output table.
          columnDescription - a list of descriptions of data in the columns.
  '''

  # Column labels of the metadata table
  columnLabels = ['cluster_id', 'local_cluster_id', 'type', 'peak_channel_index', 'peak_channel_id', 'local_peak_channel_id',
                  'rel_horz_pos', 'rel_vert_pos', 'isi_violations', 'isolation_distance', 'area', 'probe_id', 'electrode_group']
  columnDescription = [
    'A unique cluster ID formed by combining session ID, probe reference number, and unit cluster ID.',
    'A unit cluster ID. This is only unique within the probe.',
    'Single unit (unit) or multi-unit (mua).',
    'A Recording channel with the highest unit peak index relative to the tip of the probe.',
    'A corresponding unnique probe channel ID formed by combining session ID, probe reference number, and channel number relative to the tip of the probe.',
    'A corresponding channel index relative to the tip of the probe. Channel indices are only unique within a probe.',
    'Relative horizontal position in um.',
    'Probe tip-relative vertical position in um.',
    'Interspike interval violations, a cluster quality measure.',
    'Cluster isolation distance, a cluster quality measure.',
    'Unit brain area location.',
    'Probe label.',
    'Channel electrode group (e.g., shank 1).'
  ]
  
  # Data series names with different brain areas
  if isinstance(animalDerivedDataFile, str):
    derivedData = h5py.File(animalDerivedDataFile,'r')
  else:
    derivedData = animalDerivedDataFile
  dataSeriesNames = []
  for iSeries in range(11):
    dataSeriesNames.append(animalID + '_s' + sessionID + str(iSeries+1))
  dataSeriesNames.append(animalID + '_s' + sessionID)

  # Load data
  spikes = []; metadataTbl = np.array([])
  for iSeries in range(len(dataSeriesNames)):
    metadata = np.array([])
    seriesDerivedData = np.array(derivedData.get('dataStruct/seriesData/' + dataSeriesNames[iSeries]))
    if len(seriesDerivedData.shape):
      srData = np.array(derivedData.get('dataStruct/seriesData/' + dataSeriesNames[iSeries] + '/conf/samplingParams/srData'))

      # Series spike array
      spikesSeries_data = np.array(derivedData.get('dataStruct/seriesData/' + dataSeriesNames[iSeries] + '/popData/spkDB/data'))
      if not (spikesSeries_data == None).all():
        spikesSeries_ir = np.array(derivedData.get('dataStruct/seriesData/' + dataSeriesNames[iSeries] + '/popData/spkDB/ir'))
        spikesSeries_jc = np.array(derivedData.get('dataStruct/seriesData/' + dataSeriesNames[iSeries] + '/popData/spkDB/jc'))
        spikesSeries = sparse.csc_matrix((spikesSeries_data, spikesSeries_ir, spikesSeries_jc)).toarray()
      else:
        spikesSeries = np.array([])
      
      nRows = spikesSeries.shape[0]
      if nRows:
        
        # Unit metadata: [local_unit_id type local_probe_channel horizontal_position vertical_position ...
        #                 isi_violations isolation_distance anterior_posterior_ccf_coordinate ...
        #                 dorsal_ventral_ccf_coordinate left_right_ccf_coordinate]
        metadata = concatenateMat(metadata, np.transpose(np.array(derivedData.get('dataStruct/seriesData/' + dataSeriesNames[iSeries] + '/popData/muaMetadata'))))
        
        # Spike times
        if nRows != len(metadata): # means that h5py did not load the sparse matrix correctly. There are redundant rows
          spikesSeries = spikesSeries[:metadata.shape[0],:]
          nRows = spikesSeries.shape[0]
        nColumns = spikesSeries.shape[1]
        timeVector = np.linspace(1,nColumns,nColumns)/srData
        for iUnit in range(nRows):
          spikes.append(timeVector[0][spikesSeries[iUnit].astype(np.bool)])

        # Unit metadata: [metadata area]
        if iSeries+1 == 1:
          areas = ['S1'] * nRows
        elif iSeries+1 == 2:
          areas = ['VB'] * nRows
        elif iSeries+1 == 3:
          areas = ['Po'] * nRows
        elif iSeries+1 == 4:
          areas = ['LP'] * nRows
        elif iSeries+1 == 5:
          areas = ['DG'] * nRows
        elif iSeries+1 == 6:
          areas = ['CA1'] * nRows
        elif iSeries+1 == 7:
          areas = ['RSC'] * nRows
        elif iSeries+1 == 8:
          areas = ['VB'] * nRows
        elif iSeries+1 == 9:
          areas = ['LP'] * nRows
        elif iSeries+1 == 10:
          areas = ['LGN'] * nRows
        elif iSeries+1 == 11:
          areas = ['CA3'] * nRows
        elif iSeries+1 == 12:
          areas = ['VB'] * nRows
        metadata = concatenateMat(metadata.astype(object), np.matrix(areas).transpose().astype(object), 'horizontal')
        
        # Unit metadata: correct unit type
        units = np.array(derivedData.get('dataStruct/seriesData/' + dataSeriesNames[iSeries] + '/shankData/shank1/units'))
        muas = np.array(derivedData.get('dataStruct/seriesData/' + dataSeriesNames[iSeries] + '/popData/spkDB_units'))
        unitInds = np.isin(muas[0], units[0])
        unitTypes = np.array(['mua'] * len(muas[0]))
        unitTypes[np.array(unitInds)] = 'unit'
        metadata = metadata.astype(object)
        metadata[:,1] = np.matrix(unitTypes).transpose()
        
        # Unit metadata: [metadata probe_id]
        if iSeries+1 == 1:
          probeLabel = ['probe1'] * nRows
        elif iSeries+1 == 2:
          probeLabel = ['probe2'] * nRows
        elif iSeries+1 == 3:
          probeLabel = ['probe2'] * nRows
        elif iSeries+1 == 4:
          probeLabel = ['probe2'] * nRows
        elif iSeries+1 == 5:
          probeLabel = ['probe2'] * nRows
        elif iSeries+1 == 6:
          probeLabel = ['probe2'] * nRows
        elif iSeries+1 == 7:
          probeLabel = ['probe2'] * nRows
        elif iSeries+1 == 8:
          probeLabel = ['probe1'] * nRows
        elif iSeries+1 == 9:
          probeLabel = ['probe1'] * nRows
        elif iSeries+1 == 10:
          probeLabel = ['probe1'] * nRows
        elif iSeries+1 == 11:
          probeLabel = ['probe1'] * nRows
        elif iSeries+1 == 12:
          probeLabel = ['probe1'] * nRows
        metadata = concatenateMat(metadata.astype(object), np.matrix(probeLabel).transpose().astype(object), 'horizontal')

        # Unit metadata: [unit_id metadata]
        unitIDs = []
        for iUnit in range(nRows):
          if metadata[iUnit, -1] == 'probe1':
            unitID = str(sessionID) + '1'
          else:
            unitID = str(sessionID) + '2'
          if metadata[iUnit, 0] < 9:
            unitID = unitID + '000' + str(metadata[iUnit, 0])
          elif metadata[iUnit, 0] < 99:
            unitID = unitID + '00' + str(metadata[iUnit, 0])
          elif metadata[iUnit, 0] < 999:
            unitID = unitID + '0' + str(metadata[iUnit, 0])
          else:
            unitID = unitID + str(metadata[iUnit, 0])
          unitIDs.append(unitID)
        metadata = concatenateMat(np.matrix(unitIDs).transpose().astype(object), metadata, 'horizontal')
        
        # Unit metadata: [metadata[:,:3] probe_channel_index probe_channel_id metadata[:,3:] electrode_group]
        channelIndices = []
        channelIDs = []
        electrodeGroups = []
        for iUnit in range(nRows):
          ind = np.logical_and(np.isin(np.array(electrodeTbl[1].data), np.array(metadata[iUnit,3])), \
            np.isin(electrodeTbl[-1].data, metadata[iUnit,-1])) # Find channel ID on a particular probe and get its index and ID
          channelIndices.append(np.where(ind)[0]+1)
          channelIDs.append(electrodeTbl[0].data[ind[0]])
          electrodeGroups.append(electrodeTbl[8].data[ind[0]])
        metadataInit = concatenateMat(metadata[:,:3], np.matrix(channelIndices).astype(object), 'horizontal')
        metadataInit = concatenateMat(metadataInit, np.matrix(channelIDs).transpose().astype(object), 'horizontal')
        metadataInit = concatenateMat(metadataInit, metadata[:,3:], 'horizontal')
        metadata = concatenateMat(metadataInit, np.matrix(electrodeGroups).transpose().astype(object), 'horizontal')
        metadataTbl = concatenateMat(metadataTbl, metadata)

  return spikes, metadataTbl, derivedData, columnLabels, columnDescription


def reshapeWaveforms(waveforms, iEl, metadata):
  '''
  waveformsMean = reshapeWaveforms(waveforms, iEl, metadata)

  Function extracts relevant waveform information.
  Input: waveforms - a strucuture loaded from the waveforms MAT file.
                    Relevant fields are waveforms (described above),
                    maxWaveforms (same as waveforms but excluding all
                    channels except for the maximum amplitude one), and
                    cluIDs (unit cluster IDs corresponding to the
                    dimension one in waveforms and maxWaveforms).
        iEl - probe reference number.
        metadata - a numpy array unit table produced by the function
                   getSpikes.
        nCh - number of recording channels with waveforms for the same
              unit.
  Output: waveformsMean - waveforms.waveforms converted into a numpy array.
                          Rows correspond to units. MUAs are NaNs.
  '''

  # Load data fields
  if len(waveforms):
    maxWaveforms = np.array(waveforms.get('maxWaveforms')).transpose()
    cluIDs = np.array(waveforms.get('cluIDs')).transpose()
  else:
    maxWaveforms = []
    cluIDs = []

  # Load waveforms
  metadataInds = np.where(np.array(np.isin(metadata[:,11], 'probe' + str(iEl))))[0]
  metadata = metadata[metadataInds]
  waveformsMean = []
  if len(maxWaveforms) and len(maxWaveforms.shape):
    nWaveformSamples = maxWaveforms.shape[1]
  else:
    nWaveformSamples = 200

  for iUnit in range(metadata.shape[0]):
    row = np.where(np.isin(cluIDs, metadata[iUnit,1]))[0]
    if row.size:
      waveformsMean.append(maxWaveforms[row])
    else:
      waveformsMean.append(np.full([1,nWaveformSamples], np.nan))
  
  return waveformsMean


def parsePeriod(acceptablePeriod, derivedData):
  '''
  acceptablePeriod = parsePeriod(acceptablePeriod, derivedData)

  Function resolves acceptable behavioural data periods. In cases when
  there are multiple periods, each of them needs to be loaded separately.
  Input: acceptablePeriod - behavioural data period loaded from a MAT
                            file. It should be a numpy array.
         derivedData - passively loaded data from a MAT file.
  Output: acceptablePeriod - a numpy array of acceptable behavioural data
                             periods.
  '''
  if len(acceptablePeriod.shape):
    if isinstance(acceptablePeriod[0], np.ndarray) and isinstance(acceptablePeriod[0][0], h5py.h5r.Reference):
        acceptablePeriodNew = []
        for iPeriod in range(len(acceptablePeriod)):
          data_of_interest_reference = acceptablePeriod[iPeriod, 0]
          acceptablePeriodNew.append(np.array(derivedData[data_of_interest_reference]))
        acceptablePeriod = np.matrix(np.array(acceptablePeriodNew).squeeze())
    else:
      acceptablePeriod = np.matrix(acceptablePeriod)

  return acceptablePeriod


def markQualitySamples(acceptablePeriod, videoFrameTimes):
  '''
  acceptableSamples = markQualitySamples(acceptablePeriod, videoFrameTimes)

  Function marks acceptable behavioural samples given the sample times and
  the range of acceptable time periods.
  Input: acceptablePeriod - a numpy array of row vectors marking the
                            beginning and end of acceptable time periods.
         videoFrameTimes - a numpy array (a row vector) with sample times.
  Ouptut: acceptableSamples - a logical vector marking acceptable samples
                              by ones.
  '''

  if not len(acceptablePeriod) or not len(videoFrameTimes):
    acceptableSamples = []
  else:
    acceptableSamples = np.full((1, videoFrameTimes.shape[1]), False)
    for iPeriod in range(acceptablePeriod.shape[0]):
      acceptableSamples[np.logical_and(videoFrameTimes >= acceptablePeriod[iPeriod][0,0],
        videoFrameTimes <= acceptablePeriod[iPeriod][0,1])] = True

  return acceptableSamples+0