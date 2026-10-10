%% estimateUOLRecordingDurations.m
%
% Loops through every animal entry and every recording_area
% (dataStruct.seriesData entry) inside allData_uol.mat, estimating each
% area's duration from the number of columns (time samples) in
% shankData.shank1's MUAs trace divided by the data sampling rate
% conf.samplingParams.srData. A single recording_area ID is a recording
% name (animalID_s + a 14-digit YYYYMMDDHHMMSS timestamp) followed
% directly by an area ID suffix (e.g. 'M190128_A_MD_s201902050857511' is
% recording 'M190128_A_MD_s20190205085751', area '1') -- see
% determineArea.m's own areaID doc ("seriesName(15:end)") and
% seriesFromEntry.m, both reused here. Since one recording covers several
% areas, only a single duration per recording is reported: the maximum
% across all of that recording's areas.
%
% allData_uol.mat is very large (tens of GB) and cannot be loaded with a
% plain load() call. Instead this script inspects it directly as an HDF5
% file (a v7.3 MAT-file is HDF5 under the hood) using
% h5_list_group_children and h5_dataset_size, which read only object
% metadata (group listings and dataset dimensions) without touching the
% underlying data arrays, so the whole sweep runs in seconds.
%
% Note: shankData.shank1 itself is a struct (group), not a data matrix, so
% its "number of columns" is taken from its MUAs field (a 1 x nSamples
% multi-unit-activity trace sampled at srData) -- set shankField below to
% use a different sub-field (e.g. 'spk') if needed.
%
% See also: h5_list_group_children, h5_dataset_size, seriesFromEntry, determineArea

dataFile = 'D:\Datasets\infraslow-dynamics\04_data_analysis\001_uol\noRun\allData_uol.mat';
shankField = 'MUAs';
timestampLength = 14; % length of the YYYYMMDDHHMMSS recording timestamp

animalIDs = h5_list_group_children(dataFile, '/allData');

for iAnimal = 1:numel(animalIDs)
  animalID = animalIDs{iAnimal};
  seriesPath = ['/allData/' animalID '/dataStruct/seriesData'];
  recordingAreas = h5_list_group_children(dataFile, seriesPath);
  nAreas = numel(recordingAreas);

  % Duration of every individual recording_area entry, plus the recording
  % name (area-ID suffix stripped) it belongs to.
  areaDurations = nan(nAreas, 1);
  recordingNames = cell(nAreas, 1);
  for iArea = 1:nAreas
    recordingArea = recordingAreas{iArea};
    basePath = [seriesPath '/' recordingArea];

    try
      [seriesDigits, animalPrefix] = seriesFromEntry(recordingArea);
      recordingNames{iArea} = [animalPrefix '_s' seriesDigits(1:timestampLength)];
    catch
      recordingNames{iArea} = recordingArea; % fallback: treat as its own recording
    end

    try
      shankSize = h5_dataset_size(dataFile, [basePath '/shankData/shank1/' shankField]);
      nSamples = max(shankSize); % robust to MUAs being stored as [1 x N] or a plain N-length vector
      srData = h5read(dataFile, [basePath '/conf/samplingParams/srData']);
      areaDurations(iArea) = nSamples/srData; % seconds
    catch ME
      fprintf('%s | %s: FAILED to estimate duration (%s)\n', ...
        animalID, recordingArea, ME.message);
    end
  end

  % Report one duration per recording: the maximum across its areas.
  uniqueRecordings = unique(recordingNames, 'stable');
  for iRec = 1:numel(uniqueRecordings)
    recordingName = uniqueRecordings{iRec};
    inRecording = strcmp(recordingNames, recordingName);
    recordingDuration = max(areaDurations(inRecording));
    fprintf('%s | %s: duration = %.2f s (%.2f min)\n', ...
      animalID, recordingName, recordingDuration, recordingDuration/60);
  end
end
