% Run this script to perform analyses comparing LFP and pupil area data.

%% LOAD PRE-PROCESSED DATA
load(dataFile);


%% INITIALISE PARAMETERS
lfpParams


%% LOAD RECORDINGS AND EXTRACT LFP MEASURES
% LOOP THROUGH DB ENTRIES
fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf
  dbEntries = 1:numel(fnsData);
end
for dbCount = dbEntries
  chOIDB = chOI{dbCount};
  if isempty(chOIDB)
    disp(['No channel of interest specified for ' fnsData{dbCount} '. Skippig to the next db entry...']);
    continue
  end
  
  % Load the contents of dbStruct
  [dbStruct, repository, dataDir, entryName, baseFilename, probe, ~,...
    ~, ~, ~, srRecording] = get_dbStruct(dataStruct, dbCount);
  chN = dbStruct.db(dbCount).chN;
  samplingParams = dbStruct.conf.samplingParams;
  
  if strcmp(repository,'uol')
    if strcmpi(probe, 'Neuropixels')
      strInd = strfind(baseFilename, 'ap');
      baseFilename = [baseFilename(1:strInd-1) 'lf.bin'];
      if ~exist(baseFilename, 'file')
        baseFilename = [baseFilename(1:strInd-1) 'lf.dat'];
      end
    else
      baseFilename = [baseFilename '.dat']; %#ok<*AGROW>
    end
    if ~exist(baseFilename, 'file')
      error(['The supplied LFP file does not exist. Please check the location folder: ' dataDir]);
    end
  else
    seriesName = seriesFromEntry(entryName);
    seriesName = seriesName(1:14);
    if dbCount > 1
      entryNamePrev = dbStruct.db(dbCount-1).entryName;
      seriesNamePrev = seriesFromEntry(entryNamePrev);
      seriesNamePrev = seriesNamePrev(1:14);
      if ~strcmpi(seriesName, seriesNamePrev)
        %loadcsv
      end
    else
      %loadcsv
    end
  end
  
  strSep = strfind(entryName,'s');
  areaStr = entryName(strSep+15:end);
  
  % ESTIMATE LFP FREQUENCY BAND POWER MEASURES
  options.bandRange = LFPbands;
  options.chunkSize = chunkSize;
  options.srInterpInit = dssrLFPinit;
  options.srInterpFinal = dssrLFPfinal;
  options.chOI = chOI{dbCount};
  options.deleteChans = deleteChans;
  if subtractMedian
    options.lfpCAR = 'subtract';
  else
    options.lfpCAR = 'none';
  end
  options.transformFunc = dbStruct.conf.syncFunc;
  options.powerCalcMethod = 'wavelet';
  options.saturationMethod = 'combined';
  if strcmpi(areaStr, '6')
    options.spectrogram = 1;
  else
    options.spectrogram = 0;
  end
  options.rippleDuration = rippleDuration;
  options.wGaussian = wGaussian;
  options.sdGaussian = sdGaussian;
  options.saturationPlot = false; %true;
  if strcmpi(probe, 'Neuropixels')
    lfpPowerDB = lfpPowers(baseFilename, chN, samplingParams.srRecordingLFP, options);
  else
    lfpPowerDB = lfpPowers(baseFilename, chN, srRecording, options);
  end
  lfpPowerDB.dssrLFPinit = dssrLFPinit;
  lfpPowerDB.dssrLFPfinal = dssrLFPfinal;
  lfpPowerDB.LFPbands = LFPbands;
  lfpPowerDB.chOIDB = chOI{dbCount};
  lfpPowerDB.lfpTimes = lfpPowerDB.time;
  fields = {'time','options'};
  lfpPowerDB = rmfield(lfpPowerDB,fields);
  
  % SAVE DATA
  dataString = ['dataStruct.seriesData.' entryName '.lfpPowerData = lfpPowerDB;'];
  eval(dataString);
  if intermediateSaving
    save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
  end
end

if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca