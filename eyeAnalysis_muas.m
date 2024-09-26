% Run this script to perform analyses of correlation, phase, and coherence
% between MUA neural spiking and pupil area.


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITILIASE PARAMETERS
params
lists
intermediateSaving = false; % Save everytime a db series is finished being analysed (bad idea since saving is a lengthy process)


%% ADJUST THE NUMBER OF PARALLEL PROCESSORS (applies only when running on NSG Portal)
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end
if poolsize > 20 && poolsize < 25
  delete(gcp('nocreate'))
  parpool(20);
end


%% CORRELATE SPIKING ACTIVITY WITH EYE MEASURES
fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf % Here you can chose to execute only certain db entries
  dbEntriesLocal = 1:numel(fnsData);
else
  dbEntriesLocal = dbEntries;
end
for dbCount = dbEntriesLocal % Loop through dbEntries
  
  % Load the contents of dbStruct
  [dbStruct, ~, ~, entryName, ~, ~, shankIDs,...
    ~, period, ~, ~, srData, optCoh, exclRad, ~,...
    MUAsAll, spkDB, spkDB_units] = get_dbStruct(dataStruct, dbCount);
  optCoh.maxFreq = maxFreq_pupil;
  optCoh.winfactor = winfactor;
  optCoh.freqfactor = freqfactor;
  optCoh.monotoneFreq = true;
  MUAsAll = sum(MUAsAll,1);
  if isempty(MUAsAll) || ~sum(sum(MUAsAll))
    disp(['No spiking data for ' fnsData{dbCount} '. Skippig to the next db entry...']);
    continue
  end
  
  % Get eye data
  [seriesName, animal] = seriesFromEntry(entryName);
  entryNameEye = [animal '_s' seriesName(1:min([numel(seriesName) 14]))];
  if ~isfield(dataStruct, 'eyeData') || ~isfield(dataStruct.eyeData, entryNameEye) ||...
      isempty(dataStruct.eyeData.(entryNameEye).pupilArea)
    disp(['No pupil data for ' entryName '. Skippig to the next db entry...']);
    continue
  end
  eyeDataDB = dataStruct.eyeData.(entryNameEye);
  
  % Interpolate and filter pupil area data
  [seriesName, animal] = seriesFromEntry(entryName);
  runningSpeedDataEntry = [animal '_s' seriesName(1:min([14 numel(seriesName)]))];
  if excludeRunning && isfield(dataStruct, 'runningSpeedData') && ~isempty(dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod)
    period = combinePeriods(period, dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod, srData);
  end
  commonPeriod = combinePeriods(period, eyeDataDB.period, srData);
  if isempty(commonPeriod)
    continue
  end
  [areaInterpFilt, interpTimes, areaInd] = pupilFilt(eyeDataDB, srData, MUAsAll, 2*srData, commonPeriod, srData);
  
  % Get unit IDs
  units = [];
  phaseCoh = {};
  for sh = 1:numel(shankIDs) % Loop through shanks
    [~, ~, unitsSh] = get_shankStruct(dbStruct, sh);
    units = [units; unitsSh]; %#ok<*AGROW>
    if ~isempty(unitsSh)
      phaseCoh = [phaseCoh dbStruct.shankData.(['shank' num2str(sh)]).pupil.unitData];
    end
  end
  
  parfor u = 1:numel(spkDB_units)
  %for u = 1:numel(spkDB_units)
    fprintf('Started processing unit %i\n',spkDB_units(u));
    
    if ismember(spkDB_units(u), units)
      unitInd = units == spkDB_units(u);
      unitData = phaseCoh{unitInd}; %#ok<*PFBNS>
    else
      unitData = struct();
      unitData.unit = spkDB_units(u);
      spkOI = full(spkDB(u,:));
      
      % Phase and coherence analysis of unit spiking relative to pupil area
      [unitData.mfr, unitData.mfr_1sthalf, unitData.mfr_2ndhalf, unitData.lfr1, unitData.lfr5] = rateCalc(spkOI(areaInd),...
        srData);
      [unitData.psd_halves_freq, unitData.psd_halves, freqPSD, unitData.psd, unitData.psd_conf,...
        unitData.psd_numelSignal] = psdCalc(spkOI(areaInd), srData, optCoh);
      [unitData.freq, unitData.coh, unitData.phase, unitData.coh_conf, unitData.phase_confU, unitData.phase_confL,...
        unitData.coh_halves_freq, unitData.coh_halves, unitData.coh_conf_halves, unitData.phase_halves,...
        unitData.phase_conf_halves] = phaseCohCalc(areaInterpFilt, spkOI(areaInd), srData, optCoh);
      [unitData.rateadjust_kappa, unitData.rateadjust_kappa_halves] = kappaCalc(unitData.mfr, unitData.mfr_1sthalf,...
        unitData.mfr_2ndhalf, unitData.psd, unitData.psd_halves, unitData.coh, unitData.coh_halves);
      [unitData.psd_halves_freq, unitData.psd_halves, unitData.freq, unitData.psd, unitData.psd_conf,...
        unitData.coh, unitData.phase, unitData.coh_conf, unitData.phase_confU, unitData.phase_confL,...
        unitData.coh_halves_freq, unitData.coh_halves, unitData.coh_conf_halves, unitData.phase_halves,...
        unitData.phase_conf_halves, unitData.rateadjust_kappa, unitData.rateadjust_kappa_halves] = correctFreq(...
        unitData.psd_halves_freq, unitData.psd_halves, freqPSD, unitData.psd, unitData.psd_conf,...
        unitData.freq, unitData.coh, unitData.phase, unitData.coh_conf, unitData.phase_confU,...
        unitData.phase_confL, unitData.coh_halves_freq, unitData.coh_halves, unitData.coh_conf_halves,...
        unitData.phase_halves, unitData.phase_conf_halves, unitData.rateadjust_kappa, unitData.rateadjust_kappa_halves);
    end
    
    % Store temporary data
    saveParfor([entryName '_' unitNamer(u) '_pupil.mat'], unitData);
    fprintf('Finished processing unit %i\n',spkDB_units(u));
  end
  
  % Retrieve temporary data
  fileList = dir([entryName '_*_pupil.mat']);
  if ~isempty(fileList)
    pupil.unitData = struct([]);
    for i = 1:size(fileList,1)
      qPupil = load(fileList(i).name);
      pupil.unitData{end+1} = qPupil.q;
    end
    
    % Save permanent data
    dataString = ['dataStruct.seriesData.' entryName '.popData.pupil.phaseCoh = pupil;'];
    eval(dataString);
    if intermediateSaving
      save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
    end
    delete([entryName '_*_pupil.mat']);
  end
end

if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct