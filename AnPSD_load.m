% A part of the old AnPSD script adapted for loading and storing full
% spiking data in sparse matrices.


%% INITIALISE PARAMETERS
syncFuncCalc = false; % Synchronise Neuropixels recordings
overwrite = false; % Creates a clean dataStruct variable
intermediateSaving = false; % Save everytime a db series is finished being analysed (bad idea since saving is a lengthy process)

optCoh.maxFreq = maxFreq;
optCoh.winfactor = winfactor;
optCoh.freqfactor = freqfactor;
optCoh.tapers = tapers;
optCoh.monotoneFreq = true;
if optCoh.maxFreq > srData/2
  error('opt.maxFreq > Nyquist?!')
end

unitHeader = {'unit_id'; 'type'; 'probe_channel'; 'horizontal_position';...
  'vertical_position'; 'isi_violations'; 'isolation_distance';...
  'anterior_posterior_ccf_coordinate'; 'dorsal_ventral_ccf_coordinate';...
  'left_right_ccf_coordinate'}; % Unit metadata header

samplingParams.srData = srData;
samplingParams.srRecording = srRecording;
samplingParams.srRecordingLFP = srRecordingLFP;


%% LOAD PRE-PROCESSED DATA IF EXISTS
if exist('dataFile', 'var') && exist(dataFile, 'file')
  load(dataFile);
end


%% ESTIMATE SYNCHRONISATION FUNCTION FOR UOL DUAL NEUROPIXELS RECORDINGS
syncFuncWrap



%% CREATE THE DATA STRUCTURE
if overwrite && exist('dataStruct', 'var')
  seriesData = dataStruct.seriesData;
  for dbCount = 1:numel(fnsData)
    if ismember(dbCount, dbEntries)
      seriesData.(fnsData{dbCount}) = struct();
    else
      dbStruct = seriesData.(fnsData{dbCount});
      dbEntry = db(dbCount).entryName;
      dbStruct.db = db;
      dbStruct.dbSeries = db(dbCount);
      dbStruct.conf.syncFunc = syncFunc{dbCount};
      seriesData = setfield(seriesData,dbEntry,dbStruct);
      dataStruct.seriesData = seriesData;
    end
  end
elseif overwrite || (~overwrite && ~exist('dataStruct', 'var'))
  seriesData = struct();
  dataFile = [db(1).animal '.mat'];
else
  seriesData = dataStruct.seriesData;
end


%% LOAD SPIKING DATA
if ~isempty(dbEntries) && dbEntries(1) == inf % Here you can chose to execute only certain db entries
  dbEntriesLocal = 1:numel(fnsData);
else
  dbEntriesLocal = dbEntries;
end
for dbCount = dbEntriesLocal % Loop through the db entries
  shankData = struct();
  shankIDs = db(dbCount).shank;
  
  if strcmp(repository, 'uol')
    if iscell(db(dbCount).series) % This section is executed for single probe multiple area recordings
      dataDir = [topDir filesep db(dbCount).animal filesep num2str(db(dbCount).series{1})];
    else % This section is executed for single probe recordings in a single brain area
      dataDir = [topDir filesep db(dbCount).animal filesep num2str(db(dbCount).series)];
    end
  else
    dataDir = [topDir filesep db(dbCount).animal];
  end
  baseFilename = [dataDir filesep db(dbCount).basefilename];
  
  % Get the probe setup
  if strcmp(repository, 'uol')
    fileList = dir([dataDir filesep 'forPRB*.mat']);
    if isempty(fileList)
      error('The data processing folder is missing the forPRB*.mat file containing probe configuration.');
    elseif numel(fileList) > 1
      error('There is more than one forPRB*.mat file in the data processing folder. Please remove conflicting files.');
    else
      probeSetup = load([dataDir filesep fileList.name]);
      probe = fileList.name(8:end-4);
      if strcmpi(probe,'mecP3opt3')
        probe = 'imecP3opt3';
      end
    end
  else
    probe = 'Neuropixels';
  end
  
  % Load spikes on all shanks
  if strcmp(repository, 'uol') % UOL MUA data here
    load_opt.selectedUnits = 'all';
    load_opt.templates = true;
    MUAsAll = loadAsMUA_noResClu(dataDir, baseFilename, [dataDir filesep fileList.name],...
      shankIDs, numel(db(dbCount).LFPch)/db(dbCount).shanks,...
      db(dbCount).chOI, srRecording, 1/srData, load_opt);
    MUAsAll = syncAPs(MUAsAll, srData, syncFunc{dbCount});
    [spkDB, spkDB_units, chanMap] = loadAsRasterSparse(dataDir, baseFilename, [dataDir filesep fileList.name], shankIDs,...
      numel(db(dbCount).LFPch)/db(dbCount).shanks, db(dbCount).chOI, srRecording, 1/srData,...
      size(MUAsAll, 2), load_opt);
    spkDB = syncAPs(spkDB, srData, syncFunc{dbCount});
    assert(sum(sum(MUAsAll)) == sum(sum(spkDB)));
    if size(spkDB_units,2) > 1
      spkDB_units = spkDB_units(:,2);
    end
    if isempty(spkDB_units)
      muaMetadata = [];
    else
      unitChan = chanMap(:,2);
      qualitySeries = [dataDir filesep db(dbCount).basefilename '.qua.1.mat'];
      [~, ~, isolationDist, isiViolations] = qualityTest(qualitySeries, spkDB_units, -inf, inf, false);
      try
        muaMetadata = [spkDB_units spkDB_units unitChan probeSetup.xcoords(unitChan)' probeSetup.ycoords(unitChan)'...
          isiViolations isolationDist];
      catch
        muaMetadata = [spkDB_units spkDB_units unitChan probeSetup.xcoords(unitChan) probeSetup.ycoords(unitChan)...
          isiViolations isolationDist];
      end
    end
  else % Allen data here (both MUAs and units)
    [MUAsAll, spkDB, spkDB_units, muaMetadata, ~, spk, units, unitMetadata] = loadSpikes_allensdk(headerUnits,...
      unitMetadataOI, headerSpikes, spikeData, series{dbCount}, refractCont, cluDist, minSpikeCount, srData);
  end
  
  for sh = shankIDs % Loop through shanks. For Neuropixels probes there is only one shank obviously
    fprintf('%s shank %d -------------------\n', db(dbCount).entryName, sh);
    
    % Get the spikes for the shank of interest
    MUAs = MUAsAll(shankIDs(sh), :);  % population rate on this shank/tetrode
    
    % Get unit spiking data for the shank of interest. We have it already if we are dealing with allensdk data
    if strcmp(repository, 'uol') % UOL unit data here (Allen data was already fully loaded earlier)
      if exist([baseFilename, '.clu.', num2str(sh)], 'file')
        clu = load([baseFilename, '.clu.', num2str(sh)]);
      else
        clu = resCluFromKilosort(dataDir, sh, numel(db(dbCount).LFPch)/db(dbCount).shanks,...
          db(dbCount).chOI, [dataDir filesep fileList.name]);
      end
      clu = clu(2:end); % first entry is the total number of units (remove it)
      if isempty(clu)
        units = [];
        unitMetadata = [];
        spk = [];
        MUAs = [];
      else
        if exist([dataDir filesep 'waveforms.mat'],'file')
          load([dataDir filesep 'waveforms.mat'],'chanMap');
          unitsOI = chanMap(logical(sum(chanMap(:,2) == db(dbCount).chOI,2)),3)';
          clu = clu(logical(sum(clu == unitsOI,2)));
        end
        
        if isempty(clu)
          units = [];
          unitMetadata = [];
          spk = [];
          MUAs = [];
        else
          units = unique(clu);
          units = units(units > 1); % remove noise and MUA spikes
          if isempty(units)
            units = [];
            unitMetadata = [];
            spk = [];
          else
            
            % Remove units with a small spike count (< 300 spikes?)
            badUnits = [];
            for u = 1:numel(units)
              if sum(clu == units(u)) < minSpikeCount
                badUnits(end+1) = units(u); %#ok<*SAGROW>
              end
            end
            if ~isempty(badUnits)
              units = setdiff(units, badUnits);
            end
            clear badUnits clu u
            
            % Load the raster matrix containing good units only
            if isempty(units)
              units = [];
              unitMetadata = [];
              spk = [];
            else
              load_opt.selectedUnits = units;
              load_opt.templates = false;
              [spk, units] = loadAsRasterSparse(dataDir, baseFilename, [dataDir filesep fileList.name], sh,...
                numel(db(dbCount).LFPch)/db(dbCount).shanks, db(dbCount).chOI, srRecording, 1/srData,...
                size(MUAsAll, 2), load_opt);
              if isempty(units)
                units = [];
                unitMetadata = [];
                spk = [];
              else
                spk = syncAPs(spk, srData, syncFunc{dbCount});
                units = units(:,2);
                unitChan = zeros(numel(units),1);
                try
                  load([baseFilename '.chanMap.mat']);
                catch
                  try
                    load([baseFilename '.imec.ap_CAR' '.chanMap.mat']);
                  catch
                    try
                      load([dataDir filesep 'waveforms.mat'], 'chanMap');
                    catch
                      chanMap = [];
                    end
                  end
                end
                if ~isempty(chanMap)
                  for u = 1:numel(units)
                    i = find(chanMap(:,1) == units(u));
                    if ~isempty(i)
                      unitChan(u) = chanMap(i,2);
                    end
                  end
                end
                
                % Get unit quality data
                qualitySeries = [dataDir filesep db(dbCount).basefilename '.qua.1.mat'];
                [~, ~, isolationDist, isiViolations] = qualityTest(qualitySeries, units, -inf, inf, false);
                try
                  unitMetadata = [units units unitChan probeSetup.xcoords(unitChan)' probeSetup.ycoords(unitChan)'...
                    isiViolations isolationDist];
                catch
                  unitMetadata = [units units unitChan probeSetup.xcoords(unitChan) probeSetup.ycoords(unitChan)...
                    isiViolations isolationDist];
                end
              end
            end
          end
        end
      end
    end
    
    % Check if units and multiunits were loaded correctly
    %unitInds = ismember(units, spkDB_units(:,end));
    %assert(sum(unitInds) == numel(units));
    %disp([num2str(size(spk,2)) ' ' num2str(size(spkDB,2))])
    
    % Save loaded data
    if ~overwrite && exist('dataStruct', 'var')
      dbEntry = db(dbCount).entryName;
      dataStruct.seriesData.(dbEntry).shankData.(['shank' num2str(sh)]).shankID = sh;
      dataStruct.seriesData.(dbEntry).shankData.(['shank' num2str(sh)]).MUAs = MUAs;
      dataStruct.seriesData.(dbEntry).shankData.(['shank' num2str(sh)]).spk = spk;
      dataStruct.seriesData.(dbEntry).shankData.(['shank' num2str(sh)]).units = units;
      dataStruct.seriesData.(dbEntry).shankData.(['shank' num2str(sh)]).unitHeader = {unitHeader};
      dataStruct.seriesData.(dbEntry).shankData.(['shank' num2str(sh)]).unitMetadata = unitMetadata;
      shankStruct = struct('shankID',sh, 'MUAs',MUAs, 'spk',spk,...
        'units',units, 'unitHeader',{unitHeader}, 'unitMetadata',unitMetadata);
      shankEntry = ['shank' num2str(sh)];
      shankData = setfield(shankData,shankEntry,shankStruct); %#ok<*SFLD>
      dbStruct.shankData = shankData;
    else
      overwrite = true;
      shankStruct = struct('shankID',sh, 'MUAs',MUAs, 'spk',spk,...
        'units',units, 'unitHeader',{unitHeader}, 'unitMetadata',unitMetadata);
      shankEntry = ['shank' num2str(sh)];
      shankData = setfield(shankData,shankEntry,shankStruct); %#ok<*SFLD>
      dbStruct.shankData = shankData;
    end
  end
  dbStruct.popData = struct('MUAsAll',MUAsAll, 'spkDB',spkDB, 'spkDB_units',spkDB_units,...
    'muaHeader',{unitHeader}, 'muaMetadata',muaMetadata);
  dbStruct.io = struct('topDir',topDir, 'dataDir',dataDir, 'baseFilename',baseFilename);
  dbStruct.conf = struct('optCoh',optCoh, 'exclRad',exclRad, 'FOI',FOI,...
    'samplingParams',samplingParams, 'probe',probe, 'syncFunc',syncFunc{dbCount});
  dbStruct.db = db;
  dbStruct.dbSeries = db(dbCount);
  dbEntry = db(dbCount).entryName;
  seriesData = setfield(seriesData,dbEntry,dbStruct);
  if overwrite || (~overwrite && ~exist('dataStruct', 'var'))
    overwrite = true;
    dataStruct.seriesData = seriesData;
  else
    for iShank = 1:numel(dataStruct.seriesData.(dbEntry).shankData)
      assert(dataStruct.seriesData.(dbEntry).shankData.(['shank' num2str(iShank)]).shankID ==...
        dbStruct.shankData.(['shank' num2str(iShank)]).shankID);
      dataStruct.seriesData.(dbEntry).shankData.(['shank' num2str(iShank)]).MUAs =...
        dbStruct.shankData.(['shank' num2str(iShank)]).MUAs;
      dataStruct.seriesData.(dbEntry).shankData.(['shank' num2str(iShank)]).spk =...
        dbStruct.shankData.(['shank' num2str(iShank)]).spk;
      dataStruct.seriesData.(dbEntry).shankData.(['shank' num2str(iShank)]).units =...
        dbStruct.shankData.(['shank' num2str(iShank)]).units;
      dataStruct.seriesData.(dbEntry).shankData.(['shank' num2str(iShank)]).unitHeader =...
        dbStruct.shankData.(['shank' num2str(iShank)]).unitHeader;
      dataStruct.seriesData.(dbEntry).shankData.(['shank' num2str(iShank)]).unitMetadata =...
        dbStruct.shankData.(['shank' num2str(iShank)]).unitMetadata;
    end
    dataStruct.seriesData.(dbEntry).popData.MUAsAll = dbStruct.popData.MUAsAll;
    dataStruct.seriesData.(dbEntry).popData.spkDB = dbStruct.popData.spkDB;
    dataStruct.seriesData.(dbEntry).popData.spkDB_units = dbStruct.popData.spkDB_units;
    dataStruct.seriesData.(dbEntry).popData.muaHeader = dbStruct.popData.muaHeader;
    dataStruct.seriesData.(dbEntry).popData.muaMetadata = dbStruct.popData.muaMetadata;
  end
  if intermediateSaving % If executed, would produce significant delays
    save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
  end
  
  % Test for consistency between series' lengths. Ideally there shouldn't
  % be any errors in this section suggestive of incorrect files being
  % loaded.
  sInd = strfind(dbEntry,'s');
  if sInd+14 <= numel(dbEntry)
    seriesName = dbEntry(sInd+1:sInd+14);
    seriesLength = size(spkDB,2);
    if dbCount == dbEntriesLocal(1)
      prevSeriesName = seriesName;
      prevSeriesLength = seriesLength;
    end
    if prevSeriesLength > 0 && seriesLength > 0
      if strcmpi(prevSeriesName, seriesName)
        assert(abs(prevSeriesLength - seriesLength) < 360000,...
          'Durations of data vectors relating to the same data series are of unequal lengths - a cause for concern');
      else
        assert(abs(prevSeriesLength - seriesLength) > 3600,...
          'Durations of data vectors relating to different data series are of equal lengths - a cause for concern');
      end
    end
    prevSeriesName = seriesName;
    prevSeriesLength = seriesLength;
  end
end

% Save loaded spiking data if it hasn't been saved already (best to save
% only once at the end of the script execution)
if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct
