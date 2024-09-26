% Load waveform data.


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITIALISE PARAMETERS
intermediateSaving = false; % Save everytime a db series is finished being analysed (bad idea since saving is a lengthy process)


%% PERFORM COHERENCE ANALSYSES
fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf % Here you can chose to execute only certain db entries
  dbEntriesLocal = 1:numel(fnsData);
else
  dbEntriesLocal = dbEntries;
end
for dbCount = dbEntriesLocal % Loop through db entries
  
  % Load the contents of dbStruct
  [dbStruct, repository, ~, ~, ~, ~, shankIDs,...
    ~, ~, ~, ~, ~, ~, ~, ~,...
    MUAsAll] = get_dbStruct(dataStruct, dbCount);
  if isempty(MUAsAll) || ~sum(sum(MUAsAll))
    continue
  end
  
  % Load waveforms
  if strcmpi(repository, 'uol')
    if exist([dbStruct.io.dataDir filesep 'waveforms.mat'], 'file')
      waveforms = load([dbStruct.io.dataDir filesep 'waveforms.mat']);
      assert(numel(waveforms.cluIDs) == numel(waveforms.amplitudes));
    else
      continue
    end
  elseif strcmpi(repository, 'allensdk')
    waveforms = load([dbStruct.io.baseFilename '.mat'],...
      'waveformsMaxCh','waveformsTimes','waveformsUnits');
  end
  
  for sh = 1:numel(shankIDs) % Loop through shanks
    fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
    
    % Load the contents of shankStruct
    [~, ~, units] = get_shankStruct(dbStruct, sh);
    if isempty(units)
      continue
    end
    
    if strcmpi(repository, 'uol')
      [~, uInds] = intersect(waveforms.cluIDs, units);
      waveformsShank.amplitudes = waveforms.amplitudes(uInds,:);
      waveformsShank.chanMap = waveforms.chanMap(uInds,:);
      waveformsShank.cluIDs = waveforms.cluIDs(uInds,:);
      waveformsShank.datFile = waveforms.datFile;
      waveformsShank.maxChan = waveforms.maxChan(uInds,:);
      waveformsShank.maxWaveforms = waveforms.maxWaveforms(uInds,:);
      waveformsShank.spikeCentreIndex = waveforms.spikeCentreIndex;
      waveformsShank.waveforms = waveforms.waveforms(uInds,:,:);
    elseif strcmpi(repository, 'allensdk')
      [~, uInds] = intersect(waveforms.waveformsUnits, units);
      waveformsShank.cluIDs = waveforms.waveformsUnits(uInds);
      waveformsShank.maxWaveforms = waveforms.waveformsMaxCh(uInds,:);
      waveformsShank.waveformsTimes = waveforms.waveformsTimes(uInds,:);
    end
    
    dbStruct.shankData.(['shank' num2str(sh)]).waveformData = waveformsShank;
    
    fprintf('Finished processing shank %i\n',sh);
  end % loop over shanks
  
  dataStruct.seriesData.(fnsData{dbCount}) = dbStruct;
  
  if intermediateSaving % If executed, would produce significant delays
    save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
  end

  fprintf('Finished processing db entry %i\n',dbCount);
end % loop over db entries


%% SAVE DATA
if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct