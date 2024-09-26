% Filter the pupil signal and estimate the phase of the filtered signal


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITIALISE PARAMETERS
intermediateSaving = false; % Save everytime a db series is finished being analysed (bad idea since saving is a lengthy process)


%% LOAD DATA
if isfield(dataStruct, 'eyeData') && ~isempty(dataStruct.eyeData)
  fnsData = fieldnames(dataStruct.eyeData);
  if ~isempty(dbEntries) && dbEntries(1) == inf
    dbEntriesEye = 1:numel(fnsData);
  end
  for dbCount = dbEntriesEye % Loop through db entries
    
    % Get eye data
    if isempty(dataStruct.eyeData.(fnsData{dbCount}).pupilArea)
      disp(['No pupil data for ' fnsData{dbCount} '. Skippig to the next db entry...']);
      continue
    end
    eyeDataDB = dataStruct.eyeData.(fnsData{dbCount});
    
    % Interpolate the pupil signal
    dt = mean(eyeDataDB.frameTimes(2:end) - eyeDataDB.frameTimes(1:end-1));
    timesFilt = eyeDataDB.frameTimes(1):dt:eyeDataDB.frameTimes(end);
    assert(numel(timesFilt) == numel(eyeDataDB.frameTimes));
    pupilAreaInterp = interp1(eyeDataDB.frameTimes, eyeDataDB.pupilArea, timesFilt);
    pupilAreaInterp(isnan(pupilAreaInterp)) = 0;
    
    % Filter the pupil signal
    passbandRipple = 0.5;
    stopbandAttenuation = 16;
    freqFiltLP0p001Hz = [0.001 0.0015];
    d1 = designFilterLP(freqFiltLP0p001Hz(1), freqFiltLP0p001Hz(2), passbandRipple, stopbandAttenuation, 1/dt);
    pupilAreaFiltLP0p001Hz = single(filtfilt(d1, pupilAreaInterp));
    
    freqFiltLP0p01Hz = [0.01 0.015];
    d2 = designFilterLP(freqFiltLP0p01Hz(1), freqFiltLP0p01Hz(2), passbandRipple, stopbandAttenuation, 1/dt);
    pupilAreaFiltLP0p01Hz = single(filtfilt(d2, pupilAreaInterp));
    
    freqFiltLP0p1Hz = [0.1 0.15];
    d3 = designFilterLP(freqFiltLP0p1Hz(1), freqFiltLP0p1Hz(2), passbandRipple, stopbandAttenuation, 1/dt);
    pupilAreaFiltLP0p1Hz = single(filtfilt(d3, pupilAreaInterp));
    
    freqFiltBP0p01to0p05Hz = [0.005 0.01 0.03 0.05 0.1];
    d4 = designFilterBP(freqFiltBP0p01to0p05Hz', stopbandAttenuation, passbandRipple, stopbandAttenuation, 1/dt);
    pupilAreaFiltBP0p01to0p05Hz = single(filtfilt(d4, pupilAreaInterp));
    
    freqFiltBP0p1to0p5Hz = [0.05 0.1 0.3 0.5 1];
    d5 = designFilterBP(freqFiltBP0p1to0p5Hz', stopbandAttenuation, passbandRipple, stopbandAttenuation, 1/dt);
    pupilAreaFiltBP0p1to0p5Hz = single(filtfilt(d5, pupilAreaInterp));
    
    freqFiltHP0p01Hz = [0.01 0.00675];
    d6 = designFilterHP(freqFiltHP0p01Hz(1), freqFiltHP0p01Hz(2), passbandRipple, stopbandAttenuation, 1/dt);
    pupilAreaFiltHP0p01Hz = single(filtfilt(d6, pupilAreaInterp));
    
    assert(filtord(d1) == filtord(d2));
    assert(filtord(d1) == filtord(d3));
    assert(filtord(d1) == filtord(d4));
    assert(filtord(d1) == filtord(d5));
    assert(filtord(d1) == filtord(d6));
    filterOrder = filtord(d1);
    
    % Hilbert Transform of the filtered signal
    [~, phaseFiltBP0p01to0p05Hz] = hilbertTransform(pupilAreaFiltBP0p01to0p05Hz);
    phaseFiltBP0p01to0p05Hz = single(phaseFiltBP0p01to0p05Hz);
    [~, phaseFiltBP0p1to0p5Hz] = hilbertTransform(pupilAreaFiltBP0p1to0p5Hz);
    phaseFiltBP0p1to0p5Hz = single(phaseFiltBP0p1to0p5Hz);
    
    % Store variables
    eyeDataDB.pupilAreaFilt.timesFiltStart = timesFilt(1);
    eyeDataDB.pupilAreaFilt.timesFiltStop = timesFilt(end);
    eyeDataDB.pupilAreaFilt.timesFiltStep = dt;
    eyeDataDB.pupilAreaFilt.passbandRipple = passbandRipple;
    eyeDataDB.pupilAreaFilt.stopbandAttenuation = stopbandAttenuation;
    eyeDataDB.pupilAreaFilt.filterOrder = filterOrder;
    eyeDataDB.pupilAreaFilt.freqFiltLP0p001Hz = freqFiltLP0p001Hz;
    eyeDataDB.pupilAreaFilt.freqFiltLP0p01Hz = freqFiltLP0p01Hz;
    eyeDataDB.pupilAreaFilt.freqFiltLP0p1Hz = freqFiltLP0p1Hz;
    eyeDataDB.pupilAreaFilt.freqFiltBP0p01to0p05Hz = freqFiltBP0p01to0p05Hz;
    eyeDataDB.pupilAreaFilt.freqFiltBP0p1to0p5Hz = freqFiltBP0p1to0p5Hz;
    eyeDataDB.pupilAreaFilt.pupilAreaFiltLP0p001Hz = pupilAreaFiltLP0p001Hz;
    eyeDataDB.pupilAreaFilt.pupilAreaFiltLP0p01Hz = pupilAreaFiltLP0p01Hz;
    eyeDataDB.pupilAreaFilt.pupilAreaFiltLP0p1Hz = pupilAreaFiltLP0p1Hz;
    eyeDataDB.pupilAreaFilt.pupilAreaFiltBP0p01to0p05Hz = pupilAreaFiltBP0p01to0p05Hz;
    eyeDataDB.pupilAreaFilt.pupilAreaFiltBP0p1to0p5Hz = pupilAreaFiltBP0p1to0p5Hz;
    eyeDataDB.pupilAreaFilt.phaseFiltBP0p01to0p05Hz = phaseFiltBP0p01to0p05Hz;
    eyeDataDB.pupilAreaFilt.phaseFiltBP0p1to0p5Hz = phaseFiltBP0p1to0p5Hz;
    eyeDataDB.pupilAreaFilt.pupilAreaFiltHP0p01Hz = pupilAreaFiltHP0p01Hz;
    dataStruct.eyeData.(fnsData{dbCount}) = eyeDataDB;
    
    if intermediateSaving % If executed, would produce significant delays
      save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
    end
    
    fprintf('Finished processing db entry %i\n',dbCount);
  end % loop over db entries
  
  
  %% SAVE DATA
  if ~intermediateSaving
    save(dataFile,'dataStruct','-v7.3');
  end
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct