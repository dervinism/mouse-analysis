% Load and divide sparse matrices of spiking data based on the sign of
% their correlation coefficient with respect to the pupil signal.


%% INITIALIZE PARAMETERS
params
lists
slowF = [0.2 1]; % Band-passed slow frequency range
infraSlowF = [0.02 0.2]; % Band-passed infra-slow frequency range
ultraSlowF = 0.01; % Low-pass ultra-slow frequency range


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITIALISE PARAMETERS
intermediateSaving = false;


%% PERFORM COHERENCE ANALSYSES
fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf
  dbEntriesLocal = 1:numel(fnsData);
else
  dbEntriesLocal = dbEntries;
end
for dbCount = dbEntriesLocal % Loop through db entries
  
  % Load the contents of dbStruct
  [dbStruct, ~, ~, entryName, ~, ~, shankIDs,...
    ~, period, ~, ~, srData, ~, ~, ~,...
    MUAsAll, spkDB, spkDB_units] = get_dbStruct(dataStruct, dbCount);
  if isempty(MUAsAll) || ~sum(sum(MUAsAll))
    continue
  end
  
  % Get eye data
  [seriesName, animal] = seriesFromEntry(entryName);
  entryNameEye = [animal '_s' seriesName(1:min([14 numel(seriesName)]))];
  if ~isfield(dataStruct, 'eyeData') || ~isfield(dataStruct.eyeData, entryNameEye) ||...
      isempty(dataStruct.eyeData.(entryNameEye).pupilArea)
    disp(['No pupil data for ' entryName '. Skippig to the next db entry...']);
    continue
  end
  eyeDataDB = dataStruct.eyeData.(entryNameEye);
  
  % Interpolate and filter pupil area data
  commonPeriod = combinePeriods(period, eyeDataDB.period, srData);
  if isempty(commonPeriod)
    continue
  end
  downsamplingFactor = 40;
  [pupilArea, frameTimes, ~, unfilteredArea] = pupilFilt(eyeDataDB, srData/downsamplingFactor,...
    MUAsAll(1:downsamplingFactor:end), 2*(srData/downsamplingFactor), commonPeriod, srData/downsamplingFactor);
  dtInterpVec = diff(frameTimes);
  dtInterp = mean(dtInterpVec(dtInterpVec < 1));
  
  % Bandpass-filter the MUA spiking matrix to obtain the slow signal
  artifact = round(100/dtInterp);
  FpassS = [slowF(1)-(slowF(1)/2) slowF(1)...
    slowF(1)+(slowF(2)-slowF(1))/2 slowF(2) slowF(2)+(slowF(2)/2)];
  d = designFilterBP(FpassS', 6, 0.5, 6, srData);
  spkDBFiltS = filtfilt(d, full(spkDB)')';
  validInds = artifact+1:size(spkDBFiltS,2)-artifact;
  times = validInds./srData;
%   figure; plot(time, spkDB(1, artifact:end-artifact)); hold on
%   plot(time, spkDBFilt(artifact:end-artifact)); hold off
  spkDBFiltS = interp1(times, spkDBFiltS(:,validInds)', frameTimes, 'linear','extrap')';
  if min(size(spkDBFiltS)) == 1
    spkDBFiltS = torow(spkDBFiltS);
  end
  d = designFilterBP(FpassS', 6, 0.5, 6, round(1/dtInterp));
  pupilAreaFiltS = filtfilt(d, pupilArea);
  
  % Bandpass-filter the MUA spiking matrix to obtain the infra-slow signal
  artifact = round(100/dtInterp);
  FpassIS = [infraSlowF(1)-(infraSlowF(1)/2) infraSlowF(1)...
    infraSlowF(1)+(infraSlowF(2)-infraSlowF(1))/2 infraSlowF(2) infraSlowF(2)+(infraSlowF(2)/2)];
  d = designFilterBP(FpassIS', 6, 0.5, 6, srData);
  spkDBFiltIS = filtfilt(d, full(spkDB)')';
  validInds = artifact+1:size(spkDBFiltIS,2)-artifact;
  times = validInds./srData;
%   figure; plot(time, spkDB(1, artifact:end-artifact)); hold on
%   plot(time, spkDBFilt(artifact:end-artifact)); hold off
  spkDBFiltIS = interp1(times, spkDBFiltIS(:,validInds)', frameTimes, 'linear','extrap')';
  if min(size(spkDBFiltIS)) == 1
    spkDBFiltIS = torow(spkDBFiltIS);
  end
  d = designFilterBP(FpassIS', 6, 0.5, 6, round(1/dtInterp));
  pupilAreaFiltIS = filtfilt(d, pupilArea);
  
  % Lowpass-filter the MUA spiking matrix to obtain the ultra-slow signal
  d = designFilterLP(ultraSlowF, ultraSlowF+(ultraSlowF/2), 0.5, 16, srData);
  spkDBFiltUS = filtfilt(d, full(spkDB)')';
%   figure; plot(time, spkDB(1, artifact:end-artifact)); hold on
%   plot(time, spkDBFilt(artifact:end-artifact));
%   d = designFilterLP(ultraSlowF, ultraSlowF+(ultraSlowF/2), 0.5, 16, srData);
%   spkDBFilt = filtfilt(d, full(spkDB(1,:)));
%   plot(time, spkDBFilt(artifact:end-artifact));
  spkDBFiltUS = interp1(times, spkDBFiltUS(:,validInds)', frameTimes, 'linear','extrap')';
  if min(size(spkDBFiltUS)) == 1
    spkDBFiltUS = torow(spkDBFiltUS);
  end
  d = designFilterLP(ultraSlowF, ultraSlowF+(ultraSlowF/2), 0.5, 16, round(1/dtInterp));
  pupilAreaFiltUS = filtfilt(d, pupilArea);
  
  % Correlate down-sampled pupil and MUA spiking signals
  [rPearsonS, pvalPearsonS] = corrMulti(pupilArea, spkDBFiltS, 'Pearson');
  [rSpearmanS, pvalSpearmanS] = corrMulti(pupilArea, spkDBFiltS, 'Spearman');
  [rPearsonIS, pvalPearsonIS] = corrMulti(pupilArea, spkDBFiltIS, 'Pearson');
  [rSpearmanIS, pvalSpearmanIS] = corrMulti(pupilArea, spkDBFiltIS, 'Spearman');
  [rPearsonUS, pvalPearsonUS] = corrMulti(pupilArea, spkDBFiltUS, 'Pearson');
  [rSpearmanUS, pvalSpearmanUS] = corrMulti(pupilArea, spkDBFiltUS, 'Spearman');
  [rPearsonSS, pvalPearsonSS] = corrMulti(pupilAreaFiltS, spkDBFiltS, 'Pearson');
  [rSpearmanSS, pvalSpearmanSS] = corrMulti(pupilAreaFiltS, spkDBFiltS, 'Spearman');
  [rPearsonISIS, pvalPearsonISIS] = corrMulti(pupilAreaFiltIS, spkDBFiltIS, 'Pearson');
  [rSpearmanISIS, pvalSpearmanISIS] = corrMulti(pupilAreaFiltIS, spkDBFiltIS, 'Spearman');
  [rPearsonUSUS, pvalPearsonUSUS] = corrMulti(pupilAreaFiltUS, spkDBFiltUS, 'Pearson');
  [rSpearmanUSUS, pvalSpearmanUSUS] = corrMulti(pupilAreaFiltUS, spkDBFiltUS, 'Spearman');
  
  % Correlate down-sampled pupil and MUA spiking signals over 10, 20, and 30-minute windows
  [rPearsonIS10minWindows, rPearsonIS20minWindows, rPearsonIS30minWindows,...
    pvalPearsonIS10minWindows, pvalPearsonIS20minWindows, pvalPearsonIS30minWindows,...
    rSpearmanIS10minWindows, rSpearmanIS20minWindows, rSpearmanIS30minWindows,...
    pvalSpearmanIS10minWindows, pvalSpearmanIS20minWindows, pvalSpearmanIS30minWindows] = corrWindows(pupilArea, spkDBFiltIS, round(1/dtInterp));
  [rPearsonUS10minWindows, rPearsonUS20minWindows, rPearsonUS30minWindows,...
    pvalPearsonUS10minWindows, pvalPearsonUS20minWindows, pvalPearsonUS30minWindows,...
    rSpearmanUS10minWindows, rSpearmanUS20minWindows, rSpearmanUS30minWindows,...
    pvalSpearmanUS10minWindows, pvalSpearmanUS20minWindows, pvalSpearmanUS30minWindows] = corrWindows(pupilArea, spkDBFiltUS, round(1/dtInterp));
  [rPearsonISIS10minWindows, rPearsonISIS20minWindows, rPearsonISIS30minWindows,...
    pvalPearsonISIS10minWindows, pvalPearsonISIS20minWindows, pvalPearsonISIS30minWindows,...
    rSpearmanISIS10minWindows, rSpearmanISIS20minWindows, rSpearmanISIS30minWindows,...
    pvalSpearmanISIS10minWindows, pvalSpearmanISIS20minWindows, pvalSpearmanISIS30minWindows] = corrWindows(pupilAreaFiltIS, spkDBFiltIS, round(1/dtInterp));
  [rPearsonUSUS10minWindows, rPearsonUSUS20minWindows, rPearsonUSUS30minWindows,...
    pvalPearsonUSUS10minWindows, pvalPearsonUSUS20minWindows, pvalPearsonUSUS30minWindows,...
    rSpearmanUSUS10minWindows, rSpearmanUSUS20minWindows, rSpearmanUSUS30minWindows,...
    pvalSpearmanUSUS10minWindows, pvalSpearmanUSUS20minWindows, pvalSpearmanUSUS30minWindows] = corrWindows(pupilAreaFiltUS, spkDBFiltUS, round(1/dtInterp));
  
  for sh = 1:numel(shankIDs) % Loop through shanks
    fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
    
    % Load the contents of shankStruct
    [~, ~, ~, ~, ~, ~, spk] = get_shankStruct(dbStruct, sh);
    if isempty(spk)
      continue
    end
    
    % Bandpass-filter the unit spiking matrix to obtain the slow signal
    artifact = round((2000/srData)/dtInterp);
    d = designFilterBP(FpassS', 6, 0.5, 6, srData);
    spkFiltS = filtfilt(d, full(spk)')';
    validInds = artifact+1:size(spkFiltS,2)-artifact;
    times = validInds./srData;
    spkFiltS = interp1(times, spkFiltS(:,validInds)', frameTimes, 'linear','extrap')';
    if min(size(spkFiltS)) == 1
      spkFiltS = torow(spkFiltS);
    end
    
    % Bandpass-filter the unit spiking matrix to obtain the infra-slow signal
    artifact = round((2000/srData)/dtInterp);
    d = designFilterBP(FpassIS', 6, 0.5, 6, srData);
    spkFiltIS = filtfilt(d, full(spk)')';
    validInds = artifact+1:size(spkFiltIS,2)-artifact;
    times = validInds./srData;
    spkFiltIS = interp1(times, spkFiltIS(:,validInds)', frameTimes, 'linear','extrap')';
    if min(size(spkFiltIS)) == 1
      spkFiltIS = torow(spkFiltIS);
    end
    
    % Lowpass-filter the unit spiking matrix to obtain the ultra-slow signal
    d = designFilterLP(ultraSlowF, ultraSlowF+(ultraSlowF/2), 0.5, 16, srData);
    spkFiltUS = filtfilt(d, full(spk)')';
    spkFiltUS = interp1(times, spkFiltUS(:,validInds)', frameTimes, 'linear','extrap')';
    if min(size(spkFiltUS)) == 1
      spkFiltUS = torow(spkFiltUS);
    end
    
    % Correlate down-sampled pupil and unit spiking signals
    [rPearsonUnitsS, pvalPearsonUnitsS] = corrMulti(pupilArea, spkFiltS, 'Pearson');
    [rSpearmanUnitsS, pvalSpearmanUnitsS] = corrMulti(pupilArea, spkFiltS, 'Spearman');
    [rPearsonUnitsIS, pvalPearsonUnitsIS] = corrMulti(pupilArea, spkFiltIS, 'Pearson');
    [rSpearmanUnitsIS, pvalSpearmanUnitsIS] = corrMulti(pupilArea, spkFiltIS, 'Spearman');
    [rPearsonUnitsUS, pvalPearsonUnitsUS] = corrMulti(pupilArea, spkFiltUS, 'Pearson');
    [rSpearmanUnitsUS, pvalSpearmanUnitsUS] = corrMulti(pupilArea, spkFiltUS, 'Spearman');
    [rPearsonUnitsSS, pvalPearsonUnitsSS] = corrMulti(pupilAreaFiltS, spkFiltS, 'Pearson');
    [rSpearmanUnitsSS, pvalSpearmanUnitsSS] = corrMulti(pupilAreaFiltS, spkFiltS, 'Spearman');
    [rPearsonUnitsISIS, pvalPearsonUnitsISIS] = corrMulti(pupilAreaFiltIS, spkFiltIS, 'Pearson');
    [rSpearmanUnitsISIS, pvalSpearmanUnitsISIS] = corrMulti(pupilAreaFiltIS, spkFiltIS, 'Spearman');
    [rPearsonUnitsUSUS, pvalPearsonUnitsUSUS] = corrMulti(pupilAreaFiltUS, spkFiltUS, 'Pearson');
    [rSpearmanUnitsUSUS, pvalSpearmanUnitsUSUS] = corrMulti(pupilAreaFiltUS, spkFiltUS, 'Spearman');
    
    % Correlate down-sampled pupil and unit spiking signals over 10, 20, and 30-minute windows
    [rPearsonUnitsIS10minWindows, rPearsonUnitsIS20minWindows, rPearsonUnitsIS30minWindows,...
      pvalPearsonUnitsIS10minWindows, pvalPearsonUnitsIS20minWindows, pvalPearsonUnitsIS30minWindows,...
      rSpearmanUnitsIS10minWindows, rSpearmanUnitsIS20minWindows, rSpearmanUnitsIS30minWindows,...
      pvalSpearmanUnitsIS10minWindows, pvalSpearmanUnitsIS20minWindows, pvalSpearmanUnitsIS30minWindows] = corrWindows(pupilArea, spkFiltIS, round(1/dtInterp));
    [rPearsonUnitsUS10minWindows, rPearsonUnitsUS20minWindows, rPearsonUnitsUS30minWindows,...
      pvalPearsonUnitsUS10minWindows, pvalPearsonUnitsUS20minWindows, pvalPearsonUnitsUS30minWindows,...
      rSpearmanUnitsUS10minWindows, rSpearmanUnitsUS20minWindows, rSpearmanUnitsUS30minWindows,...
      pvalSpearmanUnitsUS10minWindows, pvalSpearmanUnitsUS20minWindows, pvalSpearmanUnitsUS30minWindows] = corrWindows(pupilArea, spkFiltUS, round(1/dtInterp));
    [rPearsonUnitsISIS10minWindows, rPearsonUnitsISIS20minWindows, rPearsonUnitsISIS30minWindows,...
      pvalPearsonUnitsISIS10minWindows, pvalPearsonUnitsISIS20minWindows, pvalPearsonUnitsISIS30minWindows,...
      rSpearmanUnitsISIS10minWindows, rSpearmanUnitsISIS20minWindows, rSpearmanUnitsISIS30minWindows,...
      pvalSpearmanUnitsISIS10minWindows, pvalSpearmanUnitsISIS20minWindows, pvalSpearmanUnitsISIS30minWindows] = corrWindows(pupilAreaFiltIS, spkFiltIS, round(1/dtInterp));
    [rPearsonUnitsUSUS10minWindows, rPearsonUnitsUSUS20minWindows, rPearsonUnitsUSUS30minWindows,...
      pvalPearsonUnitsUSUS10minWindows, pvalPearsonUnitsUSUS20minWindows, pvalPearsonUnitsUSUS30minWindows,...
      rSpearmanUnitsUSUS10minWindows, rSpearmanUnitsUSUS20minWindows, rSpearmanUnitsUSUS30minWindows,...
      pvalSpearmanUnitsUSUS10minWindows, pvalSpearmanUnitsUSUS20minWindows, pvalSpearmanUnitsUSUS30minWindows] = corrWindows(pupilAreaFiltUS, spkFiltUS, round(1/dtInterp));
    
    % Store shank data
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonS = rPearsonUnitsS;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonS = pvalPearsonUnitsS;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanS = rSpearmanUnitsS;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanS = pvalSpearmanUnitsS;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonIS = rPearsonUnitsIS;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonIS = pvalPearsonUnitsIS;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanIS = rSpearmanUnitsIS;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanIS = pvalSpearmanUnitsIS;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonUS = rPearsonUnitsUS;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonUS = pvalPearsonUnitsUS;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanUS = rSpearmanUnitsUS;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanUS = pvalSpearmanUnitsUS;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonSS = rPearsonUnitsSS;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonSS = pvalPearsonUnitsSS;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanSS = rSpearmanUnitsSS;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanSS = pvalSpearmanUnitsSS;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonISIS = rPearsonUnitsISIS;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonISIS = pvalPearsonUnitsISIS;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanISIS = rSpearmanUnitsISIS;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanISIS = pvalSpearmanUnitsISIS;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonUSUS = rPearsonUnitsUSUS;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonUSUS = pvalPearsonUnitsUSUS;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanUSUS = rSpearmanUnitsUSUS;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanUSUS = pvalSpearmanUnitsUSUS;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonIS10minWindows = rPearsonUnitsIS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonIS10minWindows = pvalPearsonUnitsIS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanIS10minWindows = rSpearmanUnitsIS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanIS10minWindows = pvalSpearmanUnitsIS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonUS10minWindows = rPearsonUnitsUS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonUS10minWindows = pvalPearsonUnitsUS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanUS10minWindows = rSpearmanUnitsUS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanUS10minWindows = pvalSpearmanUnitsUS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonISIS10minWindows = rPearsonUnitsISIS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonISIS10minWindows = pvalPearsonUnitsISIS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanISIS10minWindows = rSpearmanUnitsISIS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanISIS10minWindows = pvalSpearmanUnitsISIS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonUSUS10minWindows = rPearsonUnitsUSUS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonUSUS10minWindows = pvalPearsonUnitsUSUS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanUSUS10minWindows = rSpearmanUnitsUSUS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanUSUS10minWindows = pvalSpearmanUnitsUSUS10minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonIS20minWindows = rPearsonUnitsIS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonIS20minWindows = pvalPearsonUnitsIS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanIS20minWindows = rSpearmanUnitsIS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanIS20minWindows = pvalSpearmanUnitsIS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonUS20minWindows = rPearsonUnitsUS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonUS20minWindows = pvalPearsonUnitsUS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanUS20minWindows = rSpearmanUnitsUS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanUS20minWindows = pvalSpearmanUnitsUS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonISIS20minWindows = rPearsonUnitsISIS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonISIS20minWindows = pvalPearsonUnitsISIS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanISIS20minWindows = rSpearmanUnitsISIS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanISIS20minWindows = pvalSpearmanUnitsISIS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonUSUS20minWindows = rPearsonUnitsUSUS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonUSUS20minWindows = pvalPearsonUnitsUSUS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanUSUS20minWindows = rSpearmanUnitsUSUS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanUSUS20minWindows = pvalSpearmanUnitsUSUS20minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonIS30minWindows = rPearsonUnitsIS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonIS30minWindows = pvalPearsonUnitsIS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanIS30minWindows = rSpearmanUnitsIS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanIS30minWindows = pvalSpearmanUnitsIS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonUS30minWindows = rPearsonUnitsUS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonUS30minWindows = pvalPearsonUnitsUS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanUS30minWindows = rSpearmanUnitsUS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanUS30minWindows = pvalSpearmanUnitsUS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonISIS30minWindows = rPearsonUnitsISIS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonISIS30minWindows = pvalPearsonUnitsISIS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanISIS30minWindows = rSpearmanUnitsISIS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanISIS30minWindows = pvalSpearmanUnitsISIS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rPearsonUSUS30minWindows = rPearsonUnitsUSUS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearsonUSUS30minWindows = pvalPearsonUnitsUSUS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanUSUS30minWindows = rSpearmanUnitsUSUS30minWindows;
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanUSUS30minWindows = pvalSpearmanUnitsUSUS30minWindows;
  end
  
  % Store db data
  dbStruct.popData.rPearsonS = rPearsonS;
  dbStruct.popData.pvalPearsonS = pvalPearsonS;
  dbStruct.popData.rSpearmanS = rSpearmanS;
  dbStruct.popData.pvalSpearmanS = pvalSpearmanS;
  dbStruct.popData.rPearsonIS = rPearsonIS;
  dbStruct.popData.pvalPearsonIS = pvalPearsonIS;
  dbStruct.popData.rSpearmanIS = rSpearmanIS;
  dbStruct.popData.pvalSpearmanIS = pvalSpearmanIS;
  dbStruct.popData.rPearsonUS = rPearsonUS;
  dbStruct.popData.pvalPearsonUS = pvalPearsonUS;
  dbStruct.popData.rSpearmanUS = rSpearmanUS;
  dbStruct.popData.pvalSpearmanUS = pvalSpearmanUS;
  dbStruct.popData.rPearsonSS = rPearsonSS;
  dbStruct.popData.pvalPearsonSS = pvalPearsonSS;
  dbStruct.popData.rSpearmanSS = rSpearmanSS;
  dbStruct.popData.pvalSpearmanSS = pvalSpearmanSS;
  dbStruct.popData.rPearsonISIS = rPearsonISIS;
  dbStruct.popData.pvalPearsonISIS = pvalPearsonISIS;
  dbStruct.popData.rSpearmanISIS = rSpearmanISIS;
  dbStruct.popData.pvalSpearmanISIS = pvalSpearmanISIS;
  dbStruct.popData.rPearsonUSUS = rPearsonUSUS;
  dbStruct.popData.pvalPearsonUSUS = pvalPearsonUSUS;
  dbStruct.popData.rSpearmanUSUS = rSpearmanUSUS;
  dbStruct.popData.pvalSpearmanUSUS = pvalSpearmanUSUS;
  dbStruct.popData.rPearsonIS10minWindows = rPearsonIS10minWindows;
  dbStruct.popData.pvalPearsonIS10minWindows = pvalPearsonIS10minWindows;
  dbStruct.popData.rSpearmanIS10minWindows = rSpearmanIS10minWindows;
  dbStruct.popData.pvalSpearmanIS10minWindows = pvalSpearmanIS10minWindows;
  dbStruct.popData.rPearsonUS10minWindows = rPearsonUS10minWindows;
  dbStruct.popData.pvalPearsonUS10minWindows = pvalPearsonUS10minWindows;
  dbStruct.popData.rSpearmanUS10minWindows = rSpearmanUS10minWindows;
  dbStruct.popData.pvalSpearmanUS10minWindows = pvalSpearmanUS10minWindows;
  dbStruct.popData.rPearsonISIS10minWindows = rPearsonISIS10minWindows;
  dbStruct.popData.pvalPearsonISIS10minWindows = pvalPearsonISIS10minWindows;
  dbStruct.popData.rSpearmanISIS10minWindows = rSpearmanISIS10minWindows;
  dbStruct.popData.pvalSpearmanISIS10minWindows = pvalSpearmanISIS10minWindows;
  dbStruct.popData.rPearsonUSUS10minWindows = rPearsonUSUS10minWindows;
  dbStruct.popData.pvalPearsonUSUS10minWindows = pvalPearsonUSUS10minWindows;
  dbStruct.popData.rSpearmanUSUS10minWindows = rSpearmanUSUS10minWindows;
  dbStruct.popData.pvalSpearmanUSUS10minWindows = pvalSpearmanUSUS10minWindows;
  dbStruct.popData.rPearsonIS20minWindows = rPearsonIS20minWindows;
  dbStruct.popData.pvalPearsonIS20minWindows = pvalPearsonIS20minWindows;
  dbStruct.popData.rSpearmanIS20minWindows = rSpearmanIS20minWindows;
  dbStruct.popData.pvalSpearmanIS20minWindows = pvalSpearmanIS20minWindows;
  dbStruct.popData.rPearsonUS20minWindows = rPearsonUS20minWindows;
  dbStruct.popData.pvalPearsonUS20minWindows = pvalPearsonUS20minWindows;
  dbStruct.popData.rSpearmanUS20minWindows = rSpearmanUS20minWindows;
  dbStruct.popData.pvalSpearmanUS20minWindows = pvalSpearmanUS20minWindows;
  dbStruct.popData.rPearsonISIS20minWindows = rPearsonISIS20minWindows;
  dbStruct.popData.pvalPearsonISIS20minWindows = pvalPearsonISIS20minWindows;
  dbStruct.popData.rSpearmanISIS20minWindows = rSpearmanISIS20minWindows;
  dbStruct.popData.pvalSpearmanISIS20minWindows = pvalSpearmanISIS20minWindows;
  dbStruct.popData.rPearsonUSUS20minWindows = rPearsonUSUS20minWindows;
  dbStruct.popData.pvalPearsonUSUS20minWindows = pvalPearsonUSUS20minWindows;
  dbStruct.popData.rSpearmanUSUS20minWindows = rSpearmanUSUS20minWindows;
  dbStruct.popData.pvalSpearmanUSUS20minWindows = pvalSpearmanUSUS20minWindows;
  dbStruct.popData.rPearsonIS30minWindows = rPearsonIS30minWindows;
  dbStruct.popData.pvalPearsonIS30minWindows = pvalPearsonIS30minWindows;
  dbStruct.popData.rSpearmanIS30minWindows = rSpearmanIS30minWindows;
  dbStruct.popData.pvalSpearmanIS30minWindows = pvalSpearmanIS30minWindows;
  dbStruct.popData.rPearsonUS30minWindows = rPearsonUS30minWindows;
  dbStruct.popData.pvalPearsonUS30minWindows = pvalPearsonUS30minWindows;
  dbStruct.popData.rSpearmanUS30minWindows = rSpearmanUS30minWindows;
  dbStruct.popData.pvalSpearmanUS30minWindows = pvalSpearmanUS30minWindows;
  dbStruct.popData.rPearsonISIS30minWindows = rPearsonISIS30minWindows;
  dbStruct.popData.pvalPearsonISIS30minWindows = pvalPearsonISIS30minWindows;
  dbStruct.popData.rSpearmanISIS30minWindows = rSpearmanISIS30minWindows;
  dbStruct.popData.pvalSpearmanISIS30minWindows = pvalSpearmanISIS30minWindows;
  dbStruct.popData.rPearsonUSUS30minWindows = rPearsonUSUS30minWindows;
  dbStruct.popData.pvalPearsonUSUS30minWindows = pvalPearsonUSUS30minWindows;
  dbStruct.popData.rSpearmanUSUS30minWindows = rSpearmanUSUS30minWindows;
  dbStruct.popData.pvalSpearmanUSUS30minWindows = pvalSpearmanUSUS30minWindows;
  
  dataStruct.seriesData.(fnsData{dbCount}) = dbStruct;
  
  if intermediateSaving
    save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
  end

  fprintf('Finished processing db entry %i\n',dbCount);
end % loop over db entries


%% SAVE DATA
if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct



function [rPearson10minWindows, rPearson20minWindows, rPearson30minWindows,...
  pvalPearson10minWindows, pvalPearson20minWindows, pvalPearson30minWindows,...
  rSpearman10minWindows, rSpearman20minWindows, rSpearman30minWindows,...
  pvalSpearman10minWindows, pvalSpearman20minWindows, pvalSpearman30minWindows] = corrWindows(pupilArea, spkDBFilt, sr)

correlationWindowSize10min = floor(10*60*sr);
correlationWindowSize20min = floor(20*60*sr);
correlationWindowSize30min = floor(30*60*sr);
n10minWindows = numel(pupilArea)/correlationWindowSize10min;
if n10minWindows > floor(n10minWindows)+0.95
  n10minWindows = ceil(n10minWindows);
else
  n10minWindows = floor(n10minWindows);
end
n20minWindows = numel(pupilArea)/correlationWindowSize20min;
if n20minWindows > floor(n20minWindows)+0.95
  n20minWindows = ceil(n20minWindows);
else
  n20minWindows = floor(n20minWindows);
end
n30minWindows = numel(pupilArea)/correlationWindowSize30min;
if n30minWindows > floor(n30minWindows)+0.95
  n30minWindows = ceil(n30minWindows);
else
  n30minWindows = floor(n30minWindows);
end
rPearson10minWindows = cell(n10minWindows,1);
rPearson20minWindows = cell(n20minWindows,1);
rPearson30minWindows = cell(n30minWindows,1);
pvalPearson10minWindows = cell(n10minWindows,1);
pvalPearson20minWindows = cell(n20minWindows,1);
pvalPearson30minWindows = cell(n30minWindows,1);
rSpearman10minWindows = cell(n10minWindows,1);
rSpearman20minWindows = cell(n20minWindows,1);
rSpearman30minWindows = cell(n30minWindows,1);
pvalSpearman10minWindows = cell(n10minWindows,1);
pvalSpearman20minWindows = cell(n20minWindows,1);
pvalSpearman30minWindows = cell(n30minWindows,1);
for iWindow = 1:n10minWindows
  [rPearson10minWindows{iWindow}, pvalPearson10minWindows{iWindow}] =...
    corrMulti(pupilArea((iWindow-1)*correlationWindowSize10min+1:min([iWindow*correlationWindowSize10min numel(pupilArea)])),...
    spkDBFilt(:,(iWindow-1)*correlationWindowSize10min+1:min([iWindow*correlationWindowSize10min size(spkDBFilt,2)])), 'Pearson');
  [rSpearman10minWindows{iWindow}, pvalSpearman10minWindows{iWindow}] =...
    corrMulti(pupilArea((iWindow-1)*correlationWindowSize10min+1:min([iWindow*correlationWindowSize10min numel(pupilArea)])),...
    spkDBFilt(:,(iWindow-1)*correlationWindowSize10min+1:min([iWindow*correlationWindowSize10min size(spkDBFilt,2)])), 'Spearman');
end
for iWindow = 1:n20minWindows
  [rPearson20minWindows{iWindow}, pvalPearson20minWindows{iWindow}] =...
    corrMulti(pupilArea((iWindow-1)*correlationWindowSize20min+1:min([iWindow*correlationWindowSize20min numel(pupilArea)])),...
    spkDBFilt(:,(iWindow-1)*correlationWindowSize20min+1:min([iWindow*correlationWindowSize20min size(spkDBFilt,2)])), 'Pearson');
  [rSpearman20minWindows{iWindow}, pvalSpearman20minWindows{iWindow}] =...
    corrMulti(pupilArea((iWindow-1)*correlationWindowSize20min+1:min([iWindow*correlationWindowSize20min numel(pupilArea)])),...
    spkDBFilt(:,(iWindow-1)*correlationWindowSize20min+1:min([iWindow*correlationWindowSize20min size(spkDBFilt,2)])), 'Spearman');
end
for iWindow = 1:n30minWindows
  [rPearson30minWindows{iWindow}, pvalPearson30minWindows{iWindow}] =...
    corrMulti(pupilArea((iWindow-1)*correlationWindowSize30min+1:min([iWindow*correlationWindowSize30min numel(pupilArea)])),...
    spkDBFilt(:,(iWindow-1)*correlationWindowSize30min+1:min([iWindow*correlationWindowSize30min size(spkDBFilt,2)])), 'Pearson');
  [rSpearman30minWindows{iWindow}, pvalSpearman30minWindows{iWindow}] =...
    corrMulti(pupilArea((iWindow-1)*correlationWindowSize30min+1:min([iWindow*correlationWindowSize30min numel(pupilArea)])),...
    spkDBFilt(:,(iWindow-1)*correlationWindowSize30min+1:min([iWindow*correlationWindowSize30min size(spkDBFilt,2)])), 'Spearman');
end
end