% Load and divide sparse matrices of spiking data based on the sign of
% their correlation coefficient (or phase) with respect to the pupil signal.


%% INITIALIZE PARAMETERS
params
lists
if excludeRunning % Applies only to Allen data when excluding running periods, because periods without running are too short
  fRef = 0.3; % Reference frequency
  fExtrapThr = 0.5; % Extrapolate to fRef if the lowest frequency falls below fExtrapThr
else
  fRef = 0.03;
  fExtrapThr = 0.05;
end
dr = 0.2; % down-sampled sampling rate (Hz)
alpha = 0.05; % signficance level
draw = false; % Produce raster graphs


%% LOAD PRE-PROCESSED DATA
if ~exist('dataStruct','var')
  load(dataFile);
end


%% INITIALISE PARAMETERS
intermediateSaving = false; % Save everytime a db series is finished being analysed (bad idea since saving is a lengthy process)


%% LOAD DATA
fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf % Here you can chose to execute only certain db entries
  dbEntriesLocal = 1:numel(fnsData);
else
  dbEntriesLocal = dbEntries;
end
for dbCount = dbEntriesLocal % Loop through db entries
  
  % Load the contents of dbStruct
  [dbStruct, repository, ~, entryName, ~, ~, shankIDs,...
    ~, period, ~, ~, srData, ~, ~, ~,...
    MUAsAll, spkDB, spkDB_units] = get_dbStruct(dataStruct, dbCount);
  
  % Vipe out existing data
  dbStruct_positive = [];
  dbStruct_negative = [];
  
  for sh = 1:numel(shankIDs)
    dbStruct.shankData.(['shank' num2str(sh)]).rPearson = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearman = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rPearson10percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson10percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearman10percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman10percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rPearson25percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson25percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearman25percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman25percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rPearson33percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson33percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearman33percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman33percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rPearson50percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson50percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearman50percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman50percent = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rPearson10minWindows = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson10minWindows = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearman10minWindows = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman10minWindows = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rPearson20minWindows = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson20minWindows = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearman20minWindows = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman20minWindows = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rPearson30minWindows = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson30minWindows = [];
    dbStruct.shankData.(['shank' num2str(sh)]).rSpearman30minWindows = [];
    dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman30minWindows = [];
  end
  
  dbStruct.popData.meanPupilArea = [];
  dbStruct.popData.meanPupilArea10minWindows = [];
  dbStruct.popData.meanPupilArea20minWindows = [];
  dbStruct.popData.meanPupilArea30minWindows = [];
  dbStruct.popData.rPearson = [];
  dbStruct.popData.pvalPearson = [];
  dbStruct.popData.rSpearman = [];
  dbStruct.popData.pvalSpearman = [];
  dbStruct.popData.rPearson10percent = [];
  dbStruct.popData.pvalPearson10percent = [];
  dbStruct.popData.rSpearman10percent = [];
  dbStruct.popData.pvalSpearman10percent = [];
  dbStruct.popData.rPearson25percent = [];
  dbStruct.popData.pvalPearson25percent = [];
  dbStruct.popData.rSpearman25percent = [];
  dbStruct.popData.pvalSpearman25percent = [];
  dbStruct.popData.rPearson33percent = [];
  dbStruct.popData.pvalPearson33percent = [];
  dbStruct.popData.rSpearman33percent = [];
  dbStruct.popData.pvalSpearman33percent = [];
  dbStruct.popData.rPearson50percent = [];
  dbStruct.popData.pvalPearson50percent = [];
  dbStruct.popData.rSpearman50percent = [];
  dbStruct.popData.pvalSpearman50percent = [];
  dbStruct.popData.rPearson10minWindows = [];
  dbStruct.popData.pvalPearson10minWindows = [];
  dbStruct.popData.rSpearman10minWindows = [];
  dbStruct.popData.pvalSpearman10minWindows = [];
  dbStruct.popData.rPearson20minWindows = [];
  dbStruct.popData.pvalPearson20minWindows = [];
  dbStruct.popData.rSpearman20minWindows = [];
  dbStruct.popData.pvalSpearman20minWindows = [];
  dbStruct.popData.rPearson30minWindows = [];
  dbStruct.popData.pvalPearson30minWindows = [];
  dbStruct.popData.rSpearman30minWindows = [];
  dbStruct.popData.pvalSpearman30minWindows = [];
  
  dataStruct.seriesData_positive.(entryName) = dbStruct_positive;
  dataStruct.seriesData_negative.(entryName) = dbStruct_negative;
  dataStruct.seriesData.(fnsData{dbCount}) = dbStruct;
  
  % Check if any spiking data exists
  if isempty(MUAsAll) || ~sum(sum(MUAsAll))
    continue
  end
  
  % Check if eye data exists
  [seriesName, animal] = seriesFromEntry(entryName);
  [breakClause, pupilCorrCond] = series2condition(awake, anaesthesia, seriesName);
  if breakClause
    continue
  end
  entryNameEye = [animal '_s' seriesName(1:min([14 numel(seriesName)]))];
  if ~isfield(dataStruct, 'eyeData') || ~isfield(dataStruct.eyeData, entryNameEye) ||...
      isempty(dataStruct.eyeData.(entryNameEye).pupilArea)
    disp(['No pupil data for ' entryName '. Skippig to the next db entry...']);
    continue
  end

  % Determine mode boundaries
  areaCode = determineAreaFromSeries(seriesName);
  areaCode = areaCode(1);
  if strcmpi(repository, 'uol')
    modeBoundaries = phaseHistoBinCentres(modesUOL{areaCode});
  elseif strcmpi(repository, 'allensdk')
    modeBoundaries = phaseHistoBinCentres(modesAllensdk{areaCode});
  end
  
  if strcmpi(pupilCorrCond, 'coh')
    % Get the mua phase
    phaseCoh = dbStruct.popData.pupil.phaseCoh.unitData;
    if isempty(phaseCoh)
      continue
    end
    muaPhase = zeros(1,numel(phaseCoh));
    for u = 1:numel(phaseCoh)
      freq = phaseCoh{u}.freq;
      freqInterp = unique([freq fRef fExtrapThr]);
      if ~isfield(phaseCoh{u}, 'phase') || sum(isnan(phaseCoh{u}.phase)) == numel(phaseCoh{u}.phase)
        muaPhase(u) = NaN;
        continue
      end
      if freq(1) <= fExtrapThr
        phase = interp1(freq, phaseCoh{u}.phase, freqInterp, 'linear', 'extrap');
      else
        phase = interp1(freq, phaseCoh{u}.phase, freqInterp, 'linear');
      end
      [~, iFreq] = min(abs(freqInterp - fRef));
      muaPhase(u) = phase(iFreq);
    end
    positivePhase = false(size(muaPhase));
    modeBoundaries = recentrePhase(modeBoundaries, modeBoundaries(1));
    muaPhase = recentrePhase(muaPhase, modeBoundaries(1));
    positivePhase(muaPhase > modeBoundaries(end) & muaPhase <= modeBoundaries(2)) = true;
    negativePhase = false(size(muaPhase));
    modeBoundaries = recentrePhase(modeBoundaries, modeBoundaries(3));
    muaPhase = recentrePhase(muaPhase, modeBoundaries(3));
    negativePhase(muaPhase > modeBoundaries(2) & muaPhase <= modeBoundaries(end)) = true;
    muaPhase = recentrePhase(muaPhase, 0);
    
    % Divide the populations
    spkDB_positive = sparse(spkDB(positivePhase, :));
    spkDB_units_positive = spkDB_units(positivePhase, :);
    spkDB_phase_positive = muaPhase(positivePhase);
    spkDB_negative = sparse(spkDB(negativePhase, :));
    spkDB_units_negative = spkDB_units(negativePhase, :);
    spkDB_phase_negative = muaPhase(negativePhase);
    
    MUAsAll_positive = zeros(size(MUAsAll));
    MUAsAll_negative = zeros(size(MUAsAll));
    
    for sh = 1:numel(shankIDs) % Loop through shanks
      fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
      
      % Divide shank MUAs
      if size(spkDB_units_positive,2) > 1
        MUAs_positive = full(sum(spkDB_positive(ismember(spkDB_units_positive(:,1),sh),:), 1));
      else
        MUAs_positive = full(sum(spkDB_positive, 1));
      end
      if ~isempty(MUAs_positive)
        MUAsAll_positive(sh,1:numel(MUAs_positive)) = MUAs_positive;
      end
      if size(spkDB_units_negative,2) > 1
        MUAs_negative = full(sum(spkDB_negative(ismember(spkDB_units_negative(:,1),sh),:), 1));
      else
        MUAs_negative = full(sum(spkDB_negative, 1));
      end
      if ~isempty(MUAs_negative)
        MUAsAll_negative(sh,1:numel(MUAs_negative)) = MUAs_negative;
      end
      
      % Load the contents of shankStruct
      [shankStruct, ~, units, unitMetadata, xcoords, ycoords, spk, ~, phaseCoh] = get_shankStruct(dbStruct, sh);
      unitHeader = shankStruct.unitHeader;
      if isempty(phaseCoh)
        shankData_positive.(['shank' num2str(sh)]) = [];
        dbStruct_positive.shankData = shankData_positive;
        shankData_negative.(['shank' num2str(sh)]) = [];
        dbStruct_negative.shankData = shankData_negative;
        continue
      end
      
      % Get the unit phase
      phaseCoh = dbStruct.shankData.(['shank' num2str(sh)]).pupil.unitData;
      unitPhase = zeros(1,numel(phaseCoh));
      for u = 1:numel(phaseCoh)
        freq = phaseCoh{u}.freq;
        freqInterp = unique([freq fRef fExtrapThr]);
        if ~isfield(phaseCoh{u}, 'phase') || sum(isnan(phaseCoh{u}.phase)) == numel(phaseCoh{u}.phase)
          unitPhase(u) = NaN;
          continue
        end
        if freq(1) <= fExtrapThr
          phase = interp1(freq, phaseCoh{u}.phase, freqInterp, 'linear', 'extrap');
        else
          phase = interp1(freq, phaseCoh{u}.phase, freqInterp, 'linear');
        end
        [~, iFreq] = min(abs(freqInterp - fRef));
        unitPhase(u) = phase(iFreq);
      end
      positivePhase = false(size(unitPhase));
      modeBoundaries = recentrePhase(modeBoundaries, modeBoundaries(1));
      unitPhase = recentrePhase(unitPhase, modeBoundaries(1));
      positivePhase(unitPhase > modeBoundaries(end) & unitPhase <= modeBoundaries(2)) = true;
      negativePhase = false(size(unitPhase));
      modeBoundaries = recentrePhase(modeBoundaries, modeBoundaries(3));
      unitPhase = recentrePhase(unitPhase, modeBoundaries(3));
      negativePhase(unitPhase > modeBoundaries(2) & unitPhase <= modeBoundaries(end)) = true;
      unitPhase = recentrePhase(unitPhase, 0);
      
      units_positive = units(positivePhase, :);
      phase_positive = unitPhase(positivePhase);
      unitMetadata_positive = unitMetadata(positivePhase,:);
      xcoords_positive = xcoords(positivePhase);
      ycoords_positive = ycoords(positivePhase);
      spk_positive = sparse(spk(positivePhase,:));
      
      units_negative = units(negativePhase, :);
      phase_negative = unitPhase(negativePhase);
      unitMetadata_negative = unitMetadata(negativePhase,:);
      xcoords_negative = xcoords(negativePhase);
      ycoords_negative = ycoords(negativePhase);
      spk_negative = sparse(spk(negativePhase,:));
      
      if ~isempty(units_positive)
        assert(numel(spkDB_units_positive(:,end)) >= numel(units_positive(:,end)));
      end
      shankStruct_positive = struct('shankID',sh, 'MUAs',full(MUAs_positive), 'spk',spk_positive,...
        'units',units_positive, 'unitHeader',{unitHeader}, 'unitMetadata',unitMetadata_positive,...
        'phase_positive',phase_positive);
      shankEntry = ['shank' num2str(sh)];
      shankData_positive.(shankEntry) = shankStruct_positive;
      dbStruct_positive.shankData = shankData_positive;
      
      if ~isempty(units_negative)
        assert(numel(spkDB_units_negative(:,end)) >= numel(units_negative(:,end)));
      end
      shankStruct_negative = struct('shankID',sh, 'MUAs',full(MUAs_negative), 'spk',spk_negative,...
        'units',units_negative, 'unitHeader',{unitHeader}, 'unitMetadata',unitMetadata_negative,...
        'phase_negative',phase_negative);
      shankEntry = ['shank' num2str(sh)];
      shankData_negative.(shankEntry) = shankStruct_negative;
      dbStruct_negative.shankData = shankData_negative;
      
      fprintf('Finished processing shank %i\n',sh);
    end % loop over shanks
    
    dbStruct_positive.popData = struct('MUAsAll',full(MUAsAll_positive), 'spkDB',spkDB_positive,...
      'spkDB_units',spkDB_units_positive, 'spkDB_phase_positive',spkDB_phase_positive);
    dbStruct_positive.io = dbStruct.io;
    dbStruct_positive.conf = dbStruct.conf;
    dbStruct_positive.db = dbStruct.db;
    dbStruct_positive.dbSeries = dbStruct.dbSeries;
    dbStruct_positive.splitType = pupilCorrCond;
    dataStruct.seriesData_positive.(entryName) = dbStruct_positive;
    
    dbStruct_negative.popData = struct('MUAsAll',full(MUAsAll_negative), 'spkDB',spkDB_negative,...
      'spkDB_units',spkDB_units_negative, 'spkDB_phase_negative',spkDB_phase_negative);
    dbStruct_negative.io = dbStruct.io;
    dbStruct_negative.conf = dbStruct.conf;
    dbStruct_negative.db = dbStruct.db;
    dbStruct_negative.dbSeries = dbStruct.dbSeries;
    dbStruct_negative.splitType = pupilCorrCond;
    dataStruct.seriesData_negative.(entryName) = dbStruct_negative;
    
    dataStruct.seriesData.(fnsData{dbCount}) = dbStruct;
    
    if intermediateSaving
      save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
    end
    
    fprintf('Finished processing db entry %i\n',dbCount);
  else
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
    %eyeDataDB.pupilArea = double(eyeDataDB.pupilAreaFilt.pupilAreaFiltHP0p01Hz);
    %eyeDataDB.frameTimes = eyeDataDB.pupilAreaFilt.timesFiltStart:eyeDataDB.pupilAreaFilt.timesFiltStep:eyeDataDB.pupilAreaFilt.timesFiltStop;
    [pupilArea, areaTimes] = pupilFilt(eyeDataDB, srData, MUAsAll, 2*srData, commonPeriod, srData);
    
    % Average the pupil area
    meanPupilArea = mean(pupilArea,'omitnan');
    averagingWindowSize10min = floor(10*60*srData);
    averagingWindowSize20min = floor(20*60*srData);
    averagingWindowSize30min = floor(30*60*srData);
    n10minWindows = numel(pupilArea)/averagingWindowSize10min;
    if n10minWindows > floor(n10minWindows)+0.95
      n10minWindows = ceil(n10minWindows);
    else
      n10minWindows = floor(n10minWindows);
    end
    n20minWindows = numel(pupilArea)/averagingWindowSize20min;
    if n20minWindows > floor(n20minWindows)+0.95
      n20minWindows = ceil(n20minWindows);
    else
      n20minWindows = floor(n20minWindows);
    end
    n30minWindows = numel(pupilArea)/averagingWindowSize30min;
    if n30minWindows > floor(n30minWindows)+0.95
      n30minWindows = ceil(n30minWindows);
    else
      n30minWindows = floor(n30minWindows);
    end
    meanPupilArea10minWindows = zeros(n10minWindows,1);
    meanPupilArea20minWindows = zeros(n20minWindows,1);
    meanPupilArea30minWindows = zeros(n30minWindows,1);
    for iWindow = 1:n10minWindows
      meanPupilArea10minWindows(iWindow) = mean(pupilArea((iWindow-1)*averagingWindowSize10min+1:...
        min([iWindow*averagingWindowSize10min numel(pupilArea)])),'omitnan');
    end
    for iWindow = 1:n20minWindows
      meanPupilArea20minWindows(iWindow) = mean(pupilArea((iWindow-1)*averagingWindowSize20min+1:...
        min([iWindow*averagingWindowSize20min numel(pupilArea)])),'omitnan');
    end
    for iWindow = 1:n30minWindows
      meanPupilArea30minWindows(iWindow) = mean(pupilArea((iWindow-1)*averagingWindowSize30min+1:...
        min([iWindow*averagingWindowSize30min numel(pupilArea)])),'omitnan');
    end
    
    % Down-sample the spiking matrix
    if iscell(commonPeriod)
      spkDB_ds = [];
      dsTimes = [];
      for iCell = 1:numel(commonPeriod)
        [~, spkDB_ds_period] = determineInds(commonPeriod{iCell}, srData, spkDB);
        %[spkDB_ds_period, dsTimes_period] = downsampleRasterMatrix(full(spkDB_ds_period), srData, dr);
        %dsTimes_period = dsTimes_period - dsTimes_period(1)/2;
        [spkDB_ds_period, dsTimes_period] = resampleSpikeCounts( ...
          full(spkDB_ds_period), stepsize=1/srData, newStepsize=1/dr);
        dsTimes_period = commonPeriod{iCell}(1) + dsTimes_period;
        if size(spkDB_ds_period,1) == 1 || size(spkDB_ds_period,2) == 1
          spkDB_ds_period = torow(spkDB_ds_period);
        end
        spkDB_ds = [spkDB_ds spkDB_ds_period]; %#ok<*AGROW>
        dsTimes = [dsTimes dsTimes_period];
      end
    else
      [~, spkDB_ds] = determineInds(commonPeriod, srData, spkDB);
      %[spkDB_ds, dsTimes] = downsampleRasterMatrix(full(spkDB_ds), srData, dr);
      %dsTimes = dsTimes - dsTimes(1)/2;
      [spkDB_ds, dsTimes] = resampleSpikeCounts( ...
        full(spkDB_ds), stepsize=1/srData, newStepsize=1/dr);
      dsTimes = commonPeriod(1) + dsTimes;
      if size(spkDB_ds,1) == 1 || size(spkDB_ds,2) == 1
        spkDB_ds = torow(spkDB_ds);
      end
    end
    assert(size(spkDB_ds,1) == size(spkDB,1));
    
    % Down-sample the pupil signal
    pupilArea = interp1(areaTimes, pupilArea, dsTimes, 'linear','extrap');
    
    % Correlate down-sampled pupil and spiking signals
    %nonemptyInds = logical(sum(spkDB_ds));
    [rPearson, pvalPearson] = corrMulti(pupilArea, spkDB_ds, 'Pearson');
    [rSpearman, pvalSpearman] = corrMulti(pupilArea, spkDB_ds, 'Spearman');
    
    % Correlate down-sampled pupil and spiking signals based on the first 10% of the total signal duration
    [rPearson10percent, pvalPearson10percent, rSpearman10percent, pvalSpearman10percent] = fractionCorrCalc(10, pupilArea, spkDB_ds);
    
    % Correlate down-sampled pupil and spiking signals based on the first 25% of the total signal duration
    [rPearson25percent, pvalPearson25percent, rSpearman25percent, pvalSpearman25percent] = fractionCorrCalc(4, pupilArea, spkDB_ds);
    
    % Correlate down-sampled pupil and spiking signals based on the first 33% of the total signal duration
    [rPearson33percent, pvalPearson33percent, rSpearman33percent, pvalSpearman33percent] = fractionCorrCalc(3, pupilArea, spkDB_ds);
    
    % Correlate down-sampled pupil and spiking signals based on the first 50% of the total signal duration
    [rPearson50percent, pvalPearson50percent, rSpearman50percent, pvalSpearman50percent] = fractionCorrCalc(2, pupilArea, spkDB_ds);
    
    % Correlate down-sampled pupil and spiking signals over 10, 20, and 30-minute windows
    correlationWindowSize10min = floor(10*60*dr);
    correlationWindowSize20min = floor(20*60*dr);
    correlationWindowSize30min = floor(30*60*dr);
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
        spkDB_ds(:,(iWindow-1)*correlationWindowSize10min+1:min([iWindow*correlationWindowSize10min size(spkDB_ds,2)])), 'Pearson');
      [rSpearman10minWindows{iWindow}, pvalSpearman10minWindows{iWindow}] =...
        corrMulti(pupilArea((iWindow-1)*correlationWindowSize10min+1:min([iWindow*correlationWindowSize10min numel(pupilArea)])),...
        spkDB_ds(:,(iWindow-1)*correlationWindowSize10min+1:min([iWindow*correlationWindowSize10min size(spkDB_ds,2)])), 'Spearman');
    end
    for iWindow = 1:n20minWindows
      [rPearson20minWindows{iWindow}, pvalPearson20minWindows{iWindow}] =...
        corrMulti(pupilArea((iWindow-1)*correlationWindowSize20min+1:min([iWindow*correlationWindowSize20min numel(pupilArea)])),...
        spkDB_ds(:,(iWindow-1)*correlationWindowSize20min+1:min([iWindow*correlationWindowSize20min size(spkDB_ds,2)])), 'Pearson');
      [rSpearman20minWindows{iWindow}, pvalSpearman20minWindows{iWindow}] =...
        corrMulti(pupilArea((iWindow-1)*correlationWindowSize20min+1:min([iWindow*correlationWindowSize20min numel(pupilArea)])),...
        spkDB_ds(:,(iWindow-1)*correlationWindowSize20min+1:min([iWindow*correlationWindowSize20min size(spkDB_ds,2)])), 'Spearman');
    end
    for iWindow = 1:n30minWindows
      [rPearson30minWindows{iWindow}, pvalPearson30minWindows{iWindow}] =...
        corrMulti(pupilArea((iWindow-1)*correlationWindowSize30min+1:min([iWindow*correlationWindowSize30min numel(pupilArea)])),...
        spkDB_ds(:,(iWindow-1)*correlationWindowSize30min+1:min([iWindow*correlationWindowSize30min size(spkDB_ds,2)])), 'Pearson');
      [rSpearman30minWindows{iWindow}, pvalSpearman30minWindows{iWindow}] =...
        corrMulti(pupilArea((iWindow-1)*correlationWindowSize30min+1:min([iWindow*correlationWindowSize30min numel(pupilArea)])),...
        spkDB_ds(:,(iWindow-1)*correlationWindowSize30min+1:min([iWindow*correlationWindowSize30min size(spkDB_ds,2)])), 'Spearman');
    end
    
    % Sort the spiking matrix based on the correlation coefficient values
    [sortedR, iSort] = sort(rSpearman, 'descend');
    iSort = iSort';
    sortedRaster = spkDB_ds(iSort,:);
    dividingLine = numel(sortedR) - (find(sortedR < 0, 1) - 2);
    
    if draw
      
      % Plot unsorted spiking raster
      figure;
      h = pcolor(logical([zeros(1,size(spkDB_ds,2)); flipud(spkDB_ds)]));
      h.EdgeColor = 'none';
      colormap(flipud(gray));
      xTicks = xticks;
      xTicks = xTicks(1:2:end);
      xTickLabel = string(xTicks);
      yTickLabel = string(fliplr(yticks));
      yLim = ylim;
      yLim = [yLim(1) yLim(2)+1];
      titleStr = 'Unit spiking prior to sorting';
      ax1 = axesProperties(titleStr, 1, 'normal', 'off', 'w', 'Calibri', 20, 1, 2, [0.01 0.025], 'out', 'off', 'k',...
        'Time (s)', xlim, xTicks, 'off', 'k', 'Unsorted units', yLim, fliplr(size(spkDB_ds,1)-yticks+1));
      ax1.XTickLabel = xTickLabel;
      ax1.YTickLabel = yTickLabel;
      set(gcf,'color','white');
      %   filename = titleStr;
      %   hgsave(gcf, filename);
      %   print(gcf, [filename '.png'],'-dpng','-r300');
      
      % Sorted spiking
      figure;
      h = pcolor(logical(flipud([zeros(1,size(spkDB_ds,2)); spkDB_ds(iSort,:)])));
      h.EdgeColor = 'none';
      colormap(flipud(gray));
      titleStr = 'Unit spiking after sorting';
      ax1 = axesProperties(titleStr, 1, 'normal', 'off', 'w', 'Calibri', 20, 1, 2, [0.01 0.025], 'out', 'off', 'k',...
        'Time (s)', xlim, xTicks, 'off', 'k', 'Sorted units', yLim, fliplr(size(spkDB_ds,1)-yticks+1));
      ax1.XTickLabel = xTickLabel;
      ax1.YTickLabel = yTickLabel;
      set(gcf,'color','white');
      %   filename = titleStr;
      %   hgsave(gcf, filename);
      %   print(gcf, [filename '.png'],'-dpng','-r300');
      
      if ~isempty(dividingLine)
        hold on
        plot(xlim, [dividingLine dividingLine], '--r', 'LineWidth',2)
        hold off
      end
    end
    
    % Divide the populations
    spkDB_positive = sparse(spkDB(rSpearman >= 0, :));
    spkDB_units_positive = spkDB_units(rSpearman >= 0, :);
    spkDB_negative = sparse(spkDB(rSpearman < 0, :));
    spkDB_units_negative = spkDB_units(rSpearman < 0, :);
    spkDB_neutral = sparse(spkDB(pvalSpearman >= alpha, :));
    spkDB_units_neutral = spkDB_units(pvalSpearman >= alpha, :);
    
    MUAsAll_positive = zeros(size(MUAsAll));
    MUAsAll_negative = zeros(size(MUAsAll));
    MUAsAll_neutral = zeros(size(MUAsAll));
    
    for sh = 1:numel(shankIDs) % Loop through shanks
      fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
      
      % Divide shank MUAs
      if size(spkDB_units_positive,2) > 1
        MUAs_positive = full(sum(spkDB_positive(ismember(spkDB_units_positive(:,1),sh),:), 1));
      else
        MUAs_positive = full(sum(spkDB_positive, 1));
      end
      if ~isempty(MUAs_positive)
        MUAsAll_positive(sh,1:numel(MUAs_positive)) = MUAs_positive;
      end
      if size(spkDB_units_negative,2) > 1
        MUAs_negative = full(sum(spkDB_negative(ismember(spkDB_units_negative(:,1),sh),:), 1));
      else
        MUAs_negative = full(sum(spkDB_negative, 1));
      end
      if ~isempty(MUAs_negative)
        MUAsAll_negative(sh,1:numel(MUAs_negative)) = MUAs_negative;
      end
      if size(spkDB_units_neutral,2) > 1
        MUAs_neutral = full(sum(spkDB_neutral(ismember(spkDB_units_neutral(:,1),sh),:), 1));
      else
        MUAs_neutral = full(sum(spkDB_neutral, 1));
      end
      if ~isempty(MUAs_neutral)
        MUAsAll_neutral(sh,1:numel(MUAs_neutral)) = MUAs_neutral;
      end
      
      % Load the contents of shankStruct
      [shankStruct, ~, units, unitMetadata, xcoords, ycoords, spk] = get_shankStruct(dbStruct, sh);
      unitHeader = shankStruct.unitHeader;
      
      % Down-sample the unit spiking matrix
      if iscell(commonPeriod)
        spk_ds = [];
        for iCell = 1:numel(commonPeriod)
          [~, spk_ds_period] = determineInds(commonPeriod{iCell}, srData, spk);
          %spk_ds_period = downsampleRasterMatrix(full(spk_ds_period), srData, dr);
          spkDB_ds_period = resampleSpikeCounts( ...
            full(spk_ds_period), stepsize=1/srData, newStepsize=1/dr);
          if size(spk_ds_period,1) == 1 || size(spk_ds_period,2) == 1
            spk_ds_period = torow(spk_ds_period);
          end
          spk_ds = [spk_ds spk_ds_period]; %#ok<*AGROW>
        end
      else
        [~, spk_ds] = determineInds(commonPeriod, srData, spk);
        %spk_ds = downsampleRasterMatrix(full(spk_ds), srData, dr);
        spkDB_ds = resampleSpikeCounts( ...
          full(spk_ds), stepsize=1/srData, newStepsize=1/dr);
        if size(spk_ds,1) == 1 || size(spk_ds,2) == 1
          spk_ds = torow(spk_ds);
        end
      end
      assert(size(spk_ds,1) == size(spk,1));
      
      % Correlate down-sampled pupil and spiking signals
      %nonemptyInds = logical(sum(spk_ds));
      [rPearsonUnits, pvalPearsonUnits] = corrMulti(pupilArea, spk_ds, 'Pearson');
      [rSpearmanUnits, pvalSpearmanUnits] = corrMulti(pupilArea, spk_ds, 'Spearman');
      assert(numel(rSpearmanUnits) == size(spk,1));
      assert(numel(rSpearmanUnits) == size(spk_ds,1));
      
      % Correlate down-sampled pupil and spiking signals based on the first 10% of the total signal duration
      [rPearsonUnits10percent, pvalPearsonUnits10percent, rSpearmanUnits10percent, pvalSpearmanUnits10percent] = fractionCorrCalc(10, pupilArea, spk_ds);
      
      % Correlate down-sampled pupil and spiking signals based on the first 25% of the total signal duration
      [rPearsonUnits25percent, pvalPearsonUnits25percent, rSpearmanUnits25percent, pvalSpearmanUnits25percent] = fractionCorrCalc(4, pupilArea, spk_ds);
      
      % Correlate down-sampled pupil and spiking signals based on the first 33% of the total signal duration
      [rPearsonUnits33percent, pvalPearsonUnits33percent, rSpearmanUnits33percent, pvalSpearmanUnits33percent] = fractionCorrCalc(3, pupilArea, spk_ds);
      
      % Correlate down-sampled pupil and spiking signals based on the first 50% of the total signal duration
      [rPearsonUnits50percent, pvalPearsonUnits50percent, rSpearmanUnits50percent, pvalSpearmanUnits50percent] = fractionCorrCalc(2, pupilArea, spk_ds);
      
      % Sort the spiking matrix based on the correlation coefficient values
      [sortedRUnits, iSortUnits] = sort(rSpearmanUnits, 'descend');
      iSortUnits = iSortUnits';
      sortedRasterUnits = spk_ds(iSortUnits,:);
      dividingLineUnits = numel(sortedRUnits) - (find(sortedRUnits < 0, 1) - 2);
      
      % Divide shank units
      if isempty(units)
        rPearsonUnits = [];
        pvalPearsonUnits = [];
        rSpearmanUnits = [];
        pvalSpearmanUnits = [];
        rPearsonUnits10minWindows = [];
        pvalPearsonUnits10minWindows = [];
        rSpearmanUnits10minWindows = [];
        pvalSpearmanUnits10minWindows = [];
        rPearsonUnits20minWindows = [];
        pvalPearsonUnits20minWindows = [];
        rSpearmanUnits20minWindows = [];
        pvalSpearmanUnits20minWindows = [];
        rPearsonUnits30minWindows = [];
        pvalPearsonUnits30minWindows = [];
        rSpearmanUnits30minWindows = [];
        pvalSpearmanUnits30minWindows = [];
        
        unitInds_positive = [];
        units_positive = [];
        unitMetadata_positive = [];
        xcoords_positive = [];
        ycoords_positive = [];
        spk_positive = [];
        
        unitInds_negative = [];
        units_negative = [];
        unitMetadata_negative = [];
        xcoords_negative = [];
        ycoords_negative = [];
        spk_negative = [];
        
%         unitInds_neutral = [];
%         units_neutral = [];
%         unitMetadata_neutral = [];
%         xcoords_neutral = [];
%         ycoords_neutral = [];
%         spk_neutral = [];
      else
        
        % Correlate down-sampled pupil and spiking signals over 10, 20, and 30-minute windows
        rPearsonUnits10minWindows = cell(n10minWindows,1);
        rPearsonUnits20minWindows = cell(n20minWindows,1);
        rPearsonUnits30minWindows = cell(n30minWindows,1);
        pvalPearsonUnits10minWindows = cell(n10minWindows,1);
        pvalPearsonUnits20minWindows = cell(n20minWindows,1);
        pvalPearsonUnits30minWindows = cell(n30minWindows,1);
        rSpearmanUnits10minWindows = cell(n10minWindows,1);
        rSpearmanUnits20minWindows = cell(n20minWindows,1);
        rSpearmanUnits30minWindows = cell(n30minWindows,1);
        pvalSpearmanUnits10minWindows = cell(n10minWindows,1);
        pvalSpearmanUnits20minWindows = cell(n20minWindows,1);
        pvalSpearmanUnits30minWindows = cell(n30minWindows,1);
        for iWindow = 1:n10minWindows
          [rPearsonUnits10minWindows{iWindow}, pvalPearsonUnits10minWindows{iWindow}] =...
            corrMulti(pupilArea((iWindow-1)*correlationWindowSize10min+1:min([iWindow*correlationWindowSize10min numel(pupilArea)])),...
            spk_ds(:,(iWindow-1)*correlationWindowSize10min+1:min([iWindow*correlationWindowSize10min size(spk_ds,2)])), 'Pearson');
          [rSpearmanUnits10minWindows{iWindow}, pvalSpearmanUnits10minWindows{iWindow}] =...
            corrMulti(pupilArea((iWindow-1)*correlationWindowSize10min+1:min([iWindow*correlationWindowSize10min numel(pupilArea)])),...
            spk_ds(:,(iWindow-1)*correlationWindowSize10min+1:min([iWindow*correlationWindowSize10min size(spk_ds,2)])), 'Spearman');
        end
        for iWindow = 1:n20minWindows
          [rPearsonUnits20minWindows{iWindow}, pvalPearsonUnits20minWindows{iWindow}] =...
            corrMulti(pupilArea((iWindow-1)*correlationWindowSize20min+1:min([iWindow*correlationWindowSize20min numel(pupilArea)])),...
            spk_ds(:,(iWindow-1)*correlationWindowSize20min+1:min([iWindow*correlationWindowSize20min size(spk_ds,2)])), 'Pearson');
          [rSpearmanUnits20minWindows{iWindow}, pvalSpearmanUnits20minWindows{iWindow}] =...
            corrMulti(pupilArea((iWindow-1)*correlationWindowSize20min+1:min([iWindow*correlationWindowSize20min numel(pupilArea)])),...
            spk_ds(:,(iWindow-1)*correlationWindowSize20min+1:min([iWindow*correlationWindowSize20min size(spk_ds,2)])), 'Spearman');
        end
        for iWindow = 1:n30minWindows
          [rPearsonUnits30minWindows{iWindow}, pvalPearsonUnits30minWindows{iWindow}] =...
            corrMulti(pupilArea((iWindow-1)*correlationWindowSize30min+1:min([iWindow*correlationWindowSize30min numel(pupilArea)])),...
            spk_ds(:,(iWindow-1)*correlationWindowSize30min+1:min([iWindow*correlationWindowSize30min size(spk_ds,2)])), 'Pearson');
          [rSpearmanUnits30minWindows{iWindow}, pvalSpearmanUnits30minWindows{iWindow}] =...
            corrMulti(pupilArea((iWindow-1)*correlationWindowSize30min+1:min([iWindow*correlationWindowSize30min numel(pupilArea)])),...
            spk_ds(:,(iWindow-1)*correlationWindowSize30min+1:min([iWindow*correlationWindowSize30min size(spk_ds,2)])), 'Spearman');
        end
        
        unitInds_positive = rSpearmanUnits >= 0;
        units_positive = units(unitInds_positive, :);
        unitHeader{numel(unitHeader)+1} = 'r_Pearson_corr_with_pupil';
        unitHeader{numel(unitHeader)+1} = 'pVal_Pearson_corr_with_pupil';
        unitHeader{numel(unitHeader)+1} = 'r_Spearman_corr_with_pupil';
        unitHeader{numel(unitHeader)+1} = 'pVal_Spearman_corr_with_pupil';
        unitMetadata_positive = [unitMetadata(unitInds_positive,:)...
          rPearsonUnits(unitInds_positive)' pvalPearsonUnits(unitInds_positive)'...
          rSpearmanUnits(unitInds_positive)' pvalSpearmanUnits(unitInds_positive)'];
        xcoords_positive = xcoords(unitInds_positive);
        ycoords_positive = ycoords(unitInds_positive);
        spk_positive = sparse(spk(unitInds_positive,:));
        
        unitInds_negative = rSpearmanUnits < 0;
        units_negative = units(unitInds_negative, :);
        unitMetadata_negative = [unitMetadata(unitInds_negative,:)...
          rPearsonUnits(unitInds_negative)' pvalPearsonUnits(unitInds_negative)'...
          rSpearmanUnits(unitInds_negative)' pvalSpearmanUnits(unitInds_negative)'];
        xcoords_negative = xcoords(unitInds_negative);
        ycoords_negative = ycoords(unitInds_negative);
        spk_negative = sparse(spk(unitInds_negative,:));
        assert(size(spk,1) == size(spk_positive,1)+size(spk_negative,1));
        
%         unitInds_neutral = pvalSpearmanUnits >= alpha;
%         units_neutral = units(unitInds_neutral, :);
%         unitMetadata_neutral = [unitMetadata(unitInds_neutral,:)...
%           rPearsonUnits(unitInds_neutral)' pvalPearsonUnits(unitInds_neutral)'...
%           rSpearmanUnits(unitInds_neutral)' pvalSpearmanUnits(unitInds_neutral)'];
%         xcoords_neutral = xcoords(unitInds_neutral);
%         ycoords_neutral = ycoords(unitInds_neutral);
%         spk_neutral = sparse(spk(unitInds_neutral,:));
      end
      
      shankStruct_positive = struct('shankID',sh, 'MUAs',full(MUAs_positive), 'spk',spk_positive,...
        'units',units_positive, 'unitHeader',{unitHeader}, 'unitMetadata',unitMetadata_positive);
      shankEntry = ['shank' num2str(sh)];
      shankData_positive.(shankEntry) = shankStruct_positive;
      dbStruct_positive.shankData = shankData_positive;
      
      shankStruct_negative = struct('shankID',sh, 'MUAs',full(MUAs_negative), 'spk',spk_negative,...
        'units',units_negative, 'unitHeader',{unitHeader}, 'unitMetadata',unitMetadata_negative);
      shankEntry = ['shank' num2str(sh)];
      shankData_negative.(shankEntry) = shankStruct_negative;
      dbStruct_negative.shankData = shankData_negative;
      
%       shankStruct_neutral = struct('shankID',sh, 'MUAs',full(MUAs_neutral), 'spk',spk_neutral,...
%         'units',units_neutral, 'unitHeader',{unitHeader}, 'unitMetadata',unitMetadata_neutral);
%       shankEntry = ['shank' num2str(sh)];
%       shankData_neutral.(shankEntry) = shankStruct_neutral;
%       dbStruct_neutral.shankData = shankData_neutral;
      
      dbStruct.shankData.(['shank' num2str(sh)]).rPearson = rPearsonUnits;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson = pvalPearsonUnits;
      dbStruct.shankData.(['shank' num2str(sh)]).rSpearman = rSpearmanUnits;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman = pvalSpearmanUnits;
      dbStruct.shankData.(['shank' num2str(sh)]).rPearson10percent = rPearsonUnits10percent;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson10percent = pvalPearsonUnits10percent;
      dbStruct.shankData.(['shank' num2str(sh)]).rSpearman10percent = rSpearmanUnits10percent;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman10percent = pvalSpearmanUnits10percent;
      dbStruct.shankData.(['shank' num2str(sh)]).rPearson25percent = rPearsonUnits25percent;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson25percent = pvalPearsonUnits25percent;
      dbStruct.shankData.(['shank' num2str(sh)]).rSpearman25percent = rSpearmanUnits25percent;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman25percent = pvalSpearmanUnits25percent;
      dbStruct.shankData.(['shank' num2str(sh)]).rPearson33percent = rPearsonUnits33percent;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson33percent = pvalPearsonUnits33percent;
      dbStruct.shankData.(['shank' num2str(sh)]).rSpearman33percent = rSpearmanUnits33percent;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman33percent = pvalSpearmanUnits33percent;
      dbStruct.shankData.(['shank' num2str(sh)]).rPearson50percent = rPearsonUnits50percent;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson50percent = pvalPearsonUnits50percent;
      dbStruct.shankData.(['shank' num2str(sh)]).rSpearman50percent = rSpearmanUnits50percent;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman50percent = pvalSpearmanUnits50percent;
      dbStruct.shankData.(['shank' num2str(sh)]).rPearson10minWindows = rPearsonUnits10minWindows;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson10minWindows = pvalPearsonUnits10minWindows;
      dbStruct.shankData.(['shank' num2str(sh)]).rSpearman10minWindows = rSpearmanUnits10minWindows;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman10minWindows = pvalSpearmanUnits10minWindows;
      dbStruct.shankData.(['shank' num2str(sh)]).rPearson20minWindows = rPearsonUnits20minWindows;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson20minWindows = pvalPearsonUnits20minWindows;
      dbStruct.shankData.(['shank' num2str(sh)]).rSpearman20minWindows = rSpearmanUnits20minWindows;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman20minWindows = pvalSpearmanUnits20minWindows;
      dbStruct.shankData.(['shank' num2str(sh)]).rPearson30minWindows = rPearsonUnits30minWindows;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalPearson30minWindows = pvalPearsonUnits30minWindows;
      dbStruct.shankData.(['shank' num2str(sh)]).rSpearman30minWindows = rSpearmanUnits30minWindows;
      dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearman30minWindows = pvalSpearmanUnits30minWindows;
      
      fprintf('Finished processing shank %i\n',sh);
    end % loop over shanks
    
    if ~isempty(units_positive)
      assert(numel(spkDB_units_positive(:,end)) >= numel(units_positive(:,end)));
    end
    dbStruct_positive.popData = struct('MUAsAll',full(MUAsAll_positive), 'spkDB',spkDB_positive,...
      'spkDB_units',spkDB_units_positive, 'rPearson',rPearson(rSpearman>=0), 'pvalPearson',pvalPearson(rSpearman>=0),...
      'rSpearman',rSpearman(rSpearman>=0), 'pvalSpearman',pvalSpearman(rSpearman>=0),...
      'rPearson10percent',rPearson10percent(rPearson10percent>=0), 'pvalPearson10percent',pvalPearson10percent(rPearson10percent>=0),...
      'rSpearman10percent',rSpearman10percent(rSpearman10percent>=0), 'pvalSpearman10percent',pvalSpearman10percent(rSpearman10percent>=0),...
      'rPearson25percent',rPearson25percent(rPearson25percent>=0), 'pvalPearson25percent',pvalPearson25percent(rPearson25percent>=0),...
      'rSpearman25percent',rSpearman25percent(rSpearman25percent>=0), 'pvalSpearman25percent',pvalSpearman25percent(rSpearman25percent>=0),...
      'rPearson33percent',rPearson33percent(rPearson33percent>=0), 'pvalPearson33percent',pvalPearson33percent(rPearson33percent>=0),...
      'rSpearman33percent',rSpearman33percent(rSpearman33percent>=0), 'pvalSpearman33percent',pvalSpearman33percent(rSpearman33percent>=0),...
      'rPearson50percent',rPearson50percent(rPearson50percent>=0), 'pvalPearson50percent',pvalPearson50percent(rPearson50percent>=0),...
      'rSpearman50percent',rSpearman50percent(rSpearman50percent>=0), 'pvalSpearman50percent',pvalSpearman50percent(rSpearman50percent>=0));
    dbStruct_positive.io = dbStruct.io;
    dbStruct_positive.conf = dbStruct.conf;
    dbStruct_positive.conf.drData = dr;
    dbStruct_positive.db = dbStruct.db;
    dbStruct_positive.dbSeries = dbStruct.dbSeries;
    dbStruct_positive.splitType = pupilCorrCond;
    dataStruct.seriesData_positive.(entryName) = dbStruct_positive;
    
    if ~isempty(units_negative)
      assert(numel(spkDB_units_negative(:,end)) >= numel(units_negative(:,end)));
    end
    dbStruct_negative.popData = struct('MUAsAll',full(MUAsAll_negative), 'spkDB',spkDB_negative,...
      'spkDB_units',spkDB_units_negative, 'rPearson',rPearson(rSpearman<0), 'pvalPearson',pvalPearson(rSpearman<0),...
      'rSpearman',rSpearman(rSpearman<0), 'pvalSpearman',pvalSpearman(rSpearman<0),...
      'rPearson10percent',rPearson10percent(rPearson10percent<0), 'pvalPearson10percent',pvalPearson10percent(rPearson10percent<0),...
      'rSpearman10percent',rSpearman10percent(rSpearman10percent<0), 'pvalSpearman10percent',pvalSpearman10percent(rSpearman10percent<0),...
      'rPearson25percent',rPearson25percent(rPearson25percent<0), 'pvalPearson25percent',pvalPearson25percent(rPearson25percent<0),...
      'rSpearman25percent',rSpearman25percent(rSpearman25percent<0), 'pvalSpearman25percent',pvalSpearman25percent(rSpearman25percent<0),...
      'rPearson33percent',rPearson33percent(rPearson33percent<0), 'pvalPearson33percent',pvalPearson33percent(rPearson33percent<0),...
      'rSpearman33percent',rSpearman33percent(rSpearman33percent<0), 'pvalSpearman33percent',pvalSpearman33percent(rSpearman33percent<0),...
      'rPearson50percent',rPearson50percent(rPearson50percent<0), 'pvalPearson50percent',pvalPearson50percent(rPearson50percent<0),...
      'rSpearman50percent',rSpearman50percent(rSpearman50percent<0), 'pvalSpearman50percent',pvalSpearman50percent(rSpearman50percent<0));
    dbStruct_negative.io = dbStruct.io;
    dbStruct_negative.conf = dbStruct.conf;
    dbStruct_negative.conf.drData = dr;
    dbStruct_negative.db = dbStruct.db;
    dbStruct_negative.dbSeries = dbStruct.dbSeries;
    dbStruct_negative.splitType = pupilCorrCond;
    dataStruct.seriesData_negative.(entryName) = dbStruct_negative;
    
%     if ~isempty(units_neutral)
%       assert(numel(spkDB_units_neutral(:,end)) >= numel(units_neutral(:,end)));
%     end
%     dbStruct_neutral.popData = struct('MUAsAll',full(MUAsAll_neutral), 'spkDB',spkDB_neutral,...
%       'spkDB_units',spkDB_units_neutral, 'rPearson',rPearson(pvalSpearman >= alpha), 'pvalPearson',pvalPearson(pvalSpearman >= alpha),...
%       'rSpearman',rSpearman(pvalSpearman >= alpha), 'pvalSpearman',pvalSpearman(pvalSpearman >= alpha),...
%       'rPearson10percent',rPearson10percent(pvalPearson10percent >= alpha), 'pvalPearson10percent',pvalPearson10percent(pvalPearson10percent >= alpha),...
%       'rSpearman10percent',rSpearman10percent(pvalSpearman10percent >= alpha), 'pvalSpearman10percent',pvalSpearman10percent(pvalSpearman10percent >= alpha),...
%       'rPearson25percent',rPearson25percent(pvalPearson25percent >= alpha), 'pvalPearson25percent',pvalPearson25percent(pvalPearson25percent >= alpha),...
%       'rSpearman25percent',rSpearman25percent(pvalSpearman25percent >= alpha), 'pvalSpearman25percent',pvalSpearman25percent(pvalSpearman25percent >= alpha),...
%       'rPearson33percent',rPearson33percent(pvalPearson33percent >= alpha), 'pvalPearson33percent',pvalPearson33percent(pvalPearson33percent >= alpha),...
%       'rSpearman33percent',rSpearman33percent(pvalSpearman33percent >= alpha), 'pvalSpearman33percent',pvalSpearman33percent(pvalSpearman33percent >= alpha),...
%       'rPearson50percent',rPearson50percent(pvalPearson50percent >= alpha), 'pvalPearson50percent',pvalPearson50percent(pvalPearson50percent >= alpha),...
%       'rSpearman50percent',rSpearman50percent(pvalSpearman50percent >= alpha), 'pvalSpearman50percent',pvalSpearman50percent(pvalSpearman50percent >= alpha));
%     dbStruct_neutral.io = dbStruct.io;
%     dbStruct_neutral.conf = dbStruct.conf;
%     dbStruct_neutral.conf.drData = dr;
%     dbStruct_neutral.db = dbStruct.db;
%     dbStruct_neutral.dbSeries = dbStruct.dbSeries;
%     dataStruct.seriesData_neutral.(entryName) = dbStruct_neutral;
    
    dbStruct.popData.meanPupilArea = meanPupilArea;
    dbStruct.popData.meanPupilArea10minWindows = meanPupilArea10minWindows;
    dbStruct.popData.meanPupilArea20minWindows = meanPupilArea20minWindows;
    dbStruct.popData.meanPupilArea30minWindows = meanPupilArea30minWindows;
    dbStruct.popData.rPearson = rPearson;
    dbStruct.popData.pvalPearson = pvalPearson;
    dbStruct.popData.rSpearman = rSpearman;
    dbStruct.popData.pvalSpearman = pvalSpearman;
    dbStruct.popData.rPearson10percent = rPearson10percent;
    dbStruct.popData.pvalPearson10percent = pvalPearson10percent;
    dbStruct.popData.rSpearman10percent = rSpearman10percent;
    dbStruct.popData.pvalSpearman10percent = pvalSpearman10percent;
    dbStruct.popData.rPearson25percent = rPearson25percent;
    dbStruct.popData.pvalPearson25percent = pvalPearson25percent;
    dbStruct.popData.rSpearman25percent = rSpearman25percent;
    dbStruct.popData.pvalSpearman25percent = pvalSpearman25percent;
    dbStruct.popData.rPearson33percent = rPearson33percent;
    dbStruct.popData.pvalPearson33percent = pvalPearson33percent;
    dbStruct.popData.rSpearman33percent = rSpearman33percent;
    dbStruct.popData.pvalSpearman33percent = pvalSpearman33percent;
    dbStruct.popData.rPearson50percent = rPearson50percent;
    dbStruct.popData.pvalPearson50percent = pvalPearson50percent;
    dbStruct.popData.rSpearman50percent = rSpearman50percent;
    dbStruct.popData.pvalSpearman50percent = pvalSpearman50percent;
    dbStruct.popData.rPearson10minWindows = rPearson10minWindows;
    dbStruct.popData.pvalPearson10minWindows = pvalPearson10minWindows;
    dbStruct.popData.rSpearman10minWindows = rSpearman10minWindows;
    dbStruct.popData.pvalSpearman10minWindows = pvalSpearman10minWindows;
    dbStruct.popData.rPearson20minWindows = rPearson20minWindows;
    dbStruct.popData.pvalPearson20minWindows = pvalPearson20minWindows;
    dbStruct.popData.rSpearman20minWindows = rSpearman20minWindows;
    dbStruct.popData.pvalSpearman20minWindows = pvalSpearman20minWindows;
    dbStruct.popData.rPearson30minWindows = rPearson30minWindows;
    dbStruct.popData.pvalPearson30minWindows = pvalPearson30minWindows;
    dbStruct.popData.rSpearman30minWindows = rSpearman30minWindows;
    dbStruct.popData.pvalSpearman30minWindows = pvalSpearman30minWindows;
    
    dataStruct.seriesData.(fnsData{dbCount}) = dbStruct;
    
    if intermediateSaving
      save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
    end
    
    fprintf('Finished processing db entry %i\n',dbCount);
  end
end % loop over db entries


%% SAVE DATA
if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct



%% Local functions
function [rPearson, pvalPearson, rSpearman, pvalSpearman] = fractionCorrCalc(nFractions, pupilArea, spk)

if isempty(pupilArea) || isempty(spk)
  rPearson = []; pvalPearson = []; rSpearman = []; pvalSpearman = [];
else
  lengthFractionalDuration = round((1/nFractions)*numel(pupilArea));
  rPearson = zeros(nFractions, size(spk,1));
  pvalPearson = zeros(nFractions, size(spk,1));
  rSpearman = zeros(nFractions, size(spk,1));
  pvalSpearman = zeros(nFractions, size(spk,1));
  for iFraction = 1:nFractions
    if iFraction == nFractions
      [rPearson(iFraction,:), pvalPearson(iFraction,:)] = corrMulti(pupilArea((iFraction-1)*lengthFractionalDuration+1:end),...
        spk(:,(iFraction-1)*lengthFractionalDuration+1:end), 'Pearson');
      [rSpearman(iFraction,:), pvalSpearman(iFraction,:)] = corrMulti(pupilArea((iFraction-1)*lengthFractionalDuration+1:end),...
        spk(:,(iFraction-1)*lengthFractionalDuration+1:end), 'Spearman');
    else
      [rPearson(iFraction,:), pvalPearson(iFraction,:)] = corrMulti(pupilArea((iFraction-1)*lengthFractionalDuration+1:iFraction*lengthFractionalDuration),...
        spk(:,(iFraction-1)*lengthFractionalDuration+1:iFraction*lengthFractionalDuration), 'Pearson');
      [rSpearman(iFraction,:), pvalSpearman(iFraction,:)] = corrMulti(pupilArea((iFraction-1)*lengthFractionalDuration+1:iFraction*lengthFractionalDuration),...
        spk(:,(iFraction-1)*lengthFractionalDuration+1:iFraction*lengthFractionalDuration), 'Spearman');
    end
  end
end
end