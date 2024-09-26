% Run this script to perform PCA analyses.

%% LOAD PRE-PROCESSED DATA
load(dataFile);


%% INITIALISE PARAMETERS
lists
pcaParams


%% LOAD RECORDINGS AND EXTRACT LFP MEASURES
% LOOP THROUGH DB ENTRIES
fnsData = fieldnames(dataStruct.seriesData);
for dbCount = dbStart:numel(fnsData)
  if dbCount > dbEnd
    break
  end
  dbStruct = dataStruct.seriesData.(fnsData{dbCount});
  
  
  % TEST FOR EXCEPTIONS
  [seriesName, animal] = seriesFromEntry(dbStruct.db(dbCount).entryName);
  if numel(seriesName) <= 14
    continue
  elseif numel(seriesName) == 15 && (strcmpi(seriesName(15), '4') || strcmpi(seriesName(15), '5'))
    continue
  elseif numel(seriesName) >= 16
    continue
  end
  
  
  % TEST FOR RIPPLES
  ripplesExist = rippleTest(fnsData{dbCount}, fullSeries);
  
  
  % IDENTIFY AREA
  area = determineAreaFromSeries(seriesName);
  if ~area
    continue
  end
  if area == 2 && size(LFPs,1) == 2
    continue
  end
  
  
  % PICK UP LFP CHANNEL
  if ~(numel(seriesName) == 16)
    chOI = dbStruct.lfpPowerData.chOIDB;
    iCh = pickChan(area, animal, chOI);
    if isempty(iCh)
      continue
    end
  end
  
  
  % LOAD THE CONTENTS OF THE DB STRUCTURE VARIABLE
  params = dbStruct.conf.params;
  probe = dbStruct.conf.probe;
  if strcmpi(probe, 'Neuropixels')
    strInd = strfind(dbStruct.io.baseFilename, 'ap');
    baseFilename = [dbStruct.io.baseFilename(1:strInd-1) 'lf.bin'];
    if ~exist(baseFilename, 'file')
      baseFilename = [dbStruct.io.baseFilename(1:strInd-1) 'lf.dat'];
    end
  else
    baseFilename = [dbStruct.io.baseFilename '.dat'];
  end
  if ~exist(baseFilename, 'file')
    error(['The supplied LFP file does not exist. Please check the location folder: ' dbStruct.io.dataDir]);
  end
  chN = dbStruct.db(dbCount).chN;
  entryName = dbStruct.db(dbCount).entryName;
  
  
  % LOAD EXISTING LFP ANALYSIS DATA
  if strcmpi(seriesName(15), '1')
    neuralActivity = [];
    LFPs = [];
    PRs = [];
    areas = [];
    entryNames = {};
    slowPowerS1 = [];
    slowPowerTh = [];
    slowPowerHp = [];
    slowPowerRSC = [];
    fastPowerS1 = [];
    fastPowerTh = [];
    fastPowerHp = [];
    fastPowerRSC = [];
    LFP_S1 = [];
    LFP_Th = [];
    LFP_Hp = [];
    LFP_RSC = [];
  end
  if isfield(dbStruct,'lfpPowerData')
    if (area == 5 || area == 10) && medianSubtracted % It's ok to use the existing ripple rate because it is based on median subtracted LFP
      rippleRate = dbStruct.lfpPowerData.rippleRate;
    else
      rippleRate = {};
      rippleRate{iCh} = []; %#ok<*SAGROW>
    end
    theta2deltaRatio = dbStruct.lfpPowerData.theta2deltaRatio;
    chOIDB = dbStruct.lfpPowerData.chOIDB;
    dssrLFPinit = dbStruct.lfpPowerData.dssrLFPinit;
    lfpTimes = dbStruct.lfpPowerData.lfpTimes;
    rippleDuration = dbStruct.lfpPowerData.rippleDuration;
    sdGaussian = dbStruct.lfpPowerData.sdGaussian;
    wGaussian = dbStruct.lfpPowerData.wGaussian;
  else
    error(['No LFP analysis data exist for series: ' fnsData{dbCount} '. Please run lfpLoad.m']);
  end
  
  
  % RESAMPLE EXISTING LFP DATA IF NECESSARY
  if params.srData ~= dssrLFPinit
    interpTimes = 1/params.srData:1/params.srData:lfpTimes(end);
    if (area == 5 || area == 10) && medianSubtracted
      rippleRate{iCh} = interp1(lfpTimes, rippleRate{iCh}, interpTimes)';
    end
    theta2deltaRatio{iCh} = interp1(lfpTimes, theta2deltaRatio{iCh}, interpTimes)';
  else
    interpTimes = lfpTimes;
  end
  
  
  % RECALCULATE RIPPLE RATE IF EXISTING ONE IS NOT BASED ON MEDIAN-SUBTRACTED LFP
  if (area == 5 || area == 10) && ~medianSubtracted && ripplesExist
    srRippleRate = 1000;
    if strcmpi(probe, 'Neuropixels')
      [LFPbandPower] = lfpFileLoad(baseFilename, 1, chunkSize, chN, chOIDB(iCh), LFPbands, params.srRecordingLFP, srRippleRate, true, deleteChans, false);
    else
      [LFPbandPower] = lfpFileLoad(baseFilename, 1, chunkSize, chN, chOIDB(iCh), LFPbands, params.srRecording, srRippleRate, true, deleteChans, false);
    end
    rippleBandPower{iCh} = LFPbandPower{1}(end,:);
    rippleRate{iCh} = rippleRateCalculator({rippleBandPower{iCh}}, rippleDuration, sdGaussian, wGaussian, chOIDB(iCh), bandNames, srRippleRate); %#ok<CCAT1>
    rippleRate{iCh} = rippleRate{iCh}{1};
    rippleRateTimes = 1/srRippleRate:1/srRippleRate:numel(rippleRate{iCh})/srRippleRate;
    rippleRate{iCh} = interp1(rippleRateTimes, rippleRate{iCh}, interpTimes)';
  elseif (area == 5 || area == 10) && ~ripplesExist
    rippleRate{iCh} = zeros(size(theta2deltaRatio{iCh}));
  end
  
  
  % LOAD, DOWN-SAMPLE LFP DATA, AND CALCULATE WAVELET TRANSFORMS
  if strcmpi(probe, 'Neuropixels')
    [LFPbandPower, ~, ~, LFP] = lfpFileLoad(...
      baseFilename, min(LFPbands{1}), chunkSize, chN, chOIDB, LFPbands, params.srRecordingLFP, params.srData, medianSubtracted, deleteChans, false);
  else
    [LFPbandPower, ~, ~, LFP] = lfpFileLoad(...
      baseFilename, min(LFPbands{1}), chunkSize, chN, chOIDB, LFPbands, params.srRecording, params.srData, medianSubtracted, deleteChans, false);
  end
  
  
  % LOAD SPIKING DATA
  PR = sum(dbStruct.popData.MUAsAll,1);
  
  
  % CONCATENATE DATA
  PR = torow(PR);
  LFP{iCh} = torow(LFP{iCh});
  rippleRate{iCh} = torow(rippleRate{iCh});
  theta2deltaRatio{iCh} = torow(theta2deltaRatio{iCh});
  dataDuration = max([size(neuralActivity,2) size(LFPbandPower{iCh},2) numel(rippleRate{iCh}) numel(theta2deltaRatio{iCh}) numel(LFP{iCh}) numel(PR)]);
  if isempty(PR)
    PR = zeros(1,dataDuration);
  end
  if ~isempty(neuralActivity) && size(neuralActivity,2) < dataDuration
    zeroPad = zeros(size(neuralActivity,1), dataDuration - size(neuralActivity,2));
    neuralActivity = [neuralActivity zeroPad];
  end
  if size(LFPbandPower{iCh},2) < dataDuration
    zeroPad = zeros(size(LFPbandPower{iCh},1), dataDuration - size(LFPbandPower{iCh},2));
    LFPbandPower{iCh} = [LFPbandPower{iCh} zeroPad];
  end
  if (area == 5 || area == 10) && size(rippleRate{iCh},2) < dataDuration
    zeroPad = zeros(size(rippleRate{iCh},1), dataDuration - size(rippleRate{iCh},2));
    rippleRate{iCh} = [rippleRate{iCh} zeroPad];
  end
  if size(theta2deltaRatio{iCh},2) < dataDuration
    zeroPad = zeros(size(theta2deltaRatio{iCh},1), dataDuration - size(theta2deltaRatio{iCh},2));
    theta2deltaRatio{iCh} = [theta2deltaRatio{iCh} zeroPad];
  end
  if size(LFP{iCh},2) < dataDuration
    zeroPad = zeros(size(LFP{iCh},1), dataDuration - size(LFP{iCh},2));
    LFP{iCh} = [LFP{iCh} zeroPad];
  end
  if size(PR,2) < dataDuration
    zeroPad = zeros(size(PR,1), dataDuration - size(PR,2));
    PR = [PR zeroPad];
  end
  if (area == 5 || area == 10)
    neuralActivity = [neuralActivity; LFPbandPower{iCh}; rippleRate{iCh}; theta2deltaRatio{iCh}; LFP{iCh}; PR];
  else
    neuralActivity = [neuralActivity; LFPbandPower{iCh}; theta2deltaRatio{iCh}; LFP{iCh}; PR];
  end
  addedVectors = size(LFPbandPower{iCh},1) + ~isempty(rippleRate{iCh}) + size(theta2deltaRatio{iCh},1) + size(LFP{iCh},1) + ~isempty(PR);
  neuralActivity(isnan(neuralActivity)) = 0;
  neuralActivity(isnan(neuralActivity)) = 0;
  
  if size(LFPs,2) < dataDuration
    zeroPad = zeros(size(LFPs,1), dataDuration - size(LFPs,2));
    LFPs = [LFPs zeroPad];
  end
  LFPs = [LFPs; LFP{iCh}];
  
  if size(PRs,2) < dataDuration
    zeroPad = zeros(size(PRs,1), dataDuration - size(PRs,2));
    PRs = [PRs zeroPad];
  end
  PRs = [PRs; PR];
  
  areas = [areas; ones(addedVectors,1)*area];
  
  entryNames{numel(entryNames) + 1} = entryName;
  
  
  % CALCULATE SLOW AND FAST POWER IN EACH AREA
  if area == 1
    slowPowerS1 = sum(LFPbandPower{iCh}(1:4,:),1);
    fastPowerS1 = sum(LFPbandPower{iCh}(5:8,:),1);
    LFP_S1 = LFP{iCh};
  elseif area == 2
    slowPowerTh = sum(LFPbandPower{iCh}(1:4,:),1);
    fastPowerTh = sum(LFPbandPower{iCh}(5:8,:),1);
    LFP_Th = LFP{iCh};
  elseif (area == 5 || area == 10)
    slowPowerHp = sum(LFPbandPower{iCh}(1:4,:),1);
    fastPowerHp = sum(LFPbandPower{iCh}(5:8,:),1);
    LFP_Hp = LFP{iCh};
  elseif area == 7
    slowPowerRSC = sum(LFPbandPower{iCh}(1:4,:),1);
    fastPowerRSC = sum(LFPbandPower{iCh}(5:8,:),1);
    LFP_RSC = LFP{iCh};
  end
  
  
  if dbCount < numel(fnsData)
    seriesNameNext = seriesFromEntry(dbStruct.db(dbCount+1).entryName);
  end
  if dbCount < numel(fnsData) - 1
    seriesNameNext2 = seriesFromEntry(dbStruct.db(dbCount+2).entryName);
  end
  if dbCount == numel(fnsData) || ~strcmpi(seriesName(1:14), seriesNameNext(1:14)) ||...
      (dbCount == numel(fnsData)-1 && (area == 5 || area == 10)) || ((area == 5 || area == 10) && ~strcmpi(seriesName(1:14), seriesNameNext2(1:14)))
    
    % SUBTRACT THE MEAN AND Z-SCORE THE DATA
    neuralActivity = neuralActivity - mean(neuralActivity,2);
    neuralActivity = zscore(neuralActivity')';
    
    
    % RUN PCA
    [pcaCoef, explained, PCs, nPCs, prob, individualVarExplained] = pcaGeneric(neuralActivity);
    
    % Calculate a number of variables needed to explain proportions of variance
    [explained25, explained50, explained75, explained95, explainedThird, explainedTwoThirds] = explainedProportions(explained);
    
    
    % CORRELATE PCs WITH MEASURES OF AROUSAL
    % Load arousal measures
    if isfield(dbStruct.lfpPowerData, 'areaInterpFilt')
      pupilArea = dbStruct.lfpPowerData.areaInterpFilt;
      pupilAreaTimes = dbStruct.lfpPowerData.areaInterpTimes;
    else
      pupilArea = [];
      pupilAreaTimes = [];
    end
    if isfield(dbStruct.lfpPowerData, 'motionInterpFilt')
      motion = dbStruct.lfpPowerData.motionInterpFilt;
      motionTimes = dbStruct.lfpPowerData.motionInterpTimes;
    else
      motion = [];
      motionTimes = [];
    end
    
    % Truncate PCs, Power vectors, and LFPs to match arousal measures
    if ~isempty(pupilArea)
      [~, iPCs2pupilAreaStart] = min(abs(interpTimes - pupilAreaTimes(1)));
      [~, iPCs2pupilAreaEnd] = min(abs(interpTimes - pupilAreaTimes(end)));
      PCs2pupilAreaTimes = interpTimes(iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea = PCs(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      slowPowerS12pupilArea = slowPowerS1(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      slowPowerTh2pupilArea = slowPowerTh(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      slowPowerHp2pupilArea = slowPowerHp(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      slowPowerRSC2pupilArea = slowPowerRSC(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      fastPowerS12pupilArea = fastPowerS1(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      fastPowerTh2pupilArea = fastPowerTh(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      fastPowerHp2pupilArea = fastPowerHp(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      fastPowerRSC2pupilArea = fastPowerRSC(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
    end
    if ~isempty(motion)
      [~, iPCs2motionStart] = min(abs(interpTimes - motionTimes(1)));
      [~, iPCs2motionEnd] = min(abs(interpTimes - motionTimes(end)));
      PCs2motionTimes = interpTimes(iPCs2motionStart:iPCs2motionEnd);
      PCs2motion = PCs(:,iPCs2motionStart:iPCs2motionEnd);
      slowPowerS12motion = slowPowerS1(:,iPCs2motionStart:iPCs2motionEnd);
      slowPowerTh2motion = slowPowerTh(:,iPCs2motionStart:iPCs2motionEnd);
      slowPowerHp2motion = slowPowerHp(:,iPCs2motionStart:iPCs2motionEnd);
      slowPowerRSC2motion = slowPowerRSC(:,iPCs2motionStart:iPCs2motionEnd);
      fastPowerS12motion = fastPowerS1(:,iPCs2motionStart:iPCs2motionEnd);
      fastPowerTh2motion = fastPowerTh(:,iPCs2motionStart:iPCs2motionEnd);
      fastPowerHp2motion = fastPowerHp(:,iPCs2motionStart:iPCs2motionEnd);
      fastPowerRSC2motion = fastPowerRSC(:,iPCs2motionStart:iPCs2motionEnd);
    end
    
    % Resample arousal measures if necessary
    if params.srData ~= dssrLFPinit
      if ~isempty(pupilArea)
        pupilArea = interp1(pupilAreaTimes, pupilArea, PCs2pupilAreaTimes)';
      end
      if ~isempty(motion)
        motion = interp1(motionTimes, motion, PCs2motionTimes)';
      end
    end
    
    % Correlate the vectors
    if ~isempty(pupilArea)
      [rPCs2pupilAreaPearson, pvalPCs2pupilAreaPearson] = corrMulti(pupilArea', PCs2pupilArea, 'Pearson');
      [rPCs2pupilAreaSpearman, pvalPCs2pupilAreaSpearman] = corrMulti(pupilArea', PCs2pupilArea, 'Spearman');
      [rSlowPowerS12pupilAreaPearson, pvalSlowPowerS12pupilAreaPearson] = corrMulti(pupilArea', slowPowerS12pupilArea, 'Pearson');
      [rSlowPowerS12pupilAreaSpearman, pvalSlowPowerS12pupilAreaSpearman] = corrMulti(pupilArea', slowPowerS12pupilArea, 'Spearman');
      [rSlowPowerTh2pupilAreaPearson, pvalSlowPowerTh2pupilAreaPearson] = corrMulti(pupilArea', slowPowerTh2pupilArea, 'Pearson');
      [rSlowPowerTh2pupilAreaSpearman, pvalSlowPowerTh2pupilAreaSpearman] = corrMulti(pupilArea', slowPowerTh2pupilArea, 'Spearman');
      [rSlowPowerHp2pupilAreaPearson, pvalSlowPowerHp2pupilAreaPearson] = corrMulti(pupilArea', slowPowerHp2pupilArea, 'Pearson');
      [rSlowPowerHp2pupilAreaSpearman, pvalSlowPowerHp2pupilAreaSpearman] = corrMulti(pupilArea', slowPowerHp2pupilArea, 'Spearman');
      [rSlowPowerRSC2pupilAreaPearson, pvalSlowPowerRSC2pupilAreaPearson] = corrMulti(pupilArea', slowPowerRSC2pupilArea, 'Pearson');
      [rSlowPowerRSC2pupilAreaSpearman, pvalSlowPowerRSC2pupilAreaSpearman] = corrMulti(pupilArea', slowPowerRSC2pupilArea, 'Spearman');
      [rFastPowerS12pupilAreaPearson, pvalFastPowerS12pupilAreaPearson] = corrMulti(pupilArea', fastPowerS12pupilArea, 'Pearson');
      [rFastPowerS12pupilAreaSpearman, pvalFastPowerS12pupilAreaSpearman] = corrMulti(pupilArea', fastPowerS12pupilArea, 'Spearman');
      [rFastPowerTh2pupilAreaPearson, pvalFastPowerTh2pupilAreaPearson] = corrMulti(pupilArea', fastPowerTh2pupilArea, 'Pearson');
      [rFastPowerTh2pupilAreaSpearman, pvalFastPowerTh2pupilAreaSpearman] = corrMulti(pupilArea', fastPowerTh2pupilArea, 'Spearman');
      [rFastPowerHp2pupilAreaPearson, pvalFastPowerHp2pupilAreaPearson] = corrMulti(pupilArea', fastPowerHp2pupilArea, 'Pearson');
      [rFastPowerHp2pupilAreaSpearman, pvalFastPowerHp2pupilAreaSpearman] = corrMulti(pupilArea', fastPowerHp2pupilArea, 'Spearman');
      [rFastPowerRSC2pupilAreaPearson, pvalFastPowerRSC2pupilAreaPearson] = corrMulti(pupilArea', fastPowerRSC2pupilArea, 'Pearson');
      [rFastPowerRSC2pupilAreaSpearman, pvalFastPowerRSC2pupilAreaSpearman] = corrMulti(pupilArea', fastPowerRSC2pupilArea, 'Spearman');
    end
    if ~isempty(motion)
      [rPCs2motionPearson, pvalPCs2motionPearson] = corrMulti(motion', PCs2motion, 'Pearson');
      [rPCs2motionSpearman, pvalPCs2motionSpearman] = corrMulti(motion', PCs2motion, 'Spearman');
      [rSlowPowerS12motionPearson, pvalSlowPowerS12motionPearson] = corrMulti(motion', slowPowerS12motion, 'Pearson');
      [rSlowPowerS12motionSpearman, pvalSlowPowerS12motionSpearman] = corrMulti(motion', slowPowerS12motion, 'Spearman');
      [rSlowPowerTh2motionPearson, pvalSlowPowerTh2motionPearson] = corrMulti(motion', slowPowerTh2motion, 'Pearson');
      [rSlowPowerTh2motionSpearman, pvalSlowPowerTh2motionSpearman] = corrMulti(motion', slowPowerTh2motion, 'Spearman');
      [rSlowPowerHp2motionPearson, pvalSlowPowerHp2motionPearson] = corrMulti(motion', slowPowerHp2motion, 'Pearson');
      [rSlowPowerHp2motionSpearman, pvalSlowPowerHp2motionSpearman] = corrMulti(motion', slowPowerHp2motion, 'Spearman');
      [rSlowPowerRSC2motionPearson, pvalSlowPowerRSC2motionPearson] = corrMulti(motion', slowPowerRSC2motion, 'Pearson');
      [rSlowPowerRSC2motionSpearman, pvalSlowPowerRSC2motionSpearman] = corrMulti(motion', slowPowerRSC2motion, 'Spearman');
      [rFastPowerS12motionPearson, pvalFastPowerS12motionPearson] = corrMulti(motion', fastPowerS12motion, 'Pearson');
      [rFastPowerS12motionSpearman, pvalFastPowerS12motionSpearman] = corrMulti(motion', fastPowerS12motion, 'Spearman');
      [rFastPowerTh2motionPearson, pvalFastPowerTh2motionPearson] = corrMulti(motion', fastPowerTh2motion, 'Pearson');
      [rFastPowerTh2motionSpearman, pvalFastPowerTh2motionSpearman] = corrMulti(motion', fastPowerTh2motion, 'Spearman');
      [rFastPowerHp2motionPearson, pvalFastPowerHp2motionPearson] = corrMulti(motion', fastPowerHp2motion, 'Pearson');
      [rFastPowerHp2motionSpearman, pvalFastPowerHp2motionSpearman] = corrMulti(motion', fastPowerHp2motion, 'Spearman');
      [rFastPowerRSC2motionPearson, pvalFastPowerRSC2motionPearson] = corrMulti(motion', fastPowerRSC2motion, 'Pearson');
      [rFastPowerRSC2motionSpearman, pvalFastPowerRSC2motionSpearman] = corrMulti(motion', fastPowerRSC2motion, 'Spearman');
    end
    [rSlowPowerS12PCsPearson, pvalSlowPowerS12PCsPearson] = corrMulti(PCs, slowPowerS1, 'Pearson');
    [rSlowPowerS12PCsSpearman, pvalSlowPowerS12PCsSpearman] = corrMulti(PCs, slowPowerS1, 'Spearman');
    [rSlowPowerTh2PCsPearson, pvalSlowPowerTh2PCsPearson] = corrMulti(PCs, slowPowerTh, 'Pearson');
    [rSlowPowerTh2PCsSpearman, pvalSlowPowerTh2PCsSpearman] = corrMulti(PCs, slowPowerTh, 'Spearman');
    [rSlowPowerHp2PCsPearson, pvalSlowPowerHp2PCsPearson] = corrMulti(PCs, slowPowerHp, 'Pearson');
    [rSlowPowerHp2PCsSpearman, pvalSlowPowerHp2PCsSpearman] = corrMulti(PCs, slowPowerHp, 'Spearman');
    [rSlowPowerRSC2PCsPearson, pvalSlowPowerRSC2PCsPearson] = corrMulti(PCs, slowPowerRSC, 'Pearson');
    [rSlowPowerRSC2PCsSpearman, pvalSlowPowerRSC2PCsSpearman] = corrMulti(PCs, slowPowerRSC, 'Spearman');
    [rFastPowerS12PCsPearson, pvalFastPowerS12PCsPearson] = corrMulti(PCs, fastPowerS1, 'Pearson');
    [rFastPowerS12PCsSpearman, pvalFastPowerS12PCsSpearman] = corrMulti(PCs, fastPowerS1, 'Spearman');
    [rFastPowerTh2PCsPearson, pvalFastPowerTh2PCsPearson] = corrMulti(PCs, fastPowerTh, 'Pearson');
    [rFastPowerTh2PCsSpearman, pvalFastPowerTh2PCsSpearman] = corrMulti(PCs, fastPowerTh, 'Spearman');
    [rFastPowerHp2PCsPearson, pvalFastPowerHp2PCsPearson] = corrMulti(PCs, fastPowerHp, 'Pearson');
    [rFastPowerHp2PCsSpearman, pvalFastPowerHp2PCsSpearman] = corrMulti(PCs, fastPowerHp, 'Spearman');
    [rFastPowerRSC2PCsPearson, pvalFastPowerRSC2PCsPearson] = corrMulti(PCs, fastPowerRSC, 'Pearson');
    [rFastPowerRSC2PCsSpearman, pvalFastPowerRSC2PCsSpearman] = corrMulti(PCs, fastPowerRSC, 'Spearman');
    
    
    % CALCULATE EXPLAINED VARIANCE FOR EACH AREA
    areasOI = unique(areas);
    for iArea = 1:numel(areasOI)
      dataOI = areas;
      dataOI(dataOI ~= areasOI(iArea)) = 0;
      dataOI = logical(dataOI);
      pcaCoef_area = pcaCoef(dataOI);
      explained_area = 100*mean(individualVarExplained(:,dataOI), 2)';
      
      % Calculate a number of variables needed to explain proportions of variance for the area
      [explained25_area, explained50_area, explained75_area, explained95_area,...
        explainedThird_area, explainedTwoThirds_area] = explainedProportions(explained_area);
      
      
      % CORRELATE PCs WITH LFP
      [rPCs2lfpPearson_area, pvalPCs2lfpPearson_area] = corrMulti(LFPs(iArea,:), PCs, 'Pearson');
      [rPCs2lfpSpearman_area, pvalPCs2lfpSpearman_area] = corrMulti(LFPs(iArea,:), PCs, 'Spearman');
    
    
      % CORRELATE PCs WITH PR
      [rPCs2prPearson_area, pvalPCs2prPearson_area] = corrMulti(PRs(iArea,:), PCs, 'Pearson');
      [rPCs2prSpearman_area, pvalPCs2prSpearman_area] = corrMulti(PRs(iArea,:), PCs, 'Spearman');
    
    
      % SAVE DATA
      % PCA data for all areas
      pcaData.pcaCoef = pcaCoef;
      pcaData.explained = explained;
      pcaData.nPCs = nPCs;
      pcaData.prob = prob;
      pcaData.explained25 = explained25;
      pcaData.explained50 = explained50;
      pcaData.explained75 = explained75;
      pcaData.explained95 = explained95;
      pcaData.explainedThird = explainedThird;
      pcaData.explainedTwoThirds = explainedTwoThirds;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson = rPCs2pupilAreaPearson;
        pcaData.pvalPCs2pupilAreaPearson = pvalPCs2pupilAreaPearson;
        pcaData.rPCs2pupilAreaSpearman = rPCs2pupilAreaSpearman;
        pcaData.pvalPCs2pupilAreaSpearman = pvalPCs2pupilAreaSpearman;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson = rPCs2motionPearson;
        pcaData.pvalPCs2motionPearson = pvalPCs2motionPearson;
        pcaData.rPCs2motionSpearman = rPCs2motionSpearman;
        pcaData.pvalPCs2motionSpearman = pvalPCs2motionSpearman;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area = pcaCoef_area;
      pcaData.explained_area = explained_area;
      pcaData.explained25_area = explained25_area;
      pcaData.explained50_area = explained50_area;
      pcaData.explained75_area = explained75_area;
      pcaData.explained95_area = explained95_area;
      pcaData.explainedThird_area = explainedThird_area;
      pcaData.explainedTwoThirds_area = explainedTwoThirds_area;
      pcaData.rPCs2lfpPearson_area = rPCs2lfpPearson_area;
      pcaData.pvalPCs2lfpPearson_area = pvalPCs2lfpPearson_area;
      pcaData.rPCs2lfpSpearman_area = rPCs2lfpSpearman_area;
      pcaData.pvalPCs2lfpSpearman_area = pvalPCs2lfpSpearman_area;
      pcaData.rPCs2prPearson_area = rPCs2prPearson_area;
      pcaData.pvalPCs2prPearson_area = pvalPCs2prPearson_area;
      pcaData.rPCs2prSpearman_area = rPCs2prSpearman_area;
      pcaData.pvalPCs2prSpearman_area = pvalPCs2prSpearman_area;
      if ~isempty(pupilArea)
        if areasOI(iArea) == 1
          pcaData.rSlowPower2pupilAreaPearson = rSlowPowerS12pupilAreaPearson;
          pcaData.pvalSlowPower2pupilAreaPearson = pvalSlowPowerS12pupilAreaPearson;
          pcaData.rFastPower2pupilAreaPearson = rFastPowerS12pupilAreaPearson;
          pcaData.pvalFastPower2pupilAreaPearson = pvalFastPowerS12pupilAreaPearson;
        elseif areasOI(iArea) == 2
          pcaData.rSlowPower2pupilAreaPearson = rSlowPowerTh2pupilAreaPearson;
          pcaData.pvalSlowPower2pupilAreaPearson = pvalSlowPowerTh2pupilAreaPearson;
          pcaData.rFastPower2pupilAreaPearson = rFastPowerTh2pupilAreaPearson;
          pcaData.pvalFastPower2pupilAreaPearson = pvalFastPowerTh2pupilAreaPearson;
        elseif areasOI(iArea) == 5
          pcaData.rSlowPower2pupilAreaPearson = rSlowPowerHp2pupilAreaPearson;
          pcaData.pvalSlowPower2pupilAreaPearson = pvalSlowPowerHp2pupilAreaPearson;
          pcaData.rFastPower2pupilAreaPearson = rFastPowerHp2pupilAreaPearson;
          pcaData.pvalFastPower2pupilAreaPearson = pvalFastPowerHp2pupilAreaPearson;
        elseif areasOI(iArea) == 7
          pcaData.rSlowPower2pupilAreaPearson = rSlowPowerRSC2pupilAreaPearson;
          pcaData.pvalSlowPower2pupilAreaPearson = pvalSlowPowerRSC2pupilAreaPearson;
          pcaData.rFastPower2pupilAreaPearson = rFastPowerRSC2pupilAreaPearson;
          pcaData.pvalFastPower2pupilAreaPearson = pvalFastPowerRSC2pupilAreaPearson; 
        end
      end
      if ~isempty(motion)
        if areasOI(iArea) == 1
          pcaData.rSlowPower2motionPearson = rSlowPowerS12motionPearson;
          pcaData.pvalSlowPower2motionPearson = pvalSlowPowerS12motionPearson;
          pcaData.rFastPower2motionPearson = rFastPowerS12motionPearson;
          pcaData.pvalFastPower2motionPearson = pvalFastPowerS12motionPearson;
        elseif areasOI(iArea) == 2
          pcaData.rSlowPower2motionPearson = rSlowPowerTh2motionPearson;
          pcaData.pvalSlowPower2motionPearson = pvalSlowPowerTh2motionPearson;
          pcaData.rFastPower2motionPearson = rFastPowerTh2motionPearson;
          pcaData.pvalFastPower2motionPearson = pvalFastPowerTh2motionPearson;
        elseif areasOI(iArea) == 5
          pcaData.rSlowPower2motionPearson = rSlowPowerHp2motionPearson;
          pcaData.pvalSlowPower2motionPearson = pvalSlowPowerHp2motionPearson;
          pcaData.rFastPower2motionPearson = rFastPowerHp2motionPearson;
          pcaData.pvalFastPower2motionPearson = pvalFastPowerHp2motionPearson;
        elseif areasOI(iArea) == 7
          pcaData.rSlowPower2motionPearson = rSlowPowerRSC2motionPearson;
          pcaData.pvalSlowPower2motionPearson = pvalSlowPowerRSC2motionPearson;
          pcaData.rFastPower2motionPearson = rFastPowerRSC2motionPearson;
          pcaData.pvalFastPower2motionPearson = pvalFastPowerRSC2motionPearson; 
        end
      end
      if areasOI(iArea) == 1
        pcaData.rSlowPower2PCsPearson = rSlowPowerS12PCsPearson;
        pcaData.pvalSlowPower2PCsPearson = pvalSlowPowerS12PCsPearson;
        pcaData.rFastPower2PCsPearson = rFastPowerS12PCsPearson;
        pcaData.pvalFastPower2PCsPearson = pvalFastPowerS12PCsPearson;
      elseif areasOI(iArea) == 2
        pcaData.rSlowPower2PCsPearson = rSlowPowerTh2PCsPearson;
        pcaData.pvalSlowPower2PCsPearson = pvalSlowPowerTh2PCsPearson;
        pcaData.rFastPower2PCsPearson = rFastPowerTh2PCsPearson;
        pcaData.pvalFastPower2PCsPearson = pvalFastPowerTh2PCsPearson;
      elseif areasOI(iArea) == 5
        pcaData.rSlowPower2PCsPearson = rSlowPowerHp2PCsPearson;
        pcaData.pvalSlowPower2PCsPearson = pvalSlowPowerHp2PCsPearson;
        pcaData.rFastPower2PCsPearson = rFastPowerHp2PCsPearson;
        pcaData.pvalFastPower2PCsPearson = pvalFastPowerHp2PCsPearson;
      elseif areasOI(iArea) == 7
        pcaData.rSlowPower2PCsPearson = rSlowPowerRSC2PCsPearson;
        pcaData.pvalSlowPower2PCsPearson = pvalSlowPowerRSC2PCsPearson;
        pcaData.rFastPower2PCsPearson = rFastPowerRSC2PCsPearson;
        pcaData.pvalFastPower2PCsPearson = pvalFastPowerRSC2PCsPearson;
      end
      
      dataString = ['dataStruct.seriesData.' entryNames{iArea} '.pcaData2 = pcaData;'];
      eval(dataString);
      if intermediateSaving
        save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
      end
      
      
      % PLOT THE DATA
      figs{iCh} = figure('Visible', 'on'); %#ok<*SAGROW>
      plot(PCs2pupilAreaTimes, pupilArea'./mean(pupilArea), 'LineWidth',2); hold on;
      plot(PCs2motionTimes, motion'./mean(motion), 'LineWidth',2); hold on;
      lfpTimes = 1/params.srData:1/params.srData:numel(LFP_S1)/params.srData;
      if areasOI(iArea) == 1
        plot(lfpTimes, LFP_S1./mean(LFP_S1))
        plot(lfpTimes, slowPowerS1./mean(slowPowerS1))
        plot(lfpTimes, fastPowerS1./mean(fastPowerS1))
      elseif areasOI(iArea) == 2
        plot(lfpTimes, LFP_Th./mean(LFP_Th))
        plot(lfpTimes, slowPowerTh./mean(slowPowerTh))
        plot(lfpTimes, fastPowerTh./mean(fastPowerTh))
      elseif areasOI(iArea) == 5
        plot(lfpTimes, LFP_Hp./mean(LFP_Hp))
        plot(lfpTimes, slowPowerHp./mean(slowPowerHp))
        plot(lfpTimes, fastPowerHp./mean(fastPowerHp))
      elseif areasOI(iArea) == 7
        plot(lfpTimes, LFP_RSC./mean(LFP_RSC))
        plot(lfpTimes, slowPowerRSC./mean(slowPowerRSC))
        plot(lfpTimes, fastPowerRSC./mean(fastPowerRSC))
      end
      plot(lfpTimes, PCs(1,:)./mean(PCs(1,:)))
      plot(lfpTimes, PCs(2,:)./mean(PCs(2,:)))
      plot(lfpTimes, PCs(3,:)./mean(PCs(3,:)))
      legend('Pupil area','total movement','LFP','Slow power','Fast power','PC1','PC2','PC3')
      xlabel('Time (s)')
      ylabel('Normalised signal')
      figsubdirname = seriesFromEntry(entryNames{iArea});
      if ~exist(figsubdirname,'dir')
          mkdir(figsubdirname)
      end
      figName = [figsubdirname filesep entryName '_LFP_variousSignals_v_PCs_' 'ch' num2str(chOIDB(iCh))];
      set(figs{iCh}, 'Name',figName);
      hgsave(figs{iCh}, figName);
      close(figs{iCh});
    end
  end
end

if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbStart dbEnd


function [LFPbandPower, wtSpectrogram, fSpectrogram, interpDat] = lfpFileLoad(...
  baseFilename, minF, chunkSize, chN, chOI, LFPbands, srRecording, srInterp, medianSubtracted, deleteChans, spectrogram)

% File loading parameters
fid = fopen(baseFilename, 'r');
d = dir(baseFilename);
nSampsTotal = d.bytes/chN/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);

% Load the file
chunkInd = 1;
LFPbandPower = {};
wtSpectrogram = {};
interpDat = {};
while 1
  fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
  dat = fread(fid, [chN chunkSize], '*int16');
  
  % Subtract median and/or delete channels if necessary
  if medianSubtracted
    if deleteChans
      chans2include = ones(1,size(dat,1));
      chans2include(deleteChans) = zeros(1,numel(deleteChans));
      chm = zeros(size(dat,1),1);
      chm(logical(chans2include)) = median(dat(logical(chans2include),:),2);
      dat = bsxfun(@minus, dat, int16(chm)); % subtract median of each channel
      tm = int16(median(dat(logical(chans2include),:),1));
    else
      chm = median(dat,2);
      dat = bsxfun(@minus, dat, chm); % subtract median of each channel
      tm = median(dat,1);
    end
    dat = bsxfun(@minus, dat, tm); % subtract median of each time point
  end
  
  if ~isempty(dat)
    
    % Interpolate data
    originalTimes = 1/srRecording:1/srRecording:size(dat,2)/srRecording;
    interpTimes = 1/srInterp:1/srInterp:size(dat,2)/srRecording;
    interpLFP = interp1(originalTimes, double(dat'), interpTimes)';
    for iCh = 1:numel(chOI)
      
      % Calculate wavelet transforms
      fb = cwtfilterbank('SignalLength',size(interpLFP,2),'SamplingFrequency',srInterp,...
        'FrequencyLimits',[minF 200],'WaveletParameters',[3 16],'VoicesPerOctave',20);
      [wt,f] = cwt(interpLFP(chOI(iCh),:),'FilterBank',fb); % Continuous wavelet transform
      LFPbandPowerChunk = zeros(numel(LFPbands),size(wt,2));
      for iBand = 1:numel(LFPbands)
        for iLim = 1:2
          [~, LFPbandLimits{iBand}(iLim)] = min(abs(f-LFPbands{iBand}(iLim))); %#ok<*AGROW>
        end
        LFPbandPowerChunk(iBand,:) = sum(abs(wt(LFPbandLimits{iBand}(2):LFPbandLimits{iBand}(1),:)).^2,1);
      end
      if chunkInd == 1
        LFPbandPower{iCh} = LFPbandPowerChunk;
        interpDat{iCh} = interpLFP(chOI(iCh),:);
      else
        LFPbandPower{iCh} = [LFPbandPower{iCh} LFPbandPowerChunk]; %#ok<*AGROW>
        interpDat{iCh} = [interpDat{iCh} interpLFP(chOI(iCh),:)];
      end
      
      % Spectrogram
      if spectrogram
        fb = cwtfilterbank('SignalLength',size(interpLFP,2),'SamplingFrequency',srInterp,...
          'FrequencyLimits',[minF 200],'WaveletParameters',[3 16],'VoicesPerOctave',10);
        [wtSpectrogramChunk,fSpectrogram] = cwt(interpLFP(chOI(iCh),:),'FilterBank',fb); % Continuous wavelet transform
        %helperCWTTimeFreqPlot(wt,interpTimes,f,'surf','Spectrogram for CA1 channel','Seconds','Hz');
        %set(gca, 'YScale', 'log')
        if chunkInd == 1
          wtSpectrogram{iCh} = abs(wtSpectrogramChunk).^2;
        else
          wtSpectrogram{iCh} = [wtSpectrogram{iCh} abs(wtSpectrogramChunk).^2];
        end
      else
        fSpectrogram = [];
      end
      
    end
  else
    break
  end
  chunkInd = chunkInd+1;
end
end

function [rippleRate, meanRippleRate] = rippleRateCalculator(rippleBandPower, rippleDuration, sdGaussian, wGaussian, chOI, bandNames, srRippleRate)

for iCh = 1:numel(chOI)
  % CALCULATE RIPPLE RATE AS IN MCGINLEY ET AL. (2015)
  % Descriptive measures
  for iName = 1:numel(bandNames)
    if strcmpi(bandNames{iName}, 'ripples/uf')
      iRipples = iName;
      break
    end
  end
  %       figure; plot(timeInit, LFPbandPower{iCh}(iRipples,:)); hold on
  medianPower = median(rippleBandPower{iCh},'omitnan');
  stdPower = std(rippleBandPower{iCh},'omitnan');
  %       stdPowerTop = medianPower+1.96*stdPower;
  %       stdPowerBottom = medianPower-1.96*stdPower;
  stdPowerTop = medianPower+5*stdPower;
  stdPowerBottom = medianPower-5*stdPower; %#ok<*NASGU>
  %plot(medianPower*ones(1,size(LFPbandPower{iCh},2)))
  %plot(stdPowerTop*ones(1,size(LFPbandPower{iCh},2)))
  %plot(stdPowerBottom*ones(1,size(LFPbandPower{iCh},2)))
  
  % Detect initial ripples
  [ripplePowerPeaks, locations] = findpeaks(rippleBandPower{iCh});
  locations = locations(ripplePowerPeaks > stdPowerTop);
  totalRipplesInit = numel(locations);
  %       ripplePowerPeaks = ripplePowerPeaks(ripplePowerPeaks > stdPowerTop);
  %       plot(timeInit(locations),ripplePowerPeaks, 'r.', 'MarkerSize',5)
  
  % Mark ripple initiation and termination
  ripplePowerPeaksInit = zeros(1,size(rippleBandPower{iCh},2));
  ripplePowerPeaksInit(locations) = 1;
  ripplePowerPeaks = ripplePowerPeaksInit;
  for dtRipple = 1:round(rippleDuration/2)
    ripplePowerPeaks = ripplePowerPeaks + [ripplePowerPeaksInit(1+dtRipple:end) zeros(1,dtRipple)];
    ripplePowerPeaks = ripplePowerPeaks + [zeros(1,dtRipple) ripplePowerPeaksInit(1:end-dtRipple)];
  end
  ripplePowerPeaks(ripplePowerPeaks > 0) = 1;
  locations = 1:size(rippleBandPower{iCh},2);
  locations = locations(logical(ripplePowerPeaks));
  ripplePowerPeaksInit = ripplePowerPeaks;
  %       ripplePowerPeaks = ripplePowerPeaks(ripplePowerPeaks > 0);
  %       plot(timeInit(locations),ripplePowerPeaks, 'k.', 'MarkerSize',10); hold on
  
  % Total ripple event count
  ripplePowerPeaks2 = findpeaks(ripplePowerPeaksInit);
  totalRipples = numel(ripplePowerPeaks2);
  
  % Convolve with Gaussian
  w = gausswin(wGaussian*sdGaussian*srRippleRate, (wGaussian*sdGaussian*srRippleRate-1)/(2*sdGaussian*srRippleRate));
  w = w/sum(w);
  ripplePowerPeaks = filtfilt(w,1,ripplePowerPeaksInit);
  ripplePowerPeaks = ripplePowerPeaks/sum(ripplePowerPeaks)*sum(rippleBandPower{iCh});
  %plot(ripplePowerPeaks)
  %rippleRateInit = ripplePowerPeaks/mean(ripplePowerPeaks)*((totalRipplesInit/size(rippleBandPower{iCh},2))*srRippleRate); % Hz
  rippleRate{iCh} = ripplePowerPeaks/mean(ripplePowerPeaks)*((totalRipples/size(rippleBandPower{iCh},2))*srRippleRate); % Hz
  meanRippleRate{iCh} = mean(rippleRate{iCh});
  %timeInit = 1/srRippleRate:1/srRippleRate:size(rippleBandPower{iCh},2)/srRippleRate;
  %       plot(timeInit,rippleRate{iCh});
  
  % Down-sample again
  %interpTimesFinal = 1/dssrLFPfinal:1/dssrLFPfinal:size(rippleBandPower{iCh},2)/srRippleRate;
  %rippleRate{iCh} = interp1(timeInit, rippleRate{iCh}, interpTimesFinal)';
  %timeFinal = 1/dssrLFPfinal:1/dssrLFPfinal:size(rippleBandPower{iCh},2)/srRippleRate;
  %plot(timeFinal,rippleRate{iCh}); hold off
end
end

function [explained25, explained50, explained75, explained95, explainedThird, explainedTwoThirds] = explainedProportions(explained)

explainedCumulitive = cumsum(explained);
explainedCumulitive25 = explainedCumulitive - 25;
explained25 = find(explainedCumulitive25 >= 0);
explainedCumulitive50 = explainedCumulitive - 50;
explained50 = find(explainedCumulitive50 >= 0);
explainedCumulitive75 = explainedCumulitive - 75;
explained75 = find(explainedCumulitive75 >= 0);
explainedCumulitive95 = explainedCumulitive - 95;
explained95 = find(explainedCumulitive95 >= 0);
explainedCumulitiveThird = explainedCumulitive - 100/3;
explainedThird = find(explainedCumulitiveThird >= 0);
explainedCumulitiveTwoThirds = explainedCumulitive - 200/3;
explainedTwoThirds = find(explainedCumulitiveTwoThirds >= 0);
end