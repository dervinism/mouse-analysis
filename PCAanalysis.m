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
  entryName = dbStruct.db(dbCount).entryName;
  chN = dbStruct.db(dbCount).chN;
  
  
  % TEST FOR RIPPLES
  ripplesExist = rippleTest(fnsData{dbCount}, fullSeries);
  
  
  % IDENTIFY AREA
  seriesName = seriesFromEntry(entryName);
  area = determineAreaFromSeries(seriesName);
  if ~area
    continue
  end
  
  
  % LOAD EXISTING LFP ANALYSIS DATA
  if isfield(dbStruct,'lfpPowerData')
    if medianSubtracted % It's ok to use the existing ripple rate because it is based on median subtracted LFP
      rippleRate = dbStruct.lfpPowerData.rippleRate;
    end
    theta2deltaRatio = dbStruct.lfpPowerData.theta2deltaRatio;
    chOIDB = dbStruct.lfpPowerData.chOIDB;
    dssrLFPinit = dbStruct.lfpPowerData.dssrLFPinit;
    lfpTimes = dbStruct.lfpPowerData.lfpTimes;
    rippleDuration = dbStruct.lfpPowerData.rippleDuration;
    sdGaussian = dbStruct.lfpPowerData.sdGaussian;
    wGaussian = dbStruct.lfpPowerData.wGaussian;
  else
    %error('No LFP analysis data exist. Please run lfpLoad.m');
    continue
  end
  
  
  % RESAMPLE EXISTING LFP DATA IF NECESSARY
  if params.srData ~= dssrLFPinit
    interpTimes = 1/params.srData:1/params.srData:lfpTimes(end);
    for iCh = 1:numel(chOIDB)
      if (area == 5 || area == 10) && medianSubtracted
        rippleRate{iCh} = interp1(lfpTimes, rippleRate{iCh}, interpTimes)';
      end
      theta2deltaRatio{iCh} = interp1(lfpTimes, theta2deltaRatio{iCh}, interpTimes)';
    end
  else
    interpTimes = lfpTimes;
  end
  

  % RECALCULATE RIPPLE RATE IF EXISTING ONE IS NOT BASED ON MEDIAN-SUBTRACTED LFP
  if (area == 5 || area == 10) && ~medianSubtracted && ripplesExist
    srRippleRate = 1000;
    if strcmpi(probe, 'Neuropixels')
      LFPbandPower = lfpFileLoad(baseFilename, 1, chunkSize, chN, chOIDB, LFPbands, params.srRecordingLFP, srRippleRate, true, deleteChans, false);
    else
      LFPbandPower = lfpFileLoad(baseFilename, 1, chunkSize, chN, chOIDB, LFPbands, params.srRecording, srRippleRate, true, deleteChans, false);
    end
    for iCh = 1:numel(chOIDB)
      rippleBandPower{iCh} = LFPbandPower{iCh}(end,:); %#ok<*SAGROW>
    end
    rippleRate = rippleRateCalculator(rippleBandPower, rippleDuration, sdGaussian, wGaussian, chOIDB, bandNames, srRippleRate);
    rippleRateTimes = 1/srRippleRate:1/srRippleRate:numel(rippleRate{iCh})/srRippleRate;
    for iCh = 1:numel(chOIDB)
      rippleRate{iCh} = interp1(rippleRateTimes, rippleRate{iCh}, interpTimes)';
    end
  else
    for iCh = 1:numel(chOIDB)
      rippleRate{iCh} = zeros(size(theta2deltaRatio{iCh}));
    end
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
  neuralActivityBasic = {};
  neuralActivityHp = {};
  for iCh = 1:numel(chOIDB)
    if isempty(PR)
      dataDuration = min([size(LFPbandPower{iCh},2) numel(rippleRate{iCh}) numel(theta2deltaRatio{iCh}) numel(LFP{iCh})]);
      PR = zeros(1,dataDuration);
    else
      dataDuration = min([size(LFPbandPower{iCh},2) numel(rippleRate{iCh}) numel(theta2deltaRatio{iCh}) numel(LFP{iCh}) numel(PR)]);
      PR = PR(1:dataDuration);
    end
    LFPbandPower{iCh} = LFPbandPower{iCh}(:,1:dataDuration);
    theta2deltaRatio{iCh} = torow(theta2deltaRatio{iCh}(1:dataDuration));
    rippleRate{iCh} = torow(rippleRate{iCh}(1:dataDuration));
    LFP{iCh} = torow(LFP{iCh}(1:dataDuration));
    interpTimes = torow(interpTimes(1:dataDuration));
    neuralActivityBasic{iCh} = [LFPbandPower{iCh}; theta2deltaRatio{iCh}; LFP{iCh}; PR];
    neuralActivityHp{iCh} = [LFPbandPower{iCh}; rippleRate{iCh}; theta2deltaRatio{iCh}; LFP{iCh}; PR];
    neuralActivityBasic{iCh}(isnan(neuralActivityBasic{iCh})) = 0;
    neuralActivityHp{iCh}(isnan(neuralActivityHp{iCh})) = 0;
  end
  
  
  % SUBTRACT THE MEAN AND Z-SCORE THE DATA
  for iCh = 1:numel(chOIDB)
    neuralActivityBasic{iCh} = neuralActivityBasic{iCh} - mean(neuralActivityBasic{iCh},2);
    neuralActivityBasic{iCh} = zscore(neuralActivityBasic{iCh}')';
    neuralActivityHp{iCh} = neuralActivityHp{iCh} - mean(neuralActivityHp{iCh},2);
    neuralActivityHp{iCh} = zscore(neuralActivityHp{iCh}')';
  end
  
  
  % RUN PCA
  pcaCoef = {};
  explained = {};
  PCs = {};
  nPCs = {};
  prob = {};
  pcaCoefHp = {};
  explainedHp = {};
  PCsHp = {};
  nPCsHp = {};
  probHp = {};
  explained25 = {};
  explained50 = {};
  explained75 = {};
  explained95 = {};
  explainedThird = {};
  explainedTwoThirds = {};
  explainedHp25 = {};
  explainedHp50 = {};
  explainedHp75 = {};
  explainedHp95 = {};
  explainedHpThird = {};
  explainedHpTwoThirds = {};
  for iCh = 1:numel(chOIDB)
    [pcaCoef{iCh}, explained{iCh}, PCs{iCh}, nPCs{iCh}, prob{iCh}] = pcaGeneric(neuralActivityBasic{iCh});
    [pcaCoefHp{iCh}, explainedHp{iCh}, PCsHp{iCh}, nPCsHp{iCh}, probHp{iCh}] = pcaGeneric(neuralActivityHp{iCh});
    
    % Calculate a number of variables needed to explain proportions of variance
    [explained25{iCh}, explained50{iCh}, explained75{iCh}, explained95{iCh}, explainedThird{iCh},...
      explainedTwoThirds{iCh}] = explainedProportions(explained{iCh});
    [explainedHp25{iCh}, explainedHp50{iCh}, explainedHp75{iCh}, explainedHp95{iCh}, explainedHpThird{iCh},...
      explainedHpTwoThirds{iCh}] = explainedProportions(explainedHp{iCh});
  end
  
  
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
  
  % Truncate PCs to match arousal measures
  PCs2pupilArea = {};
  PCsHp2pupilArea = {};
  if ~isempty(pupilArea)
    [~, iPCs2pupilAreaStart] = min(abs(interpTimes - pupilAreaTimes(1)));
    [~, iPCs2pupilAreaEnd] = min(abs(interpTimes - pupilAreaTimes(end)));
    PCs2pupilAreaTimes = interpTimes(iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
    for iCh = 1:numel(chOIDB)
      PCs2pupilArea{iCh} = PCs{iCh}(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCsHp2pupilArea{iCh} = PCsHp{iCh}(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
    end
  end
  PCs2motion = {};
  PCsHp2motion = {};
  if ~isempty(motion)
    [~, iPCs2motionStart] = min(abs(interpTimes - motionTimes(1)));
    [~, iPCs2motionEnd] = min(abs(interpTimes - motionTimes(end)));
    PCs2motionTimes = interpTimes(iPCs2motionStart:iPCs2motionEnd);
    for iCh = 1:numel(chOIDB)
      PCs2motion{iCh} = PCs{iCh}(:,iPCs2motionStart:iPCs2motionEnd);
      PCsHp2motion{iCh} = PCsHp{iCh}(:,iPCs2motionStart:iPCs2motionEnd);
    end
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
  rPCs2pupilAreaPearson = {};
  pvalPCs2pupilAreaPearson = {};
  rPCsHp2pupilAreaPearson = {};
  pvalPCsHp2pupilAreaPearson = {};
  rPCs2pupilAreaSpearman = {};
  pvalPCs2pupilAreaSpearman = {};
  rPCsHp2pupilAreaSpearman = {};
  pvalPCsHp2pupilAreaSpearman = {};
  rPCs2motionPearson = {};
  pvalPCs2motionPearson = {};
  rPCsHp2motionPearson = {};
  pvalPCsHp2motionPearson = {};
  rPCs2motionSpearman = {};
  pvalPCs2motionSpearman = {};
  rPCsHp2motionSpearman = {};
  pvalPCsHp2motionSpearman = {};
  for iCh = 1:numel(chOIDB)
    if ~isempty(pupilArea)
      [rPCs2pupilAreaPearson{iCh}, pvalPCs2pupilAreaPearson{iCh}] = corrMulti(pupilArea', PCs2pupilArea{iCh}, 'Pearson');
      [rPCsHp2pupilAreaPearson{iCh}, pvalPCsHp2pupilAreaPearson{iCh}] = corrMulti(pupilArea', PCs2pupilArea{iCh}, 'Pearson');
      [rPCs2pupilAreaSpearman{iCh}, pvalPCs2pupilAreaSpearman{iCh}] = corrMulti(pupilArea', PCs2pupilArea{iCh}, 'Spearman');
      [rPCsHp2pupilAreaSpearman{iCh}, pvalPCsHp2pupilAreaSpearman{iCh}] = corrMulti(pupilArea', PCs2pupilArea{iCh}, 'Spearman');
    end
    if ~isempty(motion)
      [rPCs2motionPearson{iCh}, pvalPCs2motionPearson{iCh}] = corrMulti(motion', PCs2motion{iCh}, 'Pearson');
      [rPCsHp2motionPearson{iCh}, pvalPCsHp2motionPearson{iCh}] = corrMulti(motion', PCs2motion{iCh}, 'Pearson');
      [rPCs2motionSpearman{iCh}, pvalPCs2motionSpearman{iCh}] = corrMulti(motion', PCs2motion{iCh}, 'Spearman');
      [rPCsHp2motionSpearman{iCh}, pvalPCsHp2motionSpearman{iCh}] = corrMulti(motion', PCs2motion{iCh}, 'Spearman');
    end
  end
  
  
  % CORRELATE PCs WITH LFP
  rPCs2lfpPearson = {};
  pvalPCs2lfpPearson = {};
  rPCsHp2lfpPearson = {};
  pvalPCsHp2lfpPearson = {};
  rPCs2lfpSpearman = {};
  pvalPCs2lfpSpearman = {};
  rPCsHp2lfpSpearman = {};
  pvalPCsHp2lfpSpearman = {};
  for iCh = 1:numel(chOIDB)
    [rPCs2lfpPearson{iCh}, pvalPCs2lfpPearson{iCh}] = corrMulti(LFP{iCh}, PCs{iCh}, 'Pearson');
    [rPCsHp2lfpPearson{iCh}, pvalPCsHp2lfpPearson{iCh}] = corrMulti(LFP{iCh}, PCsHp{iCh}, 'Pearson');
    [rPCs2lfpSpearman{iCh}, pvalPCs2lfpSpearman{iCh}] = corrMulti(LFP{iCh}, PCs{iCh}, 'Spearman');
    [rPCsHp2lfpSpearman{iCh}, pvalPCsHp2lfpSpearman{iCh}] = corrMulti(LFP{iCh}, PCsHp{iCh}, 'Spearman');
  end
  
  
  % CORRELATE PCs WITH PR
  rPCs2prPearson = {};
  pvalPCs2prPearson = {};
  rPCsHp2prPearson = {};
  pvalPCsHp2prPearson = {};
  rPCs2prSpearman = {};
  pvalPCs2prSpearman = {};
  rPCsHp2prSpearman = {};
  pvalPCsHp2prSpearman = {};
  for iCh = 1:numel(chOIDB)
    [rPCs2prPearson{iCh}, pvalPCs2prPearson{iCh}] = corrMulti(PR, PCs{iCh}, 'Pearson');
    [rPCsHp2prPearson{iCh}, pvalPCsHp2prPearson{iCh}] = corrMulti(PR, PCsHp{iCh}, 'Pearson');
    [rPCs2prSpearman{iCh}, pvalPCs2prSpearman{iCh}] = corrMulti(PR, PCs{iCh}, 'Spearman');
    [rPCsHp2prSpearman{iCh}, pvalPCsHp2prSpearman{iCh}] = corrMulti(PR, PCsHp{iCh}, 'Spearman');
  end
  
  
  % SAVE DATA
  % PCA data for non-hippocampal areas
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
  else
    pcaData.rPCs2pupilAreaPearson = {};
    pcaData.pvalPCs2pupilAreaPearson = {};
    pcaData.rPCs2pupilAreaSpearman = {};
    pcaData.pvalPCs2pupilAreaSpearman = {};
  end
  if ~isempty(motion)
    pcaData.rPCs2motionPearson = rPCs2motionPearson;
    pcaData.pvalPCs2motionPearson = pvalPCs2motionPearson;
    pcaData.rPCs2motionSpearman = rPCs2motionSpearman;
    pcaData.pvalPCs2motionSpearman = pvalPCs2motionSpearman;
  else
    pcaData.rPCs2motionPearson = {};
    pcaData.pvalPCs2motionPearson = {};
    pcaData.rPCs2motionSpearman = {};
    pcaData.pvalPCs2motionSpearman = {};
  end
  pcaData.rPCs2lfpPearson = rPCs2lfpPearson;
  pcaData.pvalPCs2lfpPearson = pvalPCs2lfpPearson;
  pcaData.rPCs2lfpSpearman = rPCs2lfpSpearman;
  pcaData.pvalPCs2lfpSpearman = pvalPCs2lfpSpearman;
  pcaData.rPCs2prPearson = rPCs2prPearson;
  pcaData.pvalPCs2prPearson = pvalPCs2prPearson;
  pcaData.rPCs2prSpearman = rPCs2prSpearman;
  pcaData.pvalPCs2prSpearman = pvalPCs2prSpearman;
  
  % Hippocampal PCA data (including ripple rate)
  pcaData.pcaCoefHp = pcaCoefHp;
  pcaData.explainedHp = explainedHp;
  pcaData.nPCsHp = nPCsHp;
  pcaData.probHp = probHp;
  pcaData.explainedHp25 = explainedHp25;
  pcaData.explainedHp50 = explainedHp50;
  pcaData.explainedHp75 = explainedHp75;
  pcaData.explainedHp95 = explainedHp95;
  pcaData.explainedHpThird = explainedHpThird;
  pcaData.explainedHpTwoThirds = explainedHpTwoThirds;
  if ~isempty(pupilArea)
    pcaData.rPCsHp2pupilAreaPearson = rPCsHp2pupilAreaPearson;
    pcaData.pvalPCsHp2pupilAreaPearson = pvalPCsHp2pupilAreaPearson;
    pcaData.rPCsHp2pupilAreaSpearman = rPCsHp2pupilAreaSpearman;
    pcaData.pvalPCsHp2pupilAreaSpearman = pvalPCsHp2pupilAreaSpearman;
  else
    pcaData.rPCsHp2pupilAreaPearson = {};
    pcaData.pvalPCsHp2pupilAreaPearson = {};
    pcaData.rPCsHp2pupilAreaSpearman = {};
    pcaData.pvalPCsHp2pupilAreaSpearman = {};
  end
  if ~isempty(motion)
    pcaData.rPCsHp2motionPearson = rPCsHp2motionPearson;
    pcaData.pvalPCsHp2motionPearson = pvalPCsHp2motionPearson;
    pcaData.rPCsHp2motionSpearman = rPCsHp2motionSpearman;
    pcaData.pvalPCsHp2motionSpearman = pvalPCsHp2motionSpearman;
  else
    pcaData.rPCsHp2motionPearson = {};
    pcaData.pvalPCsHp2motionPearson = {};
    pcaData.rPCsHp2motionSpearman = {};
    pcaData.pvalPCsHp2motionSpearman = {};
  end
  pcaData.rPCsHp2lfpPearson = rPCsHp2lfpPearson;
  pcaData.pvalPCsHp2lfpPearson = pvalPCsHp2lfpPearson;
  pcaData.rPCsHp2lfpSpearman = rPCsHp2lfpSpearman;
  pcaData.pvalPCsHp2lfpSpearman = pvalPCsHp2lfpSpearman;
  pcaData.rPCsHp2prPearson = rPCsHp2prPearson;
  pcaData.pvalPCsHp2prPearson = pvalPCsHp2prPearson;
  pcaData.rPCsHp2prSpearman = rPCsHp2prSpearman;
  pcaData.pvalPCsHp2prSpearman = pvalPCsHp2prSpearman;
  
  dataString = ['dataStruct.seriesData.' entryName '.pcaData = pcaData;'];
  eval(dataString);
  if intermediateSaving
    save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
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