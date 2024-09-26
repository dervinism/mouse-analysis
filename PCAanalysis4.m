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
  elseif numel(seriesName) == 15 &&...
      (strcmpi(seriesName(15), '4') || strcmpi(seriesName(15), '5'))
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
    narrowSlowBands = [];
    deltaBands = [];
    alphaBands = [];
    thetaBands = [];
    highSlowBands = [];
    slowBands = [];
    betaBands = [];
    slowGammaBands = [];
    fastGammaBands = [];
    gammaBands = [];
    ultraFastBands = [];
    ultraFastBandsRipples = [];
    fastBands = [];
    LFPs = [];
    PRs = [];
    areasNarrow = [];
    areasHighSlow = [];
    areasSlow = [];
    areasFast = [];
    areasGamma = [];
    areasUltraFastBandsRipples = [];
    entryNames = {};
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
  dataDuration = max([size(LFPbandPower{iCh},2) numel(rippleRate{iCh}) numel(theta2deltaRatio{iCh}) numel(LFP{iCh})...
    numel(PR) size(narrowSlowBands,2) size(slowBands,2) size(fastBands,2) size(LFPs,2) size(PRs,2)]);
  if isempty(PR)
    PR = zeros(1,dataDuration);
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
  
  if size(narrowSlowBands,2) < dataDuration
    zeroPad = zeros(size(narrowSlowBands,1), dataDuration - size(narrowSlowBands,2));
    narrowSlowBands = [narrowSlowBands zeroPad];
  end
  narrowSlowBands = [narrowSlowBands; LFPbandPower{iCh}(1,:)];
  narrowSlowBands(isnan(narrowSlowBands)) = 0;
  
  if size(deltaBands,2) < dataDuration
    zeroPad = zeros(size(deltaBands,1), dataDuration - size(deltaBands,2));
    deltaBands = [deltaBands zeroPad];
  end
  deltaBands = [deltaBands; LFPbandPower{iCh}(2,:)];
  deltaBands(isnan(deltaBands)) = 0;
  
  if size(thetaBands,2) < dataDuration
    zeroPad = zeros(size(thetaBands,1), dataDuration - size(thetaBands,2));
    thetaBands = [thetaBands zeroPad];
  end
  thetaBands = [thetaBands; LFPbandPower{iCh}(3,:)];
  thetaBands(isnan(thetaBands)) = 0;
  
  if size(alphaBands,2) < dataDuration
    zeroPad = zeros(size(alphaBands,1), dataDuration - size(alphaBands,2));
    alphaBands = [alphaBands zeroPad];
  end
  alphaBands = [alphaBands; LFPbandPower{iCh}(4,:)];
  alphaBands(isnan(alphaBands)) = 0;
  
  if size(highSlowBands,2) < dataDuration
    zeroPad = zeros(size(highSlowBands,1), dataDuration - size(highSlowBands,2));
    highSlowBands = [highSlowBands zeroPad];
  end
  highSlowBands = [highSlowBands; LFPbandPower{iCh}(2:4,:)];
  highSlowBands(isnan(highSlowBands)) = 0;
  
  if size(slowBands,2) < dataDuration
    zeroPad = zeros(size(slowBands,1), dataDuration - size(slowBands,2));
    slowBands = [slowBands zeroPad];
  end
  slowBands = [slowBands; LFPbandPower{iCh}(1:4,:); theta2deltaRatio{iCh}];
  slowBands(isnan(slowBands)) = 0;
  
  if size(betaBands,2) < dataDuration
    zeroPad = zeros(size(betaBands,1), dataDuration - size(betaBands,2));
    betaBands = [betaBands zeroPad];
  end
  betaBands = [betaBands; LFPbandPower{iCh}(5,:)];
  betaBands(isnan(betaBands)) = 0;
  
  if size(slowGammaBands,2) < dataDuration
    zeroPad = zeros(size(slowGammaBands,1), dataDuration - size(slowGammaBands,2));
    slowGammaBands = [slowGammaBands zeroPad];
  end
  slowGammaBands = [slowGammaBands; LFPbandPower{iCh}(6,:)];
  slowGammaBands(isnan(slowGammaBands)) = 0;
  
  if size(fastGammaBands,2) < dataDuration
    zeroPad = zeros(size(fastGammaBands,1), dataDuration - size(fastGammaBands,2));
    fastGammaBands = [fastGammaBands zeroPad];
  end
  fastGammaBands = [fastGammaBands; LFPbandPower{iCh}(7,:)];
  fastGammaBands(isnan(fastGammaBands)) = 0;
  
  if size(gammaBands,2) < dataDuration
    zeroPad = zeros(size(gammaBands,1), dataDuration - size(gammaBands,2));
    gammaBands = [gammaBands zeroPad];
  end
  gammaBands = [gammaBands; LFPbandPower{iCh}(6:7,:)];
  gammaBands(isnan(gammaBands)) = 0;
  
  if size(ultraFastBands,2) < dataDuration
    zeroPad = zeros(size(ultraFastBands,1), dataDuration - size(ultraFastBands,2));
    ultraFastBands = [ultraFastBands zeroPad];
  end
  ultraFastBands = [ultraFastBands; LFPbandPower{iCh}(8,:)];
  ultraFastBands(isnan(ultraFastBands)) = 0;
  
  if size(ultraFastBandsRipples,2) < dataDuration
    zeroPad = zeros(size(ultraFastBandsRipples,1), dataDuration - size(ultraFastBandsRipples,2));
    ultraFastBandsRipples = [ultraFastBandsRipples zeroPad];
  end
  ultraFastBandsRipples = [ultraFastBandsRipples; LFPbandPower{iCh}(8,:)];
  if (area == 5 || area == 10)
    ultraFastBandsRipples = [ultraFastBandsRipples; rippleRate{iCh}];
  end
  ultraFastBandsRipples(isnan(ultraFastBandsRipples)) = 0;
    
  if size(fastBands,2) < dataDuration
    zeroPad = zeros(size(fastBands,1), dataDuration - size(fastBands,2));
    fastBands = [fastBands zeroPad];
  end
  fastBands = [fastBands; LFPbandPower{iCh}(5:8,:)];
  if (area == 5 || area == 10)
    fastBands = [fastBands; rippleRate{iCh}];
  end
  fastBands(isnan(fastBands)) = 0;
  
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
  
  areasNarrow = [areasNarrow; area];
  areasHighSlow = [areasHighSlow; ones(3,1)*area];
  areasSlow = [areasSlow; ones(5,1)*area];
  areasGamma = [areasGamma; ones(2,1)*area];
  areasUltraFastBandsRipples = [areasUltraFastBandsRipples; area];
  areasFast = [areasFast; ones(4,1)*area];
  if (area == 5 || area == 10)
    areasUltraFastBandsRipples = [areasUltraFastBandsRipples; area];
    areasFast = [areasFast; area];
  end
  
  entryNames{numel(entryNames) + 1} = entryName;
  
  
  if dbCount < numel(fnsData)
    seriesNameNext = seriesFromEntry(dbStruct.db(dbCount+1).entryName);
  end
  if dbCount < numel(fnsData) - 1
    seriesNameNext2 = seriesFromEntry(dbStruct.db(dbCount+2).entryName);
  end
  if dbCount == numel(fnsData) || ~strcmpi(seriesName(1:14), seriesNameNext(1:14)) ||...
      (dbCount == numel(fnsData)-1 && (area == 5 || area == 10)) || ((area == 5 || area == 10) && ~strcmpi(seriesName(1:14), seriesNameNext2(1:14)))
    
    % SUBTRACT THE MEAN AND Z-SCORE THE DATA
    narrowSlowBands = narrowSlowBands - mean(narrowSlowBands,2);
    narrowSlowBands = zscore(narrowSlowBands')';
    deltaBands = deltaBands - mean(deltaBands,2);
    deltaBands = zscore(deltaBands')';
    thetaBands = thetaBands - mean(thetaBands,2);
    thetaBands = zscore(thetaBands')';
    alphaBands = alphaBands - mean(alphaBands,2);
    alphaBands = zscore(alphaBands')';
    highSlowBands = highSlowBands - mean(highSlowBands,2);
    highSlowBands = zscore(highSlowBands')';
    slowBands = slowBands - mean(slowBands,2);
    slowBands = zscore(slowBands')';
    betaBands = betaBands - mean(betaBands,2);
    betaBands = zscore(betaBands')';
    slowGammaBands = slowGammaBands - mean(slowGammaBands,2);
    slowGammaBands = zscore(slowGammaBands')';
    fastGammaBands = fastGammaBands - mean(fastGammaBands,2);
    fastGammaBands = zscore(fastGammaBands')';
    gammaBands = gammaBands - mean(gammaBands,2);
    gammaBands = zscore(gammaBands')';
    ultraFastBands = ultraFastBands - mean(ultraFastBands,2);
    ultraFastBands = zscore(ultraFastBands')';
    ultraFastBandsRipples = ultraFastBandsRipples - mean(ultraFastBandsRipples,2);
    ultraFastBandsRipples = zscore(ultraFastBandsRipples')';
    fastBands = fastBands - mean(fastBands,2);
    fastBands = zscore(fastBands')';
    LFPs = LFPs - mean(LFPs,2);
    LFPs = zscore(LFPs')';
    PRs = PRs - mean(PRs,2);
    PRs = zscore(PRs')';
    
    
    % RUN PCA
    [pcaCoef_NSB, explained_NSB, PCs_NSB, nPCs_NSB, prob_NSB, individualVarExplained_NSB] = pcaGeneric(narrowSlowBands);
    [pcaCoef_delta, explained_delta, PCs_delta, nPCs_delta, prob_delta, individualVarExplained_delta] = pcaGeneric(deltaBands);
    [pcaCoef_theta, explained_theta, PCs_theta, nPCs_theta, prob_theta, individualVarExplained_theta] = pcaGeneric(thetaBands);
    [pcaCoef_alpha, explained_alpha, PCs_alpha, nPCs_alpha, prob_alpha, individualVarExplained_alpha] = pcaGeneric(alphaBands);
    [pcaCoef_HSB, explained_HSB, PCs_HSB, nPCs_HSB, prob_HSB, individualVarExplained_HSB] = pcaGeneric(highSlowBands);
    [pcaCoef_SB, explained_SB, PCs_SB, nPCs_SB, prob_SB, individualVarExplained_SB] = pcaGeneric(slowBands);
    [pcaCoef_beta, explained_beta, PCs_beta, nPCs_beta, prob_beta, individualVarExplained_beta] = pcaGeneric(betaBands);
    [pcaCoef_SG, explained_SG, PCs_SG, nPCs_SG, prob_SG, individualVarExplained_SG] = pcaGeneric(slowGammaBands);
    [pcaCoef_FG, explained_FG, PCs_FG, nPCs_FG, prob_FG, individualVarExplained_FG] = pcaGeneric(fastGammaBands);
    [pcaCoef_gamma, explained_gamma, PCs_gamma, nPCs_gamma, prob_gamma, individualVarExplained_gamma] = pcaGeneric(gammaBands);
    [pcaCoef_UF, explained_UF, PCs_UF, nPCs_UF, prob_UF, individualVarExplained_UF] = pcaGeneric(ultraFastBands);
    [pcaCoef_UFR, explained_UFR, PCs_UFR, nPCs_UFR, prob_UFR, individualVarExplained_UFR] = pcaGeneric(ultraFastBandsRipples);
    [pcaCoef_FB, explained_FB, PCs_FB, nPCs_FB, prob_FB, individualVarExplained_FB] = pcaGeneric(fastBands);
    [pcaCoef_LFP, explained_LFP, PCs_LFP, nPCs_LFP, prob_LFP, individualVarExplained_LFP] = pcaGeneric(LFPs);
    [pcaCoef_PR, explained_PR, PCs_PR, nPCs_PR, prob_PR, individualVarExplained_PR] = pcaGeneric(PRs);
    
    % Calculate a number of variables needed to explain proportions of variance
    [explained25_NSB, explained50_NSB, explained75_NSB, explained95_NSB, explainedThird_NSB, explainedTwoThirds_NSB] = explainedProportions(explained_NSB);
    [explained25_delta, explained50_delta, explained75_delta, explained95_delta, explainedThird_delta, explainedTwoThirds_delta] = explainedProportions(explained_delta);
    [explained25_theta, explained50_theta, explained75_theta, explained95_theta, explainedThird_theta, explainedTwoThirds_theta] = explainedProportions(explained_theta);
    [explained25_alpha, explained50_alpha, explained75_alpha, explained95_alpha, explainedThird_alpha, explainedTwoThirds_alpha] = explainedProportions(explained_alpha);
    [explained25_HSB, explained50_HSB, explained75_HSB, explained95_HSB, explainedThird_HSB, explainedTwoThirds_HSB] = explainedProportions(explained_HSB);
    [explained25_SB, explained50_SB, explained75_SB, explained95_SB, explainedThird_SB, explainedTwoThirds_SB] = explainedProportions(explained_SB);
    [explained25_beta, explained50_beta, explained75_beta, explained95_beta, explainedThird_beta, explainedTwoThirds_beta] = explainedProportions(explained_beta);
    [explained25_SG, explained50_SG, explained75_SG, explained95_SG, explainedThird_SG, explainedTwoThirds_SG] = explainedProportions(explained_SG);
    [explained25_FG, explained50_FG, explained75_FG, explained95_FG, explainedThird_FG, explainedTwoThirds_FG] = explainedProportions(explained_FG);
    [explained25_gamma, explained50_gamma, explained75_gamma, explained95_gamma, explainedThird_gamma, explainedTwoThirds_gamma] = explainedProportions(explained_gamma);
    [explained25_UF, explained50_UF, explained75_UF, explained95_UF, explainedThird_UF, explainedTwoThirds_UF] = explainedProportions(explained_UF);
    [explained25_UFR, explained50_UFR, explained75_UFR, explained95_UFR, explainedThird_UFR, explainedTwoThirds_UFR] = explainedProportions(explained_UFR);
    [explained25_FB, explained50_FB, explained75_FB, explained95_FB, explainedThird_FB, explainedTwoThirds_FB] = explainedProportions(explained_FB);
    [explained25_LFP, explained50_LFP, explained75_LFP, explained95_LFP, explainedThird_LFP, explainedTwoThirds_LFP] = explainedProportions(explained_LFP);
    [explained25_PR, explained50_PR, explained75_PR, explained95_PR, explainedThird_PR, explainedTwoThirds_PR] = explainedProportions(explained_PR);
    
    
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
    if ~isempty(pupilArea)
      [~, iPCs2pupilAreaStart] = min(abs(interpTimes - pupilAreaTimes(1)));
      [~, iPCs2pupilAreaEnd] = min(abs(interpTimes - pupilAreaTimes(end)));
      PCs2pupilAreaTimes = interpTimes(iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_NSB = PCs_NSB(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_delta = PCs_delta(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_theta = PCs_theta(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_alpha = PCs_alpha(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_HSB = PCs_HSB(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_SB = PCs_SB(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_beta = PCs_beta(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_SG = PCs_SG(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_FG = PCs_FG(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_gamma = PCs_gamma(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_UF = PCs_UF(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_UFR = PCs_UFR(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_FB = PCs_FB(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_LFP = PCs_LFP(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
      PCs2pupilArea_PR = PCs_PR(:,iPCs2pupilAreaStart:iPCs2pupilAreaEnd);
    end
    if ~isempty(motion)
      [~, iPCs2motionStart] = min(abs(interpTimes - motionTimes(1)));
      [~, iPCs2motionEnd] = min(abs(interpTimes - motionTimes(end)));
      PCs2motionTimes = interpTimes(iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_NSB = PCs_NSB(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_delta = PCs_delta(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_theta = PCs_theta(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_alpha = PCs_alpha(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_HSB = PCs_HSB(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_SB = PCs_SB(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_beta = PCs_beta(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_SG = PCs_SG(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_FG = PCs_FG(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_gamma = PCs_gamma(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_UF = PCs_UF(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_UFR = PCs_UFR(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_FB = PCs_FB(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_LFP = PCs_LFP(:,iPCs2motionStart:iPCs2motionEnd);
      PCs2motion_PR = PCs_PR(:,iPCs2motionStart:iPCs2motionEnd);
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
      [rPCs2pupilAreaPearson_NSB, pvalPCs2pupilAreaPearson_NSB] = corrMulti(pupilArea', PCs2pupilArea_NSB, 'Pearson');
      [rPCs2pupilAreaSpearman_NSB, pvalPCs2pupilAreaSpearman_NSB] = corrMulti(pupilArea', PCs2pupilArea_NSB, 'Spearman');
      [rPCs2pupilAreaPearson_delta, pvalPCs2pupilAreaPearson_delta] = corrMulti(pupilArea', PCs2pupilArea_delta, 'Pearson');
      [rPCs2pupilAreaSpearman_delta, pvalPCs2pupilAreaSpearman_delta] = corrMulti(pupilArea', PCs2pupilArea_delta, 'Spearman');
      [rPCs2pupilAreaPearson_theta, pvalPCs2pupilAreaPearson_theta] = corrMulti(pupilArea', PCs2pupilArea_theta, 'Pearson');
      [rPCs2pupilAreaSpearman_theta, pvalPCs2pupilAreaSpearman_theta] = corrMulti(pupilArea', PCs2pupilArea_theta, 'Spearman');
      [rPCs2pupilAreaPearson_alpha, pvalPCs2pupilAreaPearson_alpha] = corrMulti(pupilArea', PCs2pupilArea_alpha, 'Pearson');
      [rPCs2pupilAreaSpearman_alpha, pvalPCs2pupilAreaSpearman_alpha] = corrMulti(pupilArea', PCs2pupilArea_alpha, 'Spearman');
      [rPCs2pupilAreaPearson_HSB, pvalPCs2pupilAreaPearson_HSB] = corrMulti(pupilArea', PCs2pupilArea_HSB, 'Pearson');
      [rPCs2pupilAreaSpearman_HSB, pvalPCs2pupilAreaSpearman_HSB] = corrMulti(pupilArea', PCs2pupilArea_HSB, 'Spearman');
      [rPCs2pupilAreaPearson_SB, pvalPCs2pupilAreaPearson_SB] = corrMulti(pupilArea', PCs2pupilArea_SB, 'Pearson');
      [rPCs2pupilAreaSpearman_SB, pvalPCs2pupilAreaSpearman_SB] = corrMulti(pupilArea', PCs2pupilArea_SB, 'Spearman');
      [rPCs2pupilAreaPearson_beta, pvalPCs2pupilAreaPearson_beta] = corrMulti(pupilArea', PCs2pupilArea_beta, 'Pearson');
      [rPCs2pupilAreaSpearman_beta, pvalPCs2pupilAreaSpearman_beta] = corrMulti(pupilArea', PCs2pupilArea_beta, 'Spearman');
      [rPCs2pupilAreaPearson_SG, pvalPCs2pupilAreaPearson_SG] = corrMulti(pupilArea', PCs2pupilArea_SG, 'Pearson');
      [rPCs2pupilAreaSpearman_SG, pvalPCs2pupilAreaSpearman_SG] = corrMulti(pupilArea', PCs2pupilArea_SG, 'Spearman');
      [rPCs2pupilAreaPearson_FG, pvalPCs2pupilAreaPearson_FG] = corrMulti(pupilArea', PCs2pupilArea_FG, 'Pearson');
      [rPCs2pupilAreaSpearman_FG, pvalPCs2pupilAreaSpearman_FG] = corrMulti(pupilArea', PCs2pupilArea_FG, 'Spearman');
      [rPCs2pupilAreaPearson_gamma, pvalPCs2pupilAreaPearson_gamma] = corrMulti(pupilArea', PCs2pupilArea_gamma, 'Pearson');
      [rPCs2pupilAreaSpearman_gamma, pvalPCs2pupilAreaSpearman_gamma] = corrMulti(pupilArea', PCs2pupilArea_gamma, 'Spearman');
      [rPCs2pupilAreaPearson_UF, pvalPCs2pupilAreaPearson_UF] = corrMulti(pupilArea', PCs2pupilArea_UF, 'Pearson');
      [rPCs2pupilAreaSpearman_UF, pvalPCs2pupilAreaSpearman_UF] = corrMulti(pupilArea', PCs2pupilArea_UF, 'Spearman');
      [rPCs2pupilAreaPearson_UFR, pvalPCs2pupilAreaPearson_UFR] = corrMulti(pupilArea', PCs2pupilArea_UFR, 'Pearson');
      [rPCs2pupilAreaSpearman_UFR, pvalPCs2pupilAreaSpearman_UFR] = corrMulti(pupilArea', PCs2pupilArea_UFR, 'Spearman');
      [rPCs2pupilAreaPearson_FB, pvalPCs2pupilAreaPearson_FB] = corrMulti(pupilArea', PCs2pupilArea_FB, 'Pearson');
      [rPCs2pupilAreaSpearman_FB, pvalPCs2pupilAreaSpearman_FB] = corrMulti(pupilArea', PCs2pupilArea_FB, 'Spearman');
      [rPCs2pupilAreaPearson_LFP, pvalPCs2pupilAreaPearson_LFP] = corrMulti(pupilArea', PCs2pupilArea_LFP, 'Pearson');
      [rPCs2pupilAreaSpearman_LFP, pvalPCs2pupilAreaSpearman_LFP] = corrMulti(pupilArea', PCs2pupilArea_LFP, 'Spearman');
      [rPCs2pupilAreaPearson_PR, pvalPCs2pupilAreaPearson_PR] = corrMulti(pupilArea', PCs2pupilArea_PR, 'Pearson');
      [rPCs2pupilAreaSpearman_PR, pvalPCs2pupilAreaSpearman_PR] = corrMulti(pupilArea', PCs2pupilArea_PR, 'Spearman');
    end
    if ~isempty(motion)
      [rPCs2motionPearson_NSB, pvalPCs2motionPearson_NSB] = corrMulti(motion', PCs2motion_NSB, 'Pearson');
      [rPCs2motionSpearman_NSB, pvalPCs2motionSpearman_NSB] = corrMulti(motion', PCs2motion_NSB, 'Spearman');
      [rPCs2motionPearson_delta, pvalPCs2motionPearson_delta] = corrMulti(motion', PCs2motion_delta, 'Pearson');
      [rPCs2motionSpearman_delta, pvalPCs2motionSpearman_delta] = corrMulti(motion', PCs2motion_delta, 'Spearman');
      [rPCs2motionPearson_theta, pvalPCs2motionPearson_theta] = corrMulti(motion', PCs2motion_theta, 'Pearson');
      [rPCs2motionSpearman_theta, pvalPCs2motionSpearman_theta] = corrMulti(motion', PCs2motion_theta, 'Spearman');
      [rPCs2motionPearson_alpha, pvalPCs2motionPearson_alpha] = corrMulti(motion', PCs2motion_alpha, 'Pearson');
      [rPCs2motionSpearman_alpha, pvalPCs2motionSpearman_alpha] = corrMulti(motion', PCs2motion_alpha, 'Spearman');
      [rPCs2motionPearson_HSB, pvalPCs2motionPearson_HSB] = corrMulti(motion', PCs2motion_HSB, 'Pearson');
      [rPCs2motionSpearman_HSB, pvalPCs2motionSpearman_HSB] = corrMulti(motion', PCs2motion_HSB, 'Spearman');
      [rPCs2motionPearson_SB, pvalPCs2motionPearson_SB] = corrMulti(motion', PCs2motion_SB, 'Pearson');
      [rPCs2motionSpearman_SB, pvalPCs2motionSpearman_SB] = corrMulti(motion', PCs2motion_SB, 'Spearman');
      [rPCs2motionPearson_beta, pvalPCs2motionPearson_beta] = corrMulti(motion', PCs2motion_beta, 'Pearson');
      [rPCs2motionSpearman_beta, pvalPCs2motionSpearman_beta] = corrMulti(motion', PCs2motion_beta, 'Spearman');
      [rPCs2motionPearson_SG, pvalPCs2motionPearson_SG] = corrMulti(motion', PCs2motion_SG, 'Pearson');
      [rPCs2motionSpearman_SG, pvalPCs2motionSpearman_SG] = corrMulti(motion', PCs2motion_SG, 'Spearman');
      [rPCs2motionPearson_FG, pvalPCs2motionPearson_FG] = corrMulti(motion', PCs2motion_FG, 'Pearson');
      [rPCs2motionSpearman_FG, pvalPCs2motionSpearman_FG] = corrMulti(motion', PCs2motion_FG, 'Spearman');
      [rPCs2motionPearson_gamma, pvalPCs2motionPearson_gamma] = corrMulti(motion', PCs2motion_gamma, 'Pearson');
      [rPCs2motionSpearman_gamma, pvalPCs2motionSpearman_gamma] = corrMulti(motion', PCs2motion_gamma, 'Spearman');
      [rPCs2motionPearson_UF, pvalPCs2motionPearson_UF] = corrMulti(motion', PCs2motion_UF, 'Pearson');
      [rPCs2motionSpearman_UF, pvalPCs2motionSpearman_UF] = corrMulti(motion', PCs2motion_UF, 'Spearman');
      [rPCs2motionPearson_UFR, pvalPCs2motionPearson_UFR] = corrMulti(motion', PCs2motion_UFR, 'Pearson');
      [rPCs2motionSpearman_UFR, pvalPCs2motionSpearman_UFR] = corrMulti(motion', PCs2motion_UFR, 'Spearman');
      [rPCs2motionPearson_FB, pvalPCs2motionPearson_FB] = corrMulti(motion', PCs2motion_FB, 'Pearson');
      [rPCs2motionSpearman_FB, pvalPCs2motionSpearman_FB] = corrMulti(motion', PCs2motion_FB, 'Spearman');
      [rPCs2motionPearson_LFP, pvalPCs2motionPearson_LFP] = corrMulti(motion', PCs2motion_LFP, 'Pearson');
      [rPCs2motionSpearman_LFP, pvalPCs2motionSpearman_LFP] = corrMulti(motion', PCs2motion_LFP, 'Spearman');
      [rPCs2motionPearson_PR, pvalPCs2motionPearson_PR] = corrMulti(motion', PCs2motion_PR, 'Pearson');
      [rPCs2motionSpearman_PR, pvalPCs2motionSpearman_PR] = corrMulti(motion', PCs2motion_PR, 'Spearman');
    end
    
    
    % CALCULATE EXPLAINED VARIANCE FOR EACH AREA
    areasOI = unique(areasNarrow);
    for iArea = 1:numel(areasOI)
      dataOInarrow = getDataOI(areasNarrow, areasOI, iArea);
      dataOIhighSlow = getDataOI(areasHighSlow, areasOI, iArea);
      dataOIslow = getDataOI(areasSlow, areasOI, iArea);
      dataOIgamma = getDataOI(areasGamma, areasOI, iArea);
      dataOIultraFastBandsRipples = getDataOI(areasUltraFastBandsRipples, areasOI, iArea);
      dataOIfast = getDataOI(areasFast, areasOI, iArea);
      pcaCoef_area_NSB = pcaCoef_NSB(dataOInarrow);
      pcaCoef_area_delta = pcaCoef_SB(dataOInarrow);
      pcaCoef_area_theta = pcaCoef_SB(dataOInarrow);
      pcaCoef_area_alpha = pcaCoef_SB(dataOInarrow);
      pcaCoef_area_HSB = pcaCoef_SB(dataOInarrow);
      pcaCoef_area_SB = pcaCoef_SB(dataOIslow);
      pcaCoef_area_beta = pcaCoef_beta(dataOInarrow);
      pcaCoef_area_SG = pcaCoef_SG(dataOInarrow);
      pcaCoef_area_FG = pcaCoef_FG(dataOInarrow);
      pcaCoef_area_gamma = pcaCoef_gamma(dataOIgamma);
      pcaCoef_area_UF = pcaCoef_UF(dataOInarrow);
      pcaCoef_area_UFR = pcaCoef_UFR(dataOIultraFastBandsRipples);
      pcaCoef_area_FB = pcaCoef_FB(dataOIfast);
      pcaCoef_area_LFP = pcaCoef_LFP(dataOInarrow);
      pcaCoef_area_PR = pcaCoef_PR(dataOInarrow);
      explained_area_NSB = 100*mean(individualVarExplained_NSB(:,dataOInarrow), 2)';
      explained_area_delta = 100*mean(individualVarExplained_delta(:,dataOInarrow), 2)';
      explained_area_theta = 100*mean(individualVarExplained_theta(:,dataOInarrow), 2)';
      explained_area_alpha = 100*mean(individualVarExplained_alpha(:,dataOInarrow), 2)';
      explained_area_HSB = 100*mean(individualVarExplained_HSB(:,dataOInarrow), 2)';
      explained_area_SB = 100*mean(individualVarExplained_SB(:,dataOIslow), 2)';
      explained_area_beta = 100*mean(individualVarExplained_beta(:,dataOInarrow), 2)';
      explained_area_SG = 100*mean(individualVarExplained_SG(:,dataOInarrow), 2)';
      explained_area_FG = 100*mean(individualVarExplained_FG(:,dataOInarrow), 2)';
      explained_area_gamma = 100*mean(individualVarExplained_gamma(:,dataOIgamma), 2)';
      explained_area_UF = 100*mean(individualVarExplained_UF(:,dataOInarrow), 2)';
      explained_area_UFR = 100*mean(individualVarExplained_UFR(:,dataOIultraFastBandsRipples), 2)';
      explained_area_FB = 100*mean(individualVarExplained_FB(:,dataOIfast), 2)';
      explained_area_LFP = 100*mean(individualVarExplained_LFP(:,dataOInarrow), 2)';
      explained_area_PR = 100*mean(individualVarExplained_PR(:,dataOInarrow), 2)';
      
      % Calculate a number of variables needed to explain proportions of variance for the area
      [explained25_area_NSB, explained50_area_NSB, explained75_area_NSB, explained95_area_NSB,...
        explainedThird_area_NSB, explainedTwoThirds_area_NSB] = explainedProportions(explained_area_NSB);
      [explained25_area_delta, explained50_area_delta, explained75_area_delta, explained95_area_delta,...
        explainedThird_area_delta, explainedTwoThirds_area_delta] = explainedProportions(explained_area_delta);
      [explained25_area_theta, explained50_area_theta, explained75_area_theta, explained95_area_theta,...
        explainedThird_area_theta, explainedTwoThirds_area_theta] = explainedProportions(explained_area_theta);
      [explained25_area_alpha, explained50_area_alpha, explained75_area_alpha, explained95_area_alpha,...
        explainedThird_area_alpha, explainedTwoThirds_area_alpha] = explainedProportions(explained_area_alpha);
      [explained25_area_HSB, explained50_area_HSB, explained75_area_HSB, explained95_area_HSB,...
        explainedThird_area_HSB, explainedTwoThirds_area_HSB] = explainedProportions(explained_area_HSB);
      [explained25_area_SB, explained50_area_SB, explained75_area_SB, explained95_area_SB,...
        explainedThird_area_SB, explainedTwoThirds_area_SB] = explainedProportions(explained_area_SB);
      [explained25_area_beta, explained50_area_beta, explained75_area_beta, explained95_area_beta,...
        explainedThird_area_beta, explainedTwoThirds_area_beta] = explainedProportions(explained_area_beta);
      [explained25_area_SG, explained50_area_SG, explained75_area_SG, explained95_area_SG,...
        explainedThird_area_SG, explainedTwoThirds_area_SG] = explainedProportions(explained_area_SG);
      [explained25_area_FG, explained50_area_FG, explained75_area_FG, explained95_area_FG,...
        explainedThird_area_FG, explainedTwoThirds_area_FG] = explainedProportions(explained_area_FG);
      [explained25_area_gamma, explained50_area_gamma, explained75_area_gamma, explained95_area_gamma,...
        explainedThird_area_gamma, explainedTwoThirds_area_gamma] = explainedProportions(explained_area_gamma);
      [explained25_area_UF, explained50_area_UF, explained75_area_UF, explained95_area_UF,...
        explainedThird_area_UF, explainedTwoThirds_area_UF] = explainedProportions(explained_area_UF);
      [explained25_area_UFR, explained50_area_UFR, explained75_area_UFR, explained95_area_UFR,...
        explainedThird_area_UFR, explainedTwoThirds_area_UFR] = explainedProportions(explained_area_UFR);
      [explained25_area_FB, explained50_area_FB, explained75_area_FB, explained95_area_FB,...
        explainedThird_area_FB, explainedTwoThirds_area_FB] = explainedProportions(explained_area_FB);
      [explained25_area_LFP, explained50_area_LFP, explained75_area_LFP, explained95_area_LFP,...
        explainedThird_area_LFP, explainedTwoThirds_area_LFP] = explainedProportions(explained_area_LFP);
      [explained25_area_PR, explained50_area_PR, explained75_area_PR, explained95_area_PR,...
        explainedThird_area_PR, explainedTwoThirds_area_PR] = explainedProportions(explained_area_PR);
    
    
      % SAVE NARROW SLOW BAND DATA
      % PCA data for all areas
      pcaData.pcaCoef_NSB = pcaCoef_NSB;
      pcaData.explained_NSB = explained_NSB;
      pcaData.nPCs_NSB = nPCs_NSB;
      pcaData.prob_NSB = prob_NSB;
      pcaData.explained25_NSB = explained25_NSB;
      pcaData.explained50_NSB = explained50_NSB;
      pcaData.explained75_NSB = explained75_NSB;
      pcaData.explained95_NSB = explained95_NSB;
      pcaData.explainedThird_NSB = explainedThird_NSB;
      pcaData.explainedTwoThirds_NSB = explainedTwoThirds_NSB;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_NSB = rPCs2pupilAreaPearson_NSB;
        pcaData.pvalPCs2pupilAreaPearson_NSB = pvalPCs2pupilAreaPearson_NSB;
        pcaData.rPCs2pupilAreaSpearman_NSB = rPCs2pupilAreaSpearman_NSB;
        pcaData.pvalPCs2pupilAreaSpearman_NSB = pvalPCs2pupilAreaSpearman_NSB;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_NSB = rPCs2motionPearson_NSB;
        pcaData.pvalPCs2motionPearson_NSB = pvalPCs2motionPearson_NSB;
        pcaData.rPCs2motionSpearman_NSB = rPCs2motionSpearman_NSB;
        pcaData.pvalPCs2motionSpearman_NSB = pvalPCs2motionSpearman_NSB;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_NSB = pcaCoef_area_NSB;
      pcaData.explained_area_NSB = explained_area_NSB;
      pcaData.explained25_area_NSB = explained25_area_NSB;
      pcaData.explained50_area_NSB = explained50_area_NSB;
      pcaData.explained75_area_NSB = explained75_area_NSB;
      pcaData.explained95_area_NSB = explained95_area_NSB;
      pcaData.explainedThird_area_NSB = explainedThird_area_NSB;
      pcaData.explainedTwoThirds_area_NSB = explainedTwoThirds_area_NSB;
      
      
      % SAVE DELTA BAND DATA
      % PCA data for all areas
      pcaData.pcaCoef_delta = pcaCoef_delta;
      pcaData.explained_delta = explained_delta;
      pcaData.nPCs_delta = nPCs_delta;
      pcaData.prob_delta = prob_delta;
      pcaData.explained25_delta = explained25_delta;
      pcaData.explained50_delta = explained50_delta;
      pcaData.explained75_delta = explained75_delta;
      pcaData.explained95_delta = explained95_delta;
      pcaData.explainedThird_delta = explainedThird_delta;
      pcaData.explainedTwoThirds_delta = explainedTwoThirds_delta;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_delta = rPCs2pupilAreaPearson_delta;
        pcaData.pvalPCs2pupilAreaPearson_delta = pvalPCs2pupilAreaPearson_delta;
        pcaData.rPCs2pupilAreaSpearman_delta = rPCs2pupilAreaSpearman_delta;
        pcaData.pvalPCs2pupilAreaSpearman_delta = pvalPCs2pupilAreaSpearman_delta;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_delta = rPCs2motionPearson_delta;
        pcaData.pvalPCs2motionPearson_delta = pvalPCs2motionPearson_delta;
        pcaData.rPCs2motionSpearman_delta = rPCs2motionSpearman_delta;
        pcaData.pvalPCs2motionSpearman_delta = pvalPCs2motionSpearman_delta;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_delta = pcaCoef_area_delta;
      pcaData.explained_area_delta = explained_area_delta;
      pcaData.explained25_area_delta = explained25_area_delta;
      pcaData.explained50_area_delta = explained50_area_delta;
      pcaData.explained75_area_delta = explained75_area_delta;
      pcaData.explained95_area_delta = explained95_area_delta;
      pcaData.explainedThird_area_delta = explainedThird_area_delta;
      pcaData.explainedTwoThirds_area_delta = explainedTwoThirds_area_delta;
      
      
      % SAVE THETA BAND DATA
      % PCA data for all areas
      pcaData.pcaCoef_theta = pcaCoef_theta;
      pcaData.explained_theta = explained_theta;
      pcaData.nPCs_theta = nPCs_theta;
      pcaData.prob_theta = prob_theta;
      pcaData.explained25_theta = explained25_theta;
      pcaData.explained50_theta = explained50_theta;
      pcaData.explained75_theta = explained75_theta;
      pcaData.explained95_theta = explained95_theta;
      pcaData.explainedThird_theta = explainedThird_theta;
      pcaData.explainedTwoThirds_theta = explainedTwoThirds_theta;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_theta = rPCs2pupilAreaPearson_theta;
        pcaData.pvalPCs2pupilAreaPearson_theta = pvalPCs2pupilAreaPearson_theta;
        pcaData.rPCs2pupilAreaSpearman_theta = rPCs2pupilAreaSpearman_theta;
        pcaData.pvalPCs2pupilAreaSpearman_theta = pvalPCs2pupilAreaSpearman_theta;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_theta = rPCs2motionPearson_theta;
        pcaData.pvalPCs2motionPearson_theta = pvalPCs2motionPearson_theta;
        pcaData.rPCs2motionSpearman_theta = rPCs2motionSpearman_theta;
        pcaData.pvalPCs2motionSpearman_theta = pvalPCs2motionSpearman_theta;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_theta = pcaCoef_area_theta;
      pcaData.explained_area_theta = explained_area_theta;
      pcaData.explained25_area_theta = explained25_area_theta;
      pcaData.explained50_area_theta = explained50_area_theta;
      pcaData.explained75_area_theta = explained75_area_theta;
      pcaData.explained95_area_theta = explained95_area_theta;
      pcaData.explainedThird_area_theta = explainedThird_area_theta;
      pcaData.explainedTwoThirds_area_theta = explainedTwoThirds_area_theta;
      
      
      % SAVE ALPHA BAND DATA
      % PCA data for all areas
      pcaData.pcaCoef_alpha = pcaCoef_alpha;
      pcaData.explained_alpha = explained_alpha;
      pcaData.nPCs_alpha = nPCs_alpha;
      pcaData.prob_alpha = prob_alpha;
      pcaData.explained25_alpha = explained25_alpha;
      pcaData.explained50_alpha = explained50_alpha;
      pcaData.explained75_alpha = explained75_alpha;
      pcaData.explained95_alpha = explained95_alpha;
      pcaData.explainedThird_alpha = explainedThird_alpha;
      pcaData.explainedTwoThirds_alpha = explainedTwoThirds_alpha;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_alpha = rPCs2pupilAreaPearson_alpha;
        pcaData.pvalPCs2pupilAreaPearson_alpha = pvalPCs2pupilAreaPearson_alpha;
        pcaData.rPCs2pupilAreaSpearman_alpha = rPCs2pupilAreaSpearman_alpha;
        pcaData.pvalPCs2pupilAreaSpearman_alpha = pvalPCs2pupilAreaSpearman_alpha;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_alpha = rPCs2motionPearson_alpha;
        pcaData.pvalPCs2motionPearson_alpha = pvalPCs2motionPearson_alpha;
        pcaData.rPCs2motionSpearman_alpha = rPCs2motionSpearman_alpha;
        pcaData.pvalPCs2motionSpearman_alpha = pvalPCs2motionSpearman_alpha;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_alpha = pcaCoef_area_alpha;
      pcaData.explained_area_alpha = explained_area_alpha;
      pcaData.explained25_area_alpha = explained25_area_alpha;
      pcaData.explained50_area_alpha = explained50_area_alpha;
      pcaData.explained75_area_alpha = explained75_area_alpha;
      pcaData.explained95_area_alpha = explained95_area_alpha;
      pcaData.explainedThird_area_alpha = explainedThird_area_alpha;
      pcaData.explainedTwoThirds_area_alpha = explainedTwoThirds_area_alpha;
      
      
      % SAVE HIG SLOW BAND DATA
      % PCA data for all areas
      pcaData.pcaCoef_HSB = pcaCoef_HSB;
      pcaData.explained_HSB = explained_HSB;
      pcaData.nPCs_HSB = nPCs_HSB;
      pcaData.prob_HSB = prob_HSB;
      pcaData.explained25_HSB = explained25_HSB;
      pcaData.explained50_HSB = explained50_HSB;
      pcaData.explained75_HSB = explained75_HSB;
      pcaData.explained95_HSB = explained95_HSB;
      pcaData.explainedThird_HSB = explainedThird_HSB;
      pcaData.explainedTwoThirds_HSB = explainedTwoThirds_HSB;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_HSB = rPCs2pupilAreaPearson_HSB;
        pcaData.pvalPCs2pupilAreaPearson_HSB = pvalPCs2pupilAreaPearson_HSB;
        pcaData.rPCs2pupilAreaSpearman_HSB = rPCs2pupilAreaSpearman_HSB;
        pcaData.pvalPCs2pupilAreaSpearman_HSB = pvalPCs2pupilAreaSpearman_HSB;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_HSB = rPCs2motionPearson_HSB;
        pcaData.pvalPCs2motionPearson_HSB = pvalPCs2motionPearson_HSB;
        pcaData.rPCs2motionSpearman_HSB = rPCs2motionSpearman_HSB;
        pcaData.pvalPCs2motionSpearman_HSB = pvalPCs2motionSpearman_HSB;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_HSB = pcaCoef_area_HSB;
      pcaData.explained_area_HSB = explained_area_HSB;
      pcaData.explained25_area_HSB = explained25_area_HSB;
      pcaData.explained50_area_HSB = explained50_area_HSB;
      pcaData.explained75_area_HSB = explained75_area_HSB;
      pcaData.explained95_area_HSB = explained95_area_HSB;
      pcaData.explainedThird_area_HSB = explainedThird_area_HSB;
      pcaData.explainedTwoThirds_area_HSB = explainedTwoThirds_area_HSB;
      
      
      % SAVE SLOW BAND DATA
      % PCA data for all areas
      pcaData.pcaCoef_SB = pcaCoef_SB;
      pcaData.explained_SB = explained_SB;
      pcaData.nPCs_SB = nPCs_SB;
      pcaData.prob_SB = prob_SB;
      pcaData.explained25_SB = explained25_SB;
      pcaData.explained50_SB = explained50_SB;
      pcaData.explained75_SB = explained75_SB;
      pcaData.explained95_SB = explained95_SB;
      pcaData.explainedThird_SB = explainedThird_SB;
      pcaData.explainedTwoThirds_SB = explainedTwoThirds_SB;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_SB = rPCs2pupilAreaPearson_SB;
        pcaData.pvalPCs2pupilAreaPearson_SB = pvalPCs2pupilAreaPearson_SB;
        pcaData.rPCs2pupilAreaSpearman_SB = rPCs2pupilAreaSpearman_SB;
        pcaData.pvalPCs2pupilAreaSpearman_SB = pvalPCs2pupilAreaSpearman_SB;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_SB = rPCs2motionPearson_SB;
        pcaData.pvalPCs2motionPearson_SB = pvalPCs2motionPearson_SB;
        pcaData.rPCs2motionSpearman_SB = rPCs2motionSpearman_SB;
        pcaData.pvalPCs2motionSpearman_SB = pvalPCs2motionSpearman_SB;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_SB = pcaCoef_area_SB;
      pcaData.explained_area_SB = explained_area_SB;
      pcaData.explained25_area_SB = explained25_area_SB;
      pcaData.explained50_area_SB = explained50_area_SB;
      pcaData.explained75_area_SB = explained75_area_SB;
      pcaData.explained95_area_SB = explained95_area_SB;
      pcaData.explainedThird_area_SB = explainedThird_area_SB;
      pcaData.explainedTwoThirds_area_SB = explainedTwoThirds_area_SB;
      
      
      % SAVE BETA BAND DATA
      % PCA data for all areas
      pcaData.pcaCoef_beta = pcaCoef_beta;
      pcaData.explained_beta = explained_beta;
      pcaData.nPCs_beta = nPCs_beta;
      pcaData.prob_beta = prob_beta;
      pcaData.explained25_beta = explained25_beta;
      pcaData.explained50_beta = explained50_beta;
      pcaData.explained75_beta = explained75_beta;
      pcaData.explained95_beta = explained95_beta;
      pcaData.explainedThird_beta = explainedThird_beta;
      pcaData.explainedTwoThirds_beta = explainedTwoThirds_beta;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_beta = rPCs2pupilAreaPearson_beta;
        pcaData.pvalPCs2pupilAreaPearson_beta = pvalPCs2pupilAreaPearson_beta;
        pcaData.rPCs2pupilAreaSpearman_beta = rPCs2pupilAreaSpearman_beta;
        pcaData.pvalPCs2pupilAreaSpearman_beta = pvalPCs2pupilAreaSpearman_beta;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_beta = rPCs2motionPearson_beta;
        pcaData.pvalPCs2motionPearson_beta = pvalPCs2motionPearson_beta;
        pcaData.rPCs2motionSpearman_beta = rPCs2motionSpearman_beta;
        pcaData.pvalPCs2motionSpearman_beta = pvalPCs2motionSpearman_beta;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_beta = pcaCoef_area_beta;
      pcaData.explained_area_beta = explained_area_beta;
      pcaData.explained25_area_beta = explained25_area_beta;
      pcaData.explained50_area_beta = explained50_area_beta;
      pcaData.explained75_area_beta = explained75_area_beta;
      pcaData.explained95_area_beta = explained95_area_beta;
      pcaData.explainedThird_area_beta = explainedThird_area_beta;
      pcaData.explainedTwoThirds_area_beta = explainedTwoThirds_area_beta;
      
      
      % SAVE SLOW GAMMA BAND DATA
      % PCA data for all areas
      pcaData.pcaCoef_SG = pcaCoef_SG;
      pcaData.explained_SG = explained_SG;
      pcaData.nPCs_SG = nPCs_SG;
      pcaData.prob_SG = prob_SG;
      pcaData.explained25_SG = explained25_SG;
      pcaData.explained50_SG = explained50_SG;
      pcaData.explained75_SG = explained75_SG;
      pcaData.explained95_SG = explained95_SG;
      pcaData.explainedThird_SG = explainedThird_SG;
      pcaData.explainedTwoThirds_SG = explainedTwoThirds_SG;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_SG = rPCs2pupilAreaPearson_SG;
        pcaData.pvalPCs2pupilAreaPearson_SG = pvalPCs2pupilAreaPearson_SG;
        pcaData.rPCs2pupilAreaSpearman_SG = rPCs2pupilAreaSpearman_SG;
        pcaData.pvalPCs2pupilAreaSpearman_SG = pvalPCs2pupilAreaSpearman_SG;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_SG = rPCs2motionPearson_SG;
        pcaData.pvalPCs2motionPearson_SG = pvalPCs2motionPearson_SG;
        pcaData.rPCs2motionSpearman_SG = rPCs2motionSpearman_SG;
        pcaData.pvalPCs2motionSpearman_SG = pvalPCs2motionSpearman_SG;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_SG = pcaCoef_area_SG;
      pcaData.explained_area_SG = explained_area_SG;
      pcaData.explained25_area_SG = explained25_area_SG;
      pcaData.explained50_area_SG = explained50_area_SG;
      pcaData.explained75_area_SG = explained75_area_SG;
      pcaData.explained95_area_SG = explained95_area_SG;
      pcaData.explainedThird_area_SG = explainedThird_area_SG;
      pcaData.explainedTwoThirds_area_SG = explainedTwoThirds_area_SG;
      
      
      % SAVE FAST GAMMA BAND DATA
      % PCA data for all areas
      pcaData.pcaCoef_FG = pcaCoef_FG;
      pcaData.explained_FG = explained_FG;
      pcaData.nPCs_FG = nPCs_FG;
      pcaData.prob_FG = prob_FG;
      pcaData.explained25_FG = explained25_FG;
      pcaData.explained50_FG = explained50_FG;
      pcaData.explained75_FG = explained75_FG;
      pcaData.explained95_FG = explained95_FG;
      pcaData.explainedThird_FG = explainedThird_FG;
      pcaData.explainedTwoThirds_FG = explainedTwoThirds_FG;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_FG = rPCs2pupilAreaPearson_FG;
        pcaData.pvalPCs2pupilAreaPearson_FG = pvalPCs2pupilAreaPearson_FG;
        pcaData.rPCs2pupilAreaSpearman_FG = rPCs2pupilAreaSpearman_FG;
        pcaData.pvalPCs2pupilAreaSpearman_FG = pvalPCs2pupilAreaSpearman_FG;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_FG = rPCs2motionPearson_FG;
        pcaData.pvalPCs2motionPearson_FG = pvalPCs2motionPearson_FG;
        pcaData.rPCs2motionSpearman_FG = rPCs2motionSpearman_FG;
        pcaData.pvalPCs2motionSpearman_FG = pvalPCs2motionSpearman_FG;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_FG = pcaCoef_area_FG;
      pcaData.explained_area_FG = explained_area_FG;
      pcaData.explained25_area_FG = explained25_area_FG;
      pcaData.explained50_area_FG = explained50_area_FG;
      pcaData.explained75_area_FG = explained75_area_FG;
      pcaData.explained95_area_FG = explained95_area_FG;
      pcaData.explainedThird_area_FG = explainedThird_area_FG;
      pcaData.explainedTwoThirds_area_FG = explainedTwoThirds_area_FG;
      
      
      % SAVE GAMMA BAND DATA
      % PCA data for all areas
      pcaData.pcaCoef_gamma = pcaCoef_gamma;
      pcaData.explained_gamma = explained_gamma;
      pcaData.nPCs_gamma = nPCs_gamma;
      pcaData.prob_gamma = prob_gamma;
      pcaData.explained25_gamma = explained25_gamma;
      pcaData.explained50_gamma = explained50_gamma;
      pcaData.explained75_gamma = explained75_gamma;
      pcaData.explained95_gamma = explained95_gamma;
      pcaData.explainedThird_gamma = explainedThird_gamma;
      pcaData.explainedTwoThirds_gamma = explainedTwoThirds_gamma;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_gamma = rPCs2pupilAreaPearson_gamma;
        pcaData.pvalPCs2pupilAreaPearson_gamma = pvalPCs2pupilAreaPearson_gamma;
        pcaData.rPCs2pupilAreaSpearman_gamma = rPCs2pupilAreaSpearman_gamma;
        pcaData.pvalPCs2pupilAreaSpearman_gamma = pvalPCs2pupilAreaSpearman_gamma;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_gamma = rPCs2motionPearson_gamma;
        pcaData.pvalPCs2motionPearson_gamma = pvalPCs2motionPearson_gamma;
        pcaData.rPCs2motionSpearman_gamma = rPCs2motionSpearman_gamma;
        pcaData.pvalPCs2motionSpearman_gamma = pvalPCs2motionSpearman_gamma;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_gamma = pcaCoef_area_gamma;
      pcaData.explained_area_gamma = explained_area_gamma;
      pcaData.explained25_area_gamma = explained25_area_gamma;
      pcaData.explained50_area_gamma = explained50_area_gamma;
      pcaData.explained75_area_gamma = explained75_area_gamma;
      pcaData.explained95_area_gamma = explained95_area_gamma;
      pcaData.explainedThird_area_gamma = explainedThird_area_gamma;
      pcaData.explainedTwoThirds_area_gamma = explainedTwoThirds_area_gamma;
      
      
      % SAVE ULTRA FAST BAND DATA
      % PCA data for all areas
      pcaData.pcaCoef_UF = pcaCoef_UF;
      pcaData.explained_UF = explained_UF;
      pcaData.nPCs_UF = nPCs_UF;
      pcaData.prob_UF = prob_UF;
      pcaData.explained25_UF = explained25_UF;
      pcaData.explained50_UF = explained50_UF;
      pcaData.explained75_UF = explained75_UF;
      pcaData.explained95_UF = explained95_UF;
      pcaData.explainedThird_UF = explainedThird_UF;
      pcaData.explainedTwoThirds_UF = explainedTwoThirds_UF;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_UF = rPCs2pupilAreaPearson_UF;
        pcaData.pvalPCs2pupilAreaPearson_UF = pvalPCs2pupilAreaPearson_UF;
        pcaData.rPCs2pupilAreaSpearman_UF = rPCs2pupilAreaSpearman_UF;
        pcaData.pvalPCs2pupilAreaSpearman_UF = pvalPCs2pupilAreaSpearman_UF;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_UF = rPCs2motionPearson_UF;
        pcaData.pvalPCs2motionPearson_UF = pvalPCs2motionPearson_UF;
        pcaData.rPCs2motionSpearman_UF = rPCs2motionSpearman_UF;
        pcaData.pvalPCs2motionSpearman_UF = pvalPCs2motionSpearman_UF;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_UF = pcaCoef_area_UF;
      pcaData.explained_area_UF = explained_area_UF;
      pcaData.explained25_area_UF = explained25_area_UF;
      pcaData.explained50_area_UF = explained50_area_UF;
      pcaData.explained75_area_UF = explained75_area_UF;
      pcaData.explained95_area_UF = explained95_area_UF;
      pcaData.explainedThird_area_UF = explainedThird_area_UF;
      pcaData.explainedTwoThirds_area_UF = explainedTwoThirds_area_UF;
      
      
      % SAVE ULTRA FAST BAND AND RIPPLES DATA
      % PCA data for all areas
      pcaData.pcaCoef_UFR = pcaCoef_UFR;
      pcaData.explained_UFR = explained_UFR;
      pcaData.nPCs_UFR = nPCs_UFR;
      pcaData.prob_UFR = prob_UFR;
      pcaData.explained25_UFR = explained25_UFR;
      pcaData.explained50_UFR = explained50_UFR;
      pcaData.explained75_UFR = explained75_UFR;
      pcaData.explained95_UFR = explained95_UFR;
      pcaData.explainedThird_UFR = explainedThird_UFR;
      pcaData.explainedTwoThirds_UFR = explainedTwoThirds_UFR;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_UFR = rPCs2pupilAreaPearson_UFR;
        pcaData.pvalPCs2pupilAreaPearson_UFR = pvalPCs2pupilAreaPearson_UFR;
        pcaData.rPCs2pupilAreaSpearman_UFR = rPCs2pupilAreaSpearman_UFR;
        pcaData.pvalPCs2pupilAreaSpearman_UFR = pvalPCs2pupilAreaSpearman_UFR;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_UFR = rPCs2motionPearson_UFR;
        pcaData.pvalPCs2motionPearson_UFR = pvalPCs2motionPearson_UFR;
        pcaData.rPCs2motionSpearman_UFR = rPCs2motionSpearman_UFR;
        pcaData.pvalPCs2motionSpearman_UFR = pvalPCs2motionSpearman_UFR;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_UFR = pcaCoef_area_UFR;
      pcaData.explained_area_UFR = explained_area_UFR;
      pcaData.explained25_area_UFR = explained25_area_UFR;
      pcaData.explained50_area_UFR = explained50_area_UFR;
      pcaData.explained75_area_UFR = explained75_area_UFR;
      pcaData.explained95_area_UFR = explained95_area_UFR;
      pcaData.explainedThird_area_UFR = explainedThird_area_UFR;
      pcaData.explainedTwoThirds_area_UFR = explainedTwoThirds_area_UFR;
      
      
      % SAVE FAST BAND DATA
      % PCA data for all areas
      pcaData.pcaCoef_FB = pcaCoef_FB;
      pcaData.explained_FB = explained_FB;
      pcaData.nPCs_FB = nPCs_FB;
      pcaData.prob_FB = prob_FB;
      pcaData.explained25_FB = explained25_FB;
      pcaData.explained50_FB = explained50_FB;
      pcaData.explained75_FB = explained75_FB;
      pcaData.explained95_FB = explained95_FB;
      pcaData.explainedThird_FB = explainedThird_FB;
      pcaData.explainedTwoThirds_FB = explainedTwoThirds_FB;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_FB = rPCs2pupilAreaPearson_FB;
        pcaData.pvalPCs2pupilAreaPearson_FB = pvalPCs2pupilAreaPearson_FB;
        pcaData.rPCs2pupilAreaSpearman_FB = rPCs2pupilAreaSpearman_FB;
        pcaData.pvalPCs2pupilAreaSpearman_FB = pvalPCs2pupilAreaSpearman_FB;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_FB = rPCs2motionPearson_FB;
        pcaData.pvalPCs2motionPearson_FB = pvalPCs2motionPearson_FB;
        pcaData.rPCs2motionSpearman_FB = rPCs2motionSpearman_FB;
        pcaData.pvalPCs2motionSpearman_FB = pvalPCs2motionSpearman_FB;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_FB = pcaCoef_area_FB;
      pcaData.explained_area_FB = explained_area_FB;
      pcaData.explained25_area_FB = explained25_area_FB;
      pcaData.explained50_area_FB = explained50_area_FB;
      pcaData.explained75_area_FB = explained75_area_FB;
      pcaData.explained95_area_FB = explained95_area_FB;
      pcaData.explainedThird_area_FB = explainedThird_area_FB;
      pcaData.explainedTwoThirds_area_FB = explainedTwoThirds_area_FB;
      
      
      % SAVE LFP DATA
      % PCA data for all areas
      pcaData.pcaCoef_LFP = pcaCoef_LFP;
      pcaData.explained_LFP = explained_LFP;
      pcaData.nPCs_LFP = nPCs_LFP;
      pcaData.prob_LFP = prob_LFP;
      pcaData.explained25_LFP = explained25_LFP;
      pcaData.explained50_LFP = explained50_LFP;
      pcaData.explained75_LFP = explained75_LFP;
      pcaData.explained95_LFP = explained95_LFP;
      pcaData.explainedThird_LFP = explainedThird_LFP;
      pcaData.explainedTwoThirds_LFP = explainedTwoThirds_LFP;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_LFP = rPCs2pupilAreaPearson_LFP;
        pcaData.pvalPCs2pupilAreaPearson_LFP = pvalPCs2pupilAreaPearson_LFP;
        pcaData.rPCs2pupilAreaSpearman_LFP = rPCs2pupilAreaSpearman_LFP;
        pcaData.pvalPCs2pupilAreaSpearman_LFP = pvalPCs2pupilAreaSpearman_LFP;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_LFP = rPCs2motionPearson_LFP;
        pcaData.pvalPCs2motionPearson_LFP = pvalPCs2motionPearson_LFP;
        pcaData.rPCs2motionSpearman_LFP = rPCs2motionSpearman_LFP;
        pcaData.pvalPCs2motionSpearman_LFP = pvalPCs2motionSpearman_LFP;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_LFP = pcaCoef_area_LFP;
      pcaData.explained_area_LFP = explained_area_LFP;
      pcaData.explained25_area_LFP = explained25_area_LFP;
      pcaData.explained50_area_LFP = explained50_area_LFP;
      pcaData.explained75_area_LFP = explained75_area_LFP;
      pcaData.explained95_area_LFP = explained95_area_LFP;
      pcaData.explainedThird_area_LFP = explainedThird_area_LFP;
      pcaData.explainedTwoThirds_area_LFP = explainedTwoThirds_area_LFP;
      
      
      % SAVE PR DATA
      % PCA data for all areas
      pcaData.pcaCoef_PR = pcaCoef_PR;
      pcaData.explained_PR = explained_PR;
      pcaData.nPCs_PR = nPCs_PR;
      pcaData.prob_PR = prob_PR;
      pcaData.explained25_PR = explained25_PR;
      pcaData.explained50_PR = explained50_PR;
      pcaData.explained75_PR = explained75_PR;
      pcaData.explained95_PR = explained95_PR;
      pcaData.explainedThird_PR = explainedThird_PR;
      pcaData.explainedTwoThirds_PR = explainedTwoThirds_PR;
      if ~isempty(pupilArea)
        pcaData.rPCs2pupilAreaPearson_PR = rPCs2pupilAreaPearson_PR;
        pcaData.pvalPCs2pupilAreaPearson_PR = pvalPCs2pupilAreaPearson_PR;
        pcaData.rPCs2pupilAreaSpearman_PR = rPCs2pupilAreaSpearman_PR;
        pcaData.pvalPCs2pupilAreaSpearman_PR = pvalPCs2pupilAreaSpearman_PR;
      end
      if ~isempty(motion)
        pcaData.rPCs2motionPearson_PR = rPCs2motionPearson_PR;
        pcaData.pvalPCs2motionPearson_PR = pvalPCs2motionPearson_PR;
        pcaData.rPCs2motionSpearman_PR = rPCs2motionSpearman_PR;
        pcaData.pvalPCs2motionSpearman_PR = pvalPCs2motionSpearman_PR;
      end
      
      % PCA data for the area of interest
      pcaData.pcaCoef_area_PR = pcaCoef_area_PR;
      pcaData.explained_area_PR = explained_area_PR;
      pcaData.explained25_area_PR = explained25_area_PR;
      pcaData.explained50_area_PR = explained50_area_PR;
      pcaData.explained75_area_PR = explained75_area_PR;
      pcaData.explained95_area_PR = explained95_area_PR;
      pcaData.explainedThird_area_PR = explainedThird_area_PR;
      pcaData.explainedTwoThirds_area_PR = explainedTwoThirds_area_PR;
      
      dataString = ['dataStruct.seriesData.' entryNames{iArea} '.pcaData4 = pcaData;'];
      eval(dataString);
      if intermediateSaving
        save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
      end
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

function dataOI = getDataOI(areas, areasOI, iArea)

dataOI = areas;
dataOI(dataOI ~= areasOI(iArea)) = 0;
dataOI = logical(dataOI);
end