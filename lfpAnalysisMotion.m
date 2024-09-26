% Run this script to perform analyses comparing LFP and total motion data.

%% LOAD PRE-PROCESSED DATA
load(dataFile);


%% INITIALISE PARAMETERS
visibility = 'on';
intermediateSaving = false;


%% CORRELATE LFP WITH MOTION MEASURES
figs = {};
fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf
  dbEntries = 1:numel(fnsData);
end
for dbCount = dbEntries % Loop though db entries
  
  % Load the contents of dbStruct
  [dbStruct, ~, ~, entryName, ~, ~, ~,...
    ~, ~, periodLFP, ~, srData] = get_dbStruct(dataStruct, dbCount);
  try
    motionDB = dbStruct.popData.motion;
  catch
    disp(['No movement data for ' fnsData{dbCount} '. Skipping to the next db entry...']);
    continue
  end
  
  % Get motion data
  [seriesName, animal] = seriesFromEntry(entryName);
  entryNameMotion = [animal '_s' seriesName(1:14)];
  if ~isfield(dataStruct, 'motionData') || ~isfield(dataStruct.eyeData, entryNameMotion)
    disp(['No motion data for ' entryName '. Skippig to the next db entry...']);
    continue
  end
  motionDataDB = dataStruct.motionData.(entryNameMotion);
  
  % Output folder
  strSep = strfind(entryName,'s');
  figsubdirname = entryName(strSep+1:end);
  if ~exist(figsubdirname,'dir')
    mkdir(figsubdirname)
  end
  
  % Interpolate and filter movement area data
  if ~isfield(dbStruct,'lfpPowerData')
      continue
  end
  dssrLFPinit = dbStruct.lfpPowerData.dssrLFPinit;
  dssrLFPfinal = dbStruct.lfpPowerData.dssrLFPfinal;
  lfpTimes = dbStruct.lfpPowerData.lfpTimes;
  [seriesName, animal] = seriesFromEntry(entryName);
  entryNameMotion = [animal '_s' seriesName(1:14)];
  commonPeriod = combinePeriods(periodLFP, dataStruct.motionData.(entryNameMotion).period{dbCount}, srData);
  if isempty(commonPeriod)
    continue
  end
  [motionInterpFilt, interpTimes, motionInd] = motionFilt(motionDataDB.sa, motionDataDB.frameTimes,...
    dssrLFPfinal, lfpTimes, 2*dssrLFPfinal, commonPeriod, dssrLFPfinal);
  
  % LOOP THROUGH LFP CHANNELS
  chOIDB = dbStruct.lfpPowerData.chOIDB;
  for iCh = 1:numel(chOIDB)
    
    % LOAD LFP MEASURES
    rippleRate = dbStruct.lfpPowerData.rippleRate{iCh}(motionInd);
    meanRippleRate = dbStruct.lfpPowerData.meanRippleRate{iCh};
    theta2deltaRatio = dbStruct.lfpPowerData.theta2deltaRatio{iCh}(motionInd);
    slowPower = dbStruct.lfpPowerData.slowPower{iCh}(motionInd);
    fastPower = dbStruct.lfpPowerData.fastPower{iCh}(motionInd);
    ultraFastPower = dbStruct.lfpPowerData.ultraFastPower{iCh}(motionInd);
    LFPsaturations = dbStruct.lfpPowerData.LFPsaturations{iCh};
    fLFPsaturations = dbStruct.lfpPowerData.fLFPsaturations{iCh};
    meanDurationLFPsaturations = dbStruct.lfpPowerData.meanDurationLFPsaturations{iCh};
    
    % CORRELATION ANALYSES
    lfpPhaseCohDB.rippleRate{iCh}.corr = zeros(2,2);
    [lfpPhaseCohDB.rippleRate{iCh}.corr(1,1), lfpPhaseCohDB.rippleRate{iCh}.corr(1,2)] = corrSimple(motionInterpFilt, rippleRate, 'Pearson');
    [lfpPhaseCohDB.rippleRate{iCh}.corr(2,1), lfpPhaseCohDB.rippleRate{iCh}.corr(2,2)] = corrSimple(motionInterpFilt, rippleRate, 'Spearman');
    lfpPhaseCohDB.t2dRatio{iCh}.corr = zeros(2,2);
    [lfpPhaseCohDB.t2dRatio{iCh}.corr(1,1), lfpPhaseCohDB.t2dRatio{iCh}.corr(1,2)] = corrSimple(motionInterpFilt, theta2deltaRatio, 'Pearson');
    [lfpPhaseCohDB.t2dRatio{iCh}.corr(2,1), lfpPhaseCohDB.t2dRatio{iCh}.corr(2,2)] = corrSimple(motionInterpFilt, theta2deltaRatio, 'Spearman');
    lfpPhaseCohDB.slowPower{iCh}.corr = zeros(2,2);
    [lfpPhaseCohDB.slowPower{iCh}.corr(1,1), lfpPhaseCohDB.slowPower{iCh}.corr(1,2)] = corrSimple(motionInterpFilt, slowPower, 'Pearson');
    [lfpPhaseCohDB.slowPower{iCh}.corr(2,1), lfpPhaseCohDB.slowPower{iCh}.corr(2,2)] = corrSimple(motionInterpFilt, slowPower, 'Spearman');
    lfpPhaseCohDB.fastPower{iCh}.corr = zeros(2,2);
    [lfpPhaseCohDB.fastPower{iCh}.corr(1,1), lfpPhaseCohDB.fastPower{iCh}.corr(1,2)] = corrSimple(motionInterpFilt, fastPower, 'Pearson');
    [lfpPhaseCohDB.fastPower{iCh}.corr(2,1), lfpPhaseCohDB.fastPower{iCh}.corr(2,2)] = corrSimple(motionInterpFilt, fastPower, 'Spearman');
    lfpPhaseCohDB.ultraFastPower{iCh}.corr = zeros(2,2);
    [lfpPhaseCohDB.ultraFastPower{iCh}.corr(1,1), lfpPhaseCohDB.ultraFastPower{iCh}.corr(1,2)] = corrSimple(motionInterpFilt, ultraFastPower, 'Pearson');
    [lfpPhaseCohDB.ultraFastPower{iCh}.corr(2,1), lfpPhaseCohDB.ultraFastPower{iCh}.corr(2,2)] = corrSimple(motionInterpFilt, ultraFastPower, 'Spearman');
    
    % Plot graphs
    %figure; plot(motionInterpFilt,rippleRate, 'r.', 'MarkerSize',1)
    %figure; plot(motionInterpFilt,theta2deltaRatio, 'r.', 'MarkerSize',1)
    %figure; plot(motionInterpFilt,slowPower, 'r.', 'MarkerSize',1)
    %figure; plot(motionInterpFilt,fastPower, 'r.', 'MarkerSize',1)
    %figure; plot(motionInterpFilt,ultraFastPower, 'r.', 'MarkerSize',1)
    figs{iCh} = figure('Visible', visibility); %#ok<*SAGROW>
    plot(interpTimes,motionInterpFilt, 'LineWidth',2); hold on;
    plot(interpTimes,rippleRate./mean(rippleRate))
    plot(interpTimes,theta2deltaRatio./mean(theta2deltaRatio))
    plot(interpTimes,slowPower./mean(slowPower))
    plot(interpTimes,fastPower./mean(fastPower))
    plot(interpTimes,ultraFastPower./mean(ultraFastPower))
    plot(lfpTimes(logical(LFPsaturations)), zeros(size(lfpTimes(logical(LFPsaturations)))), 'r.', 'MarkerSize',10)
    legend('Total movement','Ripple rate','Theta2delta ratio','Slow power','Fast power','Ultra fast power','LFP saturations')
    xlabel('Time (s)')
    ylabel('Z-score/Normalised signal')
    title(['LFP measures vs pupil: LFP saturation rate of ' num2str(fLFPsaturations) ' min^-^1 and mean duration of '...
      num2str(meanDurationLFPsaturations) ' s']);
    figName = [figsubdirname filesep entryName '_LFP_variousSignals_v_motion_' 'ch' num2str(chOIDB(iCh))];
    set(figs{iCh}, 'Name',figName);
    hgsave(figs{iCh}, figName);
    close(figs{iCh});
    
    % CALCULATE PHASE AND COHERENCE MEASURES
    optLFP = dbStruct.conf.optCoh;
    optLFP.winfactor = winfactor;
    optLFP.freqfactor = freqfactor;
    optLFP.monotoneFreq = true;
    optLFP.typespk1 = 'c';
    optLFP.typespk2 = 'c';
    optLFP.maxFreq = 0.5/(1/dssrLFPfinal);
    % Calculations
    [lfpPhaseCohDB.rippleRate{iCh}.freq, lfpPhaseCohDB.rippleRate{iCh}.coh, lfpPhaseCohDB.rippleRate{iCh}.phase,...
      lfpPhaseCohDB.rippleRate{iCh}.coh_conf, lfpPhaseCohDB.rippleRate{iCh}.phase_confU, lfpPhaseCohDB.rippleRate{iCh}.phase_confL,...
      lfpPhaseCohDB.rippleRate{iCh}.coh_halves_freq, lfpPhaseCohDB.rippleRate{iCh}.coh_halves, lfpPhaseCohDB.rippleRate{iCh}.coh_conf_halves,...
      lfpPhaseCohDB.rippleRate{iCh}.phase_halves, lfpPhaseCohDB.rippleRate{iCh}.phase_conf_halves, lfpPhaseCohDB.rippleRate{iCh}.stPRgr0,...
      lfpPhaseCohDB.rippleRate{iCh}.stPR, lfpPhaseCohDB.rippleRate{iCh}.m, lfpPhaseCohDB.rippleRate{iCh}.mv1, lfpPhaseCohDB.rippleRate{iCh}.mv2,...
      lfpPhaseCohDB.rippleRate{iCh}.mean_var, lfpPhaseCohDB.rippleRate{iCh}.mean_var_timeBin] = phaseCohLFP(motionInterpFilt, rippleRate,...
      dssrLFPfinal, optLFP, 1e2);
    [lfpPhaseCohDB.t2dRatio{iCh}.freq, lfpPhaseCohDB.t2dRatio{iCh}.coh, lfpPhaseCohDB.t2dRatio{iCh}.phase,...
      lfpPhaseCohDB.t2dRatio{iCh}.coh_conf, lfpPhaseCohDB.t2dRatio{iCh}.phase_confU, lfpPhaseCohDB.t2dRatio{iCh}.phase_confL,...
      lfpPhaseCohDB.t2dRatio{iCh}.coh_halves_freq, lfpPhaseCohDB.t2dRatio{iCh}.coh_halves, lfpPhaseCohDB.t2dRatio{iCh}.coh_conf_halves,...
      lfpPhaseCohDB.t2dRatio{iCh}.phase_halves, lfpPhaseCohDB.t2dRatio{iCh}.phase_conf_halves, lfpPhaseCohDB.t2dRatio{iCh}.stPRgr0,...
      lfpPhaseCohDB.t2dRatio{iCh}.stPR, lfpPhaseCohDB.t2dRatio{iCh}.m, lfpPhaseCohDB.t2dRatio{iCh}.mv1, lfpPhaseCohDB.t2dRatio{iCh}.mv2,...
      lfpPhaseCohDB.t2dRatio{iCh}.mean_var, lfpPhaseCohDB.t2dRatio{iCh}.mean_var_timeBin] = phaseCohLFP(motionInterpFilt, theta2deltaRatio,...
      dssrLFPfinal, optLFP, 1e2);
    [lfpPhaseCohDB.slowPower{iCh}.freq, lfpPhaseCohDB.slowPower{iCh}.coh, lfpPhaseCohDB.slowPower{iCh}.phase,...
      lfpPhaseCohDB.slowPower{iCh}.coh_conf, lfpPhaseCohDB.slowPower{iCh}.phase_confU, lfpPhaseCohDB.slowPower{iCh}.phase_confL,...
      lfpPhaseCohDB.slowPower{iCh}.coh_halves_freq, lfpPhaseCohDB.slowPower{iCh}.coh_halves, lfpPhaseCohDB.slowPower{iCh}.coh_conf_halves,...
      lfpPhaseCohDB.slowPower{iCh}.phase_halves, lfpPhaseCohDB.slowPower{iCh}.phase_conf_halves, lfpPhaseCohDB.slowPower{iCh}.stPRgr0,...
      lfpPhaseCohDB.slowPower{iCh}.stPR, lfpPhaseCohDB.slowPower{iCh}.m, lfpPhaseCohDB.slowPower{iCh}.mv1, lfpPhaseCohDB.slowPower{iCh}.mv2,...
      lfpPhaseCohDB.slowPower{iCh}.mean_var, lfpPhaseCohDB.slowPower{iCh}.mean_var_timeBin] = phaseCohLFP(motionInterpFilt, slowPower,...
      dssrLFPfinal, optLFP, 1e2);
    [lfpPhaseCohDB.fastPower{iCh}.freq, lfpPhaseCohDB.fastPower{iCh}.coh, lfpPhaseCohDB.fastPower{iCh}.phase,...
      lfpPhaseCohDB.fastPower{iCh}.coh_conf, lfpPhaseCohDB.fastPower{iCh}.phase_confU, lfpPhaseCohDB.fastPower{iCh}.phase_confL,...
      lfpPhaseCohDB.fastPower{iCh}.coh_halves_freq, lfpPhaseCohDB.fastPower{iCh}.coh_halves, lfpPhaseCohDB.fastPower{iCh}.coh_conf_halves,...
      lfpPhaseCohDB.fastPower{iCh}.phase_halves, lfpPhaseCohDB.fastPower{iCh}.phase_conf_halves, lfpPhaseCohDB.fastPower{iCh}.stPRgr0,...
      lfpPhaseCohDB.fastPower{iCh}.stPR, lfpPhaseCohDB.fastPower{iCh}.m, lfpPhaseCohDB.fastPower{iCh}.mv1, lfpPhaseCohDB.fastPower{iCh}.mv2,...
      lfpPhaseCohDB.fastPower{iCh}.mean_var, lfpPhaseCohDB.fastPower{iCh}.mean_var_timeBin] = phaseCohLFP(motionInterpFilt, fastPower,...
      dssrLFPfinal, optLFP, 1e2);
    [lfpPhaseCohDB.ultraFastPower{iCh}.freq, lfpPhaseCohDB.ultraFastPower{iCh}.coh, lfpPhaseCohDB.ultraFastPower{iCh}.phase,...
      lfpPhaseCohDB.ultraFastPower{iCh}.coh_conf, lfpPhaseCohDB.ultraFastPower{iCh}.phase_confU, lfpPhaseCohDB.ultraFastPower{iCh}.phase_confL,...
      lfpPhaseCohDB.ultraFastPower{iCh}.coh_halves_freq, lfpPhaseCohDB.ultraFastPower{iCh}.coh_halves, lfpPhaseCohDB.ultraFastPower{iCh}.coh_conf_halves,...
      lfpPhaseCohDB.ultraFastPower{iCh}.phase_halves, lfpPhaseCohDB.ultraFastPower{iCh}.phase_conf_halves, lfpPhaseCohDB.ultraFastPower{iCh}.stPRgr0,...
      lfpPhaseCohDB.ultraFastPower{iCh}.stPR, lfpPhaseCohDB.ultraFastPower{iCh}.m, lfpPhaseCohDB.ultraFastPower{iCh}.mv1, lfpPhaseCohDB.ultraFastPower{iCh}.mv2,...
      lfpPhaseCohDB.ultraFastPower{iCh}.mean_var, lfpPhaseCohDB.ultraFastPower{iCh}.mean_var_timeBin] = phaseCohLFP(motionInterpFilt, ultraFastPower,...
      dssrLFPfinal, optLFP, 1e2);
    
    % Graphs
    if sum(~isnan(lfpPhaseCohDB.rippleRate{iCh}.freq))
      figTmp = plotPhaseSpectrum(lfpPhaseCohDB.rippleRate{iCh}.freq', lfpPhaseCohDB.rippleRate{iCh}.coh, lfpPhaseCohDB.rippleRate{iCh}.phase,...
        lfpPhaseCohDB.rippleRate{iCh}.phase_confU, lfpPhaseCohDB.rippleRate{iCh}.phase_confL, lfpPhaseCohDB.rippleRate{iCh}.coh_conf, visibility);
      figName = [figsubdirname filesep entryName '_LFP_rippleRate_v_motion' '_PSDpolar_' 'ch' num2str(chOIDB(iCh))];
      set(figTmp, 'Name',figName);
      hgsave(figTmp, figName);
      close(figTmp);
      halfPhaseCohFig(lfpPhaseCohDB.rippleRate{iCh},...
        [figsubdirname filesep entryName '_LFP_rippleRate_v_motion' '_phaseCohHalves_' 'ch' num2str(chOIDB(iCh))], 'Phase', 'Coherence', visibility);
    end
    
    if sum(~isnan(lfpPhaseCohDB.t2dRatio{iCh}.freq))
      figTmp = plotPhaseSpectrum(lfpPhaseCohDB.t2dRatio{iCh}.freq', lfpPhaseCohDB.t2dRatio{iCh}.coh, lfpPhaseCohDB.t2dRatio{iCh}.phase,...
        lfpPhaseCohDB.t2dRatio{iCh}.phase_confU, lfpPhaseCohDB.t2dRatio{iCh}.phase_confL, lfpPhaseCohDB.t2dRatio{iCh}.coh_conf, visibility);
      figName = [figsubdirname filesep entryName '_LFP_t2dRatio_v_motion' '_PSDpolar_' 'ch' num2str(chOIDB(iCh))];
      set(figTmp, 'Name',figName);
      hgsave(figTmp, figName);
      close(figTmp);
      halfPhaseCohFig(lfpPhaseCohDB.t2dRatio{iCh},...
        [figsubdirname filesep entryName '_LFP_t2dRatio_v_motion' '_phaseCohHalves_' 'ch' num2str(chOIDB(iCh))], 'Phase', 'Coherence', visibility);
    end
    
    if sum(~isnan(lfpPhaseCohDB.slowPower{iCh}.freq))
      figTmp = plotPhaseSpectrum(lfpPhaseCohDB.slowPower{iCh}.freq', lfpPhaseCohDB.slowPower{iCh}.coh, lfpPhaseCohDB.slowPower{iCh}.phase,...
        lfpPhaseCohDB.slowPower{iCh}.phase_confU, lfpPhaseCohDB.slowPower{iCh}.phase_confL, lfpPhaseCohDB.slowPower{iCh}.coh_conf, visibility);
      figName = [figsubdirname filesep entryName '_LFP_slowPower_v_motion' '_PSDpolar_' 'ch' num2str(chOIDB(iCh))];
      set(figTmp, 'Name',figName);
      hgsave(figTmp, figName);
      close(figTmp);
      halfPhaseCohFig(lfpPhaseCohDB.slowPower{iCh},...
        [figsubdirname filesep entryName '_LFP_slowPower_v_motion' '_phaseCohHalves_' 'ch' num2str(chOIDB(iCh))], 'Phase', 'Coherence', visibility);
    end
    
    if sum(~isnan(lfpPhaseCohDB.fastPower{iCh}.freq))
      figTmp = plotPhaseSpectrum(lfpPhaseCohDB.fastPower{iCh}.freq', lfpPhaseCohDB.fastPower{iCh}.coh, lfpPhaseCohDB.fastPower{iCh}.phase,...
        lfpPhaseCohDB.fastPower{iCh}.phase_confU, lfpPhaseCohDB.fastPower{iCh}.phase_confL, lfpPhaseCohDB.fastPower{iCh}.coh_conf, visibility);
      figName = [figsubdirname filesep entryName '_LFP_fastPower_v_motion' '_PSDpolar_' 'ch' num2str(chOIDB(iCh))];
      set(figTmp, 'Name',figName);
      hgsave(figTmp, figName);
      close(figTmp);
      halfPhaseCohFig(lfpPhaseCohDB.fastPower{iCh},...
        [figsubdirname filesep entryName '_LFP_fastPower_v_motion' '_phaseCohHalves_' 'ch' num2str(chOIDB(iCh))], 'Phase', 'Coherence', visibility);
    end
    
    if sum(~isnan(lfpPhaseCohDB.ultraFastPower{iCh}.freq))
      figTmp = plotPhaseSpectrum(lfpPhaseCohDB.ultraFastPower{iCh}.freq', lfpPhaseCohDB.ultraFastPower{iCh}.coh, lfpPhaseCohDB.ultraFastPower{iCh}.phase,...
        lfpPhaseCohDB.ultraFastPower{iCh}.phase_confU, lfpPhaseCohDB.ultraFastPower{iCh}.phase_confL, lfpPhaseCohDB.ultraFastPower{iCh}.coh_conf, visibility);
      figName = [figsubdirname filesep entryName '_LFP_ultraFastPower_v_motion' '_PSDpolar_' 'ch' num2str(chOIDB(iCh))];
      set(figTmp, 'Name',figName);
      hgsave(figTmp, figName);
      close(figTmp);
      halfPhaseCohFig(lfpPhaseCohDB.ultraFastPower{iCh},...
        [figsubdirname filesep entryName '_LFP_ultraFastPower_v_motion' '_phaseCohHalves_' 'ch' num2str(chOIDB(iCh))], 'Phase', 'Coherence', visibility);
    end
  end
  
  % SAVE DATA
  dataString = ['dataStruct.seriesData.' entryName '.lfpPowerData.motionInterpFilt = motionInterpFilt;'];
  eval(dataString);
  dataString = ['dataStruct.seriesData.' entryName '.lfpPowerData.motionInterpTimes = interpTimes;'];
  eval(dataString);
  dataString = ['dataStruct.seriesData.' entryName '.lfpphaseCohDataMotion = lfpPhaseCohDB;'];
  eval(dataString);
  if intermediateSaving
    save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
  end
end

if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca



function [freq, coh, phase, coh_conf, phase_confU, phase_confL, coh_halves_freq, coh_halves, coh_conf_halves, phase_halves,...
  phase_conf_halves, stPRgr0, stPR, m, mv1, mv2, mean_var, mean_var_timeBin] = phaseCohLFP(motionSignal, LFPsignal,...
  dssrLFPfinal, optLFP, stPRwindow)

% PHASE AND COHERENCE ANALYSIS OF POPULATION SPIKING RELATIVE TO TOTAL MOVEMENT
[freq, coh, phase, coh_conf, phase_confU, phase_confL, coh_halves_freq, coh_halves, coh_conf_halves, phase_halves,...
  phase_conf_halves] = phaseCohCalc(motionSignal, LFPsignal, dssrLFPfinal, optLFP);

% CROSSCORRELATE POPULATION SPIKING AND TOTAL MOVEMENT
stPRgr0 = stprCalc(LFPsignal, motionSignal, stPRwindow);
stPR = stprHalfCalc(LFPsignal, motionSignal, stPRwindow);

% COMPARE MEAN VS VARIANCE USING DIFFERENT BIN SIZES
[m, mv1, mv2, mean_var, mean_var_timeBin] = meanVarCalc(LFPsignal, dssrLFPfinal);
end

function halfPhaseCohFig(phaseCohData, figName, titleStr1, titleStr2, visibility)

% Half phase and coherence comparisons
figTmp = figure('units', 'normalized', 'position', [0.002, .04, 1/2, .88/1.04], 'Visible', visibility);
subplot(2, 1, 1);
phaseCohData.phase_halves(1,isnan(phaseCohData.phase_conf_halves(1,:))) = NaN;
phaseCohData.phase_halves(1,isnan(phaseCohData.phase_conf_halves(2,:))) = NaN;
phaseCohData.phase_halves(2,isnan(phaseCohData.phase_conf_halves(3,:))) = NaN;
phaseCohData.phase_halves(2,isnan(phaseCohData.phase_conf_halves(4,:))) = NaN;
semilogx(phaseCohData.coh_halves_freq, phaseCohData.phase_halves(1,:), 'r');
hold on, semilogx(phaseCohData.coh_halves_freq, phaseCohData.phase_halves(2,:), 'm');
hold on, semilogx(phaseCohData.coh_halves_freq, phaseCohData.phase_conf_halves(1,:), '--r');
hold on, semilogx(phaseCohData.coh_halves_freq, phaseCohData.phase_conf_halves(2,:), '--r');
hold on, semilogx(phaseCohData.coh_halves_freq, phaseCohData.phase_conf_halves(3,:), '--m');
hold on, semilogx(phaseCohData.coh_halves_freq, phaseCohData.phase_conf_halves(4,:), '--m');
title(titleStr1, 'Interpreter', 'none')
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
legend('1st half', '2nd half');

subplot(2, 1, 2);
semilogx(phaseCohData.coh_halves_freq, phaseCohData.coh_halves(1,:), 'r');
hold on, semilogx(phaseCohData.coh_halves_freq, phaseCohData.coh_halves(2,:), 'm');
hold on, semilogx(phaseCohData.coh_halves_freq, (phaseCohData.coh_halves(1,:) + phaseCohData.coh_conf_halves(1,:)), '--r');
hold on, semilogx(phaseCohData.coh_halves_freq, (phaseCohData.coh_halves(1,:) - phaseCohData.coh_conf_halves(1,:)), '--r');
hold on, semilogx(phaseCohData.coh_halves_freq, (phaseCohData.coh_halves(2,:) + phaseCohData.coh_conf_halves(2,:)), '--m');
hold on, semilogx(phaseCohData.coh_halves_freq, (phaseCohData.coh_halves(2,:) - phaseCohData.coh_conf_halves(2,:)), '--m');
title(titleStr2, 'Interpreter', 'none')
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend('1st half', '2nd half');

set(figTmp, 'Name',figName);
hgsave(figTmp, figName);
close(figTmp);
end