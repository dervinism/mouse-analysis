% Run this script to perform analyses of correlation, phase, and coherence
% between unit neural spiking and total motion.

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


%% CORRELATE SPIKING ACTIVITY WITH MOTION MEASURES
fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf % Here you can chose to execute only certain db entries
  dbEntriesLocal = 1:numel(fnsData);
else
  dbEntriesLocal = dbEntries;
end
for dbCount = dbEntriesLocal % Loop through dbEntries
  
  % Load the contents of dbStruct
  [dbStruct, ~, ~, entryName, ~, ~, shankIDs,...
    ~, period, ~, ~, srData, optCoh, ~, FOI,...
    MUAsAll] = get_dbStruct(dataStruct, dbCount);
  optCoh.maxFreq = maxFreq_motion;
  optCoh.winfactor = winfactor;
  optCoh.freqfactor = freqfactor;
  optCoh.monotoneFreq = true;
  MUAsAll = sum(MUAsAll,1);
  if isempty(MUAsAll) || ~sum(sum(MUAsAll))
    disp(['No spiking data for ' fnsData{dbCount} '. Skippig to the next db entry...']);
    continue
  end
  
  % Get motion data
  [seriesName, animal] = seriesFromEntry(entryName);
  entryNameMotion = [animal '_s' seriesName(1:14)];
  if ~isfield(dataStruct, 'motionData') || ~isfield(dataStruct.motionData, entryNameMotion) ||...
      isempty(dataStruct.motionData.(entryNameMotion).s)
    disp(['No motion data for ' entryName '. Skippig to the next db entry...']);
    continue
  end
  motionDataDB = dataStruct.motionData.(entryNameMotion);
  
% INTERPOLATE AND FILTER MOTION DATA
  [seriesName, animal] = seriesFromEntry(entryName);
  runningSpeedDataEntry = [animal '_s' seriesName(1:min([14 numel(seriesName)]))];
  if excludeRunning && isfield(dataStruct, 'runningSpeedData') && ~isempty(dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod)
    period = combinePeriods(period, dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod, srData);
  end
  commonPeriod = combinePeriods(period, motionDataDB.period, srData);
  if isempty(commonPeriod)
    continue
  end
  [motionInterpFilt, interpTimes, motionInd] = motionFilt(motionDataDB.sa, motionDataDB.frameTimes, srData, MUAsAll, 2*srData, commonPeriod, srData);

% CORRELATE TOTAL MOTION AND UNFILTERED POPULATION SPIKING
  [motionDB.popData.rMot, motionDB.popData.pvalMot] = corrSimple(motionInterpFilt, MUAsAll(motionInd));

% CORRELATE TOTAL MOTION AND FILTERED POPULATION SPIKING
  filt.filt = true;
  filt.FOI = FOI;
  filt.sr = srData;
  filtData = filterFOI(MUAsAll, filt);
  [motionDB.popData.rMotFOI, motionDB.popData.pvalMotFOI] = corrMulti(motionInterpFilt, filtData(:,motionInd));
    
% CORRELATE TOTAL MOTION AND HILBERT AMPLITUDE OF FILTERED POPULATION SPIKING
  hAmp = zeros(size(filtData));
  for iF = 1:numel(FOI)
    [~, ~, hAmp(iF,:)] = hilbertTransform(filtData(iF,:));
  end
  [motionDB.popData.rMotHilbert, motionDB.popData.pvalMotHilbert] = corrMulti(motionInterpFilt, hAmp(:,motionInd));
%   zhAmp = zscore(hAmp')';
%   [rzHilbert, pvalzHilbert] = corrMulti(motionInterpFilt, zhAmp(:,motionInd));
%   for iF = 1:numel(FOI)
%     plot([hAmp(iF,motionInd)' zhAmp(iF,motionInd)'])
%   end

% PHASE AND COHERENCE ANALYSIS OF POPULATION SPIKING RELATIVE TO TOTAL MOTION
  optMot = optCoh;
  optMot.typespk1 = 'c';
  [motionDB.popData.mfr, motionDB.popData.mfr_1sthalf, motionDB.popData.mfr_2ndhalf, motionDB.popData.lfr1,...
    motionDB.popData.lfr5] = rateCalc(MUAsAll(motionInd), srData);
  [motionDB.motData.psd_halves_freq, motionDB.motData.psd_halves, motionDB.motData.freq, motionDB.motData.psd,...
    motionDB.motData.psd_conf, motionDB.motData.psd_numelSignal] = psdCalc(motionInterpFilt, srData, optMot);
  [motionDB.popData.psd_halves_freq, motionDB.popData.psd_halves, freqPSD, motionDB.popData.psd,...
    motionDB.popData.psd_conf, motionDB.popData.psd_numelSignal] = psdCalc(MUAsAll(motionInd), srData, optCoh);
  [motionDB.popData.freq, motionDB.popData.coh, motionDB.popData.phase, motionDB.popData.coh_conf, motionDB.popData.phase_confU,...
    motionDB.popData.phase_confL, motionDB.popData.coh_halves_freq, motionDB.popData.coh_halves, motionDB.popData.coh_conf_halves,...
    motionDB.popData.phase_halves, motionDB.popData.phase_conf_halves] = phaseCohCalc(motionInterpFilt,...
    MUAsAll(motionInd), srData, optMot);
  [motionDB.popData.rateadjust_kappa, motionDB.popData.rateadjust_kappa_halves] = kappaCalc(motionDB.popData.mfr,...
    motionDB.popData.mfr_1sthalf, motionDB.popData.mfr_2ndhalf, motionDB.popData.psd, motionDB.popData.psd_halves,...
    motionDB.popData.coh, motionDB.popData.coh_halves);
  [motionDB.popData.psd_halves_freq, motionDB.popData.psd_halves, motionDB.popData.freq, motionDB.popData.psd, motionDB.popData.psd_conf,...
    motionDB.popData.coh, motionDB.popData.phase, motionDB.popData.coh_conf, motionDB.popData.phase_confU, motionDB.popData.phase_confL,...
    motionDB.popData.coh_halves_freq, motionDB.popData.coh_halves, motionDB.popData.coh_conf_halves, motionDB.popData.phase_halves,...
    motionDB.popData.phase_conf_halves, motionDB.popData.rateadjust_kappa, motionDB.popData.rateadjust_kappa_halves] = correctFreq(...
    motionDB.popData.psd_halves_freq, motionDB.popData.psd_halves, freqPSD, motionDB.popData.psd, motionDB.popData.psd_conf,...
    motionDB.popData.freq, motionDB.popData.coh, motionDB.popData.phase, motionDB.popData.coh_conf, motionDB.popData.phase_confU,...
    motionDB.popData.phase_confL, motionDB.popData.coh_halves_freq, motionDB.popData.coh_halves, motionDB.popData.coh_conf_halves,...
    motionDB.popData.phase_halves, motionDB.popData.phase_conf_halves, motionDB.popData.rateadjust_kappa, motionDB.popData.rateadjust_kappa_halves);
    
% CROSSCORRELATE POPULATION SPIKING AND TOTAL MOTION
  stPRwindow = 1e2;
  motionDB.popData.stPRglobal = stprCalc(MUAsAll(motionInd), motionInterpFilt, stPRwindow);
  motionDB.popData.stPR = stprHalfCalc(MUAsAll(motionInd), motionInterpFilt, stPRwindow);

% COMPARE MEAN VS VARIANCE USING DIFFERENT BIN SIZES
  [motionDB.popData.m, motionDB.popData.mv1, motionDB.popData.mv2, motionDB.popData.mean_var,...
    motionDB.popData.mean_var_timeBin] = meanVarCalc(MUAsAll(motionInd), srData);
  
  motionDB.unitData.phaseFOI = []; 
  motionDB.unitData.cohFOI = [];
  motionDB.unitData.coh_confFOI = [];
  motionDB.unitData.phase_1sthalfFOI = []; 
  motionDB.unitData.coh_1sthalfFOI = [];
  motionDB.unitData.coh_1sthalf_confFOI = [];
  motionDB.unitData.phase_2ndhalfFOI = []; 
  motionDB.unitData.coh_2ndhalfFOI = [];
  motionDB.unitData.coh_2ndhalf_confFOI = [];
  
% LOOP THROUGH SHANKS
  for sh = 1:numel(shankIDs)
    shankStruct = eval(['dbStruct.shankData.' shankIDs{sh}]);
    fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
    
% LOAD THE CONTENTS OF THE SHANK STRUCTURE VARIABLE
    units = shankStruct.units;
    if isempty(units)
      dataString = ['dataStruct.seriesData.' entryName '.popData.motion = motionDB;'];
      eval(dataString);
      continue
    end
    spk = shankStruct.spk;
    MUAs = shankStruct.MUAs;

% CORRELATE TOTAL MOTION AND UNFILTERED POPULATION SPIKING
    [motion.popData.rMot, motion.popData.pvalMot] = corrSimple(motionInterpFilt, MUAs(motionInd));

% CORRELATE TOTAL MOTION AND FILTERED POPULATION SPIKING
    filt.filt = true;
    filt.FOI = FOI;
    filt.sr = srData;
    filtData = filterFOI(MUAs, filt);
    [motion.popData.rMotFOI, motion.popData.pvalMotFOI] = corrMulti(motionInterpFilt, filtData(:,motionInd));
    
% CORRELATE TOTAL MOTION AND HILBERT AMPLITUDE OF FILTERED POPULATION SPIKING
    hAmp = zeros(size(filtData));
    for iF = 1:numel(FOI)
      [~, ~, hAmp(iF,:)] = hilbertTransform(filtData(iF,:));
    end
    [motion.popData.rMotHilbert, motion.popData.pvalMotHilbert] = corrMulti(motionInterpFilt, hAmp(:,motionInd));
%     zhAmp = zscore(hAmp')';
%     [rzHilbert, pvalzHilbert] = corrMulti(motionInterpFilt, zhAmp(:,motionInd));
%     for iF = 1:numel(FOI)
%       plot([hAmp(iF,motionInd)' zhAmp(iF,motionInd)'])
%     end

% PHASE AND COHERENCE ANALYSIS OF POPULATION SPIKING RELATIVE TO TOTAL MOTION
    [motion.popData.mfr, motion.popData.mfr_1sthalf, motion.popData.mfr_2ndhalf, motion.popData.lfr1,...
      motion.popData.lfr5] = rateCalc(MUAs(motionInd), srData);
    [motion.motData.psd_halves_freq, motion.motData.psd_halves, motion.motData.freq, motion.motData.psd, motion.motData.psd_conf,...
      motion.motData.psd_numelSignal] = psdCalc(motionInterpFilt, srData, optMot);
    [motion.popData.psd_halves_freq, motion.popData.psd_halves, motion.popData.freq, motion.popData.psd, motion.popData.psd_conf,...
      motion.popData.psd_numelSignal] = psdCalc(MUAs(motionInd), srData, optCoh);
    [motion.popData.freq, motion.popData.coh, motion.popData.phase, motion.popData.coh_conf, motion.popData.phase_confU,...
      motion.popData.phase_confL, motion.popData.coh_halves_freq, motion.popData.coh_halves, motion.popData.coh_conf_halves,...
      motion.popData.phase_halves, motion.popData.phase_conf_halves] = phaseCohCalc(motionInterpFilt, MUAs(motionInd),...
      srData, optMot);
    [motion.popData.rateadjust_kappa, motion.popData.rateadjust_kappa_halves] = kappaCalc(motion.popData.mfr,...
      motion.popData.mfr_1sthalf, motion.popData.mfr_2ndhalf, motion.popData.psd, motion.popData.psd_halves, motion.popData.coh,...
      motion.popData.coh_halves);
    [motion.popData.psd_halves_freq, motion.popData.psd_halves, motion.popData.freq, motion.popData.psd, motion.popData.psd_conf,...
      motion.popData.coh, motion.popData.phase, motion.popData.coh_conf, motion.popData.phase_confU, motion.popData.phase_confL,...
      motion.popData.coh_halves_freq, motion.popData.coh_halves, motion.popData.coh_conf_halves, motion.popData.phase_halves,...
      motion.popData.phase_conf_halves, motion.popData.rateadjust_kappa, motion.popData.rateadjust_kappa_halves] = correctFreq(...
      motion.popData.psd_halves_freq, motion.popData.psd_halves, freqPSD, motion.popData.psd, motion.popData.psd_conf,...
      motion.popData.freq, motion.popData.coh, motion.popData.phase, motion.popData.coh_conf, motion.popData.phase_confU,...
      motion.popData.phase_confL, motion.popData.coh_halves_freq, motion.popData.coh_halves, motion.popData.coh_conf_halves,...
      motion.popData.phase_halves, motion.popData.phase_conf_halves, motion.popData.rateadjust_kappa, motion.popData.rateadjust_kappa_halves);
    
% CROSSCORRELATE POPULATION SPIKING AND TOTAL MOTION
    stPRwindow = 1e2;
    motion.popData.stPRglobal = stprCalc(MUAs(motionInd), motionInterpFilt, stPRwindow);
    motion.popData.stPR = stprHalfCalc(MUAs(motionInd), motionInterpFilt, stPRwindow);

% COMPARE MEAN VS VARIANCE USING DIFFERENT BIN SIZES
    [motion.popData.m, motion.popData.mv1, motion.popData.mv2, motion.popData.mean_var,...
      motion.popData.mean_var_timeBin] = meanVarCalc(MUAs(motionInd), srData);

% LOOP THROUGH UNITS
    motion.unitData = {};
    parfor u = 1:numel(units)
    %for u = 1:numel(units)
      fprintf('Started processing unit %i\n',units(u));
      
      unitData = struct();
      unitData.unit = units(u);
      spkOI = full(spk(u,:));
      
% CORRELATE TOTAL MOTION AND UNFILTERED UNIT SPIKING
      [unitData.rMot, unitData.pvalMot] = corrSimple(motionInterpFilt, spkOI(motionInd));

% CORRELATE TOTAL MOTION AND FILTERED UNIT SPIKING
      filtData = filterFOI(spkOI, filt);
      [unitData.rMotFOI, unitData.pvalMotFOI] = corrMulti(motionInterpFilt, filtData(:,motionInd));
    
% CORRELATE TOTAL MOTION AND HILBERT AMPLITUDE OF FILTERED UNIT SPIKING
      hAmp = zeros(size(filtData));
      for iF = 1:numel(FOI)
        [~, ~, hAmp(iF,:)] = hilbertTransform(filtData(iF,:));
      end
      [unitData.rMotHilbert, unitData.pvalMotHilbert] = corrMulti(motionInterpFilt, hAmp(:,motionInd));
%       zhAmp = zscore(hAmp')';
%       [rzHilbert, pvalzHilbert] = corrMulti(motionInterpFilt, zhAmp(:,motionInd));
%       for iF = 1:numel(FOI)
%         plot([hAmp(iF,motionInd)' zhAmp(iF,motionInd)'])
%       end

% PHASE AND COHERENCE ANALYSIS OF UNIT SPIKING RELATIVE TO TOTAL MOTION
      [unitData.mfr, unitData.mfr_1sthalf, unitData.mfr_2ndhalf, unitData.lfr1, unitData.lfr5] = rateCalc(spkOI(motionInd),...
        srData);
      [unitData.psd_halves_freq, unitData.psd_halves, unitData.freq, unitData.psd, unitData.psd_conf,...
        unitData.psd_numelSignal] = psdCalc(spkOI(motionInd), srData, optCoh);
      [unitData.freq, unitData.coh, unitData.phase, unitData.coh_conf, unitData.phase_confU, unitData.phase_confL,...
        unitData.coh_halves_freq, unitData.coh_halves, unitData.coh_conf_halves, unitData.phase_halves,...
        unitData.phase_conf_halves] = phaseCohCalc(motionInterpFilt, spkOI(motionInd), srData, optCoh);
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

% CROSSCORRELATE UNIT SPIKING AND TOTAL MOTION
      unitData.stPRglobal = stprCalc(spkOI(motionInd), motionInterpFilt, stPRwindow);
      unitData.stPR = stprHalfCalc(spkOI(motionInd), motionInterpFilt, stPRwindow);
      
% COMPARE MEAN VS VARIANCE USING DIFFERENT BIN SIZES
      [unitData.m, unitData.mv1, unitData.mv2, unitData.mean_var, unitData.mean_var_timeBin] = meanVarCalc(spkOI(motionInd),...
        srData);
      
% EXTRACT PHASE AND COHERENCE FOR EVERY FOI
      [unitData.phaseFOI, unitData.cohFOI, unitData.coh_confFOI, unitData.actualFOI, unitData.fInds] = phaseCohFOI(FOI,...
        unitData.freq, unitData.phase, unitData.coh, unitData.coh_conf, unitData.rateadjust_kappa);
      if isempty(unitData.coh_halves_freq)
        unitData.phase_1sthalfFOI = [];
        unitData.coh_1sthalfFOI = [];
        unitData.coh_1sthalf_confFOI = [];
        unitData.actualFOI_1sthalf = [];
        unitData.fInds_1sthalf = [];
        unitData.phase_2ndhalfFOI = [];
        unitData.coh_2ndhalfFOI = [];
        unitData.coh_2ndhalf_confFOI = [];
      else
        [unitData.phase_1sthalfFOI, unitData.coh_1sthalfFOI, unitData.coh_1sthalf_confFOI, unitData.actualFOI_1sthalf,...
          unitData.fInds_1sthalf] = phaseCohFOI(FOI, unitData.coh_halves_freq, unitData.phase_halves(1,:), unitData.coh_halves(1,:),...
          unitData.coh_conf_halves(1,:), unitData.rateadjust_kappa_halves(1,:));
        [unitData.phase_2ndhalfFOI, unitData.coh_2ndhalfFOI, unitData.coh_2ndhalf_confFOI] = phaseCohFOI(FOI,...
          unitData.coh_halves_freq, unitData.phase_halves(2,:), unitData.coh_halves(2,:), unitData.coh_conf_halves(2,:),...
          unitData.rateadjust_kappa_halves(2,:));
      end

% STORE TEMPORARY DATA
      saveParfor([entryName '_' unitNamer(u) '_motion.mat'], unitData);
      fprintf('Finished processing unit %i\n',units(u));
    end
    
% RETRIEVE TEMPORARY DATA
    fileList = dir([entryName '_*_motion.mat']);
    if ~isempty(fileList)
      for i = 1:size(fileList,1)
        qmotion = load(fileList(i).name);
        motion.unitData{end+1} = qmotion.q;
        motionDB.unitData.phaseFOI = [motionDB.unitData.phaseFOI; motion.unitData{i}.phaseFOI];
        motionDB.unitData.cohFOI = [motionDB.unitData.cohFOI; motion.unitData{i}.cohFOI];
        motionDB.unitData.coh_confFOI = [motionDB.unitData.coh_confFOI; motion.unitData{i}.coh_confFOI];
        motionDB.unitData.phase_1sthalfFOI = [motionDB.unitData.phase_1sthalfFOI; motion.unitData{i}.phase_1sthalfFOI];
        motionDB.unitData.coh_1sthalfFOI = [motionDB.unitData.coh_1sthalfFOI; motion.unitData{i}.coh_1sthalfFOI];
        motionDB.unitData.coh_1sthalf_confFOI = [motionDB.unitData.coh_1sthalf_confFOI; motion.unitData{i}.coh_1sthalf_confFOI];
        motionDB.unitData.phase_2ndhalfFOI = [motionDB.unitData.phase_2ndhalfFOI; motion.unitData{i}.phase_2ndhalfFOI];
        motionDB.unitData.coh_2ndhalfFOI = [motionDB.unitData.coh_2ndhalfFOI; motion.unitData{i}.coh_2ndhalfFOI];
        motionDB.unitData.coh_2ndhalf_confFOI = [motionDB.unitData.coh_2ndhalf_confFOI; motion.unitData{i}.coh_2ndhalf_confFOI];
      end

% SAVE PERMANENT DATA
      dataString = ['dataStruct.seriesData.' entryName '.shankData.' shankIDs{sh} '.motion = motion;'];
      eval(dataString);
      dataString = ['dataStruct.seriesData.' entryName '.popData.motion = motionDB;'];
      eval(dataString);
      if intermediateSaving
        save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
      end
      delete([entryName '_*_motion.mat']);
    end
  end

% CORRELATE UNIT PHASE RELATIVE TO POPULATION SPIKING AND TOTAL MOTION
  if isfield(dbStruct.popData, 'phaseCoh') && ~isempty(motionDB.unitData.phaseFOI)
    [rPhaseMotVsPop, pvalPhaseMotVsPop] = corrMulti(motionDB.unitData.phaseFOI, dbStruct.popData.phaseCoh.phaseFOI, 'circular');
    correctedCohMot = motionDB.unitData.cohFOI;
    correctedCohPR = dbStruct.popData.phaseCoh.cohFOI;
    for u = 1:size(correctedCohMot,1)
      correctedCohMot(u, correctedCohMot(u,:) - motionDB.unitData.coh_confFOI(u,:) <= 0) = NaN;
      correctedCohPR(u, correctedCohPR(u,:) - dbStruct.popData.phaseCoh.coh_confFOI(u,:) <= 0) = NaN;
    end
    [rCohMotVsPop, pvalCohMotVsPop] = corrMulti(correctedCohMot, correctedCohPR);
    dataString = ['dataStruct.seriesData.' entryName '.popData.rPhaseMotVsPop = rPhaseMotVsPop;'];
    eval(dataString);
    dataString = ['dataStruct.seriesData.' entryName '.popData.pvalPhaseMotVsPop = pvalPhaseMotVsPop;'];
    eval(dataString);
    dataString = ['dataStruct.seriesData.' entryName '.popData.rCohMotVsPop = rCohMotVsPop;'];
    eval(dataString);
    dataString = ['dataStruct.seriesData.' entryName '.popData.pvalCohMotVsPop = pvalCohMotVsPop;'];
    eval(dataString);
    if intermediateSaving
      save(dataFile,'dataStruct','-v7.3');
    end
  end
end

if ~intermediateSaving
  save(dataFile,'dataStruct','-v7.3');
end
clearvars -except dataFile dbEntries dbEntries_c dbEntries_ca dataStruct
