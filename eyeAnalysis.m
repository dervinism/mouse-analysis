% Run this script to perform analyses of correlation, phase, and coherence
% between unit neural spiking and pupil area.


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


%% CORRELATE SPIKING ACTIVITY WITH EYE MEASURES
fnsData = fieldnames(dataStruct.seriesData);
if ~isempty(dbEntries) && dbEntries(1) == inf % Here you can chose to execute only certain db entries
  dbEntriesLocal = 1:numel(fnsData);
else
  dbEntriesLocal = dbEntries;
end
for dbCount = dbEntriesLocal % Loop through dbEntries
  
  % Load the contents of dbStruct
  [dbStruct, repository, ~, entryName, ~, ~, shankIDs,...
    ~, period, ~, ~, srData, optCoh, ~, FOI,...
    MUAsAll] = get_dbStruct(dataStruct, dbCount);
  optCoh.maxFreq = maxFreq_pupil;
  optCoh.winfactor = winfactor;
  optCoh.freqfactor = freqfactor;
  optCoh.monotoneFreq = true;
  MUAsAll = sum(MUAsAll,1);
  if isempty(MUAsAll) || ~sum(sum(MUAsAll))
    disp(['No spiking data for ' fnsData{dbCount} '. Skippig to the next db entry...']);
    continue
  end
  
  % Get eye data
  [seriesName, animal] = seriesFromEntry(entryName);
  entryNameEye = [animal '_s' seriesName(1:min([numel(seriesName) 14]))];
  if ~isfield(dataStruct, 'eyeData') || ~isfield(dataStruct.eyeData, entryNameEye) ||...
      isempty(dataStruct.eyeData.(entryNameEye).pupilArea)
    disp(['No pupil data for ' entryName '. Skippig to the next db entry...']);
    continue
  end
  eyeDataDB = dataStruct.eyeData.(entryNameEye);
  
% INTERPOLATE AND FILTER PUPIL AREA DATA
  [seriesName, animal] = seriesFromEntry(entryName);
  runningSpeedDataEntry = [animal '_s' seriesName(1:min([14 numel(seriesName)]))];
  if excludeRunning && isfield(dataStruct, 'runningSpeedData') && ~isempty(dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod)
    period = combinePeriods(period, dataStruct.runningSpeedData.(runningSpeedDataEntry).maxQuietPeriod, srData);
  end
  commonPeriod = combinePeriods(period, eyeDataDB.period, srData);
  if isempty(commonPeriod)
    continue
  end
  [areaInterpFilt, interpTimes, areaInd] = pupilFilt(eyeDataDB, srData, MUAsAll, 2*srData, commonPeriod, srData);

% CORRELATE PUPIL AREA AND UNFILTERED POPULATION SPIKING
  [pupilDB.popData.rPup, pupilDB.popData.pvalPup] = corrSimple(areaInterpFilt, MUAsAll(areaInd));

% CORRELATE PUPIL AREA AND FILTERED POPULATION SPIKING
  filt.filt = true;
  filt.FOI = FOI;
  filt.sr = srData;
  filtData = filterFOI(MUAsAll, filt);
  [pupilDB.popData.rPupFOI, pupilDB.popData.pvalPupFOI] = corrMulti(areaInterpFilt, filtData(:,areaInd));
    
% CORRELATE PUPIL AREA AND HILBERT AMPLITUDE OF FILTERED POPULATION SPIKING
  hAmp = zeros(size(filtData));
  for iF = 1:numel(FOI)
    [~, ~, hAmp(iF,:)] = hilbertTransform(filtData(iF,:));
  end
  [pupilDB.popData.rPupHilbert, pupilDB.popData.pvalPupHilbert] = corrMulti(areaInterpFilt, hAmp(:,areaInd));
%   zhAmp = zscore(hAmp')';
%   [rzHilbert, pvalzHilbert] = corrMulti(areaInterpFilt, zhAmp(:,areaInd));
%   for iF = 1:numel(FOI)
%     plot([hAmp(iF,areaInd)' zhAmp(iF,areaInd)'])
%   end

% PHASE AND COHERENCE ANALYSIS OF POPULATION SPIKING RELATIVE TO PUPIL AREA
  optPup = optCoh;
  optPup.typespk1 = 'c';
  [pupilDB.popData.mfr, pupilDB.popData.mfr_1sthalf, pupilDB.popData.mfr_2ndhalf, pupilDB.popData.lfr1,...
    pupilDB.popData.lfr5] = rateCalc(MUAsAll(areaInd), srData);
  [pupilDB.pupData.psd_halves_freq, pupilDB.pupData.psd_halves, pupilDB.pupData.freq, pupilDB.pupData.psd,...
    pupilDB.pupData.psd_conf, pupilDB.pupData.psd_numelSignal] = psdCalc(areaInterpFilt, srData, optPup);
  [pupilDB.popData.psd_halves_freq, pupilDB.popData.psd_halves, freqPSD, pupilDB.popData.psd,...
    pupilDB.popData.psd_conf, pupilDB.popData.psd_numelSignal] = psdCalc(MUAsAll(areaInd), srData, optCoh);
  [pupilDB.popData.freq, pupilDB.popData.coh, pupilDB.popData.phase, pupilDB.popData.coh_conf, pupilDB.popData.phase_confU,...
    pupilDB.popData.phase_confL, pupilDB.popData.coh_halves_freq, pupilDB.popData.coh_halves, pupilDB.popData.coh_conf_halves,...
    pupilDB.popData.phase_halves, pupilDB.popData.phase_conf_halves] = phaseCohCalc(areaInterpFilt,...
    MUAsAll(areaInd), srData, optPup);
  [pupilDB.popData.rateadjust_kappa, pupilDB.popData.rateadjust_kappa_halves] = kappaCalc(pupilDB.popData.mfr,...
    pupilDB.popData.mfr_1sthalf, pupilDB.popData.mfr_2ndhalf, pupilDB.popData.psd, pupilDB.popData.psd_halves,...
    pupilDB.popData.coh, pupilDB.popData.coh_halves);
  [pupilDB.popData.psd_halves_freq, pupilDB.popData.psd_halves, pupilDB.popData.freq, pupilDB.popData.psd, pupilDB.popData.psd_conf,...
    pupilDB.popData.coh, pupilDB.popData.phase, pupilDB.popData.coh_conf, pupilDB.popData.phase_confU, pupilDB.popData.phase_confL,...
    pupilDB.popData.coh_halves_freq, pupilDB.popData.coh_halves, pupilDB.popData.coh_conf_halves, pupilDB.popData.phase_halves,...
    pupilDB.popData.phase_conf_halves, pupilDB.popData.rateadjust_kappa, pupilDB.popData.rateadjust_kappa_halves] = correctFreq(...
    pupilDB.popData.psd_halves_freq, pupilDB.popData.psd_halves, freqPSD, pupilDB.popData.psd, pupilDB.popData.psd_conf,...
    pupilDB.popData.freq, pupilDB.popData.coh, pupilDB.popData.phase, pupilDB.popData.coh_conf, pupilDB.popData.phase_confU,...
    pupilDB.popData.phase_confL, pupilDB.popData.coh_halves_freq, pupilDB.popData.coh_halves, pupilDB.popData.coh_conf_halves,...
    pupilDB.popData.phase_halves, pupilDB.popData.phase_conf_halves, pupilDB.popData.rateadjust_kappa, pupilDB.popData.rateadjust_kappa_halves);
    
% CROSSCORRELATE POPULATION SPIKING AND PUPIL AREA
  stPRwindow = 1e2;
  pupilDB.popData.stPRglobal = stprCalc(MUAsAll(areaInd), areaInterpFilt, stPRwindow);
  pupilDB.popData.stPR = stprHalfCalc(MUAsAll(areaInd), areaInterpFilt, stPRwindow);

% COMPARE MEAN VS VARIANCE USING DIFFERENT BIN SIZES
  [pupilDB.popData.m, pupilDB.popData.mv1, pupilDB.popData.mv2, pupilDB.popData.mean_var,...
    pupilDB.popData.mean_var_timeBin] = meanVarCalc(MUAsAll(areaInd), srData);
  
  pupilDB.unitData.phaseFOI = []; 
  pupilDB.unitData.cohFOI = [];
  pupilDB.unitData.coh_confFOI = [];
  pupilDB.unitData.phase_1sthalfFOI = []; 
  pupilDB.unitData.coh_1sthalfFOI = [];
  pupilDB.unitData.coh_1sthalf_confFOI = [];
  pupilDB.unitData.phase_2ndhalfFOI = []; 
  pupilDB.unitData.coh_2ndhalfFOI = [];
  pupilDB.unitData.coh_2ndhalf_confFOI = [];
  
% LOOP THROUGH SHANKS
  for sh = 1:numel(shankIDs)
    shankStruct = eval(['dbStruct.shankData.' shankIDs{sh}]);
    fprintf('%s shank %d -------------------\n', fnsData{dbCount}, sh);
    
% LOAD THE CONTENTS OF THE SHANK STRUCTURE VARIABLE
    units = shankStruct.units;
    if isempty(units)
      dataString = ['dataStruct.seriesData.' entryName '.popData.pupil = pupilDB;'];
      eval(dataString);
      continue
    end
    spk = shankStruct.spk;
    MUAs = shankStruct.MUAs;

% CORRELATE PUPIL AREA AND UNFILTERED POPULATION SPIKING
    [pupil.popData.rPup, pupil.popData.pvalPup] = corrSimple(areaInterpFilt, MUAs(areaInd));

% CORRELATE PUPIL AREA AND FILTERED POPULATION SPIKING
    filt.filt = true;
    filt.FOI = FOI;
    filt.sr = srData;
    filtData = filterFOI(MUAs, filt);
    [pupil.popData.rPupFOI, pupil.popData.pvalPupFOI] = corrMulti(areaInterpFilt, filtData(:,areaInd));
    
% CORRELATE PUPIL AREA AND HILBERT AMPLITUDE OF FILTERED POPULATION SPIKING
    hAmp = zeros(size(filtData));
    for iF = 1:numel(FOI)
      [~, ~, hAmp(iF,:)] = hilbertTransform(filtData(iF,:));
    end
    [pupil.popData.rPupHilbert, pupil.popData.pvalPupHilbert] = corrMulti(areaInterpFilt, hAmp(:,areaInd));
%     zhAmp = zscore(hAmp')';
%     [rzHilbert, pvalzHilbert] = corrMulti(areaInterpFilt, zhAmp(:,areaInd));
%     for iF = 1:numel(FOI)
%       plot([hAmp(iF,areaInd)' zhAmp(iF,areaInd)'])
%     end

% PHASE AND COHERENCE ANALYSIS OF POPULATION SPIKING RELATIVE TO PUPIL AREA
    [pupil.popData.mfr, pupil.popData.mfr_1sthalf, pupil.popData.mfr_2ndhalf, pupil.popData.lfr1,...
      pupil.popData.lfr5] = rateCalc(MUAs(areaInd), srData);
    [pupil.pupData.psd_halves_freq, pupil.pupData.psd_halves, pupil.pupData.freq, pupil.pupData.psd, pupil.pupData.psd_conf,...
      pupil.pupData.psd_numelSignal] = psdCalc(areaInterpFilt, srData, optPup);
    [pupil.popData.psd_halves_freq, pupil.popData.psd_halves, pupil.popData.freq, pupil.popData.psd, pupil.popData.psd_conf,...
      pupil.popData.psd_numelSignal] = psdCalc(MUAs(areaInd), srData, optCoh);
    [pupil.popData.freq, pupil.popData.coh, pupil.popData.phase, pupil.popData.coh_conf, pupil.popData.phase_confU,...
      pupil.popData.phase_confL, pupil.popData.coh_halves_freq, pupil.popData.coh_halves, pupil.popData.coh_conf_halves,...
      pupil.popData.phase_halves, pupil.popData.phase_conf_halves] = phaseCohCalc(areaInterpFilt, MUAs(areaInd),...
      srData, optPup);
    [pupil.popData.rateadjust_kappa, pupil.popData.rateadjust_kappa_halves] = kappaCalc(pupil.popData.mfr,...
      pupil.popData.mfr_1sthalf, pupil.popData.mfr_2ndhalf, pupil.popData.psd, pupil.popData.psd_halves, pupil.popData.coh,...
      pupil.popData.coh_halves);
    [pupil.popData.psd_halves_freq, pupil.popData.psd_halves, pupil.popData.freq, pupil.popData.psd, pupil.popData.psd_conf,...
      pupil.popData.coh, pupil.popData.phase, pupil.popData.coh_conf, pupil.popData.phase_confU, pupil.popData.phase_confL,...
      pupil.popData.coh_halves_freq, pupil.popData.coh_halves, pupil.popData.coh_conf_halves, pupil.popData.phase_halves,...
      pupil.popData.phase_conf_halves, pupil.popData.rateadjust_kappa, pupil.popData.rateadjust_kappa_halves] = correctFreq(...
      pupil.popData.psd_halves_freq, pupil.popData.psd_halves, freqPSD, pupil.popData.psd, pupil.popData.psd_conf,...
      pupil.popData.freq, pupil.popData.coh, pupil.popData.phase, pupil.popData.coh_conf, pupil.popData.phase_confU,...
      pupil.popData.phase_confL, pupil.popData.coh_halves_freq, pupil.popData.coh_halves, pupil.popData.coh_conf_halves,...
      pupil.popData.phase_halves, pupil.popData.phase_conf_halves, pupil.popData.rateadjust_kappa, pupil.popData.rateadjust_kappa_halves);
    
% CROSSCORRELATE POPULATION SPIKING AND PUPIL AREA
    stPRwindow = 1e2;
    pupil.popData.stPRglobal = stprCalc(MUAs(areaInd), areaInterpFilt, stPRwindow);
    pupil.popData.stPR = stprHalfCalc(MUAs(areaInd), areaInterpFilt, stPRwindow);

% COMPARE MEAN VS VARIANCE USING DIFFERENT BIN SIZES
    [pupil.popData.m, pupil.popData.mv1, pupil.popData.mv2, pupil.popData.mean_var,...
      pupil.popData.mean_var_timeBin] = meanVarCalc(MUAs(areaInd), srData);

% LOOP THROUGH UNITS
    pupil.unitData = {};
    parfor u = 1:numel(units)
    %for u = 1:numel(units)
      fprintf('Started processing unit %i\n',units(u));
      
      unitData = struct();
      unitData.unit = units(u);
      spkOI = full(spk(u,:));
      
% CORRELATE PUPIL AREA AND UNFILTERED UNIT SPIKING
      [unitData.rPup, unitData.pvalPup] = corrSimple(areaInterpFilt, spkOI(areaInd));

% CORRELATE PUPIL AREA AND FILTERED UNIT SPIKING
      filtData = filterFOI(spkOI, filt);
      [unitData.rPupFOI, unitData.pvalPupFOI] = corrMulti(areaInterpFilt, filtData(:,areaInd));
    
% CORRELATE PUPIL AREA AND HILBERT AMPLITUDE OF FILTERED UNIT SPIKING
      hAmp = zeros(size(filtData));
      for iF = 1:numel(FOI)
        [~, ~, hAmp(iF,:)] = hilbertTransform(filtData(iF,:));
      end
      [unitData.rPupHilbert, unitData.pvalPupHilbert] = corrMulti(areaInterpFilt, hAmp(:,areaInd));
%       zhAmp = zscore(hAmp')';
%       [rzHilbert, pvalzHilbert] = corrMulti(areaInterpFilt, zhAmp(:,areaInd));
%       for iF = 1:numel(FOI)
%         plot([hAmp(iF,areaInd)' zhAmp(iF,areaInd)'])
%       end

% PHASE AND COHERENCE ANALYSIS OF UNIT SPIKING RELATIVE TO PUPIL AREA
      [unitData.mfr, unitData.mfr_1sthalf, unitData.mfr_2ndhalf, unitData.lfr1, unitData.lfr5] = rateCalc(spkOI(areaInd),...
        srData);
      [unitData.psd_halves_freq, unitData.psd_halves, unitData.freq, unitData.psd, unitData.psd_conf,...
        unitData.psd_numelSignal] = psdCalc(spkOI(areaInd), srData, optCoh);
      [unitData.freq, unitData.coh, unitData.phase, unitData.coh_conf, unitData.phase_confU, unitData.phase_confL,...
        unitData.coh_halves_freq, unitData.coh_halves, unitData.coh_conf_halves, unitData.phase_halves,...
        unitData.phase_conf_halves] = phaseCohCalc(areaInterpFilt, spkOI(areaInd), srData, optCoh);
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

% CROSSCORRELATE UNIT SPIKING AND PUPIL AREA
      unitData.stPRglobal = stprCalc(spkOI(areaInd), areaInterpFilt, stPRwindow);
      unitData.stPR = stprHalfCalc(spkOI(areaInd), areaInterpFilt, stPRwindow);
      
% COMPARE MEAN VS VARIANCE USING DIFFERENT BIN SIZES
      [unitData.m, unitData.mv1, unitData.mv2, unitData.mean_var, unitData.mean_var_timeBin] = meanVarCalc(spkOI(areaInd),...
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
      saveParfor([entryName '_' unitNamer(u) '_pupil.mat'], unitData);
      fprintf('Finished processing unit %i\n',units(u));
    end
    
% RETRIEVE TEMPORARY DATA
    fileList = dir([entryName '_*_pupil.mat']);
    if ~isempty(fileList)
      for i = 1:size(fileList,1)
        qPupil = load(fileList(i).name);
        pupil.unitData{end+1} = qPupil.q;
        pupilDB.unitData.phaseFOI = [pupilDB.unitData.phaseFOI; pupil.unitData{i}.phaseFOI];
        pupilDB.unitData.cohFOI = [pupilDB.unitData.cohFOI; pupil.unitData{i}.cohFOI];
        pupilDB.unitData.coh_confFOI = [pupilDB.unitData.coh_confFOI; pupil.unitData{i}.coh_confFOI];
        pupilDB.unitData.phase_1sthalfFOI = [pupilDB.unitData.phase_1sthalfFOI; pupil.unitData{i}.phase_1sthalfFOI];
        pupilDB.unitData.coh_1sthalfFOI = [pupilDB.unitData.coh_1sthalfFOI; pupil.unitData{i}.coh_1sthalfFOI];
        pupilDB.unitData.coh_1sthalf_confFOI = [pupilDB.unitData.coh_1sthalf_confFOI; pupil.unitData{i}.coh_1sthalf_confFOI];
        pupilDB.unitData.phase_2ndhalfFOI = [pupilDB.unitData.phase_2ndhalfFOI; pupil.unitData{i}.phase_2ndhalfFOI];
        pupilDB.unitData.coh_2ndhalfFOI = [pupilDB.unitData.coh_2ndhalfFOI; pupil.unitData{i}.coh_2ndhalfFOI];
        pupilDB.unitData.coh_2ndhalf_confFOI = [pupilDB.unitData.coh_2ndhalf_confFOI; pupil.unitData{i}.coh_2ndhalf_confFOI];
      end

% SAVE PERMANENT DATA
      dataString = ['dataStruct.seriesData.' entryName '.shankData.' shankIDs{sh} '.pupil = pupil;'];
      eval(dataString);
      dataString = ['dataStruct.seriesData.' entryName '.popData.pupil = pupilDB;'];
      eval(dataString);
      if intermediateSaving
        save(dataFile,'dataStruct','-v7.3'); %#ok<*UNRCH>
      end
      delete([entryName '_*_pupil.mat']);
    end
  end

% CORRELATE UNIT PHASE RELATIVE TO POPULATION SPIKING AND PUPIL AREA
  if isfield(dbStruct.popData, 'phaseCoh') && ~isempty(pupilDB.unitData.phaseFOI)
    [rPhasePupVsPop, pvalPhasePupVsPop] = corrMulti(pupilDB.unitData.phaseFOI, dbStruct.popData.phaseCoh.phaseFOI, 'circular');
    correctedCohPup = pupilDB.unitData.cohFOI;
    correctedCohPR = dbStruct.popData.phaseCoh.cohFOI;
    for u = 1:size(correctedCohPup,1)
      correctedCohPup(u, correctedCohPup(u,:) - pupilDB.unitData.coh_confFOI(u,:) <= 0) = NaN;
      correctedCohPR(u, correctedCohPR(u,:) - dbStruct.popData.phaseCoh.coh_confFOI(u,:) <= 0) = NaN;
    end
    [rCohPupVsPop, pvalCohPupVsPop] = corrMulti(correctedCohPup, correctedCohPR);
    dataString = ['dataStruct.seriesData.' entryName '.popData.rPhasePupVsPop = rPhasePupVsPop;'];
    eval(dataString);
    dataString = ['dataStruct.seriesData.' entryName '.popData.pvalPhasePupVsPop = pvalPhasePupVsPop;'];
    eval(dataString);
    dataString = ['dataStruct.seriesData.' entryName '.popData.rCohPupVsPop = rCohPupVsPop;'];
    eval(dataString);
    dataString = ['dataStruct.seriesData.' entryName '.popData.pvalCohPupVsPop = pvalCohPupVsPop;'];
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