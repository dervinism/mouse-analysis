clearvars -except repository subpop reverse qualityCheck allData fullRun


%% INITIALISE PARAMETERS
params
lists

if ~exist('repository', 'var')
  repository = 'uol';
end
if ~exist('subpop', 'var')
  subpop = 'all';
end
if ~exist('fullRun', 'var')
  fullRun = true;
end
if ~exist('qualityCheck', 'var')
  qualityCheck = false;
end

if strcmp(repository,'all')
  error('globalConditionsAnalysis is only available for uol repository')
elseif strcmp(repository,'uol')
  if strcmp(subpop, 'all')
    rootFolder = [outputDir filesep laDir_uol];
  elseif strcmp(subpop, 'positive')
    rootFolder = [outputDir filesep laDir_uol_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [outputDir filesep laDir_uol_negative];
  end
  animals = animalsUOLOI;
elseif strcmp(repository,'allensdk')
  error('globalConditionsAnalysis is only available for uol repository')
end

if qualityCheck
  rootFolder = [rootFolder filesep qualityUnitsFolder];
  mainFolder = [rootFolder filesep conditionsSubfolder]; %#ok<*UNRCH>
else
  rootFolder = [rootFolder filesep unitsFolder];
  mainFolder = [rootFolder filesep conditionsSubfolder];
end


%% COMPUTE VARIABLES AND INITIALISE FIGURES NEEDED FOR COMPARING BETWEEN VIGILANCE CONDITIONS
if fullRun
  for animal = 1:numel(animals) % Loop through animals
    disp(['processing animal ' animals{animal}]);
    if exist('allData', 'var')
      dataStruct = allData.(animals{animal}).dataStruct;
    else
      load([dataDir_local filesep animals{animal} filesep animals{animal} '.mat'])
    end
    fnsData = fieldnames(dataStruct.seriesData);
    
    % Initialise variables
    srData = dataStruct.seriesData.(fnsData{1}).conf.samplingParams.srData;
    if animal == 1 || ~exist('areaFRPercentileIndividual_50', 'var')
      arearSpearman5secIndividual = {};
      areapvalSpearman5secIndividual = {};
      areaFRPupilIndividual = {};
      areaFRPercentileIndividual_50 = {};
      areaFRPercentileIndividual_25 = {};
      areaFRPercentileIndividual_12p5 = {};
      areaFRPercentileIndividual_third = {};
      areaFRRiseDecayIndividual = {};
      areaCohFullPupilIndividual = {};
      areaCohConfUFullPupilIndividual = {};
      areaCohConfLFullPupilIndividual = {};
      areaPhaseFullPupilIndividual = {};
      areaPhaseConfUFullPupilIndividual = {};
      areaPhaseConfLFullPupilIndividual = {};
      areaFreqFullPupilIndividual = {};
      areaCohFullInterpPupilIndividual = {};
      areaCohConfUFullInterpPupilIndividual = {};
      areaCohConfLFullInterpPupilIndividual = {};
      areaPhaseFullInterpPupilIndividual = {};
      areaPhaseConfUFullInterpPupilIndividual = {};
      areaPhaseConfLFullInterpPupilIndividual = {};
      areaFreqFullInterpPupilIndividual = {};
      for iCond = 1:numel(conditions)
        for iArea = 1:numel(areas)
          arearSpearman5secIndividual{iCond}{iArea} = []; %#ok<*SAGROW>
          areapvalSpearman5secIndividual{iCond}{iArea} = [];
          areaFRPupilIndividual{iCond}{iArea} = [];
          areaFRPercentileIndividual_50{iCond}{iArea} = [];
          areaFRPercentileIndividual_25{iCond}{iArea} = [];
          areaFRPercentileIndividual_12p5{iCond}{iArea} = [];
          areaFRPercentileIndividual_third{iCond}{iArea} = [];
          areaFRRiseDecayIndividual{iCond}{iArea} = [];
          areaCohFullPupilIndividual{iCond}{iArea} = {};
          areaCohConfUFullPupilIndividual{iCond}{iArea} = {};
          areaCohConfLFullPupilIndividual{iCond}{iArea} = {};
          areaPhaseFullPupilIndividual{iCond}{iArea} = {};
          areaPhaseConfUFullPupilIndividual{iCond}{iArea} = {};
          areaPhaseConfLFullPupilIndividual{iCond}{iArea} = {};
          areaFreqFullPupilIndividual{iCond}{iArea} = {};
          areaCohFullInterpPupilIndividual{iCond}{iArea} = [];
          areaCohConfUFullInterpPupilIndividual{iCond}{iArea} = [];
          areaCohConfLFullInterpPupilIndividual{iCond}{iArea} = [];
          areaPhaseFullInterpPupilIndividual{iCond}{iArea} = [];
          areaPhaseConfUFullInterpPupilIndividual{iCond}{iArea} = [];
          areaPhaseConfLFullInterpPupilIndividual{iCond}{iArea} = [];
          areaFreqFullInterpPupilIndividual{iCond}{iArea} = [];
          areaFRIndividual{iCond}{iArea} = [];
          areaFRFiltLP0p001HzIndividual_50{iCond}{iArea} = [];
          areaFRFiltLP0p001HzIndividual_25{iCond}{iArea} = [];
          areaFRFiltLP0p001HzIndividual_12p5{iCond}{iArea} = [];
          areaFRFiltLP0p001HzIndividual_third{iCond}{iArea} = [];
          areaFRFiltLP0p01HzIndividual_50{iCond}{iArea} = [];
          areaFRFiltLP0p1HzIndividual_50{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_pi{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver2{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_2piOver3{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver3{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver4{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver6{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver8{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver12{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver16{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piShifted{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_2piOver3Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver2Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver3Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver4Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver6Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver8Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver12Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver16Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_pi{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_2piOver3{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver2{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver3{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver4{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver6{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver8{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver12{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver16{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piShifted{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver2Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_2piOver3Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver3Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver4Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver6Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver8Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver12Shifted{iCond}{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver16Shifted{iCond}{iArea} = [];
        end
      end
    end
    
    for dbCount = 1:numel(fnsData) % Loop through db entries
      dbStruct = dataStruct.seriesData.(fnsData{dbCount});
      seriesName = seriesFromEntry(fnsData{dbCount});
      disp(['             series ' seriesName]);
      
      % Determine if series pupil data exist
      if isempty(dbStruct.popData)
        continue
      end
      if ~isfield(dbStruct.popData, 'pupil') || (isfield(dbStruct.popData, 'pupil') && isempty(dbStruct.popData.pupil))
        continue
      end
      if isempty(dbStruct.popData.pupil.popData) || isempty(dbStruct.popData.pupil.popData.phase)
        continue
      end
      
      % Test for exceptions
      if exceptionTest(except, seriesName)
        continue
      end
      
      % Determine if population rate > 0
      if firingRateTest(sum(dbStruct.popData.MUAsAll,1), dbStruct.conf.samplingParams.srData)
        continue
      end
      
      % Determine recording area
      if strcmp(repository,'all')
        error('Only allensdk and uol repositories are supported currently.');
      elseif strcmp(repository,'uol')
        [~, ~, ~, ~, ~, area] = determineArea(seriesName);
      elseif strcmp(repository,'allensdk')
        area = determineArea(seriesName);
      end
      
      % Determine recording condition (i.e., awake or anaesthesia)
      [breakClause, iCond] = series2condition(awake, anaesthesia, seriesName);
      if breakClause
        continue
      end
      
      % Disqualify low quality units if needed
      units = [];
      if (strcmpi(subpop, 'positive') || strcmpi(subpop, 'negative')) &&...
          ~isfield(dbStruct.shankData.(['shank' num2str(1)]), 'pupil')
        continue
      else
        phase = [];
        rSpearman = [];
      end
      for sh = 1:numel(fieldnames(dbStruct.shankData))
        units = [units; dbStruct.shankData.(['shank' num2str(sh)]).units]; %#ok<*AGROW>
        if strcmpi(subpop, 'positive') || strcmpi(subpop, 'negative')
          if pupilCorrCond == 1
            if ~isfield(dbStruct.shankData.(['shank' num2str(sh)]), 'pupil')
              continue
            else
              phase = [phase; spkPhase(dbStruct.shankData.(['shank' num2str(sh)]).pupil.unitData, fRef)'];
            end
          elseif pupilCorrCond == 2
            if ~isfield(dbStruct.shankData.(['shank' num2str(sh)]), 'rSpearman')
              continue
            else
              rSpearman = [rSpearman; dbStruct.shankData.(['shank' num2str(sh)]).rSpearman'];
            end
          end
        end
      end
      if strcmp(subpop, 'all')
        correlatedInd = 1:numel(units);
      elseif strcmp(subpop, 'positive')
        if pupilCorrCond == 1
          correlatedInd = false(size(phase));
          correlatedInd(recentrePhase(phase, 0) > -pi/2 & recentrePhase(phase, 0) <= pi/2) = true;
        elseif pupilCorrCond == 2
          correlatedInd = find(rSpearman >= 0);
        end
      elseif strcmp(subpop, 'negative')
        if pupilCorrCond == 1
          correlatedInd = false(size(phase));
          correlatedInd(recentrePhase(phase, pi) > pi/2 & recentrePhase(phase, pi) <= 3*pi/2) = true;
        elseif pupilCorrCond == 2
          correlatedInd = find(rSpearman < 0);
        end
      end
      if qualityCheck
        unitMetadata = [];
        for sh = 1:numel(fieldnames(dbStruct.shankData))
          unitMetadata = [unitMetadata; dbStruct.shankData.(['shank' num2str(sh)]).unitMetadata];
        end
        if isempty(unitMetadata)
          continue
        end
        [~, qualityUnitInd] = qualityTest2(unitMetadata, cluDist, refractCont, false);
      else
        qualityUnitInd = 1:numel(units);
      end
      qualityUnitInd = intersect(qualityUnitInd, find(correlatedInd));
      if isempty(qualityUnitInd)
        continue
      end
      
      if numel(conditions) > 1
        condLoop = [iCond numel(conditions)];
      else
        condLoop = iCond;
      end
      for iCondPlusAll = condLoop % Loop through the main and pooled conditions
        for iAreaPlusAll = area % Loop through the main and pooled areas
          close all
          
          % Get unit spiking pupil correlation data
          shankIDs = fieldnames(dbStruct.shankData);
          spk = [];
          rSpearman = [];
          pvalSpearman = [];
          for sh = 1:numel(shankIDs)
            if isempty(spk)
              spk = full(dbStruct.shankData.(shankIDs{sh}).spk);
              rSpearman = dbStruct.shankData.(shankIDs{sh}).rSpearman';
              pvalSpearman = dbStruct.shankData.(shankIDs{sh}).pvalSpearman';
            else
              spk = concatenateMat(spk, full(dbStruct.shankData.(shankIDs{sh}).spk));
              rSpearman = concatenateMat(rSpearman, dbStruct.shankData.(shankIDs{sh}).rSpearman');
              pvalSpearman = concatenateMat(pvalSpearman, dbStruct.shankData.(shankIDs{sh}).pvalSpearman');
            end
          end
          arearSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll} = [arearSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll}; rSpearman(qualityUnitInd)];
          areapvalSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll} = [areapvalSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll}; pvalSpearman(qualityUnitInd)];
          period = dbStruct.db(dbCount).period;
          
          % Get eye data
          eyeData = dataStruct.eyeData.([animals{animal} '_s' seriesName(1:min([14 numel(seriesName)]))]);
          
          % Combine periods
          commonPeriod = combinePeriods(period, eyeData.period, srData);
          
          % Interpolate and filter pupil area data
          if isempty(commonPeriod)
            continue
          end
          [pupilArea, interpTimes, interpInds] = pupilFilt(eyeData, srData, sum(spk,1), 0, commonPeriod, srData);
          pupilAreaDiff = [diff(pupilArea) 0];
          pupilAreaFiltTimes = eyeData.pupilAreaFilt.timesFiltStart:eyeData.pupilAreaFilt.timesFiltStep:eyeData.pupilAreaFilt.timesFiltStop;
          pupilAreaFiltLP0p001Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.pupilAreaFiltLP0p001Hz, interpTimes);
          pupilAreaFiltLP0p01Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.pupilAreaFiltLP0p01Hz, interpTimes);
          pupilAreaFiltLP0p1Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.pupilAreaFiltLP0p1Hz, interpTimes);
          pupilAreaFiltBP0p01to0p05Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.pupilAreaFiltBP0p01to0p05Hz, interpTimes);
          pupilAreaFiltBP0p1to0p5Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.pupilAreaFiltBP0p1to0p5Hz, interpTimes);
          phaseFiltBP0p01to0p05Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.phaseFiltBP0p01to0p05Hz, interpTimes);
          phaseFiltBP0p1to0p5Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.phaseFiltBP0p1to0p5Hz, interpTimes);
          
          % Truncate spiking data
          spk = spk(:, interpInds);
          
          % Calculate firing rates
          % Full rates
          areaFRIndividual{iCondPlusAll}{iAreaPlusAll} = [areaFRIndividual{iCondPlusAll}{iAreaPlusAll}; (sum(spk(qualityUnitInd,:),2).*srData)./size(spk(qualityUnitInd,:),2)];
          
          % Percentiles
          percentiles = [12.5 25 100/3 37.5 50 62.5 200/3 75 87.5];
          percentileValues = prctile(pupilArea,percentiles);
          
          areaFRPercentileIndividual_50DB = zeros(numel(qualityUnitInd),2);
          areaFRPercentileIndividual_50DB(:,1) = (sum(spk(qualityUnitInd,pupilArea<=percentileValues(5)),2).*srData)./size(spk(qualityUnitInd,pupilArea<=percentileValues(5)),2);
          areaFRPercentileIndividual_50DB(:,2) = (sum(spk(qualityUnitInd,pupilArea>percentileValues(5)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValues(5)),2);
          areaFRPercentileIndividual_50{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_50{iCondPlusAll}{iAreaPlusAll}; areaFRPercentileIndividual_50DB];
          
          areaFRPercentileIndividual_thirdDB = zeros(numel(qualityUnitInd),3);
          areaFRPercentileIndividual_thirdDB(:,1) = (sum(spk(qualityUnitInd,pupilArea<=percentileValues(3)),2).*srData)./size(spk(qualityUnitInd,pupilArea<=percentileValues(3)),2);
          areaFRPercentileIndividual_thirdDB(:,2) = (sum(spk(qualityUnitInd,pupilArea>percentileValues(3) & pupilArea<=percentileValues(7)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValues(3) & pupilArea<=percentileValues(7)),2);
          areaFRPercentileIndividual_thirdDB(:,3) = (sum(spk(qualityUnitInd,pupilArea>percentileValues(7)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValues(7)),2);
          areaFRPercentileIndividual_third{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_third{iCondPlusAll}{iAreaPlusAll}; areaFRPercentileIndividual_thirdDB];
          
          areaFRPercentileIndividual_25DB = zeros(numel(qualityUnitInd),4);
          areaFRPercentileIndividual_25DB(:,1) = (sum(spk(qualityUnitInd,pupilArea<=percentileValues(2)),2).*srData)./size(spk(qualityUnitInd,pupilArea<=percentileValues(2)),2);
          areaFRPercentileIndividual_25DB(:,2) = (sum(spk(qualityUnitInd,pupilArea>percentileValues(2) & pupilArea<=percentileValues(5)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValues(2) & pupilArea<=percentileValues(5)),2);
          areaFRPercentileIndividual_25DB(:,3) = (sum(spk(qualityUnitInd,pupilArea>percentileValues(5) & pupilArea<=percentileValues(8)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValues(5) & pupilArea<=percentileValues(8)),2);
          areaFRPercentileIndividual_25DB(:,4) = (sum(spk(qualityUnitInd,pupilArea>percentileValues(8)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValues(8)),2);
          areaFRPercentileIndividual_25{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_25{iCondPlusAll}{iAreaPlusAll}; areaFRPercentileIndividual_25DB];
          
          areaFRPercentileIndividual_12p5DB = zeros(numel(qualityUnitInd),8);
          areaFRPercentileIndividual_12p5DB(:,1) = (sum(spk(qualityUnitInd,pupilArea<=percentileValues(1)),2).*srData)./size(spk(qualityUnitInd,pupilArea<=percentileValues(1)),2);
          areaFRPercentileIndividual_12p5DB(:,2) = (sum(spk(qualityUnitInd,pupilArea>percentileValues(1) & pupilArea<=percentileValues(2)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValues(1) & pupilArea<=percentileValues(2)),2);
          areaFRPercentileIndividual_12p5DB(:,3) = (sum(spk(qualityUnitInd,pupilArea>percentileValues(2) & pupilArea<=percentileValues(4)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValues(2) & pupilArea<=percentileValues(4)),2);
          areaFRPercentileIndividual_12p5DB(:,4) = (sum(spk(qualityUnitInd,pupilArea>percentileValues(4) & pupilArea<=percentileValues(5)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValues(4) & pupilArea<=percentileValues(5)),2);
          areaFRPercentileIndividual_12p5DB(:,5) = (sum(spk(qualityUnitInd,pupilArea>percentileValues(5) & pupilArea<=percentileValues(6)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValues(5) & pupilArea<=percentileValues(6)),2);
          areaFRPercentileIndividual_12p5DB(:,6) = (sum(spk(qualityUnitInd,pupilArea>percentileValues(6) & pupilArea<=percentileValues(8)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValues(6) & pupilArea<=percentileValues(8)),2);
          areaFRPercentileIndividual_12p5DB(:,7) = (sum(spk(qualityUnitInd,pupilArea>percentileValues(8) & pupilArea<=percentileValues(9)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValues(8) & pupilArea<=percentileValues(9)),2);
          areaFRPercentileIndividual_12p5DB(:,8) = (sum(spk(qualityUnitInd,pupilArea>percentileValues(9)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValues(9)),2);
          areaFRPercentileIndividual_12p5{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_12p5{iCondPlusAll}{iAreaPlusAll}; areaFRPercentileIndividual_12p5DB];
          
          % Rise and decay
          areaFRRiseDecayIndividualDB = zeros(numel(qualityUnitInd),2);
          areaFRRiseDecayIndividualDB(:,1) = (sum(spk(qualityUnitInd,pupilAreaDiff>0),2).*srData)./size(spk(qualityUnitInd,pupilAreaDiff>0),2);
          areaFRRiseDecayIndividualDB(:,2) = (sum(spk(qualityUnitInd,pupilAreaDiff<0),2).*srData)./size(spk(qualityUnitInd,pupilAreaDiff<0),2);
          areaFRRiseDecayIndividual{iCondPlusAll}{iAreaPlusAll} = [areaFRRiseDecayIndividual{iCondPlusAll}{iAreaPlusAll}; areaFRRiseDecayIndividualDB];
          
          % Filtered percentiles
          percentileValues = percentileValues - percentileValues(5);
          percentileValuesFilt = zeros(numel(percentileValues),numel(pupilAreaFiltLP0p001Hz));
          for iCent = 1:numel(percentileValues)
            percentileValuesFilt(iCent,:) = pupilAreaFiltLP0p001Hz + percentileValues(iCent);
          end
          
          areaFRFiltLP0p001HzIndividual_50DB = zeros(numel(qualityUnitInd),2);
          areaFRFiltLP0p001HzIndividual_50DB(:,1) = (sum(spk(qualityUnitInd,pupilArea<=percentileValuesFilt(5,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea<=percentileValuesFilt(5,:)),2);
          areaFRFiltLP0p001HzIndividual_50DB(:,2) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt(5,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt(5,:)),2);
          areaFRFiltLP0p001HzIndividual_50{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltLP0p001HzIndividual_50{iCondPlusAll}{iAreaPlusAll}; areaFRFiltLP0p001HzIndividual_50DB];
          
          areaFRFiltLP0p001HzIndividual_thirdDB = zeros(numel(qualityUnitInd),3);
          areaFRFiltLP0p001HzIndividual_thirdDB(:,1) = (sum(spk(qualityUnitInd,pupilArea<=percentileValuesFilt(3,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea<=percentileValuesFilt(3,:)),2);
          areaFRFiltLP0p001HzIndividual_thirdDB(:,2) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt(3,:) & pupilArea<=percentileValuesFilt(7,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt(3,:) & pupilArea<=percentileValuesFilt(7,:)),2);
          areaFRFiltLP0p001HzIndividual_thirdDB(:,3) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt(7,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt(7,:)),2);
          areaFRFiltLP0p001HzIndividual_third{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltLP0p001HzIndividual_third{iCondPlusAll}{iAreaPlusAll}; areaFRFiltLP0p001HzIndividual_thirdDB];
          
          areaFRFiltLP0p001HzIndividual_25DB = zeros(numel(qualityUnitInd),4);
          areaFRFiltLP0p001HzIndividual_25DB(:,1) = (sum(spk(qualityUnitInd,pupilArea<=percentileValuesFilt(2,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea<=percentileValuesFilt(2,:)),2);
          areaFRFiltLP0p001HzIndividual_25DB(:,2) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt(2,:) & pupilArea<=percentileValuesFilt(5,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt(2,:) & pupilArea<=percentileValuesFilt(5,:)),2);
          areaFRFiltLP0p001HzIndividual_25DB(:,3) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt(5,:) & pupilArea<=percentileValuesFilt(8,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt(5,:) & pupilArea<=percentileValuesFilt(8,:)),2);
          areaFRFiltLP0p001HzIndividual_25DB(:,4) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt(8,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt(8,:)),2);
          areaFRFiltLP0p001HzIndividual_25{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltLP0p001HzIndividual_25{iCondPlusAll}{iAreaPlusAll}; areaFRFiltLP0p001HzIndividual_25DB];
          
          areaFRFiltLP0p001HzIndividual_12p5DB = zeros(numel(qualityUnitInd),8);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,1) = (sum(spk(qualityUnitInd,pupilArea<=percentileValuesFilt(1,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea<=percentileValuesFilt(1,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,2) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt(1,:) & pupilArea<=percentileValuesFilt(2,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt(1,:) & pupilArea<=percentileValuesFilt(2,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,3) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt(2,:) & pupilArea<=percentileValuesFilt(4,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt(2,:) & pupilArea<=percentileValuesFilt(4,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,4) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt(4,:) & pupilArea<=percentileValuesFilt(5,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt(4,:) & pupilArea<=percentileValuesFilt(5,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,5) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt(5,:) & pupilArea<=percentileValuesFilt(6,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt(5,:) & pupilArea<=percentileValuesFilt(6,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,6) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt(6,:) & pupilArea<=percentileValuesFilt(8,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt(6,:) & pupilArea<=percentileValuesFilt(8,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,7) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt(8,:) & pupilArea<=percentileValuesFilt(9,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt(8,:) & pupilArea<=percentileValuesFilt(9,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,8) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt(9,:)),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt(9,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltLP0p001HzIndividual_12p5{iCondPlusAll}{iAreaPlusAll}; areaFRFiltLP0p001HzIndividual_12p5DB];
          
          percentileValuesFilt = pupilAreaFiltLP0p01Hz;
          areaFRFiltLP0p01HzIndividual_50DB = zeros(numel(qualityUnitInd),2);
          areaFRFiltLP0p01HzIndividual_50DB(:,1) = (sum(spk(qualityUnitInd,pupilArea<=percentileValuesFilt),2).*srData)./size(spk(qualityUnitInd,pupilArea<=percentileValuesFilt),2);
          areaFRFiltLP0p01HzIndividual_50DB(:,2) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt),2);
          areaFRFiltLP0p01HzIndividual_50{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltLP0p01HzIndividual_50{iCondPlusAll}{iAreaPlusAll}; areaFRFiltLP0p01HzIndividual_50DB];
          
          percentileValuesFilt = pupilAreaFiltLP0p1Hz;
          areaFRFiltLP0p1HzIndividual_50DB = zeros(numel(qualityUnitInd),2);
          areaFRFiltLP0p1HzIndividual_50DB(:,1) = (sum(spk(qualityUnitInd,pupilArea<=percentileValuesFilt),2).*srData)./size(spk(qualityUnitInd,pupilArea<=percentileValuesFilt),2);
          areaFRFiltLP0p1HzIndividual_50DB(:,2) = (sum(spk(qualityUnitInd,pupilArea>percentileValuesFilt),2).*srData)./size(spk(qualityUnitInd,pupilArea>percentileValuesFilt),2);
          areaFRFiltLP0p1HzIndividual_50{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltLP0p1HzIndividual_50{iCondPlusAll}{iAreaPlusAll}; areaFRFiltLP0p1HzIndividual_50DB];
          
          % Phases
          areaFRFiltBP0p01to0p05HzIndividual_pi{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_pi{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi, 0, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_2piOver3{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_2piOver3{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, 2*pi/3, 0, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver2{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver2{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/2, 0, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver3{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver3{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/3, 0, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver4{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver4{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/4, 0, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver6{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver6{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/6, 0, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver8{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver8{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/8, 0, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver12{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver12{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/12, 0, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver16{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver16{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/16, 0, qualityUnitInd)];
          
          areaFRFiltBP0p1to0p5HzIndividual_pi{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_pi{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi, 0, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_2piOver3{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_2piOver3{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, 2*pi/3, 0, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver2{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver2{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/2, 0, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver3{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver3{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/3, 0, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver4{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver4{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/4, 0, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver6{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver6{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/6, 0, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver8{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver8{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/8, 0, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver12{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver12{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/12, 0, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver16{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver16{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/16, 0, qualityUnitInd)];
          
          % Shifted phases
          areaFRFiltBP0p01to0p05HzIndividual_piShifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piShifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi, pi/2, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_2piOver3Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_2piOver3Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, 2*pi/3, pi/3, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver2Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver2Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/2, pi/4, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver3Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver3Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/3, pi/6, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver4Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver4Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/4, pi/8, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver6Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver6Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/6, pi/12, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver8Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver8Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/8, pi/16, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver12Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver12Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/12, pi/24, qualityUnitInd)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver16Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver16Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/16, pi/32, qualityUnitInd)];
          
          areaFRFiltBP0p1to0p5HzIndividual_piShifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piShifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi, pi/2, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_2piOver3Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_2piOver3Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, 2*pi/3, pi/3, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver2Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver2Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/2, pi/4, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver3Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver3Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/3, pi/6, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver4Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver4Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/4, pi/8, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver6Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver6Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/6, pi/12, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver8Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver8Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/8, pi/16, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver12Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver12Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/12, pi/24, qualityUnitInd)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver16Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver16Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/16, pi/32, qualityUnitInd)];
        end
      end
    end
  end
end

% Determine the file name and either save or load the data
if qualityCheck
  filename = [rootFolder filesep 'globalUnits_quality.mat'];
else
  filename = [rootFolder filesep 'globalUnits.mat'];
end
if fullRun
  if ~exist(rootFolder, 'file')
    mkdir(rootFolder);
  end
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  if ~exist(filename, 'file')
    save(filename, 'conditions','areas','arearSpearman5secIndividual','areapvalSpearman5secIndividual','areaFRIndividual',...
      'areaFRPercentileIndividual_50','areaFRPercentileIndividual_25','areaFRPercentileIndividual_12p5','areaFRPercentileIndividual_third','areaFRRiseDecayIndividual',...
      'areaFRFiltLP0p001HzIndividual_50','areaFRFiltLP0p001HzIndividual_25','areaFRFiltLP0p001HzIndividual_12p5','areaFRFiltLP0p001HzIndividual_third',...
      'areaFRFiltLP0p01HzIndividual_50','areaFRFiltLP0p1HzIndividual_50','areaFRFiltBP0p01to0p05HzIndividual_pi','areaFRFiltBP0p01to0p05HzIndividual_piOver2','areaFRFiltBP0p01to0p05HzIndividual_2piOver3',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver3','areaFRFiltBP0p01to0p05HzIndividual_piOver4','areaFRFiltBP0p01to0p05HzIndividual_piOver6',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver8','areaFRFiltBP0p01to0p05HzIndividual_piOver12','areaFRFiltBP0p01to0p05HzIndividual_piOver16',...
      'areaFRFiltBP0p01to0p05HzIndividual_piShifted','areaFRFiltBP0p01to0p05HzIndividual_2piOver3Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver2Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver3Shifted',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver4Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver6Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver8Shifted',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver12Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver16Shifted','areaFRFiltBP0p1to0p5HzIndividual_pi','areaFRFiltBP0p1to0p5HzIndividual_2piOver3',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver2','areaFRFiltBP0p1to0p5HzIndividual_piOver3','areaFRFiltBP0p1to0p5HzIndividual_piOver4',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver6','areaFRFiltBP0p1to0p5HzIndividual_piOver8','areaFRFiltBP0p1to0p5HzIndividual_piOver12',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver16','areaFRFiltBP0p1to0p5HzIndividual_piShifted','areaFRFiltBP0p1to0p5HzIndividual_piOver2Shifted','areaFRFiltBP0p1to0p5HzIndividual_2piOver3Shifted',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver3Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver4Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver6Shifted',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver8Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver12Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver16Shifted');
  else
    save(filename, 'conditions','areas','arearSpearman5secIndividual','areapvalSpearman5secIndividual','areaFRIndividual',...
      'areaFRPercentileIndividual_50','areaFRPercentileIndividual_25','areaFRPercentileIndividual_12p5','areaFRPercentileIndividual_third','areaFRRiseDecayIndividual',...
      'areaFRFiltLP0p001HzIndividual_50','areaFRFiltLP0p001HzIndividual_25','areaFRFiltLP0p001HzIndividual_12p5','areaFRFiltLP0p001HzIndividual_third',...
      'areaFRFiltLP0p01HzIndividual_50','areaFRFiltLP0p1HzIndividual_50','areaFRFiltBP0p01to0p05HzIndividual_pi','areaFRFiltBP0p01to0p05HzIndividual_piOver2','areaFRFiltBP0p01to0p05HzIndividual_2piOver3',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver3','areaFRFiltBP0p01to0p05HzIndividual_piOver4','areaFRFiltBP0p01to0p05HzIndividual_piOver6',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver8','areaFRFiltBP0p01to0p05HzIndividual_piOver12','areaFRFiltBP0p01to0p05HzIndividual_piOver16',...
      'areaFRFiltBP0p01to0p05HzIndividual_piShifted','areaFRFiltBP0p01to0p05HzIndividual_2piOver3Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver2Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver3Shifted',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver4Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver6Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver8Shifted',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver12Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver16Shifted','areaFRFiltBP0p1to0p5HzIndividual_pi','areaFRFiltBP0p1to0p5HzIndividual_2piOver3',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver2','areaFRFiltBP0p1to0p5HzIndividual_piOver3','areaFRFiltBP0p1to0p5HzIndividual_piOver4',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver6','areaFRFiltBP0p1to0p5HzIndividual_piOver8','areaFRFiltBP0p1to0p5HzIndividual_piOver12',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver16','areaFRFiltBP0p1to0p5HzIndividual_piShifted','areaFRFiltBP0p1to0p5HzIndividual_piOver2Shifted','areaFRFiltBP0p1to0p5HzIndividual_2piOver3Shifted',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver3Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver4Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver6Shifted',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver8Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver12Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver16Shifted', '-append');
  end
else
  load(filename);
end



%% Local functions
function unitFiringRates = firingRateByPhase(spk, srSpk, phase, phaseStep, phaseShift, qualityUnitInd)

phaseValues = (pi:-phaseStep:-pi) + phaseShift;
phase = recentrePhase(phase, phaseShift);
unitFiringRates = zeros(numel(qualityUnitInd),numel(phaseValues)-1);
for iPhase = 1:numel(phaseValues)-1
  unitFiringRates(:,iPhase) = (sum(spk(qualityUnitInd,phase<phaseValues(iPhase) & phase>=phaseValues(iPhase+1)),2).*srSpk)./size(spk(qualityUnitInd,phase<phaseValues(iPhase) & phase>=phaseValues(iPhase+1)),2);
end
end