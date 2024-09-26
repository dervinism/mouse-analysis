clearvars -except repository subpop reverse qualityCheck allData fullRun includeRuns


%% INITIALISE PARAMETERS
params
lists

if ~exist('repository', 'var')
  repository = 'all';
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

dataDir = [dataDir filesep includeRuns];
if strcmp(repository,'all')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep laDir];
    rootFolderEye = [dataDir filesep area2pupilDir];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep laDir_positive];
    rootFolderEye = [dataDir filesep area2pupilDir_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep laDir_negative];
    rootFolderEye = [dataDir filesep area2pupilDir_negative];
  end
  animals = animalsOI;
elseif strcmp(repository,'uol')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep laDir_uol];
    rootFolderEye = [dataDir filesep area2pupilDir_uol];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep laDir_uol_positive];
    rootFolderEye = [dataDir filesep area2pupilDir_uol_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep laDir_uol_negative];
    rootFolderEye = [dataDir filesep area2pupilDir_uol_negative];
  end
  animals = animalsUOLOI;
elseif strcmp(repository,'allensdk')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep laDir_allensdk];
    rootFolderEye = [dataDir filesep area2pupilDir_allensdk];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep laDir_allensdk_positive];
    rootFolderEye = [dataDir filesep area2pupilDir_allensdk_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep laDir_allensdk_negative];
    rootFolderEye = [dataDir filesep area2pupilDir_allensdk_negative];
  end
  animals = animalsAllensdk;
  conditions = {'awake'};
end

significantUnits = false;
phaseProfiles = [true true true true];
drawCorrelationsLog = [true true true true];
drawFractions = true;

if qualityCheck
  rootFolder = [rootFolder filesep qualityUnitsFolder]; %#ok<*UNRCH>
  rootFolderEye = [rootFolderEye filesep qualityUnitsFolder];
else
  rootFolder = [rootFolder filesep unitsFolder];
  rootFolderEye = [rootFolderEye filesep unitsFolder];
end
mainFolder = [rootFolder filesep firingRatesSubfolder2]; %#ok<*NASGU>
if significantUnits
  rootFolder = [mainFolder filesep significantUnitsFolder];
  mainFolder = rootFolder;
end


%% COMPUTE VARIABLES AND INITIALISE FIGURES NEEDED FOR DISPLAYING FIRING RATE PLOTS
if fullRun
  for animal = 1:numel(animals) % Loop through animals
    disp(['processing animal ' animals{animal}]);
    if exist('allData', 'var')
      dataStruct = allData.(animals{animal}).dataStruct;
    else
      load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
    end
    fnsData = fieldnames(dataStruct.seriesData);
    
    % Initialise variables
    srData = dataStruct.seriesData.(fnsData{1}).conf.samplingParams.srData;
    if animal == 1 || ~exist('areaFRPercentileIndividual_50', 'var')
      arearSpearman5secIndividual = {};
      areapvalSpearman5secIndividual = {};
      areaFRIndividual = {};
      areaFRPercentileIndividual_50 = {};
      areaFRPercentileIndividual_25 = {};
      areaFRPercentileIndividual_12p5 = {};
      areaFRPercentileIndividual_third = {};
      areaPhaseEyeFullIndividual = {};
      areaFreqEyeFullIndividual = {};
      areaPhaseEyeFullInterpIndividual = {};
      areaFreqEyeFullInterpIndividual = {};
      areaPhaseLocalPRFullIndividual = {};
      areaFreqLocalPRFullIndividual = {};
      areaPhaseLocalPRFullInterpIndividual = {};
      areaFreqLocalPRFullInterpIndividual = {};
      areaEyeDataExistsIndividual = {};
      for iCond = 1:numel(conditions)
        arearSpearman5secIndividualCond = {};
        areapvalSpearman5secIndividualCond = {};
        areaFRIndividualCond = {};
        areaFRPercentileIndividual_50Cond = {};
        areaFRPercentileIndividual_25Cond = {};
        areaFRPercentileIndividual_12p5Cond = {};
        areaFRPercentileIndividual_thirdCond = {};
        areaPhaseEyeFullIndividualCond = {};
        areaFreqEyeFullIndividualCond = {};
        areaPhaseEyeFullInterpIndividualCond = {};
        areaFreqEyeFullInterpIndividualCond = {};
        areaPhaseLocalPRFullIndividualCond = {};
        areaFreqLocalPRFullIndividualCond = {};
        areaPhaseLocalPRFullInterpIndividualCond = {};
        areaFreqLocalPRFullInterpIndividualCond = {};
        areaEyeDataExistsIndividualCond = {};
        for iArea = 1:numel(areas)
          arearSpearman5secIndividualCond{iArea} = [];
          areapvalSpearman5secIndividualCond{iArea} = [];
          areaFRIndividualCond{iArea} = []; %#ok<*SAGROW>
          areaFRPercentileIndividual_50Cond{iArea} = [];
          areaFRPercentileIndividual_25Cond{iArea} = [];
          areaFRPercentileIndividual_12p5Cond{iArea} = [];
          areaFRPercentileIndividual_thirdCond{iArea} = [];
          areaPhaseEyeFullIndividualCond{iArea} = {};
          areaFreqEyeFullIndividualCond{iArea} = {};
          areaPhaseEyeFullInterpIndividualCond{iArea} = [];
          areaFreqEyeFullInterpIndividualCond{iArea} = [];
          areaPhaseLocalPRFullIndividualCond{iArea} = {};
          areaFreqLocalPRFullIndividualCond{iArea} = {};
          areaPhaseLocalPRFullInterpIndividualCond{iArea} = [];
          areaFreqLocalPRFullInterpIndividualCond{iArea} = [];
          areaEyeDataExistsIndividualCond{iArea} = [];
        end
        arearSpearman5secIndividual{iCond} = arearSpearman5secIndividualCond;
        areapvalSpearman5secIndividual{iCond} = areapvalSpearman5secIndividualCond;
        areaFRIndividual{iCond} = areaFRIndividualCond;
        areaFRPercentileIndividual_50{iCond} = areaFRPercentileIndividual_50Cond;
        areaFRPercentileIndividual_25{iCond} = areaFRPercentileIndividual_25Cond;
        areaFRPercentileIndividual_12p5{iCond} = areaFRPercentileIndividual_12p5Cond;
        areaFRPercentileIndividual_third{iCond} = areaFRPercentileIndividual_thirdCond;
        areaPhaseEyeFullIndividual{iCond} = areaPhaseEyeFullIndividualCond;
        areaFreqEyeFullIndividual{iCond} = areaFreqEyeFullIndividualCond;
        areaPhaseEyeFullInterpIndividual{iCond} = areaPhaseEyeFullInterpIndividualCond;
        areaFreqEyeFullInterpIndividual{iCond} = areaFreqEyeFullInterpIndividualCond;
        areaPhaseLocalPRFullIndividual{iCond} = areaPhaseLocalPRFullIndividualCond;
        areaFreqLocalPRFullIndividual{iCond} = areaFreqLocalPRFullIndividualCond;
        areaPhaseLocalPRFullInterpIndividual{iCond} = areaPhaseLocalPRFullInterpIndividualCond;
        areaFreqLocalPRFullInterpIndividual{iCond} = areaFreqLocalPRFullInterpIndividualCond;
        areaEyeDataExistsIndividual{iCond} = areaEyeDataExistsIndividualCond;
      end
    end
    
    for dbCount = 1:numel(fnsData) % Loop through db entries
      dbStruct = dataStruct.seriesData.(fnsData{dbCount});
      if isempty(dbStruct)
        continue
      end
      seriesName = seriesFromEntry(fnsData{dbCount});
      seriesNameEyeData = [animals{animal} '_s' seriesName(1:min([14 numel(seriesName)]))];
      disp(['             series ' seriesName]);
      
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
      for sh = 1:numel(fieldnames(dbStruct.shankData))
        units = [units; dbStruct.shankData.(['shank' num2str(sh)]).units]; %#ok<*AGROW>
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
          
          % Get unit spiking data
          shankIDs = fieldnames(dbStruct.shankData);
          spk = [];
          for sh = 1:numel(shankIDs)
            if isempty(spk)
              spk = full(dbStruct.shankData.(shankIDs{sh}).spk);
            else
              spk = concatenateMat(spk, full(dbStruct.shankData.(shankIDs{sh}).spk));
            end
          end
          period = dbStruct.db(dbCount).period;
          
          if isfield(dbStruct.popData, 'pupil') && ~isempty(dbStruct.popData.pupil)
            
            % Get unit spiking pupil correlation data
            rSpearman = [];
            pvalSpearman = [];
            for sh = 1:numel(shankIDs)
              if isempty(spk)
                rSpearman = dbStruct.shankData.(shankIDs{sh}).rSpearman';
                pvalSpearman = dbStruct.shankData.(shankIDs{sh}).pvalSpearman';
              else
                rSpearman = concatenateMat(rSpearman, dbStruct.shankData.(shankIDs{sh}).rSpearman');
                pvalSpearman = concatenateMat(pvalSpearman, dbStruct.shankData.(shankIDs{sh}).pvalSpearman');
              end
            end
            arearSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll} = [arearSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll}; rSpearman(qualityUnitInd)];
            areapvalSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll} = [areapvalSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll}; pvalSpearman(qualityUnitInd)];
            period = dbStruct.db(dbCount).period;
            
            % Get eye data
            eyeData = dataStruct.eyeData.([animals{animal} '_s' seriesName(1:min([14 numel(seriesName)]))]);
            %eyeData.pupilArea = double(eyeData.pupilAreaFilt.pupilAreaFiltHP0p01Hz);
            %eyeData.frameTimes = eyeData.pupilAreaFilt.timesFiltStart:eyeData.pupilAreaFilt.timesFiltStep:eyeData.pupilAreaFilt.timesFiltStop;
            
            % Combine periods
            commonPeriod = combinePeriods(period, eyeData.period, srData);
            if isempty(commonPeriod)
              
              % Truncate spiking data
              spk = spk(qualityUnitInd, max([ceil(period(1)*srData) 1]):min([floor(period(2)*srData) size(spk,2)]));
              
              % Calculate firing rates
              % Full rates
              areaFRIndividual{iCondPlusAll}{iAreaPlusAll} = [areaFRIndividual{iCondPlusAll}{iAreaPlusAll}; (sum(spk,2).*srData)./size(spk,2)];
              
              % Phase
              for sh = 1:numel(shankIDs)
                shankStruct = dbStruct.shankData.(shankIDs{sh});
                for u = 1:numel(units)
                  iU = find(shankStruct.units == units(u));
                  if isempty(iU)
                    continue
                  end
                  areaPhaseEyeFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseEyeFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = NaN;
                  areaFreqEyeFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaFreqEyeFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = NaN;
                  if isfield(shankStruct.phaseCoh{iU}, 'phase') && ~isempty(shankStruct.phaseCoh{iU}.phase)
                    areaPhaseLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = shankStruct.phaseCoh{iU}.phase;
                  else
                    areaPhaseLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = NaN(size(shankStruct.phaseCoh{iU}.freq));
                  end
                  areaFreqLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaFreqLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = shankStruct.phaseCoh{iU}.freq;
                end
              end
              
              % Mark if eye data exists
              arearSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll} = [arearSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll}; NaN(size(spk,1),1)];
              areapvalSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll} = [areapvalSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll}; NaN(size(spk,1),1)];
              areaFRPercentileIndividual_50{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_50{iCondPlusAll}{iAreaPlusAll}; NaN(size(spk,1),2)];
              areaFRPercentileIndividual_25{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_25{iCondPlusAll}{iAreaPlusAll}; NaN(size(spk,1),4)];
              areaFRPercentileIndividual_12p5{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_12p5{iCondPlusAll}{iAreaPlusAll}; NaN(size(spk,1),8)];
              areaFRPercentileIndividual_third{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_third{iCondPlusAll}{iAreaPlusAll}; NaN(size(spk,1),3)];
              areaEyeDataExistsIndividual{iCondPlusAll}{iAreaPlusAll} = [areaEyeDataExistsIndividual{iCondPlusAll}{iAreaPlusAll}; false(size(spk,1),1)];
              continue
            end
            
            % Interpolate and filter pupil area data
            [pupilArea, interpTimes, interpInds] = pupilFilt(eyeData, srData, sum(spk,1), 0, commonPeriod, srData);
            
            % Truncate spiking data
            spk = spk(qualityUnitInd, interpInds);
            
            % Calculate firing rates
            % Full rates
            areaFRIndividual{iCondPlusAll}{iAreaPlusAll} = [areaFRIndividual{iCondPlusAll}{iAreaPlusAll}; (sum(spk,2).*srData)./size(spk,2)];
            
            % Percentiles
            percentiles = [12.5 25 100/3 37.5 50 62.5 200/3 75 87.5];
            percentileValues = prctile(pupilArea,percentiles);
            
            areaFRPercentileIndividual_50DB = zeros(numel(qualityUnitInd),2);
            areaFRPercentileIndividual_50DB(:,1) = (sum(spk(:,pupilArea<=percentileValues(5)),2).*srData)./size(spk(:,pupilArea<=percentileValues(5)),2);
            areaFRPercentileIndividual_50DB(:,2) = (sum(spk(:,pupilArea>percentileValues(5)),2).*srData)./size(spk(:,pupilArea>percentileValues(5)),2);
            areaFRPercentileIndividual_50{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_50{iCondPlusAll}{iAreaPlusAll}; areaFRPercentileIndividual_50DB];
            
            areaFRPercentileIndividual_thirdDB = zeros(numel(qualityUnitInd),3);
            areaFRPercentileIndividual_thirdDB(:,1) = (sum(spk(:,pupilArea<=percentileValues(3)),2).*srData)./size(spk(:,pupilArea<=percentileValues(3)),2);
            areaFRPercentileIndividual_thirdDB(:,2) = (sum(spk(:,pupilArea>percentileValues(3) & pupilArea<=percentileValues(7)),2).*srData)./size(spk(:,pupilArea>percentileValues(3) & pupilArea<=percentileValues(7)),2);
            areaFRPercentileIndividual_thirdDB(:,3) = (sum(spk(:,pupilArea>percentileValues(7)),2).*srData)./size(spk(:,pupilArea>percentileValues(7)),2);
            areaFRPercentileIndividual_third{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_third{iCondPlusAll}{iAreaPlusAll}; areaFRPercentileIndividual_thirdDB];
            
            areaFRPercentileIndividual_25DB = zeros(numel(qualityUnitInd),4);
            areaFRPercentileIndividual_25DB(:,1) = (sum(spk(:,pupilArea<=percentileValues(2)),2).*srData)./size(spk(:,pupilArea<=percentileValues(2)),2);
            areaFRPercentileIndividual_25DB(:,2) = (sum(spk(:,pupilArea>percentileValues(2) & pupilArea<=percentileValues(5)),2).*srData)./size(spk(:,pupilArea>percentileValues(2) & pupilArea<=percentileValues(5)),2);
            areaFRPercentileIndividual_25DB(:,3) = (sum(spk(:,pupilArea>percentileValues(5) & pupilArea<=percentileValues(8)),2).*srData)./size(spk(:,pupilArea>percentileValues(5) & pupilArea<=percentileValues(8)),2);
            areaFRPercentileIndividual_25DB(:,4) = (sum(spk(:,pupilArea>percentileValues(8)),2).*srData)./size(spk(:,pupilArea>percentileValues(8)),2);
            areaFRPercentileIndividual_25{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_25{iCondPlusAll}{iAreaPlusAll}; areaFRPercentileIndividual_25DB];
            
            areaFRPercentileIndividual_12p5DB = zeros(numel(qualityUnitInd),8);
            areaFRPercentileIndividual_12p5DB(:,1) = (sum(spk(:,pupilArea<=percentileValues(1)),2).*srData)./size(spk(:,pupilArea<=percentileValues(1)),2);
            areaFRPercentileIndividual_12p5DB(:,2) = (sum(spk(:,pupilArea>percentileValues(1) & pupilArea<=percentileValues(2)),2).*srData)./size(spk(:,pupilArea>percentileValues(1) & pupilArea<=percentileValues(2)),2);
            areaFRPercentileIndividual_12p5DB(:,3) = (sum(spk(:,pupilArea>percentileValues(2) & pupilArea<=percentileValues(4)),2).*srData)./size(spk(:,pupilArea>percentileValues(2) & pupilArea<=percentileValues(4)),2);
            areaFRPercentileIndividual_12p5DB(:,4) = (sum(spk(:,pupilArea>percentileValues(4) & pupilArea<=percentileValues(5)),2).*srData)./size(spk(:,pupilArea>percentileValues(4) & pupilArea<=percentileValues(5)),2);
            areaFRPercentileIndividual_12p5DB(:,5) = (sum(spk(:,pupilArea>percentileValues(5) & pupilArea<=percentileValues(6)),2).*srData)./size(spk(:,pupilArea>percentileValues(5) & pupilArea<=percentileValues(6)),2);
            areaFRPercentileIndividual_12p5DB(:,6) = (sum(spk(:,pupilArea>percentileValues(6) & pupilArea<=percentileValues(8)),2).*srData)./size(spk(:,pupilArea>percentileValues(6) & pupilArea<=percentileValues(8)),2);
            areaFRPercentileIndividual_12p5DB(:,7) = (sum(spk(:,pupilArea>percentileValues(8) & pupilArea<=percentileValues(9)),2).*srData)./size(spk(:,pupilArea>percentileValues(8) & pupilArea<=percentileValues(9)),2);
            areaFRPercentileIndividual_12p5DB(:,8) = (sum(spk(:,pupilArea>percentileValues(9)),2).*srData)./size(spk(:,pupilArea>percentileValues(9)),2);
            areaFRPercentileIndividual_12p5{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_12p5{iCondPlusAll}{iAreaPlusAll}; areaFRPercentileIndividual_12p5DB];
            
            % Phase
            for sh = 1:numel(shankIDs)
              shankStruct = dbStruct.shankData.(shankIDs{sh});
              for u = 1:numel(units)
                iU = find(shankStruct.units == units(u));
                if isempty(iU)
                  continue
                end
                areaPhaseEyeFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseEyeFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = shankStruct.pupil.unitData{iU}.phase;
                areaFreqEyeFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaFreqEyeFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = dbStruct.popData.pupil.popData.freq;
                if isfield(shankStruct.phaseCoh{iU}, 'phase') && ~isempty(shankStruct.phaseCoh{iU}.phase)
                  areaPhaseLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = shankStruct.phaseCoh{iU}.phase;
                else
                  areaPhaseLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = NaN(size(shankStruct.phaseCoh{iU}.freq));
                end
                areaFreqLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaFreqLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = shankStruct.phaseCoh{iU}.freq;
              end
            end
            
            % Mark if eye data exists
            areaEyeDataExistsIndividual{iCondPlusAll}{iAreaPlusAll} = [areaEyeDataExistsIndividual{iCondPlusAll}{iAreaPlusAll}; true(size(spk,1),1)];
          else
            % Truncate spiking data
            spk = spk(qualityUnitInd, max([ceil(period(1)*srData) 1]):min([floor(period(2)*srData) size(spk,2)]));
            
            % Calculate firing rates
            % Full rates
            areaFRIndividual{iCondPlusAll}{iAreaPlusAll} = [areaFRIndividual{iCondPlusAll}{iAreaPlusAll}; (sum(spk,2).*srData)./size(spk,2)];
            
            % Phase
            for sh = 1:numel(shankIDs)
              shankStruct = dbStruct.shankData.(shankIDs{sh});
              for u = 1:numel(units)
                iU = find(shankStruct.units == units(u));
                if isempty(iU)
                  continue
                end
                areaPhaseEyeFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseEyeFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = NaN;
                areaFreqEyeFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaFreqEyeFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = NaN;
                if isfield(shankStruct.phaseCoh{iU}, 'phase') && ~isempty(shankStruct.phaseCoh{iU}.phase)
                  areaPhaseLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = shankStruct.phaseCoh{iU}.phase;
                else
                  areaPhaseLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = NaN(size(shankStruct.phaseCoh{iU}.freq));
                end
                areaFreqLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaFreqLocalPRFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = shankStruct.phaseCoh{iU}.freq;
              end
            end
            
            % Mark if eye data exists
            arearSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll} = [arearSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll}; NaN(size(spk,1),1)];
            areapvalSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll} = [areapvalSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll}; NaN(size(spk,1),1)];
            areaFRPercentileIndividual_50{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_50{iCondPlusAll}{iAreaPlusAll}; NaN(size(spk,1),2)];
            areaFRPercentileIndividual_25{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_25{iCondPlusAll}{iAreaPlusAll}; NaN(size(spk,1),4)];
            areaFRPercentileIndividual_12p5{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_12p5{iCondPlusAll}{iAreaPlusAll}; NaN(size(spk,1),8)];
            areaFRPercentileIndividual_third{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_third{iCondPlusAll}{iAreaPlusAll}; NaN(size(spk,1),3)];
            areaEyeDataExistsIndividual{iCondPlusAll}{iAreaPlusAll} = [areaEyeDataExistsIndividual{iCondPlusAll}{iAreaPlusAll}; false(size(spk,1),1)];
          end
        end
      end
    end
  end
  
  % Interpolate and store full phases wrt pupil size
  freqCombined = FOI;
  for iCond = 1:numel(conditions)
    for iArea = 1:numel(areas)
      nRec = numel(areaFreqEyeFullIndividual{iCond}{iArea});
      for iRec = 1:nRec
        freqCombined = unique([freqCombined areaFreqEyeFullIndividual{iCond}{iArea}{iRec}]);
      end
    end
  end
  freqCombined = unique(freqCombined(~isnan(freqCombined)));
  for iCond = 1:numel(conditions)
    for iArea = 1:numel(areas)
      nRec = numel(areaFreqEyeFullIndividual{iCond}{iArea});
      areaPhaseEyeFullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      for iRec = 1:nRec
        if ~isempty(areaPhaseEyeFullIndividual{iCond}{iArea}) && sum(~isnan(areaPhaseEyeFullIndividual{iCond}{iArea}{iRec}))
          areaPhaseEyeFullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqEyeFullIndividual{iCond}{iArea}{iRec},...
            areaPhaseEyeFullIndividual{iCond}{iArea}{iRec}, freqCombined);
        end
      end
    end
  end
  areaFreqEyeFullInterpIndividual = freqCombined;
  
  % Interpolate and store full phases wrt local PR
  freqCombined = FOI;
  for iCond = 1:numel(conditions)
    for iArea = 1:numel(areas)
      nRec = numel(areaFreqLocalPRFullIndividual{iCond}{iArea});
      for iRec = 1:nRec
        freqCombined = unique([freqCombined areaFreqLocalPRFullIndividual{iCond}{iArea}{iRec}]);
      end
    end
  end
  freqCombined = unique(freqCombined(~isnan(freqCombined)));
  for iCond = 1:numel(conditions)
    for iArea = 1:numel(areas)
      nRec = numel(areaFreqLocalPRFullIndividual{iCond}{iArea});
      areaPhaseLocalPRFullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      for iRec = 1:nRec
        if ~isempty(areaPhaseLocalPRFullIndividual{iCond}{iArea}) && sum(~isnan(areaPhaseLocalPRFullIndividual{iCond}{iArea}{iRec}))
          areaPhaseLocalPRFullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqLocalPRFullIndividual{iCond}{iArea}{iRec},...
            areaPhaseLocalPRFullIndividual{iCond}{iArea}{iRec}, freqCombined);
        end
      end
    end
  end
  areaFreqLocalPRFullInterpIndividual = freqCombined;
end

% Determine the file name and either save or load the data
if qualityCheck
  filename = [mainFolder filesep 'globalUnits_quality.mat'];
  filename2 = [rootFolder filesep 'globalUnits_quality.mat'];
  filename3 = [rootFolderEye filesep 'globalUnits_area2pupil_quality.mat'];
else
  filename = [mainFolder filesep 'globalUnits.mat'];
  filename2 = [rootFolder filesep 'globalUnits.mat'];
  filename3 = [rootFolderEye filesep 'globalUnits_area2pupil.mat'];
end
if fullRun
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  if ~exist(filename, 'file')
    save(filename, 'conditions','areas','arearSpearman5secIndividual','areapvalSpearman5secIndividual','areaFRIndividual',...
      'areaFRPercentileIndividual_50','areaFRPercentileIndividual_25','areaFRPercentileIndividual_12p5','areaFRPercentileIndividual_third',...
      'areaPhaseEyeFullIndividual','areaFreqEyeFullIndividual','areaPhaseEyeFullInterpIndividual','areaFreqEyeFullInterpIndividual',...
      'areaPhaseLocalPRFullIndividual','areaFreqLocalPRFullIndividual','areaPhaseLocalPRFullInterpIndividual',...
      'areaFreqLocalPRFullInterpIndividual','areaEyeDataExistsIndividual');
  else
    save(filename, 'conditions','areas','arearSpearman5secIndividual','areapvalSpearman5secIndividual','areaFRIndividual',...
      'areaFRPercentileIndividual_50','areaFRPercentileIndividual_25','areaFRPercentileIndividual_12p5','areaFRPercentileIndividual_third',...
      'areaPhaseEyeFullIndividual','areaFreqEyeFullIndividual','areaPhaseEyeFullInterpIndividual','areaFreqEyeFullInterpIndividual',...
      'areaPhaseLocalPRFullIndividual','areaFreqLocalPRFullIndividual','areaPhaseLocalPRFullInterpIndividual',...
      'areaFreqLocalPRFullInterpIndividual','areaEyeDataExistsIndividual', '-append');
  end
end
fullData = load(filename);
significantData = load(filename2);
eyeData = load(filename3);


%% SIGNIFICANT LOCAL PHASE VS 5-SEC CORRELATION
if phaseProfiles(1)
  options = struct();
  options.corrType = 'circlinearnp';
  options.figFolder = mainFolder;
  options.diagonal = false;
  options.fitLine = false;
  options.fitLineDisplay = 'off';
  options.xLim = [];
  options.xTicks = [];
  options.yLim = [];
  options.yTicks = [];
  options.figSize = figSize;
  options.saveFig = true;
  significantPhase = significantData.areaPhaseFullInterpIndividual; % Data with eye signal only
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if ~isempty(significantPhase{iCond}{iAreasOI(iArea)})
        freq = repmat(significantData.areaFreqFullInterpIndividual, size(significantPhase{iCond}{iAreasOI(iArea)},1), 1);
        freq = reshape(freq, numel(freq), 1);
        phase = reshape(recentrePhase(significantPhase{iCond}{iAreasOI(iArea)}, phaseCentre), numel(significantPhase{iCond}{iAreasOI(iArea)}), 1);
        rSpearman = repmat(fullData.arearSpearman5secIndividual{iCond}{iAreasOI(iArea)}, size(significantPhase{iCond}{iAreasOI(iArea)},2), 1);
        pvalSpearman = repmat(fullData.areapvalSpearman5secIndividual{iCond}{iAreasOI(iArea)}, size(significantPhase{iCond}{iAreasOI(iArea)},2), 1);
        
        options2 = options;
        options2.colourVector = ones(size(rSpearman));
        options2.colourVector(rSpearman >= 0 & pvalSpearman < 0.05) = 2;
        options2.colourVector(rSpearman >= 0 & pvalSpearman >= 0.05) = 3;
        options2.colourVector(rSpearman < 0 & pvalSpearman < 0.05) = 4;
        options2.colourVector(rSpearman < 0 & pvalSpearman >= 0.05) = 5;
        options2.colourCodes = zeros(numel(unique(options2.colourVector)), 3);
        options2.colourCodes(1,:) = [200 200 200]./255; % Grey
        options2.colourCodes(2,:) = matlabColours(10); % Red
        options2.colourCodes(3,:) = [255 128 128]./255; % Light red
        options2.colourCodes(4,:) = matlabColours(8); % Blue
        options2.colourCodes(5,:) = [128 128 255]./255; % Light blue
        
        options2.figTitle = ['Pupil corr frequency local PR phase profile in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['pupilCorrFreqLocalPhaseSignificantProfile in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Frequency (Hz)';
        options2.yLabel = 'Phase (rad)';
        options2.axesType = 'semilogx';
        options2.xLim = [min(freq(~isnan(phase))) max(freq(~isnan(phase)))];
        options2.yLim = [-pi/2 3*pi/2];
        corrPlot(freq, phase, options2);
      end
    end
  end
end


%% FULL LOCAL PHASE VS 5-SEC CORRELATION
if phaseProfiles(2)
  options = struct();
  options.corrType = 'circlinearnp';
  options.figFolder = mainFolder;
  options.diagonal = false;
  options.fitLine = false;
  options.fitLineDisplay = 'off';
  options.xLim = [];
  options.xTicks = [];
  options.yLim = [];
  options.yTicks = [];
  options.figSize = figSize;
  options.saveFig = true;
  fullPhase = fullData.areaPhaseLocalPRFullInterpIndividual; % Data with and without the eye signal
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if ~isempty(fullPhase{iCond}{iAreasOI(iArea)})
        freq = repmat(fullData.areaFreqLocalPRFullInterpIndividual, size(fullPhase{iCond}{iAreasOI(iArea)},1), 1);
        freq = reshape(freq, numel(freq), 1);
        phase = reshape(recentrePhase(fullPhase{iCond}{iAreasOI(iArea)}, phaseCentre), numel(fullPhase{iCond}{iAreasOI(iArea)}), 1);
        rSpearman = repmat(fullData.arearSpearman5secIndividual{iCond}{iAreasOI(iArea)}, size(fullPhase{iCond}{iAreasOI(iArea)},2), 1);
        pvalSpearman = repmat(fullData.areapvalSpearman5secIndividual{iCond}{iAreasOI(iArea)}, size(fullPhase{iCond}{iAreasOI(iArea)},2), 1);
        
        options2 = options;
        options2.colourVector = ones(size(rSpearman));
        options2.colourVector(rSpearman >= 0 & pvalSpearman < 0.05) = 2;
        options2.colourVector(rSpearman >= 0 & pvalSpearman >= 0.05) = 3;
        options2.colourVector(rSpearman < 0 & pvalSpearman < 0.05) = 4;
        options2.colourVector(rSpearman < 0 & pvalSpearman >= 0.05) = 5;
        options2.colourCodes = zeros(numel(unique(options2.colourVector)), 3);
        options2.colourCodes(1,:) = [200 200 200]./255; % Grey
        options2.colourCodes(2,:) = matlabColours(10); % Red
        options2.colourCodes(3,:) = [255 128 128]./255; % Light red
        options2.colourCodes(4,:) = matlabColours(8); % Blue
        options2.colourCodes(5,:) = [128 128 255]./255; % Light blue
        
        options2.figTitle = ['Pupil corr frequency local PR phase profile in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['pupilCorrFreqLocalPhaseFullProfile in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Frequency (Hz)';
        options2.yLabel = 'Phase (rad)';
        options2.axesType = 'semilogx';
        options2.xLim = [min(freq(~isnan(phase))) max(freq(~isnan(phase)))];
        options2.yLim = [-pi/2 3*pi/2];
        corrPlot(freq, phase, options2);
      end
    end
  end
end


%% SIGNIFICANT PUPIL PHASE VS 5-SEC CORRELATION
if phaseProfiles(3)
  options = struct();
  options.corrType = 'circlinearnp';
  options.figFolder = mainFolder;
  options.diagonal = false;
  options.fitLine = false;
  options.fitLineDisplay = 'off';
  options.xLim = [];
  options.xTicks = [];
  options.yLim = [];
  options.yTicks = [];
  options.figSize = figSize;
  options.saveFig = true;
  pupilPhase = eyeData.areaPhaseFullInterpIndividual;
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if ~isempty(pupilPhase{iCond}{iAreasOI(iArea)})
        freq = repmat(eyeData.areaFreqFullInterpIndividual, size(pupilPhase{iCond}{iAreasOI(iArea)},1), 1);
        freq = reshape(freq, numel(freq), 1);
        phase = reshape(recentrePhase(pupilPhase{iCond}{iAreasOI(iArea)}, phaseCentre), numel(pupilPhase{iCond}{iAreasOI(iArea)}), 1);
        rSpearman = repmat(significantData.arearSpearman5secIndividual{iCond}{iAreasOI(iArea)}, size(significantPhase{iCond}{iAreasOI(iArea)},2), 1);
        pvalSpearman = repmat(significantData.areapvalSpearman5secIndividual{iCond}{iAreasOI(iArea)}, size(significantPhase{iCond}{iAreasOI(iArea)},2), 1);
        
        options2 = options;
        options2.colourVector = zeros(size(rSpearman));
        options2.colourVector(rSpearman >= 0 & pvalSpearman < 0.05) = 1;
        options2.colourVector(rSpearman >= 0 & pvalSpearman >= 0.05) = 2;
        options2.colourVector(rSpearman < 0 & pvalSpearman < 0.05) = 3;
        options2.colourVector(rSpearman < 0 & pvalSpearman >= 0.05) = 4;
        assert(numel(unique(options2.colourVector)) <= 4);
        options2.colourCodes = zeros(numel(unique(options2.colourVector)), 3);
        options2.colourCodes(1,:) = matlabColours(10); % Red
        options2.colourCodes(2,:) = [255 128 128]./255; % Light red
        options2.colourCodes(3,:) = matlabColours(8); % Blue
        options2.colourCodes(4,:) = [128 128 255]./255; % Light blue
        
        options2.figTitle = ['Pupil corr frequency pupil phase profile in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['pupilCorrFreqPupilPhaseSignificantProfile in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Frequency (Hz)';
        options2.yLabel = 'Phase (rad)';
        options2.axesType = 'semilogx';
        options2.xLim = [min(freq(~isnan(phase))) 2];
        options2.yLim = [-pi/2 3*pi/2];
        corrPlot(freq, phase, options2);
      end
    end
  end
end


%% FULL PUPIL PHASE VS 5-SEC CORRELATION
if phaseProfiles(4)
  options = struct();
  options.corrType = 'circlinearnp';
  options.figFolder = mainFolder;
  options.diagonal = false;
  options.fitLine = false;
  options.fitLineDisplay = 'off';
  options.xLim = [];
  options.xTicks = [];
  options.yLim = [];
  options.yTicks = [];
  options.figSize = figSize;
  options.saveFig = true;
  pupilPhase = fullData.areaPhaseEyeFullInterpIndividual;
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if ~isempty(pupilPhase{iCond}{iAreasOI(iArea)})
        freq = repmat(fullData.areaFreqEyeFullInterpIndividual, size(pupilPhase{iCond}{iAreasOI(iArea)},1), 1);
        freq = reshape(freq, numel(freq), 1);
        phase = reshape(recentrePhase(pupilPhase{iCond}{iAreasOI(iArea)}, phaseCentre), numel(pupilPhase{iCond}{iAreasOI(iArea)}), 1);
        rSpearman = repmat(fullData.arearSpearman5secIndividual{iCond}{iAreasOI(iArea)}, size(fullPhase{iCond}{iAreasOI(iArea)},2), 1);
        pvalSpearman = repmat(fullData.areapvalSpearman5secIndividual{iCond}{iAreasOI(iArea)}, size(fullPhase{iCond}{iAreasOI(iArea)},2), 1);
        
        options2 = options;
        options2.colourVector = ones(size(rSpearman));
        options2.colourVector(rSpearman >= 0 & pvalSpearman < 0.05) = 2;
        options2.colourVector(rSpearman >= 0 & pvalSpearman >= 0.05) = 3;
        options2.colourVector(rSpearman < 0 & pvalSpearman < 0.05) = 4;
        options2.colourVector(rSpearman < 0 & pvalSpearman >= 0.05) = 5;
        options2.colourCodes = zeros(numel(unique(options2.colourVector)), 3);
        options2.colourCodes(1,:) = [200 200 200]./255; % Grey
        options2.colourCodes(2,:) = matlabColours(10); % Red
        options2.colourCodes(3,:) = [255 128 128]./255; % Light red
        options2.colourCodes(4,:) = matlabColours(8); % Blue
        options2.colourCodes(5,:) = [128 128 255]./255; % Light blue
        
        options2.figTitle = ['Pupil corr frequency pupil phase profile in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['pupilCorrFreqPupilPhaseFullProfile in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Frequency (Hz)';
        options2.yLabel = 'Phase (rad)';
        options2.axesType = 'semilogx';
        options2.xLim = [min(freq(~isnan(phase))) 2];
        options2.yLim = [-pi/2 3*pi/2];
        corrPlot(freq, phase, options2);
      end
    end
  end
end


%% LOG FIRING RATE CORRELATIONS: SIGNIFICANT PHASE ONLY
if drawCorrelationsLog(1)
  options = struct();
  options.corrType = 'Spearman';
  options.figFolder = mainFolder;
  options.diagonal = true;
  options.fitLine = true;
  options.fitLineDisplay = 'on';
  options.xLim = [];
  options.xTicks = [];
  options.yLim = [];
  options.yTicks = [];
  options.figSize = figSize;
  options.saveFig = true;
  pupilPhase = eyeData.areaPhaseFullInterpIndividual;
  freq = eyeData.areaFreqFullInterpIndividual;
  reducedFOI = FOI(FOI <= 2);
  [~, indsFOI] = ismember(reducedFOI, freq);
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if ~isempty(pupilPhase{iCond}{iAreasOI(iArea)})
        phase = pupilPhase{iCond}{iAreasOI(iArea)}(:,indsFOI);
        for f = 1:numel(reducedFOI)
          fPhase = recentrePhase(phase(:,f), pi);
          indsNaN = isnan(fPhase);
          if sum(indsNaN) < numel(fPhase)
            options2 = options;
            options2.faceColourVector = hsv2rgb([fPhase./(2*pi) ones(size(fPhase)) ones(size(fPhase))]);
            options2.faceColourVector = options2.faceColourVector(~indsNaN,:);
            
            rSpearman = significantData.arearSpearman5secIndividual{iCond}{iAreasOI(iArea)}(~indsNaN);
            pvalSpearman = significantData.areapvalSpearman5secIndividual{iCond}{iAreasOI(iArea)}(~indsNaN);
            options2.edgeColourVector = zeros(size(rSpearman));
            options2.edgeColourVector(rSpearman >= 0 & pvalSpearman < 0.05) = 1;
            options2.edgeColourVector(rSpearman >= 0 & pvalSpearman >= 0.05) = 2;
            options2.edgeColourVector(rSpearman < 0 & pvalSpearman < 0.05) = 3;
            options2.edgeColourVector(rSpearman < 0 & pvalSpearman >= 0.05) = 4;
            options2.edgeColourCodes = zeros(numel(unique(options2.edgeColourVector)), 3);
            options2.edgeColourCodes(1,:) = matlabColours(10); % Red
            options2.edgeColourCodes(2,:) = [255 128 128]./255; % Light red
            options2.edgeColourCodes(3,:) = matlabColours(8); % Blue
            options2.edgeColourCodes(4,:) = [128 128 255]./255; % Light blue
            
            % 50 cent
            if ~isempty(significantData.areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)})
              options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
              options2.figName = ['frLogConstrictedVsfrLogDilated50Significant in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond} ' at frequency ' num2str(reducedFOI(f)) 'Hz'];
              options2.xLabel = 'Constricted log_{10}(firing rate)';
              options2.yLabel = 'Dilated log_{10}(firing rate)';
              options2.axesType = 'regular';
              corrPlot(getLog(significantData.areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(~indsNaN,1)), getLog(significantData.areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(~indsNaN,end)), options2);
            end
            
            % 25 cent
            if ~isempty(significantData.areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)})
              options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
              options2.figName = ['frLogConstrictedVsfrLogDilated25Significant in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond} ' at frequency ' num2str(reducedFOI(f)) 'Hz'];
              options2.xLabel = 'Constricted log_{10}(firing rate)';
              options2.yLabel = 'Dilated log_{10}(firing rate)';
              options2.axesType = 'regular';
              corrPlot(getLog(significantData.areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(~indsNaN,1)), getLog(significantData.areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(~indsNaN,end)), options2);
            end
            
            % 12.5 cent
            if ~isempty(significantData.areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)})
              options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
              options2.figName = ['frLogConstrictedVsfrLogDilated12p5Significant in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond} ' at frequency ' num2str(reducedFOI(f)) 'Hz'];
              options2.xLabel = 'Constricted log_{10}(firing rate)';
              options2.yLabel = 'Dilated log_{10}(firing rate)';
              options2.axesType = 'regular';
              corrPlot(getLog(significantData.areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(~indsNaN,1)), getLog(significantData.areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(~indsNaN,end)), options2);
            end
            
            % 33.33 cent
            if ~isempty(significantData.areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)})
              options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
              options2.figName = ['frLogConstrictedVsfrLogDilatedThirdSignificant in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond} ' at frequency ' num2str(reducedFOI(f)) 'Hz'];
              options2.xLabel = 'Constricted log_{10}(firing rate)';
              options2.yLabel = 'Dilated log_{10}(firing rate)';
              options2.axesType = 'regular';
              corrPlot(getLog(significantData.areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(~indsNaN,1)), getLog(significantData.areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(~indsNaN,end)), options2);
            end
          end
        end
      end
    end
  end
end


%% LOG FIRING RATE CORRELATIONS: SIGNIFICANT PHASE EXTRA
if drawCorrelationsLog(2)
  options = struct();
  options.corrType = 'Spearman';
  options.figFolder = mainFolder;
  options.diagonal = true;
  options.fitLine = true;
  options.fitLineDisplay = 'on';
  options.xLim = [];
  options.xTicks = [];
  options.yLim = [];
  options.yTicks = [];
  options.figSize = figSize;
  options.saveFig = true;
  pupilPhase = eyeData.areaPhaseFullInterpIndividual;
  freq = eyeData.areaFreqFullInterpIndividual;
  reducedFOI = FOI(FOI <= 2);
  [~, indsFOI] = ismember(reducedFOI, freq);
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if ~isempty(pupilPhase{iCond}{iAreasOI(iArea)})
        rSpearman = significantData.arearSpearman5secIndividual{iCond}{iAreasOI(iArea)};
        pvalSpearman = significantData.areapvalSpearman5secIndividual{iCond}{iAreasOI(iArea)};
        options2 = options;
        options2.edgeColourVector = zeros(size(rSpearman));
        options2.edgeColourVector(rSpearman >= 0 & pvalSpearman < 0.05) = 1;
        options2.edgeColourVector(rSpearman >= 0 & pvalSpearman >= 0.05) = 2;
        options2.edgeColourVector(rSpearman < 0 & pvalSpearman < 0.05) = 3;
        options2.edgeColourVector(rSpearman < 0 & pvalSpearman >= 0.05) = 4;
        options2.edgeColourCodes = zeros(numel(unique(options2.edgeColourVector)), 3);
        options2.edgeColourCodes(1,:) = matlabColours(10); % Red
        options2.edgeColourCodes(2,:) = [255 128 128]./255; % Light red
        options2.edgeColourCodes(3,:) = matlabColours(8); % Blue
        options2.edgeColourCodes(4,:) = [128 128 255]./255; % Light blue
        
        phase = pupilPhase{iCond}{iAreasOI(iArea)}(:,indsFOI);
        for f = 1:numel(reducedFOI)
          fPhase = recentrePhase(phase(:,f), pi);
          indsNaN = isnan(fPhase);
          if sum(indsNaN) < numel(fPhase)
            options2.faceColourVector = hsv2rgb([fPhase./(2*pi) ones(size(fPhase)) ones(size(fPhase))]);
            options2.faceColourVector(isnan(options2.faceColourVector)) = 1;
            
            % 50 cent
            if ~isempty(significantData.areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)})
              options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
              options2.figName = ['frLogConstrictedVsfrLogDilated50Extra in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond} ' at frequency ' num2str(reducedFOI(f)) 'Hz'];
              options2.xLabel = 'Constricted log_{10}(firing rate)';
              options2.yLabel = 'Dilated log_{10}(firing rate)';
              options2.axesType = 'regular';
              corrPlot(getLog(significantData.areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,1)), getLog(significantData.areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,end)), options2);
            end
            
            % 25 cent
            if ~isempty(significantData.areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)})
              options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
              options2.figName = ['frLogConstrictedVsfrLogDilated25Extra in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond} ' at frequency ' num2str(reducedFOI(f)) 'Hz'];
              options2.xLabel = 'Constricted log_{10}(firing rate)';
              options2.yLabel = 'Dilated log_{10}(firing rate)';
              options2.axesType = 'regular';
              corrPlot(getLog(significantData.areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,1)), getLog(significantData.areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,end)), options2);
            end
            
            % 12.5 cent
            if ~isempty(significantData.areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)})
              options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
              options2.figName = ['frLogConstrictedVsfrLogDilated12p5Extra in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond} ' at frequency ' num2str(reducedFOI(f)) 'Hz'];
              options2.xLabel = 'Constricted log_{10}(firing rate)';
              options2.yLabel = 'Dilated log_{10}(firing rate)';
              options2.axesType = 'regular';
              corrPlot(getLog(significantData.areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,1)), getLog(significantData.areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,end)), options2);
            end
            
            % 33.33 cent
            if ~isempty(significantData.areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)})
              options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
              options2.figName = ['frLogConstrictedVsfrLogDilatedThirdExtra in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond} ' at frequency ' num2str(reducedFOI(f)) 'Hz'];
              options2.xLabel = 'Constricted log_{10}(firing rate)';
              options2.yLabel = 'Dilated log_{10}(firing rate)';
              options2.axesType = 'regular';
              corrPlot(getLog(significantData.areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,1)), getLog(significantData.areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,end)), options2);
            end
          end
        end
      end
    end
  end
end


% LOG FIRING RATE CORRELATIONS: FULL PHASE
if drawCorrelationsLog(3)
  options = struct();
  options.corrType = 'Spearman';
  options.figFolder = mainFolder;
  options.diagonal = true;
  options.fitLine = true;
  options.fitLineDisplay = 'on';
  options.xLim = [];
  options.xTicks = [];
  options.yLim = [];
  options.yTicks = [];
  options.figSize = figSize;
  options.saveFig = true;
  pupilPhaseSignificant = eyeData.areaPhaseFullInterpIndividual;
  pupilPhase = fullData.areaPhaseEyeFullInterpIndividual;
  freq = fullData.areaFreqEyeFullInterpIndividual;
  eyeDataExists = fullData.areaEyeDataExistsIndividual;
  reducedFOI = FOI(FOI <= 2);
  [~, indsFOI] = ismember(reducedFOI, freq);
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if ~isempty(pupilPhase{iCond}{iAreasOI(iArea)})
        rSpearman = fullData.arearSpearman5secIndividual{iCond}{iAreasOI(iArea)}(logical(eyeDataExists{iCond}{iAreasOI(iArea)}));
        pvalSpearman = fullData.areapvalSpearman5secIndividual{iCond}{iAreasOI(iArea)}(logical(eyeDataExists{iCond}{iAreasOI(iArea)}));
        options2 = options;
        options2.edgeColourVector = zeros(size(rSpearman));
        options2.edgeColourVector(rSpearman >= 0 & pvalSpearman < 0.05) = 1;
        options2.edgeColourVector(rSpearman >= 0 & pvalSpearman >= 0.05) = 2;
        options2.edgeColourVector(rSpearman < 0 & pvalSpearman < 0.05) = 3;
        options2.edgeColourVector(rSpearman < 0 & pvalSpearman >= 0.05) = 4;
        options2.edgeColourCodes = zeros(numel(unique(options2.edgeColourVector)), 3);
        options2.edgeColourCodes(1,:) = matlabColours(10); % Red
        options2.edgeColourCodes(2,:) = [255 128 128]./255; % Light red
        options2.edgeColourCodes(3,:) = matlabColours(8); % Blue
        options2.edgeColourCodes(4,:) = [128 128 255]./255; % Light blue
        
        phaseSignificant = pupilPhaseSignificant{iCond}{iAreasOI(iArea)}(:,indsFOI);
        phase = pupilPhase{iCond}{iAreasOI(iArea)}(logical(eyeDataExists{iCond}{iAreasOI(iArea)}),indsFOI);
        for f = 1:numel(reducedFOI)
          fPhaseSignificant = recentrePhase(phaseSignificant(:,f), pi);
          fPhase = recentrePhase(phase(:,f), pi);
          indsNaN = isnan(fPhase);
          if sum(indsNaN) < numel(fPhase)
            saturation = ones(size(fPhase));
            saturation(isnan(fPhaseSignificant)) = 0.5;
            options2.faceColourVector = hsv2rgb([fPhase./(2*pi) saturation ones(size(fPhase))]);
            options2.faceColourVector(isnan(options2.faceColourVector)) = 1;
            
            % 50 cent
            if ~isempty(significantData.areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)})
              options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
              options2.figName = ['frLogConstrictedVsfrLogDilated50Full in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond} ' at frequency ' num2str(reducedFOI(f)) 'Hz'];
              options2.xLabel = 'Constricted log_{10}(firing rate)';
              options2.yLabel = 'Dilated log_{10}(firing rate)';
              options2.axesType = 'regular';
              corrPlot(getLog(significantData.areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,1)), getLog(significantData.areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,end)), options2);
            end
            
            % 25 cent
            if ~isempty(significantData.areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)})
              options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
              options2.figName = ['frLogConstrictedVsfrLogDilated25Full in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond} ' at frequency ' num2str(reducedFOI(f)) 'Hz'];
              options2.xLabel = 'Constricted log_{10}(firing rate)';
              options2.yLabel = 'Dilated log_{10}(firing rate)';
              options2.axesType = 'regular';
              corrPlot(getLog(significantData.areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,1)), getLog(significantData.areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,end)), options2);
            end
            
            % 12.5 cent
            if ~isempty(significantData.areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)})
              options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
              options2.figName = ['frLogConstrictedVsfrLogDilated12p5Full in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond} ' at frequency ' num2str(reducedFOI(f)) 'Hz'];
              options2.xLabel = 'Constricted log_{10}(firing rate)';
              options2.yLabel = 'Dilated log_{10}(firing rate)';
              options2.axesType = 'regular';
              corrPlot(getLog(significantData.areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,1)), getLog(significantData.areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,end)), options2);
            end
            
            % 33.33 cent
            if ~isempty(significantData.areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)})
              options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
              options2.figName = ['frLogConstrictedVsfrLogDilatedThirdFull in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond} ' at frequency ' num2str(reducedFOI(f)) 'Hz'];
              options2.xLabel = 'Constricted log_{10}(firing rate)';
              options2.yLabel = 'Dilated log_{10}(firing rate)';
              options2.axesType = 'regular';
              corrPlot(getLog(significantData.areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,1)), getLog(significantData.areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,end)), options2);
            end
          end
        end
      end
    end
  end
end


%% LOG FIRING RATE VS SIGNIFICANT PHASE
if drawCorrelationsLog(4)
  options = struct();
  options.corrType = 'circlinearnp';
  options.figFolder = mainFolder;
  options.diagonal = true;
  options.fitLine = true;
  options.fitLineDisplay = 'on';
  options.xLim = [];
  options.xTicks = [];
  options.yLim = [];
  options.yTicks = [];
  options.figSize = figSize;
  options.saveFig = true;
  pupilPhase = fullData.areaPhaseEyeFullInterpIndividual;
  freq = eyeData.areaFreqFullInterpIndividual;
  eyeDataExists = fullData.areaEyeDataExistsIndividual;
  reducedFOI = FOI(FOI <= 2);
  [~, indsFOI] = ismember(reducedFOI, freq);
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if ~isempty(pupilPhase{iCond}{iAreasOI(iArea)})
        fr = fullData.areaFRIndividual{iCond}{iAreasOI(iArea)}(logical(eyeDataExists{iCond}{iAreasOI(iArea)}),:);
        phase = pupilPhase{iCond}{iAreasOI(iArea)}(logical(eyeDataExists{iCond}{iAreasOI(iArea)}),indsFOI);
        for f = 1:numel(reducedFOI)
          fPhaseSignificant = recentrePhase(phase(:,f), pi);
          indsNaN = isnan(fPhaseSignificant);
          if sum(indsNaN) < numel(fPhaseSignificant)
            
            options.figTitle = ['Unit log10(firing rate) vs unit firing rate phase wrt pupil in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
            options.figName = ['frLogVsPhaseSignificant in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond} ' at frequency ' num2str(reducedFOI(f)) 'Hz'];
            options.xLabel = 'log_{10}(firing rate)';
            options.yLabel = 'Phase (rad)';
            options.axesType = 'regular';
            corrPlot(fr, fPhaseSignificant, options);
            
          end
        end
      end
    end
  end
end


if drawFractions
  
  %% Violin plots for units:
  if strcmp(repository, 'uol')
    area1 = 'Th'; area2 = 'VB'; area3 = 'lS1'; area4 = 'lRSC'; area5 = 'CA'; area6 = 'DG';
    areasOIFull = {area1; area2; area3; area4; area5; area6};
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    iCond = 1;
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} ' awake'];
      violinData{(iArea*2)-1} = data.positiveUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(violinData{(iArea*2)-1});
      violinAreas{(iArea*2)} = [areasOIFull{iArea} ' anaest'];
      violinData{(iArea*2)} = data.negativeUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(violinData{(iArea*2)});
      statsFractionsUnits{iArea} = meanTest([violinData{(iArea*2)-1} violinData{(iArea*2)}], 'ANOVARM', 'off');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 10.75]);
    end
    set(fH, 'Name','Proportion of units positively and negatively correlated to pupil size in different brain areas');
    filename = [mainFolder filesep 'pupilCorrFractionsUnitsViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
  end
  
  
  %% Violin plots for MUAs:
  dataMean = zeros(1, numel(areasOIFull)*2);
  dataCI95 = zeros(2, numel(areasOIFull)*2);
  iCond = 1;
  for iArea = 1:numel(areasOIFull)
    area = determineArea(areasOIFull{iArea});
    violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Pos'];
    violinData{(iArea*2)-1} = data.positiveMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
    [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(violinData{(iArea*2)-1});
    violinAreas{(iArea*2)} = [areasOIFull{iArea} 'Neg'];
    violinData{(iArea*2)} = data.negativeMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
    [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(violinData{(iArea*2)});
    statsFractionsMUAs{iArea} = meanTest([violinData{(iArea*2)-1} violinData{(iArea*2)}], 'ANOVARM', 'off');
  end
  options.yLim = [0 1];
  options.yLabel = 'Fraction';
  options.showNotches = false;
  options.medianPlot = false;
  fH2 = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
  if strcmp(repository, 'uol')
    xlim([0.25 12.75]);
  elseif strcmp(repository, 'allensdk')
    xlim([0.25 10.75]);
  end
  set(fH2, 'Name','Proportion of units and MUAs positively and negatively correlated to pupil size in different brain areas');
  filename = [mainFolder filesep 'pupilCorrFractionsMUAsViolins'];
  savefig(fH2, filename, 'compact');
  print(fH2, [filename '.png'],'-dpng','-r300');
end