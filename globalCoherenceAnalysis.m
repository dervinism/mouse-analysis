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
if (exist('fullRun', 'var') && ~islogical(fullRun) && strcmp(fullRun, 'custom')) || ~exist('fullRun', 'var')
  fullRun = true;
end
if ~exist('qualityCheck', 'var')
  qualityCheck = false;
end

outputDir = [outputDir filesep includeRuns];
if strcmp(repository,'uol')
  dataDir = [dataDir_local filesep '001_uol'];
elseif strcmp(repository,'allensdk')
  dataDir = [dataDir_local filesep '002_allen'];
end
if strcmp(repository,'all')
  if strcmp(subpop, 'all')
    rootFolder = [outputDir filesep laDir];
    rootFolder_ca = [outputDir filesep caDir];
  elseif strcmp(subpop, 'positive')
    rootFolder = [outputDir filesep laDir_positive];
    rootFolder_ca = [outputDir filesep caDir_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [outputDir filesep laDir_negative];
    rootFolder_ca = [outputDir filesep caDir_negative];
  end
  animals = animalsOI;
  xLim = [FOI(end) FOI(1)];
elseif strcmp(repository,'uol')
  if strcmp(subpop, 'all')
    rootFolder = [outputDir filesep laDir_uol];
    rootFolder_ca = [outputDir filesep caDir_uol];
  elseif strcmp(subpop, 'positive')
    rootFolder = [outputDir filesep laDir_uol_positive];
    rootFolder_ca = [outputDir filesep caDir_uol_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [outputDir filesep laDir_uol_negative];
    rootFolder_ca = [outputDir filesep caDir_uol_negative];
  end
  animals = animalsUOLOI;
  xLim = [FOI(end) FOI(1)];
elseif strcmp(repository,'allensdk')
  if strcmp(subpop, 'all')
    rootFolder = [outputDir filesep laDir_allensdk];
    rootFolder_ca = [outputDir filesep caDir_allensdk];
  elseif strcmp(subpop, 'positive')
    rootFolder = [outputDir filesep laDir_allensdk_positive];
    rootFolder_ca = [outputDir filesep caDir_allensdk_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [outputDir filesep laDir_allensdk_negative];
    rootFolder_ca = [outputDir filesep caDir_allensdk_negative];
  end
  animals = animalsAllensdk;
  xLim = [0.045 FOI(1)];
  conditions = {'awake'};
end
if qualityCheck
  mainFolder = [rootFolder filesep qualityUnitsFolder filesep coherenceCorrelationsSubfolder]; %#ok<*UNRCH>
  mainFolder_ca = [rootFolder_ca filesep qualityUnitsFolder filesep coherenceCorrelationsSubfolder];
else
  mainFolder = [rootFolder filesep unitsFolder filesep coherenceCorrelationsSubfolder];
  mainFolder_ca = [rootFolder_ca filesep unitsFolder filesep coherenceCorrelationsSubfolder];
end
areas = areas2compare;

individualGraphs = true;
summaryGraphs = false;


%% COMPUTE VARIABLES
if fullRun
  for animal = 1:numel(animals) % Loop through animals
    disp(['processing animal ' animals{animal}]);
    if exist('allData', 'var')
      dataStruct = allData.(animals{animal}).dataStruct;
    else
      load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
    end
    if strcmp(subpop, 'all')
      if ~isfield(dataStruct, 'seriesData_ca') || isempty(dataStruct.seriesData_ca)
        disp(['Animal ' animals{animal} ' is missing series comparison data. Skipping to the next animal...'])
        continue
      end
      fnsData_ca = fieldnames(dataStruct.seriesData_ca);
    elseif strcmp(subpop, 'positive')
      if ~isfield(dataStruct, 'seriesData_ca_positive') || isempty(dataStruct.seriesData_ca_positive)
        disp(['Animal ' animals{animal} ' is missing series comparison data. Skipping to the next animal...'])
        continue
      end
      fnsData_ca = fieldnames(dataStruct.seriesData_ca_positive);
    elseif strcmp(subpop, 'negative')
      if ~isfield(dataStruct, 'seriesData_ca_negative') || isempty(dataStruct.seriesData_ca_negative)
        disp(['Animal ' animals{animal} ' is missing series comparison data. Skipping to the next animal...'])
        continue
      end
      fnsData_ca = fieldnames(dataStruct.seriesData_ca_negative);
    end
    
    % Initialise variables
    try
      if strcmp(subpop, 'all')
        fnsData = fieldnames(dataStruct.seriesData);
        FOI = dataStruct.seriesData.(fnsData{1}).conf.FOI;
      elseif strcmp(subpop, 'positive')
        fnsData = fieldnames(dataStruct.seriesData_positive);
        FOI = dataStruct.seriesData_positive.(fnsData{1}).conf.FOI;
      elseif strcmp(subpop, 'negative')
        fnsData = fieldnames(dataStruct.seriesData_negative);
        FOI = dataStruct.seriesData_negative.(fnsData{1}).conf.FOI;
      end
    catch
      % do nothing
    end
    if animal == 1 || ~exist('areaCohFullInterpIndividual_la2ca', 'var')
      areaCohFullInterpIndividual_la2ca = {};
      areaCohFOIindividual_la2ca = {};
      areaCohFullIndividual = {};
      areaFreqFullIndividual = {};
      for iCond = 1:numel(conditions)
        areaCohFullInterpIndividual_la2caCond = {};
        areaCohFOIindividual_la2caCond = {};
        areaCohFullIndividualCond = {};
        areaFreqFullIndividualCond = {};
        for iArea = 1:numel(areas)
          areaCohFullInterpIndividual_la2caCond{iArea} = []; %#ok<*SAGROW>
          areaCohFOIindividual_la2caCond{iArea} = [];
          areaCohFullIndividualCond{iArea} = {};
          areaFreqFullIndividualCond{iArea} = {};
        end
        areaCohFullInterpIndividual_la2ca{iCond} = areaCohFullInterpIndividual_la2caCond;
        areaCohFOIindividual_la2ca{iCond} = areaCohFOIindividual_la2caCond;
        areaCohFullIndividual{iCond} = areaCohFullIndividualCond;
        areaFreqFullIndividual{iCond} = areaFreqFullIndividualCond;
      end
    end
    
    for dbCount = 1:numel(fnsData_ca) % Loop through db entries
      if strcmp(subpop, 'all')
        dbStruct_ca = dataStruct.seriesData_ca.(fnsData_ca{dbCount});
      elseif strcmp(subpop, 'positive')
        dbStruct_ca = dataStruct.seriesData_ca_positive.(fnsData_ca{dbCount});
      elseif strcmp(subpop, 'negative')
        dbStruct_ca = dataStruct.seriesData_ca_negative.(fnsData_ca{dbCount});
      end
      
      % Determine if series coherence data exist
      if exist('seriesName1', 'var')
        prevRec = seriesName1(1:14);
      else
        prevRec = '';
      end
      [seriesName1, seriesName2] = seriesNames(fnsData_ca{dbCount});
      recording = seriesName1(1:14);
      if ~isfield(dbStruct_ca, 'popData') || ~isfield(dbStruct_ca.popData, 'phaseCoh') || isempty(dbStruct_ca.popData.phaseCoh)
        continue
      end
      
      % Test for exceptions
      if exceptionTest(except, seriesName1, seriesName2)
        continue
      end
      
      % Determine if population rate > 0
      if strcmp(subpop, 'all')
        [breakClause, spkDBmfr] = firingRateTest(sum(dataStruct.seriesData.([animals{animal} '_s' seriesName1]).popData.MUAsAll,1),...
          dataStruct.seriesData.([animals{animal} '_s' seriesName1]).conf.samplingParams.srData);
      elseif strcmp(subpop, 'positive')
        [breakClause, spkDBmfr] = firingRateTest(sum(dataStruct.seriesData_positive.([animals{animal} '_s' seriesName1]).popData.MUAsAll,1),...
          dataStruct.seriesData_positive.([animals{animal} '_s' seriesName1]).conf.samplingParams.srData);
      elseif strcmp(subpop, 'negative')
        [breakClause, spkDBmfr] = firingRateTest(sum(dataStruct.seriesData_negative.([animals{animal} '_s' seriesName1]).popData.MUAsAll,1),...
          dataStruct.seriesData_negative.([animals{animal} '_s' seriesName1]).conf.samplingParams.srData);
      end
      if breakClause
        continue
      end
      if strcmp(subpop, 'all')
        [breakClause, PRmfr] = firingRateTest(sum(dataStruct.seriesData.([animals{animal} '_s' seriesName2]).popData.MUAsAll,1),...
          dataStruct.seriesData.([animals{animal} '_s' seriesName2]).conf.samplingParams.srData);
      elseif strcmp(subpop, 'positive')
        [breakClause, PRmfr] = firingRateTest(sum(dataStruct.seriesData_positive.([animals{animal} '_s' seriesName2]).popData.MUAsAll,1),...
          dataStruct.seriesData_positive.([animals{animal} '_s' seriesName2]).conf.samplingParams.srData);
      elseif strcmp(subpop, 'negative')
        [breakClause, PRmfr] = firingRateTest(sum(dataStruct.seriesData_negative.([animals{animal} '_s' seriesName2]).popData.MUAsAll,1),...
          dataStruct.seriesData_negative.([animals{animal} '_s' seriesName2]).conf.samplingParams.srData);
      end
      if breakClause
        continue
      end
      
      % Determine area comparison and any grouped area comparisons
      if strcmp(repository, 'all')
        breakClause = true;
      elseif strcmp(repository, 'uol')
        [breakClause, comp, compNames, areasReverse] = series2Comparison(areas, seriesName1, seriesName2, reverse, true);
      elseif strcmp(repository, 'allensdk')
        [breakClause, comp, compNames, areasReverse] = series2Comparison(areas, seriesName1, seriesName2, reverse);
      end
      if breakClause
        continue
      end
      
      % Determine recording condition (i.e., awake or anaesthesia)
      [breakClause, iCond] = series2condition(awake, anaesthesia, seriesName1, seriesName2);
      if breakClause
        continue
      end
      
      % Disqualify low quality units if needed
      if strcmp(subpop, 'all')
        dbStruct = dataStruct.seriesData.([animals{animal} '_s' seriesName1]);
      elseif strcmp(subpop, 'positive')
        dbStruct = dataStruct.seriesData_positive.([animals{animal} '_s' seriesName1]);
      elseif strcmp(subpop, 'negative')
        dbStruct = dataStruct.seriesData_negative.([animals{animal} '_s' seriesName1]);
      end
      units_ca = [];
      units = [];
      for sh = 1:numel(dbStruct_ca.shankData)
        units_ca = [units_ca; dbStruct_ca.shankData.(['shank' num2str(sh)]).units]; %#ok<*AGROW>
        units = [units; dbStruct.shankData.(['shank' num2str(sh)]).units];
      end
      assert(numel(units_ca) == numel(units));
      if qualityCheck
        fullSeriesName1 = [animals{animal} '_s' seriesName1];
        unitMetadata = [];
        for sh = 1:numel(dbStruct_ca.shankData)
          if strcmp(subpop, 'all')
            unitMetadata = [unitMetadata; dataStruct.seriesData.(fullSeriesName1).shankData.(['shank' num2str(sh)]).unitMetadata];
          elseif strcmp(subpop, 'positive')
            unitMetadata = [unitMetadata; dataStruct.seriesData_positive.(fullSeriesName1).shankData.(['shank' num2str(sh)]).unitMetadata];
          elseif strcmp(subpop, 'negative')
            unitMetadata = [unitMetadata; dataStruct.seriesData_negative.(fullSeriesName1).shankData.(['shank' num2str(sh)]).unitMetadata];
          end
        end
        [~, qualityUnitInd] = qualityTest2(unitMetadata, cluDist, refractCont, false);
      else
        qualityUnitInd = 1:numel(units_ca);
      end
      if isempty(qualityUnitInd)
        continue
      end
      units_ca = units_ca(qualityUnitInd);
      
      for iComp = 1:numel(comp) % Loop through non-grouped and grouped area comparisons
        area = comp(iComp);
        
        if numel(conditions) > 1
          condLoop = [iCond numel(conditions)];
        else
          condLoop = iCond;
        end
        for iCondPlusAll = condLoop % Loop through the main and pooled conditions
          
          shankIDs = fieldnames(dbStruct.shankData);
          for sh = 1:numel(shankIDs)
            shankStruct = dbStruct.shankData.(shankIDs{sh});
            for u = 1:numel(units_ca)
              iU = find(shankStruct.units == units_ca(u));
              if isempty(iU)
                continue
              end
              
              freq = shankStruct.phaseCoh{iU}.freq;  
              
              % Get coherence values
              if isfield(shankStruct.phaseCoh{iU}, 'coh')
                coh = shankStruct.phaseCoh{iU}.coh;
                cohConf = shankStruct.phaseCoh{iU}.coh_conf;
                rateadjust_kappa = shankStruct.phaseCoh{iU}.rateadjust_kappa;
                phase = shankStruct.phaseCoh{iU}.phase;
                phaseConfU = shankStruct.phaseCoh{iU}.phase_confU;
                phaseConfL = shankStruct.phaseCoh{iU}.phase_confL;
                [~,~,~, coh] = correctPhaseCoh(phase, phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa);
                areaCohFullIndividual{iCondPlusAll}{area}{numel(areaCohFullIndividual{iCondPlusAll}{area})+1} = coh;
              else
                areaCohFullIndividual{iCondPlusAll}{area}{numel(areaCohFullIndividual{iCondPlusAll}{area})+1} = NaN(size(freq));
              end
              areaFreqFullIndividual{iCondPlusAll}{area}{numel(areaFreqFullIndividual{iCondPlusAll}{area})+1} = freq;
            end
          end
        end
      end
    end
  end
  
  areas_la2ca = areasReverse;
  
  % Interpolate and store full coherences
  freqCombined = FOI;
  for iCond = 1:numel(conditions)
    for iArea = 1:numel(areas_la2ca)
      nRec = numel(areaFreqFullIndividual{iCond}{iArea});
      for iRec = 1:nRec
        freqCombined = unique([freqCombined areaFreqFullIndividual{iCond}{iArea}{iRec}]);
      end
    end
  end
  freqCombined = unique(freqCombined(~isnan(freqCombined)));
  for iCond = 1:numel(conditions)
    for iArea = 1:numel(areas_la2ca)
      nRec = numel(areaFreqFullIndividual{iCond}{iArea});
      areaCohFullInterpIndividual_la2ca{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      for iRec = 1:nRec
        if ~isempty(areaCohFullIndividual{iCond}{iArea}) && sum(~isnan(areaCohFullIndividual{iCond}{iArea}{iRec}))
          areaCohFullInterpIndividual_la2ca{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaCohFullIndividual{iCond}{iArea}{iRec}, freqCombined);
        end
      end
    end
  end
  areaFreqFullInterpIndividual_la2ca = freqCombined;
  
  % Obtain coherence values for FOI
  FOIInds_la2ca = zeros(size(FOI));
  for f = 1:numel(FOI)
    FOIInds_la2ca(f) = find(FOI(f) == areaFreqFullInterpIndividual_la2ca);
  end
  for iCond = 1:numel(conditions)
    for iArea = 1:numel(areas_la2ca)
      if ~isempty(areaCohFullInterpIndividual_la2ca{iCond}{iArea})
        areaCohFOIindividual_la2ca{iCond}{iArea} = areaCohFullInterpIndividual_la2ca{iCond}{iArea}(:,FOIInds_la2ca);
      else
        areaCohFOIindividual_la2ca{iCond}{iArea} = [];
      end
    end
  end
end


%% SAVE OR LOAD THE DATA
clearvars -except individualGraphs summaryGraphs mainFolder mainFolder_ca FOI conditions areas_la2ca areaCohFullInterpIndividual_la2ca areaFreqFullInterpIndividual_la2ca areaCohFOIindividual_la2ca repository subpop reverse qualityCheck allData fullRun areas2compareCritical includeRuns
if reverse
  if qualityCheck
    filename = [fileparts(mainFolder) filesep 'globalUnits_quality.mat'];
    filename_ca = [fileparts(mainFolder_ca) filesep 'globalUnits_ca_reverse_quality.mat'];
  else
    filename = [fileparts(mainFolder) filesep 'globalUnits.mat'];
    filename_ca = [fileparts(mainFolder_ca) filesep 'globalUnits_ca_reverse.mat'];
  end
else
  if qualityCheck
    filename = [fileparts(mainFolder) filesep 'globalUnits_quality.mat'];
    filename_ca = [fileparts(mainFolder_ca) filesep 'globalUnits_ca_quality.mat'];
  else
    filename = [fileparts(mainFolder) filesep 'globalUnits.mat'];
    filename_ca = [fileparts(mainFolder_ca) filesep 'globalUnits_ca.mat'];
  end
end
if fullRun
  if ~exist(fileparts(mainFolder), 'file')
    mkdir(fileparts(mainFolder));
    save(filename, 'areas_la2ca','areaCohFullInterpIndividual_la2ca','areaFreqFullInterpIndividual_la2ca','areaCohFOIindividual_la2ca');
  else
    save(filename, 'areas_la2ca','areaCohFullInterpIndividual_la2ca','areaFreqFullInterpIndividual_la2ca','areaCohFOIindividual_la2ca', '-append');
  end
else
  load(filename);
end
iAreas2compareOI = find(ismember(areas_la2ca,areas2compareCritical));

globalVariables = load(filename_ca);
FOIInds_ca = zeros(size(FOI));
for f = 1:numel(FOI)
  FOIInds_ca(f) = find(FOI(f) == globalVariables.areaFreqFullInterpIndividual);
end



%% COHERENCE CORRELATIONS
if ~exist(mainFolder, 'file')
  mkdir(mainFolder);
end
for iCond = 1:min([2 numel(conditions)])
  for iArea = 1:numel(iAreas2compareOI)
    coh1 = areaCohFOIindividual_la2ca{iCond}{iAreas2compareOI(iArea)};
    if ~isempty(globalVariables.areaCohFullInterpIndividual{iCond}{iAreas2compareOI(iArea)})
      coh2 = globalVariables.areaCohFullInterpIndividual{iCond}{iAreas2compareOI(iArea)}(:,FOIInds_ca);
    else
      coh2 = [];
    end
    disp(['Processing unit coherence data for ' conditions{iCond} ' ' areas_la2ca{iAreas2compareOI(iArea)}...
      ' (comparison # ' num2str((iCond-1)*numel(areas_la2ca) + iAreas2compareOI(iArea)) '/' num2str(numel(conditions)*numel(areas_la2ca)) ')']);
    if ~isempty(coh1) && ~isempty(coh2)
      assert(numel(coh1) == numel(coh2));
      
      % Correlations and individual graphs
      area1 = areas_la2ca{iAreas2compareOI(iArea)}(1:strfind(areas_la2ca{iAreas2compareOI(iArea)},'Vs')-1);
      opt.xLabel = [area1 ' unit coherence with local PR'];
      opt.yLabel = [areas_la2ca{iAreas2compareOI(iArea)} ' unit coherence with outside PR'];
      [figCoh, rFOI_laVca{iCond}{iAreas2compareOI(iArea)}, pvalFOI_laVca{iCond}{iAreas2compareOI(iArea)}, nFOI_laVca{iCond}{iAreas2compareOI(iArea)}] = halfCorrPlot_coh(...
        coh1, coh2, FOI, conditions{iCond}, area1, areas_la2ca{iAreas2compareOI(iArea)}, mainFolder,...
        'Local and cross-area coherence correlations for ', individualGraphs, 'SPIKING', '', opt); %#ok<*SAGROW>
      close all
      
      % Summary subplots
      if summaryGraphs
        sbCoh = halfCorrSummary_coh(coh1, coh2, FOI, rFOI_laVca{iCond}{iAreas2compareOI(iArea)}, pvalFOI_laVca{iCond}{iAreas2compareOI(iArea)},...
          nFOI_laVca{iCond}{iAreas2compareOI(iArea)}, conditions{iCond}, area1, areas_la2ca{iAreas2compareOI(iArea)}, mainFolder,...
          '_local_and_crossarea_correlations_at_p01_p03_p05_p1_p3_p5_1_4_10_Hz', 'SPIKING', '');
        close(sbCoh)
      end
    end
  end
end

% Save the results of correlation analyses
save(filename, 'rFOI_laVca','pvalFOI_laVca','nFOI_laVca', '-append');