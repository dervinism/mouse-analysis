clearvars -except repository subpop reverse qualityCheck allData fullRun includeRuns


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

dataDir = [dataDir filesep includeRuns];
if strcmp(repository,'all')
  error('globalConditionsAnalysis is only available for uol repository')
elseif strcmp(repository,'uol')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep laDir_uol];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep laDir_uol_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep laDir_uol_negative];
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

drawFRCorrelations = true;
individualPhaseGraphs = true;
summaryPhaseGraphs = true;
individualCohGraphs = true;
summaryCohGraphs = true;


%% COMPUTE VARIABLES AND INITIALISE FIGURES NEEDED FOR COMPARING BETWEEN VIGILANCE CONDITIONS
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
    try
      if strcmp(subpop, 'all')
        fnsData = fieldnames(dataStruct.seriesData);
        FOI = dataStruct.seriesData.(fnsData{1}).conf.FOI;
        srData = dataStruct.seriesData.(fnsData{1}).conf.samplingParams.srData;
      elseif strcmp(subpop, 'positive')
        fnsData = fieldnames(dataStruct.seriesData_positive);
        FOI = dataStruct.seriesData_positive.(fnsData{1}).conf.FOI;
        srData = dataStruct.seriesData_positive.(fnsData{1}).conf.samplingParams.srData;
      elseif strcmp(subpop, 'negative')
        fnsData = fieldnames(dataStruct.seriesData_negative);
        FOI = dataStruct.seriesData_negative.(fnsData{1}).conf.FOI;
        srData = dataStruct.seriesData_negative.(fnsData{1}).conf.samplingParams.srData;
      end
    catch
      % do nothing
    end
    if animal == 1 || ~exist('areaFRIndividual', 'var')
      areaFRIndividual = {};
      areaCohFullIndividual = {};
      areaCohConfUFullIndividual = {};
      areaCohConfLFullIndividual = {};
      areaPhaseFullIndividual = {};
      areaPhaseConfUFullIndividual = {};
      areaPhaseConfLFullIndividual = {};
      areaFreqFullIndividual = {};
      areaCohFullInterpIndividual = {};
      areaCohConfUFullInterpIndividual = {};
      areaCohConfLFullInterpIndividual = {};
      areaPhaseFullInterpIndividual = {};
      areaPhaseConfUFullInterpIndividual = {};
      areaPhaseConfLFullInterpIndividual = {};
      for iCond = 1:2
        for iArea = 1:numel(areas)
          areaFRIndividual{iCond}{iArea} = []; %#ok<*SAGROW>
          areaCohFullIndividual{iCond}{iArea} = {};
          areaCohConfUFullIndividual{iCond}{iArea} = {};
          areaCohConfLFullIndividual{iCond}{iArea} = {};
          areaPhaseFullIndividual{iCond}{iArea} = {};
          areaPhaseConfUFullIndividual{iCond}{iArea} = {};
          areaPhaseConfLFullIndividual{iCond}{iArea} = {};
          areaFreqFullIndividual{iCond}{iArea} = {};
          areaCohFullInterpIndividual{iCond}{iArea} = [];
          areaCohConfUFullInterpIndividual{iCond}{iArea} = [];
          areaCohConfLFullInterpIndividual{iCond}{iArea} = [];
          areaPhaseFullInterpIndividual{iCond}{iArea} = [];
          areaPhaseConfUFullInterpIndividual{iCond}{iArea} = [];
          areaPhaseConfLFullInterpIndividual{iCond}{iArea} = [];
        end
      end
    end
    
    for dbCount = 1:numel(fnsData) % Loop through db entries
      dbStruct = dataStruct.seriesData.(fnsData{dbCount});
      
      % Test the wakefulness recording for inclusion
      seriesName = seriesFromEntry(fnsData{dbCount});
      [continueClause, dbCount2] = exceptionTest(awakePaired, seriesName(1:min([numel(seriesName) 14])));
      if ~continueClause
        continue
      end
      
      % Determine if series phase and coherence data exist
      if ~isfield(dbStruct, 'popData') || ~isfield(dbStruct.popData, 'phaseCoh') || isempty(dbStruct.popData.phaseCoh)
        continue
      end
      
      % Test for exceptions
      if exceptionTest(except, seriesName)
        continue
      end
      
      % Determine if population rate > 0
      if firingRateTest(sum(dbStruct.popData.MUAsAll,1), srData)
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
          if iCond == 1
            if ~isfield(dbStruct.shankData.(['shank' num2str(sh)]), 'pupil')
              continue
            else
              phase = [phase; spkPhase(dbStruct.shankData.(['shank' num2str(sh)]).pupil.unitData, fRef)'];
            end
          elseif iCond == 2
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
        if iCond == 1
          correlatedInd = false(size(phase));
          correlatedInd(recentrePhase(phase, 0) > -pi/2 & recentrePhase(phase, 0) <= pi/2) = true;
        elseif iCond == 2
          correlatedInd = find(rSpearman >= 0);
        end
      elseif strcmp(subpop, 'negative')
        if iCond == 1
          correlatedInd = false(size(phase));
          correlatedInd(recentrePhase(phase, pi) > pi/2 & recentrePhase(phase, pi) <= 3*pi/2) = true;
        elseif iCond == 2
          correlatedInd = find(rSpearman < 0);
        end
      end
      if qualityCheck
        unitMetadata = [];
        for sh = 1:numel(fieldnames(dbStruct.shankData))
          unitMetadata = [unitMetadata; dbStruct.shankData.(['shank' num2str(sh)]).unitMetadata];
        end
        [~, qualityUnitInd] = qualityTest2(unitMetadata, cluDist, refractCont, false);
      else
        qualityUnitInd = 1:numel(units);
      end
      qualityUnitInd = intersect(qualityUnitInd, find(correlatedInd));
      if isempty(qualityUnitInd)
        continue
      end
      units = units(qualityUnitInd);
      
      % Identify anaesthesia recording
      if numel(seriesName) > 14
        seriesName2 = [anaesthesiaPaired{dbCount2} seriesName(15:end)];
      else
        seriesName2 = anaesthesiaPaired{dbCount2};
      end
      dbCount2 = find(ismember(fnsData, [animals{animal} '_s' seriesName2]));
      assert(numel(dbCount2) == 1);
      dbStruct2 = dataStruct.seriesData.(fnsData{dbCount2});
      
      % Determine if series phase and coherence data exist
      if ~isfield(dbStruct2, 'popData') || ~isfield(dbStruct2.popData, 'phaseCoh') || isempty(dbStruct2.popData.phaseCoh)
        continue
      end
      
      % Test for exceptions
      if exceptionTest(except, seriesName)
        continue
      end
      
      % Determine if population rate > 0
      if firingRateTest(sum(dbStruct2.popData.MUAsAll,1), srData)
        continue
      end
      
      % Identify units present during both wakefulness and anaesthesia recording conditions
      units2 = [];
      for sh = 1:numel(dbStruct2.shankData)
        units2 = [units2; dbStruct2.shankData.(['shank' num2str(sh)]).units]; %#ok<*AGROW>
      end
      [unitsOverlap, iUnits1, iUnits2] = intersect(units,units2,'stable');
      
      % Get unit spiking data
      shankIDs = fieldnames(dbStruct.shankData);
      spk = [];
      spk2 = [];
      for sh = 1:numel(shankIDs)
        if isempty(spk)
          spk = full(dbStruct.shankData.(shankIDs{sh}).spk);
          spk2 = full(dbStruct2.shankData.(shankIDs{sh}).spk);
        else
          spk = concatenateMat(spk, full(dbStruct.shankData.(shankIDs{sh}).spk));
          spk2 = concatenateMat(spk2, full(dbStruct2.shankData.(shankIDs{sh}).spk));
        end
      end
      spk = spk(iUnits1,:);
      spk2 = spk2(iUnits2,:);

      disp(['           series ' seriesName ' and ' seriesName2 ': ' num2str(numel(unitsOverlap)) ' units']);
      
      for iAreaPlusAll = area % Loop through the main and pooled areas
        close all
        
        % Calculate full firing rates
        areaFRIndividual{1}{iAreaPlusAll} = [areaFRIndividual{1}{iAreaPlusAll}; (sum(spk,2).*srData)./size(spk,2)];
        areaFRIndividual{2}{iAreaPlusAll} = [areaFRIndividual{2}{iAreaPlusAll}; (sum(spk2,2).*srData)./size(spk2,2)];
        
        for sh = 1:numel(shankIDs)
          shankStruct = dbStruct.shankData.(shankIDs{sh});
          shankStruct2 = dbStruct2.shankData.(shankIDs{sh});
          for u = 1:numel(unitsOverlap)
            iU = find(shankStruct.units == unitsOverlap(u));
            if isempty(iU)
              continue
            end
            iU2 = find(shankStruct2.units == unitsOverlap(u));
            if isempty(iU2)
              continue
            end

            freq = shankStruct.phaseCoh{iU}.freq;
            freq2 = shankStruct2.phaseCoh{iU2}.freq;
            
            % Get phase and coherence values
            if isfield(shankStruct.phaseCoh{iU}, 'coh') && isfield(shankStruct2.phaseCoh{iU2}, 'coh')

              % Wakefulness
              coh = shankStruct.phaseCoh{iU}.coh;
              cohConf = shankStruct.phaseCoh{iU}.coh_conf;
              rateadjust_kappa = shankStruct.phaseCoh{iU}.rateadjust_kappa;
              phase = shankStruct.phaseCoh{iU}.phase;
              phaseConfU = shankStruct.phaseCoh{iU}.phase_confU;
              phaseConfL = shankStruct.phaseCoh{iU}.phase_confL;
              [phase, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase,...
                phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa);
              
              areaCohFullIndividual{1}{iAreaPlusAll}{numel(areaCohFullIndividual{1}{iAreaPlusAll})+1} = coh;
              areaCohConfUFullIndividual{1}{iAreaPlusAll}{numel(areaCohConfUFullIndividual{1}{iAreaPlusAll})+1} = cohConfU;
              areaCohConfLFullIndividual{1}{iAreaPlusAll}{numel(areaCohConfLFullIndividual{1}{iAreaPlusAll})+1} = cohConfL;
              areaPhaseFullIndividual{1}{iAreaPlusAll}{numel(areaPhaseFullIndividual{1}{iAreaPlusAll})+1} = phase;
              areaPhaseConfUFullIndividual{1}{iAreaPlusAll}{numel(areaPhaseConfUFullIndividual{1}{iAreaPlusAll})+1} = phaseConfU;
              areaPhaseConfLFullIndividual{1}{iAreaPlusAll}{numel(areaPhaseConfLFullIndividual{1}{iAreaPlusAll})+1} = phaseConfL;

              % Anaesthesia
              coh = shankStruct2.phaseCoh{iU2}.coh;
              cohConf = shankStruct2.phaseCoh{iU2}.coh_conf;
              rateadjust_kappa = shankStruct2.phaseCoh{iU2}.rateadjust_kappa;
              phase = shankStruct2.phaseCoh{iU2}.phase;
              phaseConfU = shankStruct2.phaseCoh{iU2}.phase_confU;
              phaseConfL = shankStruct2.phaseCoh{iU2}.phase_confL;
              [phase, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase,...
                phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa);
              
              areaCohFullIndividual{2}{iAreaPlusAll}{numel(areaCohFullIndividual{2}{iAreaPlusAll})+1} = coh;
              areaCohConfUFullIndividual{2}{iAreaPlusAll}{numel(areaCohConfUFullIndividual{2}{iAreaPlusAll})+1} = cohConfU;
              areaCohConfLFullIndividual{2}{iAreaPlusAll}{numel(areaCohConfLFullIndividual{2}{iAreaPlusAll})+1} = cohConfL;
              areaPhaseFullIndividual{2}{iAreaPlusAll}{numel(areaPhaseFullIndividual{2}{iAreaPlusAll})+1} = phase;
              areaPhaseConfUFullIndividual{2}{iAreaPlusAll}{numel(areaPhaseConfUFullIndividual{2}{iAreaPlusAll})+1} = phaseConfU;
              areaPhaseConfLFullIndividual{2}{iAreaPlusAll}{numel(areaPhaseConfLFullIndividual{2}{iAreaPlusAll})+1} = phaseConfL;
            else
              areaCohFullIndividual{1}{iAreaPlusAll}{numel(areaCohFullIndividual{1}{iAreaPlusAll})+1} = NaN(size(freq));
              areaCohConfUFullIndividual{1}{iAreaPlusAll}{numel(areaCohConfUFullIndividual{1}{iAreaPlusAll})+1} = NaN(size(freq));
              areaCohConfLFullIndividual{1}{iAreaPlusAll}{numel(areaCohConfLFullIndividual{1}{iAreaPlusAll})+1} = NaN(size(freq));
              areaPhaseFullIndividual{1}{iAreaPlusAll}{numel(areaPhaseFullIndividual{1}{iAreaPlusAll})+1} = NaN(size(freq));
              areaPhaseConfUFullIndividual{1}{iAreaPlusAll}{numel(areaPhaseConfUFullIndividual{1}{iAreaPlusAll})+1} = NaN(size(freq));
              areaPhaseConfLFullIndividual{1}{iAreaPlusAll}{numel(areaPhaseConfLFullIndividual{1}{iAreaPlusAll})+1} = NaN(size(freq));
              areaCohFullIndividual{2}{iAreaPlusAll}{numel(areaCohFullIndividual{2}{iAreaPlusAll})+1} = NaN(size(freq2));
              areaCohConfUFullIndividual{2}{iAreaPlusAll}{numel(areaCohConfUFullIndividual{2}{iAreaPlusAll})+1} = NaN(size(freq2));
              areaCohConfLFullIndividual{2}{iAreaPlusAll}{numel(areaCohConfLFullIndividual{2}{iAreaPlusAll})+1} = NaN(size(freq2));
              areaPhaseFullIndividual{2}{iAreaPlusAll}{numel(areaPhaseFullIndividual{2}{iAreaPlusAll})+1} = NaN(size(freq2));
              areaPhaseConfUFullIndividual{2}{iAreaPlusAll}{numel(areaPhaseConfUFullIndividual{2}{iAreaPlusAll})+1} = NaN(size(freq2));
              areaPhaseConfLFullIndividual{2}{iAreaPlusAll}{numel(areaPhaseConfLFullIndividual{2}{iAreaPlusAll})+1} = NaN(size(freq2));
            end
            areaFreqFullIndividual{1}{iAreaPlusAll}{numel(areaFreqFullIndividual{1}{iAreaPlusAll})+1} = freq;
            areaFreqFullIndividual{2}{iAreaPlusAll}{numel(areaFreqFullIndividual{2}{iAreaPlusAll})+1} = freq2;
          end
        end
      end
    end
  end

  % Interpolate and store full phases and coherences
  freqCombined = FOI;
  for iCond = 1:2
    for iArea = 1:numel(areas)
      nRec = numel(areaFreqFullIndividual{iCond}{iArea});
      for iRec = 1:nRec
        freqCombined = unique([freqCombined areaFreqFullIndividual{iCond}{iArea}{iRec}]);
      end
    end
  end
  freqCombined = unique(freqCombined(~isnan(freqCombined)));
  for iCond = 1:2
    for iArea = 1:numel(areas)
      nRec = numel(areaFreqFullIndividual{iCond}{iArea});
      areaPhaseFullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaPhaseConfUFullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaPhaseConfLFullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaCohFullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaCohConfUFullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaCohConfLFullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      for iRec = 1:nRec
        if ~isempty(areaPhaseFullIndividual{iCond}{iArea}) && sum(~isnan(areaPhaseFullIndividual{iCond}{iArea}{iRec}))
          areaPhaseFullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaPhaseFullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaPhaseConfUFullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaPhaseConfUFullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaPhaseConfLFullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaPhaseConfLFullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaCohFullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaCohFullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaCohConfUFullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaCohConfUFullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaCohConfLFullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaCohConfLFullIndividual{iCond}{iArea}{iRec}, freqCombined);
        end
      end
    end
  end
  areaFreqFullInterpIndividual = freqCombined;
end

% Determine the file name and either save or load the data
if qualityCheck
  filename = [rootFolder filesep 'globalUnitsConditions_quality.mat'];
else
  filename = [rootFolder filesep 'globalUnitsConditions.mat'];
end
if fullRun
  if ~exist(rootFolder, 'file')
    mkdir(rootFolder);
  end
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  if ~exist(filename, 'file')
    save(filename, 'conditions','areas','areaFRIndividual','areaCohFullIndividual','areaCohConfUFullIndividual','areaCohConfLFullIndividual',...
      'areaPhaseFullIndividual','areaPhaseConfUFullIndividual','areaPhaseConfLFullIndividual','areaFreqFullIndividual',...
      'areaCohFullInterpIndividual','areaCohConfUFullInterpIndividual','areaCohConfLFullInterpIndividual','areaPhaseFullInterpIndividual',...
      'areaPhaseConfUFullInterpIndividual','areaPhaseConfLFullInterpIndividual','areaFreqFullInterpIndividual');
  else
    save(filename, 'conditions','areas','areaFRIndividual','areaCohFullIndividual','areaCohConfUFullIndividual','areaCohConfLFullIndividual',...
      'areaPhaseFullIndividual','areaPhaseConfUFullIndividual','areaPhaseConfLFullIndividual','areaFreqFullIndividual',...
      'areaCohFullInterpIndividual','areaCohConfUFullInterpIndividual','areaCohConfLFullInterpIndividual','areaPhaseFullInterpIndividual',...
      'areaPhaseConfUFullInterpIndividual','areaPhaseConfLFullInterpIndividual','areaFreqFullInterpIndividual', '-append');
  end
else
  load(filename);
end


%% FIRING RATE CORRELATIONS
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
options.saveFig = true; %#ok<*STRNU>
if drawFRCorrelations
  for iArea = 1:numel(iAreasOI)

    % Firing rate
    fr1 = areaFRIndividual{1}{iAreasOI(iArea)};
    fr2 = areaFRIndividual{2}{iAreasOI(iArea)};
    options.figTitle = ['Firing rate during wakefulness vs firing rate during anaesthesia in ' areas{iAreasOI(iArea)}];
    options.figName = ['SPIKING fr in ' areas{iAreasOI(iArea)}];
    options.xLabel = 'Unit firing rate during wakefulness (APs/s)';
    options.yLabel = 'Unit firing rate during anaesthesia (APs/s)';
    options.axesType = 'loglog';
    [~, rFR{iAreasOI(iArea)}, pvalFR{iAreasOI(iArea)}, modelFR{iAreasOI(iArea)}] = corrPlot(fr1, fr2, options);

    % Log firing rate
    fr1 = getLog(areaFRIndividual{1}{iAreasOI(iArea)});
    fr2 = getLog(areaFRIndividual{2}{iAreasOI(iArea)});
    options.figTitle = ['Log firing rate during wakefulness vs log firing rate during anaesthesia in ' areas{iAreasOI(iArea)}];
    options.figName = ['SPIKING frLog in ' areas{iAreasOI(iArea)}];
    options.xLabel = 'Unit log firing rate during wakefulness';
    options.yLabel = 'Unit log firing rate during anaesthesia';
    options.axesType = 'regular';
    [~, rFRLog{iAreasOI(iArea)}, pvalFRLog{iAreasOI(iArea)}, modelFRLog{iAreasOI(iArea)}] = corrPlot(fr1, fr2, options);
  end

  % Save variables
  save(filename, 'rFR','pvalFR','modelFR','rFRLog','pvalFRLog','modelFRLog', '-append');
end


%% PHASE CORRELATIONS
options = struct();
options.figFolder = mainFolder;
options.figSize = figSize;
options.figTitle = 'SPIKING';
options.xLabel = 'Wakefulness phase (rad)';
options.yLabel = 'Anaesthesia phase (rad)';
options.phaseLim = phaseLim;
options.individualGraphs = individualPhaseGraphs;
options.summaryGraphs = summaryPhaseGraphs;
options.iAreasOI = iAreasOI;
[rPhase, pvalPhase, nPhase] = halfPhaseCorr(areas, {'all'}, areaFreqFullInterpIndividual,...
  areaPhaseFullInterpIndividual(1), areaPhaseFullInterpIndividual(2), options);


%% COHERENCE CORRELATIONS
options = struct();
options.figFolder = mainFolder;
options.figSize = figSize;
options.figTitle = 'SPIKING';
options.xLabel = 'Wakefulness coherence';
options.yLabel = 'Anaesthesia coherence';
options.cohLim = [0 1];
options.individualGraphs = individualCohGraphs;
options.summaryGraphs = summaryCohGraphs;
options.iAreasOI = iAreasOI;
[rCoh, pvalCoh, nCoh] = halfCohCorr(areas, {'all'}, areaFreqFullInterpIndividual,...
  areaCohFullInterpIndividual(1), areaCohFullInterpIndividual(2), options);


%% SAVE THE RESULTS OF CORRELATION ANALYSES
save(filename, 'rPhase','pvalPhase','nPhase','rCoh','pvalCoh','nCoh', '-append');