% Run this script to produce unit phase and coherence frequency profiles
% and phase frequency histograms for local comparisons.
% 
% The following data files are produced:
% dataDir\laDir\unitsFolder\histosSubfolder\phaseHistos_units.mat or
%   dataDir\laDir\unitsFolder\histosSubfolder\phaseHistos_units_reverse.mat or
%   dataDir\laDir\unitsFolder\histosSubfolder\phaseHistos_units_quality.mat or
%   dataDir\laDir\unitsFolder\histosSubfolder\phaseHistos_units_reverse_quality.mat. All
%   of these files contain unit phase frequency profiles, phase frequency
%   histograms, and statistical test results.
% 
% Unit phase frequency profile figures are saved in dataDir\laDir\unitsFolder\phaseFrequencyProfilesSubfolder
%   or dataDir\laDir\qualityUnitsFolder\phaseFrequencyProfilesSubfolder.
% Individual unit phase frequency histogram figures and summary suplots are
%   saved in dataDir\laDir\unitsFolder\histosSubfolder or dataDir\laDir\qualityUnitsFolder\histosSubfolder.
% Phase frequency maps are saved in dataDir\laDir\unitsFolder\mapsSubfolder or
%   dataDir\laDir\qualityUnitsFolder\mapsSubfolder.

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
if strcmpi(includeRuns, 'noRun')
  fRef = 0.3;
else
  fRef = 0.03;
end

dataDir = [dataDir filesep includeRuns];
if strcmp(repository,'all')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep laDir];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep laDir_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep laDir_negative];
  end
  animals = animalsOI;
  xLim = freqLimUOL;
elseif strcmp(repository,'uol')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep laDir_uol];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep laDir_uol_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep laDir_uol_negative];
  end
  animals = animalsUOLOI;
  xLim = freqLimUOL;
elseif strcmp(repository,'allensdk')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep laDir_allensdk];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep laDir_allensdk_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep laDir_allensdk_negative];
  end
  animals = animalsAllensdk;
  conditions = {'awake'};
  xLim = freqLimAllen;
end
if qualityCheck
  mainFolder = [rootFolder filesep qualityUnitsFolder];
else
  mainFolder = [rootFolder filesep unitsFolder];
end

drawPhaseProfiles = true;
drawCohProfiles = true;
drawPhaseHistos = [false true true];
drawACGDs = false;
drawBetas = false;


%% COMPUTE VARIABLES AND INITIALISE FIGURES NEEDED FOR DISPLAYING UNIT PHASE AND COHERENCE FREQUENCY PROFILES
if fullRun
  for animal = 1:numel(animals) % Loop through animals
    disp(['processing animal ' animals{animal}]);
    if exist('allData', 'var')
      dataStruct = allData.(animals{animal}).dataStruct;
    else
      load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
    end
    if strcmp(subpop, 'all')
      if ~isfield(dataStruct, 'seriesData') && (~isfield(dataStruct, 'seriesData') || isempty(dataStruct.seriesData))
        disp(['Animal ' animals{animal} ' is missing series data. Skipping to the next animal...'])
        continue
      end
    elseif strcmp(subpop, 'positive') && (~isfield(dataStruct, 'seriesData_positive') || isempty(dataStruct.seriesData_positive))
      if ~isfield(dataStruct, 'seriesData_positive')
        disp(['Animal ' animals{animal} ' is missing postive series data. Skipping to the next animal...'])
        continue
      end
    elseif strcmp(subpop, 'negative') && (~isfield(dataStruct, 'seriesData_negative') || isempty(dataStruct.seriesData_negative))
      if ~isfield(dataStruct, 'seriesData_negative')
        disp(['Animal ' animals{animal} ' is missing negative series data. Skipping to the next animal...'])
        continue
      end
    end
    
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
    if animal == 1 || ~exist('areaCohFOIindividual', 'var')
      areaCohFOIindividual = {};
      areaPhaseFOIindividual = {};
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
      areaFreqFullInterpIndividual = {};
      areaAutoCorrsIndividual = {};
      for iCond = 1:numel(conditions)
        areaCohFOIindividualCond = {};
        areaPhaseFOIindividualCond = {};
        areaCohFullIndividualCond = {};
        areaCohConfUFullIndividualCond = {};
        areaCohConfLFullIndividualCond = {};
        areaPhaseFullIndividualCond = {};
        areaPhaseConfUFullIndividualCond = {};
        areaPhaseConfLFullIndividualCond = {};
        areaFreqFullIndividualCond = {};
        areaCohFullInterpIndividualCond = {};
        areaCohConfUFullInterpIndividualCond = {};
        areaCohConfLFullInterpIndividualCond = {};
        areaPhaseFullInterpIndividualCond = {};
        areaPhaseConfUFullInterpIndividualCond = {};
        areaPhaseConfLFullInterpIndividualCond = {};
        areaFreqFullInterpIndividualCond = {};
        areaAutoCorrsIndividualCond = {};
        for iArea = 1:numel(areas)
          areaCohFOIindividualCond{iArea} = []; %#ok<*SAGROW>
          areaPhaseFOIindividualCond{iArea} = [];
          areaCohFullIndividualCond{iArea} = {};
          areaCohConfUFullIndividualCond{iArea} = {};
          areaCohConfLFullIndividualCond{iArea} = {};
          areaPhaseFullIndividualCond{iArea} = {};
          areaPhaseConfUFullIndividualCond{iArea} = {};
          areaPhaseConfLFullIndividualCond{iArea} = {};
          areaFreqFullIndividualCond{iArea} = {};
          areaCohFullInterpIndividualCond{iArea} = [];
          areaCohConfUFullInterpIndividualCond{iArea} = [];
          areaCohConfLFullInterpIndividualCond{iArea} = [];
          areaPhaseFullInterpIndividualCond{iArea} = [];
          areaPhaseConfUFullInterpIndividualCond{iArea} = [];
          areaPhaseConfLFullInterpIndividualCond{iArea} = [];
          areaFreqFullInterpIndividualCond{iArea} = [];
          areaAutoCorrsIndividualCond{iArea} = [];
        end
        areaCohFOIindividual{iCond} = areaCohFOIindividualCond;
        areaPhaseFOIindividual{iCond} = areaPhaseFOIindividualCond;
        areaCohFullIndividual{iCond} = areaCohFullIndividualCond;
        areaCohConfUFullIndividual{iCond} = areaCohConfUFullIndividualCond;
        areaCohConfLFullIndividual{iCond} = areaCohConfLFullIndividualCond;
        areaPhaseFullIndividual{iCond} = areaPhaseFullIndividualCond;
        areaPhaseConfUFullIndividual{iCond} = areaPhaseConfUFullIndividualCond;
        areaPhaseConfLFullIndividual{iCond} = areaPhaseConfLFullIndividualCond;
        areaFreqFullIndividual{iCond} = areaFreqFullIndividualCond;
        areaCohFullInterpIndividual{iCond} = areaCohFullInterpIndividualCond;
        areaCohConfUFullInterpIndividual{iCond} = areaCohConfUFullInterpIndividualCond;
        areaCohConfLFullInterpIndividual{iCond} = areaCohConfLFullInterpIndividualCond;
        areaPhaseFullInterpIndividual{iCond} = areaPhaseFullInterpIndividualCond;
        areaPhaseConfUFullInterpIndividual{iCond} = areaPhaseConfUFullInterpIndividualCond;
        areaPhaseConfLFullInterpIndividual{iCond} = areaPhaseConfLFullInterpIndividualCond;
        areaFreqFullInterpIndividual{iCond} = areaFreqFullInterpIndividualCond;
        areaAutoCorrsIndividual{iCond} = areaAutoCorrsIndividualCond;
      end
    end
    
    for dbCount = 1:numel(fnsData) % Loop through db entries
      dbStruct = dataStruct.seriesData.(fnsData{dbCount});
      if isempty(dbStruct)
        continue
      end
            
      % Determine if series phase and coherence data exist
      if exist('seriesName', 'var')
        prevRec = seriesName(1:min([14 numel(seriesName)]));
      else
        prevRec = '';
      end
      seriesName = seriesFromEntry(fnsData{dbCount});
      recording = seriesName(1:min([14 numel(seriesName)]));
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
      
      % Determine recording condition (i.e., awake or anaesthesia)
      [breakClause, iCond] = series2condition(awake, anaesthesia, seriesName);
      if breakClause
        continue
      end
      
      % Determine mode boundaries
      if strcmpi(repository, 'uol')
        modeBoundaries = phaseHistoBinCentres(modesUOL{area(1)});
      elseif strcmpi(repository, 'allensdk')
        modeBoundaries = phaseHistoBinCentres(modesAllensdk{area(1)});
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
          modeBoundaries = recentrePhase(modeBoundaries, modeBoundaries(1));
          phase = recentrePhase(phase, modeBoundaries(1));
          correlatedInd(phase > modeBoundaries(end) & phase <= modeBoundaries(2)) = true;
        elseif iCond == 2
          correlatedInd = find(rSpearman >= 0);
        end
      elseif strcmp(subpop, 'negative')
        if iCond == 1
          correlatedInd = false(size(phase));
          modeBoundaries = recentrePhase(modeBoundaries, modeBoundaries(3));
          phase = recentrePhase(phase, modeBoundaries(3));
          correlatedInd(phase > modeBoundaries(2) & phase <= modeBoundaries(end)) = true;
        elseif iCond == 2
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
      
      for iAreaPlusAll = area % Loop through the main and pooled areas
        if numel(conditions) > 1
          condLoop = [iCond numel(conditions)];
        else
          condLoop = iCond;
        end
        for iCondPlusAll = condLoop % Loop through the main and pooled conditions
          
          % Load and store phase and coherence data for units
          spk = [];
          shankIDs = fieldnames(dbStruct.shankData);
          for sh = 1:numel(shankIDs)
            shankStruct = dbStruct.shankData.(shankIDs{sh});
            if isempty(spk)
              spk = shankStruct.spk;
            else
              spk = concatenateMat(spk, shankStruct.spk);
            end
            for u = 1:numel(units)
              iU = find(shankStruct.units == units(u));
              if isempty(iU)
                continue
              end
              
              freq = shankStruct.phaseCoh{iU}.freq;  
              
              % Get phase and coherence values
              if isfield(shankStruct.phaseCoh{iU}, 'coh')
                coh = shankStruct.phaseCoh{iU}.coh;
                cohConf = shankStruct.phaseCoh{iU}.coh_conf;
                rateadjust_kappa = shankStruct.phaseCoh{iU}.rateadjust_kappa;
                phase = shankStruct.phaseCoh{iU}.phase;
                phaseConfU = shankStruct.phaseCoh{iU}.phase_confU;
                phaseConfL = shankStruct.phaseCoh{iU}.phase_confL;
                [phase, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase,...
                  phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa);
%                 [~, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase,...
%                   phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa);
                
                % Obtain and store phase and coherence values for FOI
                if sum(~isnan(phase)) && sum(~isnan(coh))
                  [phaseFOI, cohFOI] = phaseCohFOI(FOI, freq, phase, coh, [cohConfU; cohConfL]);
                else
                  phaseFOI = NaN(size(FOI));
                  cohFOI = NaN(size(FOI));
                end
                
                % Store values
                areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll}; cohFOI];
                areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll}; phaseFOI];
                areaCohFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCohFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = coh;
                areaCohConfUFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCohConfUFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = cohConfU;
                areaCohConfLFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCohConfLFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = cohConfL;
                areaPhaseFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = phase;
                areaPhaseConfUFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseConfUFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = phaseConfU;
                areaPhaseConfLFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseConfLFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = phaseConfL;
              else
                areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll}; NaN(size(FOI))];
                areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll}; NaN(size(FOI))];
                areaCohFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCohFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = NaN(size(freq));
                areaCohConfUFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCohConfUFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = NaN(size(freq));
                areaCohConfLFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCohConfLFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = NaN(size(freq));
                areaPhaseFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = NaN(size(freq));
                areaPhaseConfUFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseConfUFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = NaN(size(freq));
                areaPhaseConfLFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseConfLFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = NaN(size(freq));
              end
              areaFreqFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaFreqFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = freq;
            end
          end
          
          % Calculate and store autocorrelations
          for u = 1:numel(qualityUnitInd)
            [xc,lags] = xcorr(full(spk(qualityUnitInd(u),:)),'coeff');
            %stem(lags,xc)
            inds = ceil(numel(lags)/2)+1:ceil(numel(lags)/2)+acgdPeriod*srData;
            areaAutoCorrsIndividual{iCondPlusAll}{iAreaPlusAll} = [areaAutoCorrsIndividual{iCondPlusAll}{iAreaPlusAll}; xc(inds)];
          end
        end
      end
    end
  end
  lags = (1:acgdPeriod*srData).*(1/srData);
  
  % Interpolate and store full phases and coherences
  freqCombined = FOI;
  for iCond = 1:numel(conditions)
    for iArea = 1:numel(areas)
      nRec = numel(areaFreqFullIndividual{iCond}{iArea});
      for iRec = 1:nRec
        freqCombined = unique([freqCombined areaFreqFullIndividual{iCond}{iArea}{iRec}]);
      end
    end
  end
  freqCombined = unique(freqCombined(~isnan(freqCombined)));
  for iCond = 1:min([2 numel(conditions)])
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


%% SAVE OR LOAD THE DATA
if qualityCheck
  filename = [mainFolder filesep 'globalUnits_quality.mat'];
else
  filename = [mainFolder filesep 'globalUnits.mat'];
end
if fullRun
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  save(filename, 'conditions','areas','FOI','srData','areaCohFOIindividual','areaPhaseFOIindividual',...
    'areaCohFullIndividual','areaCohConfUFullIndividual','areaCohConfLFullIndividual',...
    'areaPhaseFullIndividual','areaPhaseConfUFullIndividual','areaPhaseConfLFullIndividual',...
    'areaFreqFullIndividual','areaCohFullInterpIndividual','areaCohConfUFullInterpIndividual',...
    'areaCohConfLFullInterpIndividual','areaPhaseFullInterpIndividual','areaPhaseConfUFullInterpIndividual',...
    'areaPhaseConfLFullInterpIndividual','areaFreqFullInterpIndividual','areaAutoCorrsIndividual','lags');
else
  load(filename); %#ok<*UNRCH>
end


%% GENERATE PHASE MEAN FREQUENCY PROFILE FIGURES AND SAVE THEM
if drawPhaseProfiles
  figFileName = 'Means_only__%s';
  options = struct();
  options.figTitle = 'Unit mean phase comparisons: %s';
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.phaseFrequencyProfilesSubfolder = phaseFrequencyProfilesSubfolder;
  options.freqLim = xLim;
  options.phaseLim = [-pi pi] + [-pi/4 pi/4];
  options.cutoffFreq = xLim(2);
  options.iAreasOI = iAreasOI;
  phaseFreqProfilePlotMeans(areas, conditions(1:min([numel(conditions) 2])),...
    areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual, [], figFileName, options);
end


%% GENERATE COHERENCE MEAN FREQUENCY PROFILE FIGURES AND SAVE THEM
if drawCohProfiles
  figFileName = 'Means_only__%s';
  options = struct();
  options.figTitle = 'Unit mean coherence comparisons: %s';
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.cohFrequencyProfilesSubfolder = coherenceFrequencyProfilesSubfolder;
  options.freqLim = xLim;
  options.cohLim = [0 1];
  options.iAreasOI = iAreasOI;
  cohFreqProfilePlotMeans(areas, conditions(1:min([numel(conditions) 2])),...
    areaFreqFullInterpIndividual, areaCohFullInterpIndividual, [], figFileName, options);
end


%% DRAW PHASE FREQUENCY HISTOGRAMS AND MAPS
options = struct();
options.mainFolder = mainFolder;
options.histosSubfolder = histosSubfolder;
options.mapsSubfolder = mapsSubfolder;
options.figSize = figSize;
options.figTitle = 'SPIKING';
options.freqLim = [xLim(1)+0.002 xLim(2)];
options.phaseLimHisto = phaseLim;
options.phaseLimMap = phaseLim;
options.xLabelHist = '# units';
options.iAreasOI = iAreasOI;
[phaseHistos, distributionStats] = phaseHistosPlotMaster(drawPhaseHistos, areas, conditions(1:min([numel(conditions) 2])),...
  areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual, edges, options);


%% STATS ON MEAN PHASE FREQUENCY PROFILES
if ~strcmp(repository, 'allensdk')
  [fPEst, fWTest, strPMethod, pEst, U2, pObs, U2Obs] = phaseComparisonStats(areas, areaFreqFullInterpIndividual, FOI, areaPhaseFullInterpIndividual, [areasCritical; areas2compareCritical]);
else
  fPEst = []; fWTest = []; strPMethod = []; pEst = []; U2 = []; pObs = []; U2Obs = [];
end
save(filename, 'phaseHistos','distributionStats','fPEst','fWTest','strPMethod','pEst','U2','pObs','U2Obs', '-append');


%% FIT AUTOCORRELATION DECAYS AND DISPLAY FIGURES
if drawACGDs
  if strcmp(repository, 'uol')
    areasOI = {'VB','S1','RSC','CA'};
    areasOIExtra = {'Th','VB','Po','lS1','lRSC','CA','DG'};
  elseif strcmp(repository, 'allensdk')
    areasOI = {'VB','LGN','V1','CA'};
    areasOIExtra = {'Th','LGN','LP','V1','V2','CA','DG'};
  end
  
  if ~exist([mainFolder filesep acgdSubfolder], 'file')
    mkdir([mainFolder filesep acgdSubfolder]);
  end
  binSize = [srData, 200, 100, 50, 20, 10, 8, 5];
  
  % Unit mean analyses
  for iSize = 1:numel(binSize)
    [fH1{iSize}, acgdMeanFits{iSize}] = acgdPlotUnitMeans(lags, areaAutoCorrsIndividual, areasOI, 1, binSize(iSize), acgdPeriod);
    set(fH1{iSize}, 'Name',['Mean unit auto-correlation decay across brain areas: ' num2str(1/binSize(iSize)) 's bin size']);
    filenameFig = [mainFolder filesep acgdSubfolder filesep 'autocorrelationDecaysAcrossAreas' num2str(1/binSize(iSize)) '.fig'];
    savefig(fH1{iSize}, filenameFig, 'compact');
    print(fH1{iSize}, [filenameFig '.png'],'-dpng','-r300');
  end
  close all
  
  for iSize = 1:numel(binSize)
    [fH2{iSize}, acgdMeanFitsExtra{iSize}] = acgdPlotUnitMeans(lags, areaAutoCorrsIndividual, areasOIExtra, 1, binSize(iSize), acgdPeriod);
    set(fH2{iSize}, 'Name',['Mean unit auto-correlation decay across brain areas: ' num2str(1/binSize(iSize)) 's bin size']);
    filenameFig = [mainFolder filesep acgdSubfolder filesep 'autocorrelationDecaysAcrossAreas' num2str(1/binSize(iSize)) 'Extra.fig'];
    savefig(fH2{iSize}, filenameFig, 'compact');
    print(fH2{iSize}, [filenameFig '.png'],'-dpng','-r300');
  end
  close all
  
  % Unit analyses
  for iSize = 1:numel(binSize)
    [fH3{iSize}, acgdUnitFits{iSize}] = acgdPlotUnits(lags, areaAutoCorrsIndividual, areasOI, 1, binSize(iSize), acgdPeriod);
    set(fH3{iSize}(1), 'Name',['Unit auto-correlation effective decay time constant across brain areas: ' num2str(1/binSize(iSize)) 's bin size']);
    filenameFig = [mainFolder filesep acgdSubfolder filesep 'autocorrelationDecayTauEffAcrossAreas' num2str(1/binSize(iSize)) '.fig'];
    savefig(fH3{iSize}(1), filenameFig, 'compact');
    print(fH3{iSize}(1), [filenameFig '.png'],'-dpng','-r300');
    set(fH3{iSize}(2), 'Name',['Unit auto-correlation fitted decay time constant across brain areas: ' num2str(1/binSize(iSize)) 's bin size']);
    filenameFig = [mainFolder filesep acgdSubfolder filesep 'autocorrelationDecayTauAcrossAreas' num2str(1/binSize(iSize)) '.fig'];
    savefig(fH3{iSize}(2), filenameFig, 'compact');
    print(fH3{iSize}(2), [filenameFig '.png'],'-dpng','-r300');
  end
  close all
  
  for iSize = 1:numel(binSize)
    [fH4{iSize}, acgdUnitFitsExtra{iSize}] = acgdPlotUnits(lags, areaAutoCorrsIndividual, areasOIExtra, 1, binSize(iSize), acgdPeriod);
    set(fH4{iSize}(1), 'Name',['Unit auto-correlation effective decay time constant across brain areas: ' num2str(1/binSize(iSize)) 's bin size']);
    filenameFig = [mainFolder filesep acgdSubfolder filesep 'autocorrelationDecayTauEffAcrossAreas' num2str(1/binSize(iSize)) '.fig'];
    savefig(fH4{iSize}(1), filenameFig, 'compact');
    print(fH4{iSize}(1), [filenameFig '.png'],'-dpng','-r300');
    set(fH4{iSize}(2), 'Name',['Unit auto-correlation fitted decay time constant across brain areas: ' num2str(1/binSize(iSize)) 's bin size']);
    filenameFig = [mainFolder filesep acgdSubfolder filesep 'autocorrelationDecayTauAcrossAreas' num2str(1/binSize(iSize)) '.fig'];
    savefig(fH4{iSize}(2), filenameFig, 'compact');
    print(fH4{iSize}(2), [filenameFig '.png'],'-dpng','-r300');
  end
  close all
  
  % Save variables
  save(filename, 'acgdMeanFits','acgdMeanFitsExtra','acgdUnitFits','acgdUnitFitsExtra','binSize', '-append');
end


%% CORRELATE BETA EXPONENTS AND AUTOCORRELATION DECAY TIME CONSTANTS
if drawBetas
  if strcmp(repository, 'uol')
    areasOI = {'Th','VB','Po','lS1','lRSC','CA','DG'};
  elseif strcmp(repository, 'allensdk')
    areasOI = {'Th','LGN','LP','V1','V2','CA','DG'};
  end
  
  if qualityCheck
    betasData = load([mainFolder filesep betaCorrelationsSubfolder filesep 'globalUnitsBeta_quality.mat']);
  else
    betasData = load([mainFolder filesep betaCorrelationsSubfolder filesep 'globalUnitsBeta.mat']);
  end
  binSize = [srData, 200, 100, 50, 20, 10, 8, 5];
  
  % Effective taus
  options.tauType = 'effective';
  options.figFolder = [mainFolder filesep betaCorrelationsSubfolder];
  options.yScale = 'regular';
  options.yLim = [0 1];
  options.figSize = figSize;
  [rBetaTauEff, pvalBetaTauEff] = betaTauCorr(betasData.areaBetaIndividual, acgdUnitFitsExtra, binSize, conditions, areasOI, options);
  
  options.yScale = 'log';
  options.yLim = [0 inf];
  betaTauCorr(betasData.areaBetaIndividual, acgdUnitFitsExtra, binSize, conditions, areasOI, options);
  
  % Fitted taus
  options.tauType = 'fitted';
  options.figFolder = [mainFolder filesep betaCorrelationsSubfolder];
  options.yScale = 'regular';
  options.yLim = [0 3];
  options.figSize = figSize;
  [rBetaTauFit, pvalBetaTauFit] = betaTauCorr(betasData.areaBetaIndividual, acgdUnitFitsExtra, binSize, conditions, areasOI, options);
  
  options.yScale = 'log';
  options.yLim = [0 inf];
  betaTauCorr(betasData.areaBetaIndividual, acgdUnitFitsExtra, binSize, conditions, areasOI, options);
  
  % Save variables
  save(filename, 'rBetaTauEff','pvalBetaTauEff','rBetaTauFit','pvalBetaTauFit','binSize', '-append');
end



%% Local functions
function [rBetaTau, pvalBetaTau] = betaTauCorr(areaBetaIndividual, acgdUnitFits, binSize, conditions, areas, options) %#ok<*DEFNU>

for iSize = 1:numel(binSize)
  for iCond = 1:numel(conditions(1))
    for iArea = 1:numel(areas)
      areaCode = determineArea(areas{iArea});
      betas = areaBetaIndividual{iCond}{areaCode};
      if strcmpi(options.tauType, 'effective')
        taus = acgdUnitFits{iSize}.effectiveTau{iArea};
      elseif strcmpi(options.tauType, 'fitted')
        taus = acgdUnitFits{iSize}.fittedTau{iArea};
      end
      if ~isempty(betas) && ~isempty(taus) && sum(~isnan(betas)) && sum(~isnan(taus))
        
        % Correlations and individual graphs
        inds = ~isnan(betas) & ~isnan(taus');
        [rBetaTau{iCond}{areaCode}{iSize}, pvalBetaTau{iCond}{areaCode}{iSize}] = corrMulti(torow(betas(inds)), torow(taus(inds)), 'Spearman');
        nBetaTau{iCond}{areaCode}{iSize}.significant = sum(inds);
        nBetaTau{iCond}{areaCode}{iSize}.total = numel(inds);
        
        figBetaTauEff = figure;
        plot(betas(inds), taus(inds), '.', 'MarkerSize',5)
        if strcmpi(options.tauType, 'effective') && strcmpi(options.yScale, 'log')
          set(gca,'Yscale','log')
          figTitle = ['Unit PSD beta and effective \tau correlations for ' areas{iArea} ' ' conditions{iCond} ' ' num2str(1/binSize(iSize))...
            ' r=' num2str(rBetaTau{iCond}{areaCode}{iSize}) ' p=' num2str(pvalBetaTau{iCond}{areaCode}{iSize})...
            ' n=' num2str(nBetaTau{iCond}{areaCode}{iSize}.significant) '/' num2str(nBetaTau{iCond}{areaCode}{iSize}.total)];
          title(figTitle);
          figName = ['PSDbetaVsTauEffLog_in_' areas{iArea} '_during_' conditions{iCond} '_binSize_' num2str(1/binSize(iSize))];
        elseif strcmpi(options.tauType, 'effective') && ~strcmpi(options.yScale, 'log')
          figTitle = ['Unit PSD beta and effective \tau correlations for ' areas{iArea} ' ' conditions{iCond} ' ' num2str(1/binSize(iSize))...
            ' r=' num2str(rBetaTau{iCond}{areaCode}{iSize}) ' p=' num2str(pvalBetaTau{iCond}{areaCode}{iSize})...
            ' n=' num2str(nBetaTau{iCond}{areaCode}{iSize}.significant) '/' num2str(nBetaTau{iCond}{areaCode}{iSize}.total)];
          title(figTitle);
          figName = ['PSDbetaVsTauEff_in_' areas{iArea} '_during_' conditions{iCond} '_binSize_' num2str(1/binSize(iSize))];
        elseif strcmpi(options.tauType, 'fitted') && strcmpi(options.yScale, 'log')
          set(gca,'Yscale','log')
          figTitle = ['Unit PSD beta and fitted \tau correlations for ' areas{iArea} ' ' conditions{iCond} ' ' num2str(1/binSize(iSize))...
            ' r=' num2str(rBetaTau{iCond}{areaCode}{iSize}) ' p=' num2str(pvalBetaTau{iCond}{areaCode}{iSize})...
            ' n=' num2str(nBetaTau{iCond}{areaCode}{iSize}.significant) '/' num2str(nBetaTau{iCond}{areaCode}{iSize}.total)];
          title(figTitle);
          figName = ['PSDbetaVsTauFitLog_in_' areas{iArea} '_during_' conditions{iCond} '_binSize_' num2str(1/binSize(iSize))];
        elseif strcmpi(options.tauType, 'fitted') && ~strcmpi(options.yScale, 'log')
          figTitle = ['Unit PSD beta and fitted \tau correlations for ' areas{iArea} ' ' conditions{iCond} ' ' num2str(1/binSize(iSize))...
            ' r=' num2str(rBetaTau{iCond}{areaCode}{iSize}) ' p=' num2str(pvalBetaTau{iCond}{areaCode}{iSize})...
            ' n=' num2str(nBetaTau{iCond}{areaCode}{iSize}.significant) '/' num2str(nBetaTau{iCond}{areaCode}{iSize}.total)];
          title(figTitle);
          figName = ['PSDbetaVsTauFit_in_' areas{iArea} '_during_' conditions{iCond} '_binSize_' num2str(1/binSize(iSize))];
        end
        set(gcf, 'Name',figName);
        figName = [options.figFolder filesep figName];
        figName = strrep(figName, ' ', '_');
        figName = strrep(figName, '.', 'p');
        
        xLim = xlim;
        xAxisLength = xLim(2)-xLim(1);
        yLim = ylim;
        yLim(1) = max([options.yLim(1) yLim(1)]);
        yLim(2) = min([options.yLim(2) yLim(2)]);
        ylim(yLim);
        yAxisLength = yLim(2)-yLim(1);
        text(xLim(2)-xAxisLength*0.35, yLim(1)+yAxisLength*0.2, ['r=' num2str(rBetaTau{iCond}{areaCode}{iSize})], 'FontSize',16);
        text(xLim(2)-xAxisLength*0.35, yLim(1)+yAxisLength*0.1, ['p=' num2str(pvalBetaTau{iCond}{areaCode}{iSize})], 'FontSize',16);
        text(xLim(1)+xAxisLength*0.035, yLim(2)-yAxisLength*0.05,...
          ['n=' num2str(nBetaTau{iCond}{areaCode}{iSize}.significant) '/' num2str(nBetaTau{iCond}{areaCode}{iSize}.total)], 'FontSize',16);
        
        if strcmpi(options.tauType, 'effective')
          ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
            'on', 'k', 'Unit PSD exponent', xLim, [0 0.5 1],...
            'on', 'k', '\tau_{effective} (s)', yLim, yticks);
        else
          ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
            'on', 'k', 'Unit PSD exponent', xLim, [0 0.5 1],...
            'on', 'k', '\tau_{fitted} (s)', yLim, yticks);
        end
        
        label = [2 1.6];
        margin = [0.3 0.55];
        width = 1*options.figSize-label(1)-margin(1);
        height = 1*options.figSize-label(2)-margin(2);
        paperSize = resizeFig(figBetaTauEff, ax1, width, height, label, margin, 0);
        exportFig(figBetaTauEff, [figName '.png'],'-dpng','-r300', paperSize);
        hgsave(figBetaTauEff, [figName '.fig']);
        close(figBetaTauEff);
      end
    end
  end
end
end