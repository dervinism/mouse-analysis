% Run this script to produce unit phase and coherence frequency profiles
% and phase frequency histograms for laminar comparisons.
% 
% The following data files are produced:
% dataDir\clDir\unitsFolder\histosSubfolder\phaseHistos_units.mat or
%   dataDir\clDir\unitsFolder\histosSubfolder\phaseHistos_units_reverse.mat or
%   dataDir\clDir\unitsFolder\histosSubfolder\phaseHistos_units_quality.mat or
%   dataDir\clDir\unitsFolder\histosSubfolder\phaseHistos_units_reverse_quality.mat. All
%   of these files contain unit phase frequency profiles, phase frequency
%   histograms, and statistical test results.
% 
% Unit phase frequency profile figures are saved in dataDir\clDir\unitsFolder\phaseFrequencyProfilesSubfolder
%   or dataDir\clDir\qualityUnitsFolder\phaseFrequencyProfilesSubfolder.
% Individual unit phase frequency histogram figures and summary suplots are
%   saved in dataDir\clDir\unitsFolder\histosSubfolder or dataDir\clDir\qualityUnitsFolder\histosSubfolder.
% Phase frequency maps are saved in dataDir\clDir\unitsFolder\mapsSubfolder or
%   dataDir\clDir\qualityUnitsFolder\mapsSubfolder.

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
dataDir = [dataDir filesep includeRuns];

if strcmp(repository,'uol')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep clDir_uol];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep clDir_uol_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep clDir_uol_negative];
  end
  animals = animalsUOLOI;
  xLim = freqLimUOL;
elseif strcmp(repository,'allensdk')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep clDir_allensdk];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep clDir_allensdk_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep clDir_allensdk_negative];
  end
  animals = animalsAllensdk;
  conditions = {'awake'};
  xLim = freqLimAllen;
end
mainFolder = [rootFolder filesep MUAsFolder];

drawPhaseProfiles = true;
drawCohProfiles = true;
drawPhaseHistos = [false true true];
type = 'Superficial'; %'Deep'; %'Superficial'
minCond = 1;
maxCond = 2;

mainFolder = [mainFolder filesep lower(type)];


%% COMPUTE VARIABLES AND INITIALISE FIGURES NEEDED FOR DISPLAYING UNIT PHASE AND COHERENCE FREQUENCY PROFILES
if fullRun
  for animal = 1:numel(animals) % Loop through animals
    disp(['processing animal ' animals{animal}]);
    if exist('allData', 'var')
      dataStruct = allData.(animals{animal}).dataStruct;
    else
      load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
    end
    if ~isfield(dataStruct, 'seriesDataLaminar') ||...
        (strcmpi(subpop,'positive') && ~isfield(dataStruct, 'seriesDataLaminar_positive')) ||...
        (strcmpi(subpop,'negative') && ~isfield(dataStruct, 'seriesDataLaminar_negative'))
      disp(['Animal ' animals{animal} ' is missing laminar series data. Skipping to the next animal...'])
      continue
    end
    
    % Initialise variables
    fnsData = fieldnames(dataStruct.seriesDataLaminar);
    FOI = dataStruct.seriesData.(fnsData{1}).conf.FOI;
    srData = dataStruct.seriesData.(fnsData{1}).conf.samplingParams.srData;
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
      if strcmpi(subpop, 'positive')
        fnsDataPositive = fieldnames(dataStruct.seriesDataLaminar_positive);
        dbCountPositive = find(ismember(fnsDataPositive, fnsData{dbCount}));
        if isempty(dbCountPositive)
          continue
        end
        dbStructLaminar = dataStruct.seriesDataLaminar_positive.(fnsDataPositive{dbCountPositive});
      elseif strcmpi(subpop, 'negative')
        fnsDataNegative = fieldnames(dataStruct.seriesDataLaminar_negative);
        dbCountNegative = find(ismember(fnsDataNegative, fnsData{dbCount}));
        if isempty(dbCountNegative)
          continue
        end
        dbStructLaminar = dataStruct.seriesDataLaminar_negative.(fnsDataNegative{dbCountNegative});
      else
        dbStructLaminar = dataStruct.seriesDataLaminar.(fnsData{dbCount});
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
        [~, ~, areaName, ~, ~, area] = determineArea(seriesName);
      elseif strcmp(repository,'allensdk')
        [area, ~, areaName] = determineArea(seriesName);
      end
      if ~contains(areaName,'S1') && ~contains(areaName,'RSC') &&...
          ~contains(areaName,'V1') && ~contains(areaName,'V2') && ~contains(areaName,'VIS')
        continue
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
      units = units(qualityUnitInd);
      
      for iAreaPlusAll = area % Loop through the main and pooled areas
        if numel(conditions) > 1
          condLoop = [iCond numel(conditions)];
        else
          condLoop = iCond;
        end
        for iCondPlusAll = condLoop % Loop through the main and pooled conditions
          
          % Load and store phase and coherence data for units
          sh = 1;
          phaseCoh = dbStructLaminar.shankData.(['shank' num2str(sh)]).(['phaseCoh' type]);
          if isempty(phaseCoh)
            continue
          end
          for u = 1:numel(phaseCoh)
            if ~sum(phaseCoh{u}.unit == units)
              continue
            end
            
            freq = phaseCoh{u}.freq;
            
            % Get phase and coherence values
            if isfield(phaseCoh{u}, 'coh')
              coh = phaseCoh{u}.coh;
              cohConf = phaseCoh{u}.coh_conf;
              rateadjust_kappa = phaseCoh{u}.rateadjust_kappa;
              phase = phaseCoh{u}.phase;
              phaseConfU = phaseCoh{u}.phase_confU;
              phaseConfL = phaseCoh{u}.phase_confL;
%                 [phase, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase,...
%                   phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa);
                [~, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase,...
                  phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa);
              
              % Obtain and store phase and coherence values for FOI
              if sum(~isnan(phase)) && sum(~isnan(coh))
                [phaseFOI, cohFOI] = phaseCohFOI(FOI, freq, phase, coh, [cohConfU; cohConfL], rateadjust_kappa);
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
      end
    end
  end
  
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
  for iCond = 1:numel(conditions)
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
    'areaPhaseConfLFullInterpIndividual','areaFreqFullInterpIndividual');
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
  options.phaseLim = [-pi/4 pi/4];
  options.cutoffFreq = xLim(2);
  options.iAreasOI = iAreasOI;
  phaseFreqProfilePlotMeans(areas, conditions(minCond:min([numel(conditions) maxCond])),...
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
  cohFreqProfilePlotMeans(areas, conditions(minCond:min([numel(conditions) maxCond])),...
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
[phaseHistos, distributionStats] = phaseHistosPlotMaster(drawPhaseHistos, areas, conditions(minCond:min([numel(conditions) maxCond])),...
  areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual, edges, options);


%% STATS ON MEAN PHASE FREQUENCY PROFILES
if ~strcmp(repository, 'allensdk')
  [fPEst, fWTest, strPMethod, pEst, U2, pObs, U2Obs] = phaseComparisonStats(areas, areaFreqFullInterpIndividual, FOI, areaPhaseFullInterpIndividual, [areasCritical; areas2compareCritical]);
else
  fPEst = []; fWTest = []; strPMethod = []; pEst = []; U2 = []; pObs = []; U2Obs = [];
end
save(filename, 'phaseHistos','distributionStats','fPEst','fWTest','strPMethod','pEst','U2','pObs','U2Obs', '-append');