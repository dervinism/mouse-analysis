% Run this script to produce MUA phase and coherence frequency profiles
% and phase frequency histograms for cross-area comparisons.
%
% The following data files are produced:
% dataDir\caDir\MUAsFolder\globalMUAs_ca.mat or
%   dataDir\caDir\MUAsFolder\globalMUAs_ca_reverse.mat. All
%   of these files contain MUA phase frequency profiles, phase frequency
%   histograms, and statistical test results.
%
% MUA phase frequency profile figures are saved in
%   dataDir\caDir\MUAsFolder\phaseFrequencyProfilesSubfolder.
% Individual MUA phase frequency histogram figures and summary suplots are
%   saved in dataDir\caDir\MUAsFolder\histosSubfolder.
% Phase frequency maps are saved in dataDir\caDir\MUAsFolder\mapsSubfolder.

clearvars -except repository subpop reverse allData fullRun includeRuns


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
if ~exist('reverse', 'var')
  reverse = false;
end
dataDir = [dataDir filesep includeRuns];

if strcmp(repository,'uol')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep caDir_uol];
    rootFolder_la = [dataDir filesep laDir_uol];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep caDir_uol_positive];
    rootFolder_la = [dataDir filesep laDir_uol_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep caDir_uol_negative];
    rootFolder_la = [dataDir filesep laDir_uol_negative];
  end
  animals = animalsUOLOI;
  xLim = freqLimUOL;
elseif strcmp(repository,'allensdk')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep caDir_allensdk];
    rootFolder_la = [dataDir filesep laDir_allensdk];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep caDir_allensdk_positive];
    rootFolder_la = [dataDir filesep laDir_allensdk_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep caDir_allensdk_negative];
    rootFolder_la = [dataDir filesep laDir_allensdk_negative];
  end
  animals = animalsAllensdk;
  conditions = {'awake'};
  xLim = freqLimAllen;
end
mainFolder = [rootFolder filesep MUAsFolder];
mainFolder_la = [rootFolder_la filesep MUAsFolder];
areas = areas2compare;

drawPhaseProfiles = false;
drawCohProfiles = false;
drawPhaseHistos = [false false false];
doStats = false;
doStats_ca = false;


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
      if ~isfield(dataStruct, 'seriesData_ca_muas') || isempty(dataStruct.seriesData_ca_muas)
        disp(['Animal ' animals{animal} ' is missing series comparison data. Skipping to the next animal...'])
        continue
      end
      fnsData_ca = fieldnames(dataStruct.seriesData_ca_muas);
    elseif strcmp(subpop, 'positive')
      if ~isfield(dataStruct, 'seriesData_ca_muas_positive') || isempty(dataStruct.seriesData_ca_muas_positive)
        disp(['Animal ' animals{animal} ' is missing series comparison data. Skipping to the next animal...'])
        continue
      end
      fnsData_ca = fieldnames(dataStruct.seriesData_ca_muas_positive);
    elseif strcmp(subpop, 'negative')
      if ~isfield(dataStruct, 'seriesData_ca_muas_negative') || isempty(dataStruct.seriesData_ca_muas_negative)
        disp(['Animal ' animals{animal} ' is missing series comparison data. Skipping to the next animal...'])
        continue
      end
      fnsData_ca = fieldnames(dataStruct.seriesData_ca_muas_negative);
    end
    
    % Initialise variables
    if strcmp(subpop, 'all')
      fnsData = fieldnames(dataStruct.seriesData);
      try
        FOI = dataStruct.seriesData.(fnsData{1}).conf.FOI;
      catch
        FOI = dataStruct.seriesData.(fnsData{end}).conf.FOI;
      end
    elseif strcmp(subpop, 'positive')
      fnsData = fieldnames(dataStruct.seriesData_positive);
      try
        FOI = dataStruct.seriesData_positive.(fnsData{1}).conf.FOI;
      catch
        FOI = dataStruct.seriesData_positive.(fnsData{end}).conf.FOI;
      end
    elseif strcmp(subpop, 'negative')
      fnsData = fieldnames(dataStruct.seriesData_negative);
      try
        FOI = dataStruct.seriesData_negative.(fnsData{1}).conf.FOI;
      catch
        FOI = dataStruct.seriesData_negative.(fnsData{end}).conf.FOI;
      end
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
      end
    end
    
    for dbCount = 1:numel(fnsData_ca) % Loop through db entries
      if strcmp(subpop, 'all')
        dbStruct_ca = dataStruct.seriesData_ca_muas.(fnsData_ca{dbCount});
      elseif strcmp(subpop, 'positive')
        dbStruct_ca = dataStruct.seriesData_ca_muas_positive.(fnsData_ca{dbCount});
      elseif strcmp(subpop, 'negative')
        dbStruct_ca = dataStruct.seriesData_ca_muas_negative.(fnsData_ca{dbCount});
      end
      if isempty(dbStruct_ca)
        continue
      end
      
      % Determine if series phase and coherence data exist
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
      
      units = dbStruct_ca.units;
      
      for iComp = 1:numel(comp) % Loop through non-grouped and grouped area comparisons
        area = comp(iComp);
        
        if numel(conditions) > 1
          condLoop = [iCond numel(conditions)];
        else
          condLoop = iCond;
        end
        for iCondPlusAll = condLoop % Loop through the main and pooled conditions
          
          % Load and store phase and coherence data for units
          for u = 1:numel(units)
            % Get phase and coherence values
            coh = dbStruct_ca.phaseCoh{u}.coh;
            cohConf = dbStruct_ca.phaseCoh{u}.coh_conf;
            rateadjust_kappa = dbStruct_ca.phaseCoh{u}.rateadjust_kappa;
            phase = dbStruct_ca.phaseCoh{u}.phase;
            phaseConfU = dbStruct_ca.phaseCoh{u}.phase_confU;
            phaseConfL = dbStruct_ca.phaseCoh{u}.phase_confL;
            %[phase, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase,...
            %  phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa);
            [~, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase,...
              phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa);
            
            % Obtain and store phase and coherence values for FOI
            freq = dbStruct_ca.phaseCoh{u}.freq;
            if sum(~isnan(phase)) && sum(~isnan(coh))
              [phaseFOI, cohFOI] = phaseCohFOI(FOI, freq, phase, coh, [cohConfU; cohConfL]);
            else
              phaseFOI = NaN(size(FOI));
              cohFOI = NaN(size(FOI));
            end
            
            % Store values
            areaCohFOIindividual{iCondPlusAll}{area} = [areaCohFOIindividual{iCondPlusAll}{area}; cohFOI];
            areaPhaseFOIindividual{iCondPlusAll}{area} = [areaPhaseFOIindividual{iCondPlusAll}{area}; phaseFOI];
            areaCohFullIndividual{iCondPlusAll}{area}{numel(areaCohFullIndividual{iCondPlusAll}{area})+1} = coh;
            areaCohConfUFullIndividual{iCondPlusAll}{area}{numel(areaCohConfUFullIndividual{iCondPlusAll}{area})+1} = cohConfU;
            areaCohConfLFullIndividual{iCondPlusAll}{area}{numel(areaCohConfLFullIndividual{iCondPlusAll}{area})+1} = cohConfL;
            areaPhaseFullIndividual{iCondPlusAll}{area}{numel(areaPhaseFullIndividual{iCondPlusAll}{area})+1} = phase;
            areaPhaseConfUFullIndividual{iCondPlusAll}{area}{numel(areaPhaseConfUFullIndividual{iCondPlusAll}{area})+1} = phaseConfU;
            areaPhaseConfLFullIndividual{iCondPlusAll}{area}{numel(areaPhaseConfLFullIndividual{iCondPlusAll}{area})+1} = phaseConfL;
            areaFreqFullIndividual{iCondPlusAll}{area}{numel(areaFreqFullIndividual{iCondPlusAll}{area})+1} = freq;
          end
        end
      end
    end
  end
  
  areas = areasReverse;
  
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
if reverse
  filename = [mainFolder filesep 'globalMUAs_ca_reverse.mat']; %#ok<*UNRCH>
else
  filename = [mainFolder filesep 'globalMUAs_ca.mat'];
end
filename_la = [mainFolder_la filesep 'globalMUAs.mat'];
if fullRun
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  if exist(filename, 'file')
    save(filename, 'conditions','areas','FOI','areaCohFOIindividual','areaPhaseFOIindividual',...
      'areaCohFullIndividual','areaCohConfUFullIndividual','areaCohConfLFullIndividual',...
      'areaPhaseFullIndividual','areaPhaseConfUFullIndividual','areaPhaseConfLFullIndividual',...
      'areaFreqFullIndividual','areaCohFullInterpIndividual','areaCohConfUFullInterpIndividual',...
      'areaCohConfLFullInterpIndividual','areaPhaseFullInterpIndividual','areaPhaseConfUFullInterpIndividual',...
      'areaPhaseConfLFullInterpIndividual','areaFreqFullInterpIndividual', '-append');
  else
    save(filename, 'conditions','areas','FOI','areaCohFOIindividual','areaPhaseFOIindividual',...
      'areaCohFullIndividual','areaCohConfUFullIndividual','areaCohConfLFullIndividual',...
      'areaPhaseFullIndividual','areaPhaseConfUFullIndividual','areaPhaseConfLFullIndividual',...
      'areaFreqFullIndividual','areaCohFullInterpIndividual','areaCohConfUFullInterpIndividual',...
      'areaCohConfLFullInterpIndividual','areaPhaseFullInterpIndividual','areaPhaseConfUFullInterpIndividual',...
      'areaPhaseConfLFullInterpIndividual','areaFreqFullInterpIndividual');
  end
else
  load(filename);
end
localData = load(filename_la);
iAreas2compareOI = find(ismember(areas,areas2compareCritical));

% Load PR data
try
  if reverse
    fullData = load([rootFolder filesep PRsFolder filesep 'globalPRs_ca_reverse.mat']);
  else
    fullData = load([rootFolder filesep PRsFolder filesep 'globalPRs_ca.mat']);
  end
  areaPhaseFullInterpIndividualPR = fullData.areaPhaseFullInterpIndividual;
  areaCohFullInterpIndividualPR = fullData.areaCohFullInterpIndividual;
catch
  areaPhaseFullInterpIndividualPR = [];
  areaCohFullInterpIndividualPR = [];
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
  options.iAreasOI = iAreas2compareOI;
  phaseFreqProfilePlotMeans(areas, conditions(1:min([numel(conditions) 2])),...
    areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual, areaPhaseFullInterpIndividualPR, figFileName, options);
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
  options.iAreasOI = iAreas2compareOI;
  cohFreqProfilePlotMeans(areas, conditions(1:min([numel(conditions) 2])),...
    areaFreqFullInterpIndividual, areaCohFullInterpIndividual, areaCohFullInterpIndividualPR, figFileName, options);
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
options.iAreasOI = iAreas2compareOI;
[phaseHistos, distributionStats] = phaseHistosPlotMaster(drawPhaseHistos, areas, conditions(1:min([numel(conditions) 2])),...
  areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual, edges, options);
if ~isempty(phaseHistos) && ~isempty(distributionStats)
  save(filename, 'phaseHistos','distributionStats', '-append');
end


%% STATS ON MEAN PHASE FREQUENCY PROFILES
if doStats
  if ~strcmp(repository, 'allensdk')
    [fPEst, fWTest, strPMethod, pEst, U2, pObs, U2Obs] = phaseComparisonStats(areas, areaFreqFullInterpIndividual, FOI, areaPhaseFullInterpIndividual, [areasCritical; areas2compareCritical]);
  else
    fPEst = []; fWTest = []; strPMethod = []; pEst = []; U2 = []; pObs = []; U2Obs = [];
  end
  save(filename, 'fPEst','fWTest','strPMethod','pEst','U2','pObs','U2Obs', '-append');
end

if doStats_ca
  [fPEst_ca, fWTest_ca, strPMethod_ca, pEst_ca, U2_ca, pObs_ca, U2Obs_ca, comparisons] = phaseComparisonStats_ca(areas, localData.areas, areaFreqFullInterpIndividual, localData.areaFreqFullInterpIndividual, [0.3 0.03],...
    areaPhaseFullInterpIndividual, localData.areaPhaseFullInterpIndividual, areas2compareCriticalReduced, areasCriticalReduced);
  save(filename, 'fPEst_ca','fWTest_ca','strPMethod_ca','pEst_ca','U2_ca','pObs_ca','U2Obs_ca','comparisons', '-append');
end