% Run this script to perform analyses comparing PR and pupil area data.

clearvars -except repository subpop reverse qualityCheck allData fullRun includeRuns 


%% INITIALISE PARAMETERS
params
lists

if ~exist('repository', 'var')
  repository = 'all';
end
if ~exist('fullRun', 'var')
  fullRun = true;
end

outputDir = [outputDir filesep includeRuns];
if strcmp(repository,'uol')
  dataDir = [dataDir_local filesep '001_uol'];
elseif strcmp(repository,'allensdk')
  dataDir = [dataDir_local filesep '002_allen'];
end
if strcmp(repository,'all')
  if strcmp(subpop, 'all')
    mainFolder = [outputDir filesep area2pupilDir filesep PRsFolder];
  elseif strcmp(subpop, 'positive')
    mainFolder = [outputDir filesep area2pupilDir_positive filesep PRsFolder];
  elseif strcmp(subpop, 'negative')
    mainFolder = [outputDir filesep area2pupilDir_negative filesep PRsFolder];
  end
  animals = animalsOI;
  xLim = freqLimUOL;
elseif strcmp(repository,'uol')
  if strcmp(subpop, 'all')
    mainFolder = [outputDir filesep area2pupilDir_uol filesep PRsFolder];
  elseif strcmp(subpop, 'positive')
    mainFolder = [outputDir filesep area2pupilDir_uol_positive filesep PRsFolder];
  elseif strcmp(subpop, 'negative')
    mainFolder = [outputDir filesep area2pupilDir_uol_negative filesep PRsFolder];
  end
  animals = animalsUOLOI;
  xLim = freqLimUOL;
elseif strcmp(repository,'allensdk')
  if strcmp(subpop, 'all')
    mainFolder = [outputDir filesep area2pupilDir_allensdk filesep PRsFolder];
  elseif strcmp(subpop, 'positive')
    mainFolder = [outputDir filesep area2pupilDir_allensdk_positive filesep PRsFolder];
  elseif strcmp(subpop, 'negative')
    mainFolder = [outputDir filesep area2pupilDir_allensdk_negative filesep PRsFolder];
  end
  animals = animalsAllensdk;
  conditions = {'awake'};
  xLim = freqLimAllen;
end

drawPhaseProfiles = [true true];
drawCohProfiles = [true true];
drawPhaseHistos = [false true true];
doStats = false;


%% COMPUTE VARIABLES AND INITIALISE FIGURES NEEDED FOR DISPLAYING UNIT PHASE AND COHERENCE FREQUENCY PROFILES
if fullRun
  for animal = 1:numel(animals) % Loop through animals
    disp(['processing animal ' animals{animal}]);
    if exist('allData', 'var')
      dataStruct = allData.(animals{animal}).dataStruct;
    else
      load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
    end
    animalColour = animalColours(animals(animal));
    fnsData = fieldnames(dataStruct.seriesData);
    
    % Initialise figures and storage variables
    FOI = dataStruct.seriesData.(fnsData{1}).conf.FOI;
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
    
    for dbCount = 1:numel(fnsData) % Loop through database entries
      dbStruct = dataStruct.seriesData.(fnsData{dbCount});
      if isempty(dbStruct)
        continue
      end
      seriesName = seriesFromEntry(fnsData{dbCount});
      
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
      
      for iAreaPlusAll = area % Loop through the main and pooled areas
        if numel(conditions) > 1
          condLoop = [iCond numel(conditions)];
        else
          condLoop = iCond;
        end
        for iCondPlusAll = condLoop % Loop through the main and pooled conditions
          
          freq = dbStruct.popData.pupil.popData.freq;
          %[~, endFreq] = min(abs(freq - 2));
          %freq = freq(1:endFreq);
          
          % Get coherence and phase values
          coh = dbStruct.popData.pupil.popData.coh;
          coh_conf = dbStruct.popData.pupil.popData.coh_conf;
          rateadjust_kappa = dbStruct.popData.pupil.popData.rateadjust_kappa;
          phase = dbStruct.popData.pupil.popData.phase;
          phase_confU = dbStruct.popData.pupil.popData.phase_confU;
          phase_confL = dbStruct.popData.pupil.popData.phase_confL;
          [phase, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase,...
            phase_confU, phase_confL, coh, coh_conf, rateadjust_kappa);
%           [~, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase,...
%             phase_confU, phase_confL, coh, coh_conf, rateadjust_kappa);
          
          % Get phase and coherence values for FOI
          if sum(~isnan(phase)) && sum(~isnan(coh))
            [phaseFOI, cohFOI] = phaseCohFOI(FOI, freq, phase, coh, [cohConfU; cohConfL]);
          else
            phaseFOI = NaN(size(FOI));
            cohFOI = NaN(size(FOI));
          end
          
          % Store the values
          areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaCohFOIindividual{iCondPlusAll}{iAreaPlusAll}; cohFOI];
          areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaPhaseFOIindividual{iCondPlusAll}{iAreaPlusAll}; phaseFOI];
          areaCohFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCohFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = coh;
          areaCohConfUFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCohConfUFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = cohConfU;
          areaCohConfLFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCohConfLFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = cohConfL;
          areaPhaseFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = phase;
          areaPhaseConfUFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseConfUFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = phaseConfU;
          areaPhaseConfLFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseConfLFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = phaseConfL;
          areaFreqFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaFreqFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = freq;
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
filename = [mainFolder filesep 'globalPRs_area2pupil.mat'];
if fullRun
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  save(filename, 'conditions','areas','FOI','areaPhaseFOIindividual','areaCohFOIindividual',...
    'areaPhaseFullIndividual','areaPhaseConfUFullIndividual','areaPhaseConfLFullIndividual',...
    'areaCohFullIndividual','areaCohConfUFullIndividual','areaCohConfLFullIndividual',...
    'areaFreqFullIndividual','areaPhaseFullInterpIndividual','areaPhaseConfUFullInterpIndividual',...
    'areaPhaseConfLFullInterpIndividual','areaCohFullInterpIndividual','areaCohConfUFullInterpIndividual',...
    'areaCohConfLFullInterpIndividual','areaFreqFullInterpIndividual', '-v7.3');
else
  load(filename); %#ok<*UNRCH>
end


%% GENERATE PHASE FREQUENCY PROFILE FIGURES WITH MEANS AND SAVE THEM
if drawPhaseProfiles(1)
  figFileName = '%sVsPupil_%s_phase';
  options = struct();
  options.figTitle = 'Population rate phase comparisons: %sVsPupil %s';
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.phaseFrequencyProfilesSubfolder = phaseFrequencyProfilesSubfolder;
  options.freqLim = [xLim(1) 2];
  options.phaseLim = [-pi pi] + [-pi/4 pi/4];
  options.cutoffFreq = 0.2;
  options.iAreasOI = iAreasOI;
  phaseFreqProfilePlotIndividual(areas, conditions(1:min([numel(conditions) 2])),...
    areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual, [], [], figFileName, options);
end

if drawPhaseProfiles(2)
  figFileName = 'Means_only__%sVsPupil';
  options = struct();
  options.figTitle = 'Population rate mean phase comparisons: %sVsPupil';
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.phaseFrequencyProfilesSubfolder = phaseFrequencyProfilesSubfolder;
  options.freqLim = [xLim(1) 2];
  options.phaseLim = [-pi pi] + [-pi/4 pi/4];
  options.cutoffFreq = 0.2;
  options.iAreasOI = iAreasOI;
  phaseFreqProfilePlotMeans(areas, conditions(1:min([numel(conditions) 2])),...
    areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual, [], figFileName, options);
end


%% GENERATE COHERENCE FREQUENCY PROFILE FIGURES WITH MEANS AND SAVE THEM
if drawCohProfiles(1)
  figFileName = '%sVsPupil_%s_coherence';
  options = struct();
  options.figTitle = 'Population rate coherence comparisons: %sVsPupil %s';
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.cohFrequencyProfilesSubfolder = coherenceFrequencyProfilesSubfolder;
  options.freqLim = [xLim(1) 2];
  options.cohLim = [0 1];
  options.iAreasOI = iAreasOI;
  cohFreqProfilePlotIndividual(areas, conditions(1:min([numel(conditions) 2])),...
    areaFreqFullInterpIndividual, areaCohFullInterpIndividual, [], [], figFileName, options);
end

if drawCohProfiles(2)
  figFileName = 'Means_only__%sVsPupil';
  options = struct();
  options.figTitle = 'Population rate mean coherence comparisons: %sVsPupil';
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.cohFrequencyProfilesSubfolder = coherenceFrequencyProfilesSubfolder;
  options.freqLim = [xLim(1) 2];
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
options.figTitle = 'PUPIL';
options.freqLim = [xLim(1)+0.002 2];
phaseShift = pi/2;
options.phaseLim = phaseLim + phaseShift;
options.xLabelHist = '# recordings';
options.iAreasOI = iAreasOI;
[phaseHistos, distributionStats] = phaseHistosPlotMaster(drawPhaseHistos, areas, conditions(1:min([numel(conditions) 2])),...
  areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual, edges+phaseShift, options);


%% STATS ON MEAN PHASE FREQUENCY PROFILES
if doStats
  if ~strcmp(repository, 'allensdk')
    [fPEst, fWTest, strPMethod, pEst, U2, pObs, U2Obs] = phaseComparisonStats(areas, areaFreqFullInterpIndividual, FOI, areaPhaseFullInterpIndividual, [areasCritical; areas2compareCritical]);
  else
    fPEst = []; fWTest = []; strPMethod = []; pEst = []; U2 = []; pObs = []; U2Obs = [];
  end
  save(filename, 'phaseHistos','distributionStats','fPEst','fWTest','strPMethod','pEst','U2','pObs','U2Obs', '-append');
end