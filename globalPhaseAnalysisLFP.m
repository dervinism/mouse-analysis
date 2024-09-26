% Run this script to produce unit phase and coherence frequency profiles
% and phase frequency histograms for cross-area LFP comparisons.
%
% The following data files are produced:
%   dataDir\lfp2lfpDir\histosSubfolder\globalLFPs_ca.mat or
%   dataDir\lfp2lfpDir\histosSubfolder\globalLFPs_ca_reverse.mat. All
%   of these files contain unit phase frequency profiles, phase frequency
%   histograms, and statistical test results.
%
% LFP phase frequency profile figures are saved in
%   dataDir\lfp2lfpDir\phaseFrequencyProfilesSubfolder.
% LFP phase frequency histogram figures and summary suplots are saved in
%   dataDir\lfp2lfpDir\histosSubfolder.
% Phase frequency maps are saved in dataDir\lfp2lfpDir\mapsSubfolder.

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
if ~exist('reverse', 'var')
  reverse = false;
end
dataDir = [dataDir filesep includeRuns];

if strcmp(repository,'all')
  mainFolder = [dataDir filesep lfp2lfpDir];
  animals = animalsOI;
  xLim = freqLimUOL;
elseif strcmp(repository,'uol')
  mainFolder = [dataDir filesep lfp2lfpDir_uol];
  animals = animalsUOLOI;
  xLim = freqLimUOL;
elseif strcmp(repository,'allensdk')
  mainFolder = [dataDir filesep lfp2lfpDir_allensdk];
  animals = animalsAllensdk;
  conditions = {'awake'};
  xLim = freqLimAllen;
end
areas = areas2compare;

drawPhaseProfiles = [true true];
drawCohProfiles = [true true];
drawPhaseHistos = [true true true];


%% COMPUTE VARIABLES NEEDED FOR DISPLAYING LFP PHASE AND COHERENCE FREQUENCY PROFILES
if fullRun
  for animal = 1:numel(animals) % Loop through animals
    disp(['processing animal ' animals{animal}]);
    if exist('allData', 'var')
      dataStruct = allData.(animals{animal}).dataStruct;
    else
      load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
    end
    animalColour = animalColours(animals(animal));
    if strcmp(subpop, 'all')
      if ~isfield(dataStruct, 'seriesData_ca')
        disp(['Animal ' animals{animal} ' is missing series comparison data. Skipping to the next animal...'])
        continue
      end
      fnsData_ca = fieldnames(dataStruct.seriesData_ca);
    end
    
    % Initialise storage variables
    fnsData = fieldnames(dataStruct.seriesData);
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
    
    for dbCount = 1:numel(fnsData_ca) % Loop through database entries
      dbStruct_ca = dataStruct.seriesData_ca.(fnsData_ca{dbCount});
      if isempty(dbStruct_ca)
        continue
      end
      [seriesName1, seriesName2] = seriesNames(fnsData_ca{dbCount});
      
      % Test for exceptions
      if exceptionTest(except, seriesName1, seriesName2)
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
      
      for iComp = 1:numel(comp) % Loop through non-grouped and grouped area comparisons
        area = comp(iComp);
        
        if numel(conditions) > 1
          condLoop = [iCond numel(conditions)];
        else
          condLoop = iCond;
        end
        for iCondPlusAll = condLoop % Loop through the main and pooled conditions
          
          lfpData = dbStruct_ca.lfpphaseCohDataLFP;
          for ch = 1:numel(lfpData.fastPower) % Loop through local channels
            for ch_ca = 1:numel(lfpData.fastPower{ch}) % Loop trough external channels
              phaseCoh = lfpData.fastPower{ch}{ch_ca};
              
              % Load coherence
              coh = phaseCoh.coh;
              coh_conf = phaseCoh.coh_conf;
              
              coh_confU = coh + coh_conf;
              coh_confL = coh - coh_conf;
              coh_confU(coh_confL <= 0) = NaN;
              coh_confL(coh_confL <= 0) = NaN;
              coh(isnan(coh_confU) | isnan(coh_confL)) = NaN;
              
              % Load phase of the population rate
              phase = bestUnwrap(phaseCoh.phase);
              phase_confU = phaseCoh.phase_confU;
              phase_confL = phaseCoh.phase_confL;
              phase(isnan(phase_confU) | isnan(phase_confL) | isnan(coh) | isnan(coh_confU) | isnan(coh_confL)) = NaN;
              freq = phaseCoh.freq;
              
              % Ascertain that significant number of entries is the same for both phase and coherence
              coh(isnan(phase) | isnan(phase_confU) | isnan(phase_confL)) = NaN;
              coh_confU(isnan(phase) | isnan(phase_confU) | isnan(phase_confL)) = NaN;
              coh_confL(isnan(phase) | isnan(phase_confU) | isnan(phase_confL)) = NaN;
              assert(sum(isnan(phase)) == sum(isnan(coh)));
              
              % Obtain and store phase and coherence values for FOI
              if sum(~isnan(phase)) && sum(~isnan(coh))
                [phaseFOI, cohFOI, cohConfFOI] = phaseCohFOI(FOI, freq, phase, coh, [coh_confU; coh_confL]);
              else
                phaseFOI = NaN(size(FOI));
                cohFOI = NaN(size(FOI));
                cohConfFOI = NaN(2,numel(FOI));
              end
              cohConfUFOI = cohConfFOI(1,:);
              cohConfLFOI = cohConfFOI(2,:);
              areaCohFOIindividual{iCondPlusAll}{area} = [areaCohFOIindividual{iCondPlusAll}{area}; cohFOI];
              areaPhaseFOIindividual{iCondPlusAll}{area} = [areaPhaseFOIindividual{iCondPlusAll}{area}; phaseFOI];
              
              % Obtain and store original full phase and coherence values
              areaCohFullIndividual{iCondPlusAll}{area}{numel(areaCohFullIndividual{iCondPlusAll}{area})+1} = coh;
              areaCohConfUFullIndividual{iCondPlusAll}{area}{numel(areaCohConfUFullIndividual{iCondPlusAll}{area})+1} = coh_confU;
              areaCohConfLFullIndividual{iCondPlusAll}{area}{numel(areaCohConfLFullIndividual{iCondPlusAll}{area})+1} = coh_confL;
              areaPhaseFullIndividual{iCondPlusAll}{area}{numel(areaPhaseFullIndividual{iCondPlusAll}{area})+1} = phase;
              areaPhaseConfUFullIndividual{iCondPlusAll}{area}{numel(areaPhaseConfUFullIndividual{iCondPlusAll}{area})+1} = phase_confU;
              areaPhaseConfLFullIndividual{iCondPlusAll}{area}{numel(areaPhaseConfLFullIndividual{iCondPlusAll}{area})+1} = phase_confL;
              areaFreqFullIndividual{iCondPlusAll}{area}{numel(areaFreqFullIndividual{iCondPlusAll}{area})+1} = freq;
              
            end
          end
        end
      end
    end
  end
  
  % Update area names in case they are reversed
  areas = areasReverse;
  areaRecCountMeaning = areas;
  
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
  filename = [mainFolder filesep 'globalLFPs_ca_reverse.mat']; %#ok<*UNRCH>
else
  filename = [mainFolder filesep 'globalLFPs_ca.mat'];
end
if fullRun
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  save(filename, 'conditions','areas','FOI','areaPhaseFOIindividual','areaCohFOIindividual','areaPhaseFullIndividual',...
    'areaPhaseConfUFullIndividual','areaPhaseConfLFullIndividual','areaCohFullIndividual',...
    'areaCohConfUFullIndividual','areaCohConfLFullIndividual','areaFreqFullIndividual','areaPhaseFullInterpIndividual',...
    'areaPhaseConfUFullInterpIndividual','areaPhaseConfLFullInterpIndividual','areaCohFullInterpIndividual',...
    'areaCohConfUFullInterpIndividual','areaCohConfLFullInterpIndividual','areaFreqFullInterpIndividual', '-v7.3');
else
  load(filename);
end
iAreas2compareOI = find(ismember(areas,areas2compareCritical));


%% GENERATE PHASE FREQUENCY PROFILE FIGURES WITH MEANS AND SAVE THEM
if drawPhaseProfiles(1)
  figFileName = '%s_%s_phase';
  options = struct();
  options.figTitle = 'LFP phase comparisons: %s %s';
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.phaseFrequencyProfilesSubfolder = phaseFrequencyProfilesSubfolder;
  options.freqLim = xLim;
  options.phaseLim = [-pi pi] + [-pi/4 pi/4];
  options.iAreasOI = iAreas2compareOI;
  phaseFreqProfilePlotIndividual(areas, conditions(1:min([numel(conditions) 2])),...
    areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual, [], [], figFileName, options);
end

if drawPhaseProfiles(2)
  figFileName = 'Means_only__%s';
  options = struct();
  options.figTitle = 'LFP mean phase comparisons: %s';
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.phaseFrequencyProfilesSubfolder = phaseFrequencyProfilesSubfolder;
  options.freqLim = xLim;
  options.phaseLim = [-pi pi] + [-pi/4 pi/4];
  options.iAreasOI = iAreas2compareOI;
  phaseFreqProfilePlotMeans(areas, conditions(1:min([numel(conditions) 2])),...
    areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual, [], figFileName, options);
end


%% GENERATE COHERENCE FREQUENCY PROFILE FIGURES WITH MEANS AND SAVE THEM
if drawCohProfiles(1)
  figFileName = '%s_%s_coherence';
  options = struct();
  options.figTitle = 'LFP coherence comparisons: %s %s';
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.cohFrequencyProfilesSubfolder = coherenceFrequencyProfilesSubfolder;
  options.freqLim = xLim;
  options.cohLim = [0 1];
  options.iAreasOI = iAreas2compareOI;
  cohFreqProfilePlotIndividual(areas, conditions(1:min([numel(conditions) 2])),...
    areaFreqFullInterpIndividual, areaCohFullInterpIndividual, [], [], figFileName, options);
end

if drawCohProfiles(2)
  figFileName = 'Means_only__%s';
  options = struct();
  options.figTitle = 'LFP mean coherence comparisons: %s';
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.cohFrequencyProfilesSubfolder = coherenceFrequencyProfilesSubfolder;
  options.freqLim = xLim;
  options.cohLim = [0 1];
  options.iAreasOI = iAreas2compareOI;
  cohFreqProfilePlotMeans(areas, conditions(1:min([numel(conditions) 2])),...
    areaFreqFullInterpIndividual, areaCohFullInterpIndividual, [], figFileName, options);
end


%% DRAW PHASE FREQUENCY HISTOGRAMS AND MAPS
options = struct();
options.mainFolder = mainFolder;
options.histosSubfolder = histosSubfolder;
options.mapsSubfolder = mapsSubfolder;
options.figSize = figSize;
options.figTitle = 'LFP';
options.freqLim = [xLim(1)+0.002 xLim(2)];
options.phaseLimHisto = phaseLim;
options.phaseLimMap = phaseLim;
options.xLabelHist = '# recordings';
options.iAreasOI = iAreas2compareOI;
[phaseHistos, distributionStats] = phaseHistosPlotMaster(drawPhaseHistos, areas, conditions(1:min([numel(conditions) 2])),...
  areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual, edges, options);


%% STATS ON MEAN PHASE FREQUENCY PROFILES
if ~strcmp(repository, 'allensdk')
  [fPEst, fWTest, strPMethod, pEst, U2, pObs, U2Obs] = phaseComparisonStats(areas, areaFreqFullInterpIndividual, FOI, areaPhaseFullInterpIndividual, [areasCritical; areas2compareCritical]);
else
  fPEst = []; fWTest = []; strPMethod = []; pEst = []; U2 = []; pObs = []; U2Obs = [];
end
save(filename, 'phaseHistos','distributionStats','fPEst','fWTest','strPMethod','pEst','U2','pObs','U2Obs', '-append');