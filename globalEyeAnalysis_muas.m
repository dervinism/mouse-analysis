% Run this script to produce unit phase and coherence frequency profiles
% and phase frequency histograms for area-to-pupil comparisons.
%
% The following data files are produced:
% dataDir\area2pupilDir\unitsFolder\histosSubfolder\phaseHistos_units.mat or
%   dataDir\area2pupilDir\unitsFolder\histosSubfolder\phaseHistos_units_quality.mat or
%   All of these files contain unit phase frequency profiles, phase frequency
%   histograms, and statistical test results.
%
% Unit phase frequency profile figures are saved in dataDir\area2pupilDir\unitsFolder\phaseFrequencyProfilesSubfolder
%   or dataDir\area2pupilDir\qualityUnitsFolder\phaseFrequencyProfilesSubfolder.
% Individual unit phase frequency histogram figures and summary suplots are
%   saved in dataDir\area2pupilDir\unitsFolder\histosSubfolder or dataDir\area2pupilDir\qualityUnitsFolder\histosSubfolder.
% Phase frequency maps are saved in dataDir\area2pupilDir\unitsFolder\mapsSubfolder or
%   dataDir\area2pupilDir\qualityUnitsFolder\mapsSubfolder.

clearvars -except repository subpop reverse qualityCheck allData fullRun includeRuns


%% INITIALISE PARAMETERS
params
lists

if ~exist('repository', 'var')
  repository = 'uol';
end
if ~exist('fullRun', 'var')
  fullRun = true;
end
if strcmpi(includeRuns, 'noRun')
  fRef = 0.3;
else
  fRef = 0.03;
end

dataDir = [dataDir filesep includeRuns];
if strcmp(repository,'uol')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep area2pupilDir_uol];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep area2pupilDir_uol_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep area2pupilDir_uol_negative];
  end
  animals = animalsUOLOI;
  xLim = freqLimUOL;
elseif strcmp(repository,'allensdk')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep area2pupilDir_allensdk];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep area2pupilDir_allensdk_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep area2pupilDir_allensdk_negative];
  end
  animals = animalsAllensdk;
  conditions = {'awake'};
  xLim = freqLimAllen;
end
mainFolder = [rootFolder filesep MUAsFolder];

drawPhaseProfiles = false; %true;
drawCohProfiles = false; %true;
drawPhaseHistos = [true false false]; %[true false true];
statsPhaseProfiles = false;
drawPhaseSpans = false;
drawBetaCorrelations = false;
drawPhase2phaseCorrelations = false; %true;
drawPhase2coherenceCorrelations = false; %true;


%% COMPUTE VARIABLES AND INITIALISE FIGURES NEEDED FOR DISPLAYING UNIT PHASE AND COHERENCE FREQUENCY PROFILES
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
      areaPhaseSpanIndividual = {};
      areaBetaIndividual = {};
      areaAnimalIDsIndividual = {};
      areaSeriesIDsIndividual = {};
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
        areaPhaseSpanIndividualCond = {};
        areaBetaIndividualCond = {};
        areaAnimalIDsIndividualCond = {};
        areaSeriesIDsIndividualCond = {};
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
          areaPhaseSpanIndividualCond{iArea} = [];
          areaBetaIndividualCond{iArea} = [];
          areaAnimalIDsIndividualCond{iArea} = {};
          areaSeriesIDsIndividualCond{iArea} = {};
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
        areaPhaseSpanIndividual{iCond} = areaPhaseSpanIndividualCond;
        areaBetaIndividual{iCond} = areaBetaIndividualCond;
        areaAnimalIDsIndividual{iCond} = areaAnimalIDsIndividualCond;
        areaSeriesIDsIndividual{iCond} = areaSeriesIDsIndividualCond;
      end
    end

    for dbCount = 1:numel(fnsData) % Loop through db entries
      dbStruct = dataStruct.seriesData.(fnsData{dbCount});
      if isempty(dbStruct)
        continue
      end
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
      
      % Determine mode boundaries
      if strcmpi(repository, 'uol')
        modeBoundaries = phaseHistoBinCentres(modesUOL{area(1)});
      elseif strcmpi(repository, 'allensdk')
        modeBoundaries = phaseHistoBinCentres(modesAllensdk{area(1)});
      end

      % Select correlated units if needed
      if (strcmpi(subpop, 'positive') || strcmpi(subpop, 'negative')) &&...
          ~isfield(dbStruct.shankData.(['shank' num2str(1)]), 'pupil')
        continue
      end
      units = dbStruct.popData.spkDB_units;
      if strcmpi(subpop, 'positive') || strcmpi(subpop, 'negative')
        if pupilCorrCond == 1
          if ~isfield(dbStruct.popData, 'pupil')
            continue
          else
            phase = spkPhase(dbStruct.popData.pupil.phaseCoh.unitData, fRef)';
          end
        elseif pupilCorrCond == 2
          if ~isfield(dbStruct.popData, 'rSpearman')
            continue
          else
            rSpearman = dbStruct.popData.rSpearman';
          end
        end
      end
      if strcmp(subpop, 'all')
        correlatedInd = 1:numel(units);
      elseif strcmp(subpop, 'positive')
        if pupilCorrCond == 1
          correlatedInd = false(size(phase));
          modeBoundaries = recentrePhase(modeBoundaries, modeBoundaries(1));
          phase = recentrePhase(phase, modeBoundaries(1));
          correlatedInd(phase > modeBoundaries(end) & phase <= modeBoundaries(2)) = true;
        elseif pupilCorrCond == 2
          correlatedInd = find(rSpearman >= 0);
        end
      elseif strcmp(subpop, 'negative')
        if pupilCorrCond == 1
          correlatedInd = false(size(phase));
          modeBoundaries = recentrePhase(modeBoundaries, modeBoundaries(3));
          phase = recentrePhase(phase, modeBoundaries(3));
          correlatedInd(phase > modeBoundaries(2) & phase <= modeBoundaries(end)) = true;
        elseif pupilCorrCond == 2
          correlatedInd = find(rSpearman < 0);
        end
      end
      units = units(correlatedInd);

      for iAreaPlusAll = area % Loop through the main and pooled areas
        if numel(conditions) > 1
          condLoop = [iCond numel(conditions)];
        else
          condLoop = iCond;
        end
        for iCondPlusAll = condLoop % Loop through the main and pooled conditions
          
          % Load and store phase and coherence data for units
          for u = 1:numel(units)
            iU = find(dbStruct.popData.spkDB_units == units(u));
            if isempty(iU)
              continue
            end
            
            % Get phase and coherence values
            freq = dbStruct.popData.pupil.phaseCoh.unitData{iU}.freq;
            coh = dbStruct.popData.pupil.phaseCoh.unitData{iU}.coh;
            cohConf = dbStruct.popData.pupil.phaseCoh.unitData{iU}.coh_conf;
            rateadjust_kappa = dbStruct.popData.pupil.phaseCoh.unitData{iU}.rateadjust_kappa;
            phase = dbStruct.popData.pupil.phaseCoh.unitData{iU}.phase;
            phaseConfU = dbStruct.popData.pupil.phaseCoh.unitData{iU}.phase_confU;
            phaseConfL = dbStruct.popData.pupil.phaseCoh.unitData{iU}.phase_confL;
            [phase, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase,...
              phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa);
            %               [~, phaseConfU, phaseConfL, coh, cohConfU, cohConfL] = correctPhaseCoh(phase,...
            %                 phaseConfU, phaseConfL, coh, cohConf, rateadjust_kappa);
            phaseSlow = phase(freq <= 1);
            phaseSpan = abs(max(phaseSlow,[],2)-min(phaseSlow,[],2));
            if isempty(phaseSpan)
              phaseSpan = NaN;
            end
            
            % Obtain and store phase and coherence values for FOI
            if sum(isnan(freq)) == numel(freq) || (~sum(~isnan(freq)) && sum(~isnan(coh)))
              freq = dbStruct.popData.pupil.phaseCoh.unitData{iU}.freq;
            end
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
            areaFreqFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaFreqFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = freq;
            areaPhaseSpanIndividual{iCondPlusAll}{iAreaPlusAll} = [areaPhaseSpanIndividual{iCondPlusAll}{iAreaPlusAll}; phaseSpan];
            %areaBetaIndividual{iCondPlusAll}{iAreaPlusAll} = [areaBetaIndividual{iCondPlusAll}{iAreaPlusAll}; dbStruct.popData.phaseCoh{iU}.beta];
            areaAnimalIDsIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaAnimalIDsIndividual{iCondPlusAll}{iAreaPlusAll})+1} = animals{animal};
            areaSeriesIDsIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaSeriesIDsIndividual{iCondPlusAll}{iAreaPlusAll})+1} = seriesName;
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
filename = [mainFolder filesep 'globalMUAs_area2pupil.mat']; %#ok<*UNRCH>
if fullRun
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  save(filename, 'conditions','areas','FOI','areaCohFOIindividual','areaPhaseFOIindividual',...
    'areaCohFullIndividual','areaCohConfUFullIndividual','areaCohConfLFullIndividual',...
    'areaPhaseFullIndividual','areaPhaseConfUFullIndividual','areaPhaseConfLFullIndividual',...
    'areaFreqFullIndividual','areaCohFullInterpIndividual','areaCohConfUFullInterpIndividual',...
    'areaCohConfLFullInterpIndividual','areaPhaseFullInterpIndividual','areaPhaseConfUFullInterpIndividual',...
    'areaPhaseConfLFullInterpIndividual','areaFreqFullInterpIndividual','areaPhaseSpanIndividual',...
    'areaPhaseSpanIndividual','areaBetaIndividual','areaAnimalIDsIndividual',...
    'areaSeriesIDsIndividual');
else
  load(filename);
end


%% GENERATE PHASE MEAN FREQUENCY PROFILE FIGURES AND SAVE THEM
if drawPhaseProfiles
  figFileName = 'Means_only__%sVsPupil';
  options = struct();
  options.figTitle = 'Unit mean phase comparisons: %sVsPupil';
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


%% GENERATE COHERENCE MEAN FREQUENCY PROFILE FIGURES AND SAVE THEM
if drawCohProfiles
  figFileName = 'Means_only__%sVsPupil';
  options = struct();
  options.figTitle = 'Unit mean coherence comparisons: %sVsPupil';
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
options.phaseLimHisto = phaseLim + phaseShift;
options.phaseLimMap = phaseLim - phaseShift;
options.xLabelHist = '# units';
options.limExpansion = pi/2;
options.mask = {[options.freqLim(1) options.phaseLimMap(end)+options.limExpansion; 2 options.phaseLimMap(end)+options.limExpansion;...
                 options.freqLim(1) pi/2; 2 options.phaseLimMap(end)+options.limExpansion];...
                [0.05 options.phaseLimMap(1)-options.limExpansion; options.freqLim(end) -pi/2;...
                 0.05 options.phaseLimMap(1)-options.limExpansion; options.freqLim(end) options.phaseLimMap(1)-options.limExpansion]};
options.iAreasOI = iAreasOI;
[phaseHistos, distributionStats] = phaseHistosPlotMaster(drawPhaseHistos, areas, conditions(1:min([numel(conditions) 2])),...
  areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual, edges+phaseShift, options);
save(filename, 'phaseHistos','distributionStats', '-append');


%% STATS ON MEAN PHASE FREQUENCY PROFILES
if statsPhaseProfiles
  if ~strcmp(repository, 'allensdk')
    [fPEst, fWTest, strPMethod, pEst, U2, pObs, U2Obs] = phaseComparisonStats(areas, areaFreqFullInterpIndividual, FOI, areaPhaseFullInterpIndividual, [areasCritical; areas2compareCritical]);
  else
    fPEst = []; fWTest = []; strPMethod = []; pEst = []; U2 = []; pObs = []; U2Obs = [];
  end
  save(filename, 'fPEst','fWTest','strPMethod','pEst','U2','pObs','U2Obs', '-append');
end


%% PLOT MEAN PHASE SPANS ACROSS AREAS
if drawPhaseSpans
  if strcmp(repository, 'uol')
    areasOI = {'VB','Po','Th','S1','RSC','Cx','CA','DG','Hp'};
  elseif strcmp(repository, 'allensdk')
    areasOI = {'LGN','LP','Th','V1','V2','VIS','CA','DG','Hp'};
  end
  [fH1, pval_ttestPhaseSpans, scatterGroups] = barPlotUnits(areaPhaseSpanIndividual, areasOI, 1);
  set(fH1(1), 'Name','Mean unit phase span across brain areas');
  filenameFig = [mainFolder filesep 'phaseSpanAcrossAreas'];
  savefig(fH1(1), filenameFig, 'compact');
  print(fH1(1), [filenameFig '.png'],'-dpng','-r300');
  
  anovaPhaseSpans = runAnova(scatterGroups, areasOI);
  save(filename, 'pval_ttestPhaseSpans','anovaPhaseSpans', '-append');
end


%% INDIVIDUAL UNIT PHASE PROFILES: AREA/ANIMAL/SERIES
if drawPhaseProfiles
  figFileName = '%s_%s_%s_%s_phase';
  options = struct();
  options.figTitle = 'Population rate phase comparisons: %s %s %s %s';
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.phaseFrequencyProfilesSubfolder = phaseFrequencyProfilesSubfolder;
  options.freqLim = [xLim(1) 2];
  options.phaseLim = [-pi pi] + [-pi/2 pi/2];
  options.cutoffFreq = 0.2;
  options.iAreasOI = iAreasOI;
  phaseFreqProfilePlotIndividual(areas, conditions(1:min([numel(conditions) 2])), areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual,...
    areaAnimalIDsIndividual, areaSeriesIDsIndividual, figFileName, options);
end


%% INDIVIDUAL UNIT PHASE PROFILES: AREA/ANIMAL
if drawPhaseProfiles
  figFileName = '%s_%s_%s_phase';
  options = struct();
  options.figTitle = 'Population rate phase comparisons: %s %s %s';
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.phaseFrequencyProfilesSubfolder = phaseFrequencyProfilesSubfolder;
  options.freqLim = [xLim(1) 2];
  options.phaseLim = [-pi pi] + [-pi/2 pi/2];
  options.cutoffFreq = 0.2;
  options.iAreasOI = iAreasOI;
  phaseFreqProfilePlotIndividual(areas, conditions(1:min([numel(conditions) 2])), areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual,...
    areaAnimalIDsIndividual, [], figFileName, options);
end


%% INDIVIDUAL UNIT PHASE PROFILES: AREA
% if drawPhaseProfiles
%   figFileName = '%s_%s_phase';
%   options = struct();
%   options.figTitle = 'Population rate phase comparisons: %s %s';
%   options.figSize = figSize;
%   options.mainFolder = mainFolder;
%   options.phaseFrequencyProfilesSubfolder = phaseFrequencyProfilesSubfolder;
%   options.freqLim = [xLim(1) 2];
%   options.phaseLim = [-pi pi] + [-pi/4 pi/4];
%   options.cutoffFreq = 0.2;
%   phaseFreqProfilePlotIndividual(areas, conditions(1:min([numel(conditions) 2])), areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual,...
%     [], [], figFileName, options);
% end


%% BETA AND COHERENCE CORRELATIONS
if drawBetaCorrelations
  if ~exist([mainFolder filesep betaCorrelationsSubfolder], 'dir')
    mkdir([mainFolder filesep betaCorrelationsSubfolder]);
  end
  for iCond = 1:numel(conditions(1:min([numel(conditions) 2])))
    for iArea = 1:numel(iAreasOI)
      beta = areaBetaIndividual{iCond}{iAreasOI(iArea)};
      coh = areaCohFullInterpIndividual{iCond}{iAreasOI(iArea)};
      if ~isempty(beta) && ~isempty(coh)
        disp(['Processing unit beta data for ' conditions{iCond} ' ' areas{iAreasOI(iArea)}...
          ' (comparison # ' num2str((iCond-1)*numel(areas) + iAreasOI(iArea)) '/' num2str(numel(conditions)*numel(areas)) ')']);
        
        % Correlations and individual graphs
        [rBeta{iCond}{iAreasOI(iArea)}, pvalBeta{iCond}{iAreasOI(iArea)}] = corrMulti(torow(beta), coh', 'Spearman');
        for f = 1:numel(areaFreqFullInterpIndividual)
          inds = ~isnan(beta) & ~isnan(coh(:,f));
          nBeta{iCond}{iAreasOI(iArea)}(f).significant = sum(inds);
          nBeta{iCond}{iAreasOI(iArea)}(f).total = numel(inds);
        end
        f = areaFreqFullInterpIndividual == 0.3;
        inds = ~isnan(beta) & ~isnan(coh(:,f));
        
        figBeta = figure;
        plot(beta(inds), coh(inds,f), '.', 'MarkerSize',5)
        figTitle = ['Unit PSD beta and pupil coherence correlations for ' areas{iAreasOI(iArea)} ' ' conditions{iCond} ' '...
          areaFreqFullInterpIndividual(f) 'Hz: r=' num2str(rBeta{iCond}{iAreasOI(iArea)}(f)) ' p=' num2str(pvalBeta{iCond}{iAreasOI(iArea)}(f))...
          ' n=' num2str(nBeta{iCond}{iAreasOI(iArea)}(f).significant) '/' num2str(nBeta{iCond}{iAreasOI(iArea)}(f).total)];
        ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
          'on', 'k', 'Unit PSD exponent \beta', xlim, xticks,...
          'on', 'k', 'Coherence with pupil', [-0.05 1.05], [0 0.2 0.4 0.6 0.8 1]);
        title(figTitle);
        figName = ['PSDbetaVsPupilCoherence_in_' areas{iAreasOI(iArea)} 'VsPupil_during_' conditions{iCond} '_at_'...
          num2str(areaFreqFullInterpIndividual(f)) '_Hz'];
        set(gcf, 'Name',figName);
        figName = [mainFolder filesep betaCorrelationsSubfolder filesep figName]; %#ok<*AGROW>
        figName = strrep(figName, ' ', '_');
        
        xLimLocal = xlim;
        xAxisLength = xLimLocal(2)-xLimLocal(1);
        yLim = ylim;
        yAxisLength = yLim(2)-yLim(1);
        text(xLimLocal(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.2, ['r=' num2str(rBeta{iCond}{iAreasOI(iArea)}(f))], 'FontSize',16);
        text(xLimLocal(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.1, ['p=' num2str(pvalBeta{iCond}{iAreasOI(iArea)}(f))], 'FontSize',16);
        text(xLimLocal(1)+xAxisLength*0.025, yLim(2)-yAxisLength*0.05,...
          ['n=' num2str(nBeta{iCond}{iAreasOI(iArea)}(f).significant) '/' num2str(nBeta{iCond}{iAreasOI(iArea)}(f).total)], 'FontSize',16);
        
        label = [2 1.6];
        margin = [0.3 0.55];
        width = 1*figSize-label(1)-margin(1);
        height = 1*figSize-label(2)-margin(2);
        paperSize = resizeFig(figBeta, ax1, width, height, label, margin, 0);
        exportFig(figBeta, [figName '.png'],'-dpng','-r300', paperSize);
        hgsave(figBeta, [figName '.fig']);
        close(figBeta);
      end
    end
  end
  
  save(filename, 'rBeta','pvalBeta','nBeta', '-append');
end


%% PHASE CORRELATIONS BETWEEN 0.03 AND 0.3 Hz
options = struct();
options.corrType = 'circularnp';
options.figFolder = [mainFolder filesep phaseCorrelationsSubfolder];
options.diagonal = false;
options.fitLine = false;
options.fitLineDisplay = 'on';
options.xLabel = 'Phase at 0.03 Hz (rad)';
options.xLim = [];
options.xTicks = [-2*pi -15*pi/8 -14*pi/8 -13*pi/8 -12*pi/8 -11*pi/8 -10*pi/8 -9*pi/8 -pi -7*pi/8 -6*pi/8 -5*pi/8 -4*pi/8 -3*pi/8 -2*pi/8 -pi/8 0 ...
  pi/8 2*pi/8 3*pi/8 4*pi/8 5*pi/8 6*pi/8 7*pi/8 pi 9*pi/8 10*pi/8 11*pi/8 12*pi/8 13*pi/8 14*pi/8 15*pi/8 2*pi];
options.yLabel = 'Phase at 0.3 Hz (rad)';
options.yLim = [];
options.yTicks = [-2*pi -15*pi/8 -14*pi/8 -13*pi/8 -12*pi/8 -11*pi/8 -10*pi/8 -9*pi/8 -pi -7*pi/8 -6*pi/8 -5*pi/8 -4*pi/8 -3*pi/8 -2*pi/8 -pi/8 0 ...
  pi/8 2*pi/8 3*pi/8 4*pi/8 5*pi/8 6*pi/8 7*pi/8 pi 9*pi/8 10*pi/8 11*pi/8 12*pi/8 13*pi/8 14*pi/8 15*pi/8 2*pi];
options.figSize = figSize;
options.saveFig = false;
options.printText = false;
if drawPhase2phaseCorrelations
  if ~exist([mainFolder filesep phaseCorrelationsSubfolder], 'dir')
    mkdir([mainFolder filesep phaseCorrelationsSubfolder]);
  end
  for iCond = 1:numel(conditions(1:min([numel(conditions) 2])))
    for iArea = 1:numel(iAreasOI)
      freq = areaFreqFullInterpIndividual;
      phase03 = recentrePhase(areaPhaseFullInterpIndividual{iCond}{iAreasOI(iArea)}(:,freq==0.03), mean(phaseLim));
      phase3 = recentrePhase(areaPhaseFullInterpIndividual{iCond}{iAreasOI(iArea)}(:,freq==0.3), mean(phaseLim));
      if ~isempty(phase03) && ~isempty(phase3)
        disp(['Processing unit phase2phase correlation data for ' conditions{iCond} ' ' areas{iAreasOI(iArea)}...
          ' (comparison # ' num2str((iCond-1)*numel(areas) + iAreasOI(iArea)) '/' num2str(numel(conditions)*numel(areas)) ')']);
        
        % Correlations and individual graphs
        inds = ~isnan(phase03) & ~isnan(phase3);
        nPhase2phase{iCond}{iAreasOI(iArea)}.oneSignificant = sum(~isnan(phase03) | ~isnan(phase3));
        nPhase2phase{iCond}{iAreasOI(iArea)}.bothSignificant = sum(inds);
        nPhase2phase{iCond}{iAreasOI(iArea)}.total = numel(inds);
        options.figName = ['phase2phaseCorrelations_in_' areas{iAreasOI(iArea)} 'VsPupil_during_' conditions{iCond} '_at_0p03_0p3_Hz'];
        
        [figPhase2phase, rPhase2phase{iCond}{iAreasOI(iArea)}, pvalPhase2phase{iCond}{iAreasOI(iArea)},...
          modelPhase2phase{iCond}{iAreasOI(iArea)}] = corrPlot(torow(phase03), torow(phase3), options);
        ax1 = gca;
        ax1.XTickLabel = {'-2\pi','-15\pi/8','-14\pi/8','-13\pi/8','-12\pi/8','-11\pi/8','-10\pi/8','-9\pi/8','-\pi','-7\pi/8','-6\pi/8','-5\pi/8','-4\pi/8','-3\pi/8','-2\pi/8','-\pi/8','0',...
                          '\pi/8','2\pi/8','3\pi/8','4\pi/8','5\pi/8','6\pi/8','7\pi/8','\pi','9\pi/8','10\pi/8','11\pi/8','12\pi/8','13\pi/8','14\pi/8','15\pi/8','2\pi'};
        ax1.YTickLabel = {'-2\pi','-15\pi/8','-14\pi/8','-13\pi/8','-12\pi/8','-11\pi/8','-10\pi/8','-9\pi/8','-\pi','-7\pi/8','-6\pi/8','-5\pi/8','-4\pi/8','-3\pi/8','-2\pi/8','-\pi/8','0',...
                          '\pi/8','2\pi/8','3\pi/8','4\pi/8','5\pi/8','6\pi/8','7\pi/8','\pi','9\pi/8','10\pi/8','11\pi/8','12\pi/8','13\pi/8','14\pi/8','15\pi/8','2\pi'};
        title(['Unit phase at 0.03 and 0.3 Hz correlations for ' areas{iAreasOI(iArea)} ' ' conditions{iCond} ' :'...
          ' r=' num2str(rPhase2phase{iCond}{iAreasOI(iArea)}) ' p=' num2str(pvalPhase2phase{iCond}{iAreasOI(iArea)})...
          ' n=' num2str(nPhase2phase{iCond}{iAreasOI(iArea)}.bothSignificant) '/' num2str(nPhase2phase{iCond}{iAreasOI(iArea)}.oneSignificant) '/'...
          num2str(nPhase2phase{iCond}{iAreasOI(iArea)}.total)]);
        figName = [mainFolder filesep phaseCorrelationsSubfolder filesep options.figName]; %#ok<*AGROW>
        figName = strrep(figName, ' ', '_');
        
        xLimLocal = xlim;
        xAxisLength = xLimLocal(2)-xLimLocal(1);
        yLim = ylim;
        yAxisLength = yLim(2)-yLim(1);
        text(xLimLocal(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.2, ['r=' num2str(rPhase2phase{iCond}{iAreasOI(iArea)})], 'FontSize',16);
        text(xLimLocal(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.1, ['p=' num2str(pvalPhase2phase{iCond}{iAreasOI(iArea)})], 'FontSize',16);
        text(xLimLocal(1)+xAxisLength*0.025, yLim(2)-yAxisLength*0.05,...
          ['n=' num2str(nPhase2phase{iCond}{iAreasOI(iArea)}.bothSignificant) '/' num2str(nPhase2phase{iCond}{iAreasOI(iArea)}.oneSignificant) '/'...
          num2str(nPhase2phase{iCond}{iAreasOI(iArea)}.total)], 'FontSize',16);
        if ~isempty(modelPhase2phase{iCond}{iAreasOI(iArea)})
          text(xLimLocal(1)+xAxisLength*0.7, yLim(2)-yAxisLength*0.05,...
            ['y=' num2str(round(modelPhase2phase{iCond}{iAreasOI(iArea)}(1),2)) 'x+' num2str(round(modelPhase2phase{iCond}{iAreasOI(iArea)}(2),2))], 'FontSize',16);
        end
        
        label = [2.1 2.1];
        margin = [0.3 0.55];
        width = 1*figSize-label(1)-margin(1);
        height = 1*figSize-label(2)-margin(2);
        paperSize = resizeFig(figPhase2phase, ax1, width, height, label, margin, 0);
        exportFig(figPhase2phase, [figName '.png'],'-dpng','-r300', paperSize);
        hgsave(figPhase2phase, [figName '.fig']);
        close(figPhase2phase);
      end
    end
  end
  
  save(filename, 'rPhase2phase','pvalPhase2phase','modelPhase2phase','nPhase2phase', '-append');
end


%% PHASE AND COHERENCE CORRELATIONS BETWEEN 0.03 AND 0.03 Hz
if drawPhase2coherenceCorrelations
  if ~exist([mainFolder filesep phaseCorrelationsSubfolder], 'dir')
    mkdir([mainFolder filesep phaseCorrelationsSubfolder]);
  end
  for iCond = 1:numel(conditions(1:min([numel(conditions) 2])))
    for iArea = 1:numel(iAreasOI)
      freq = areaFreqFullInterpIndividual;
      f = [0.03 0.3];
      for iF = 1:numel(f)
        phase = recentrePhase(areaPhaseFullInterpIndividual{iCond}{iAreasOI(iArea)}(:,freq==f(iF)), mean(phaseLim));
        coh = areaCohFullInterpIndividual{iCond}{iAreasOI(iArea)}(:,freq==f(iF));
        if ~isempty(phase) && ~isempty(coh)
          disp(['Processing unit phase2phase correlation data for ' conditions{iCond} ' ' areas{iAreasOI(iArea)}...
            ' (comparison # ' num2str((iCond-1)*numel(areas) + iAreasOI(iArea)) '/' num2str(numel(conditions)*numel(areas)) ')']);
          
          % Correlations and individual graphs
          [rPhase2coherence{iCond}{iAreasOI(iArea)}(iF), pvalPhase2coherence{iCond}{iAreasOI(iArea)}(iF)] = corrMulti(torow(phase), torow(coh), 'circlinear');
          inds = ~isnan(phase) & ~isnan(coh);
          nPhase2coherence{iCond}{iAreasOI(iArea)}(iF).significant = sum(inds);
          nPhase2coherence{iCond}{iAreasOI(iArea)}(iF).total = numel(inds);
          
          figPhase2coherence = figure;
          plot(phase(inds), coh(inds), '.', 'MarkerSize',10)
          figTitle = ['Unit phase and coherence correlation at ' num2str(f(iF)) ' Hz for ' areas{iAreasOI(iArea)} ' ' conditions{iCond} ' :'...
            ' r=' num2str(rPhase2coherence{iCond}{iAreasOI(iArea)}(iF)) ' p=' num2str(rPhase2coherence{iCond}{iAreasOI(iArea)}(iF))...
            ' n=' num2str(nPhase2coherence{iCond}{iAreasOI(iArea)}(iF).significant) '/' num2str(nPhase2coherence{iCond}{iAreasOI(iArea)}(iF).total)];
          ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
            'on', 'k', ['Phase at ' num2str(f(iF)) ' Hz (rad)'], xlim, [-2*pi -15*pi/8 -14*pi/8 -13*pi/8 -12*pi/8 -11*pi/8 -10*pi/8 -9*pi/8 -pi -7*pi/8 -6*pi/8 -5*pi/8 -4*pi/8 -3*pi/8 -2*pi/8 -pi/8 0 ...
            pi/8 2*pi/8 3*pi/8 4*pi/8 5*pi/8 6*pi/8 7*pi/8 pi 9*pi/8 10*pi/8 11*pi/8 12*pi/8 13*pi/8 14*pi/8 15*pi/8 2*pi],...
            'on', 'k', ['Coherence at ' num2str(f(iF)) ' Hz (rad)'], ylim, yticks);
          ax1.XTickLabel = {'-2\pi','-15\pi/8','-14\pi/8','-13\pi/8','-12\pi/8','-11\pi/8','-10\pi/8','-9\pi/8','-\pi','-7\pi/8','-6\pi/8','-5\pi/8','-4\pi/8','-3\pi/8','-2\pi/8','-\pi/8','0',...
            '\pi/8','2\pi/8','3\pi/8','4\pi/8','5\pi/8','6\pi/8','7\pi/8','\pi','9\pi/8','10\pi/8','11\pi/8','12\pi/8','13\pi/8','14\pi/8','15\pi/8','2\pi'};
          title(figTitle);
          figName = ['phase2coherenceCorrelations_in_' areas{iAreasOI(iArea)} 'VsPupil_during_' conditions{iCond} '_at_' num2str(f(iF)) '_Hz'];
          set(figPhase2coherence, 'Name',figName);
          figName = [mainFolder filesep phaseCorrelationsSubfolder filesep figName]; %#ok<*AGROW>
          figName = strrep(figName, ' ', '_');
          figName = strrep(figName, '.', 'p');
          
          xLimLocal = xlim;
          xAxisLength = xLimLocal(2)-xLimLocal(1);
          yLim = ylim;
          yAxisLength = yLim(2)-yLim(1);
          text(xLimLocal(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.2, ['r=' num2str(rPhase2coherence{iCond}{iAreasOI(iArea)}(iF))], 'FontSize',16);
          text(xLimLocal(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.1, ['p=' num2str(pvalPhase2coherence{iCond}{iAreasOI(iArea)}(iF))], 'FontSize',16);
          text(xLimLocal(1)+xAxisLength*0.025, yLim(2)-yAxisLength*0.05,...
            ['n=' num2str(nPhase2coherence{iCond}{iAreasOI(iArea)}(iF).significant) '/' num2str(nPhase2coherence{iCond}{iAreasOI(iArea)}(iF).total)], 'FontSize',16);
          
          label = [2.1 2.1];
          margin = [0.3 0.55];
          width = 1*figSize-label(1)-margin(1);
          height = 1*figSize-label(2)-margin(2);
          paperSize = resizeFig(figPhase2coherence, ax1, width, height, label, margin, 0);
          exportFig(figPhase2coherence, [figName '.png'],'-dpng','-r300', paperSize);
          hgsave(figPhase2coherence, [figName '.fig']);
          close(figPhase2coherence);
        end
      end
    end
  end
  
  save(filename, 'rPhase2coherence','pvalPhase2coherence','nPhase2coherence', '-append');
end



%% Local functions
function [fH, stats, dataScatter, areas, colourGroups] = barPlotUnits(data, areas, condition, yLim, scaleType) %#ok<*DEFNU>

if nargin < 8
  scaleType = 'normal';
end
if nargin < 7
  yLim = [];
end
if nargin < 6
  condition = 1;
end
nAreas = numel(areas);

[stats.pval_ttest, stats.area1, stats.area2] = ttestGroup(areas, data, 1);


%% Draw combined bar graphs
fH = figProperties('Bar plot for units', 'normalized', [0, .005, .97, .90], 'w', 'on');
hold on


% Data bars
gap = 1;
nBarsPerGroup = 1;
dataMean = zeros(size(areas));
dataCI95 = zeros(2,nAreas);
bars = zeros(size(areas));
dataScatter = cell(size(areas));
colourGroups = cell(size(areas));
for iArea = 1:nAreas
  areaCode = determineArea(areas{iArea});
  if strcmp(areas{iArea},'VB') || strcmp(areas{iArea},'lVB') || strcmp(areas{iArea},'rVB')
    colour = matlabColours(1);
  elseif strcmp(areas{iArea},'LGN') || strcmp(areas{iArea},'lLGN') || strcmp(areas{iArea},'rLGN')
    colour = matlabColours(2);
  elseif strcmp(areas{iArea},'S1') || strcmp(areas{iArea},'lS1') || strcmp(areas{iArea},'rS1')
    colour = matlabColours(13);
  elseif strcmp(areas{iArea},'V1') || strcmp(areas{iArea},'lV1') || strcmp(areas{iArea},'rV1')
    colour = matlabColours(11);
  elseif strcmp(areas{iArea},'RSC') || strcmp(areas{iArea},'lRSC') || strcmp(areas{iArea},'rRSC')
    colour = matlabColours(3);
  elseif strcmp(areas{iArea},'CA') || strcmp(areas{iArea},'lCA') || strcmp(areas{iArea},'rCA')
    colour = matlabColours(4);
  end
  bars(iArea) = 1+iArea*(nBarsPerGroup+gap);
  [dataMean(iArea), dataCI95(:,iArea)] = datamean(data{condition}{areaCode});
  bar(bars(iArea), dataMean(iArea), 'FaceColor',colour, 'EdgeColor',colour, 'FaceAlpha',0.2);
  
  dataScatter{iArea} = data{condition}{areaCode};
  colourGroups{iArea} = colour;
end


% Scatter
for iBar = 1:numel(bars)
  scatter(bars(iBar)*ones(size(dataScatter{iBar}))', dataScatter{iBar}',...
    'MarkerEdgeColor',colourGroups{iBar}, 'jitter','on'); %colourGroups(iBar,:));
end


% Error bars
er = errorbar(bars,dataMean,dataCI95(2,:),dataCI95(1,:));
er.Color = [0 0 0];
er.LineStyle = 'none';


% Graph adjustments
xTickPos = bars;
if isempty(yLim)
  yLim = ylim;
end
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Arial', 25, 4/3, 2, [0.003 0.003], 'out', 'on', 'k', {}, [],...
  xTickPos, 'on', 'k', {'Mean phase span (rad)'}, yLim, yticks);
xTickLabel = areas;
for iLabel = 1:numel(xTickLabel)
  xTickLabel{iLabel} = strrep(xTickLabel{iLabel}, 'l','');
  xTickLabel{iLabel} = strrep(xTickLabel{iLabel}, 'r','');
end
ax1.XTickLabel = xTickLabel;

xLim = xlim;
xAxisLength = xLim(2)-xLim(1);
yAxisLength = yLim(2)-yLim(1);
textStr = ['t-test 1v2v3v4v1 p-val:   ' num2str(stats.pval_ttest(1,ismember(stats.area1,areas{1}) & ismember(stats.area2,areas{2}))) '   '...
  num2str(stats.pval_ttest(1,ismember(stats.area1,areas{2}) & ismember(stats.area2,areas{3}))) '   '...
  num2str(stats.pval_ttest(1,ismember(stats.area1,areas{3}) & ismember(stats.area2,areas{4}))) '   '...
  num2str(stats.pval_ttest(1,ismember(stats.area1,areas{1}) & ismember(stats.area2,areas{4})))];
text(xLim(2)-xAxisLength*0.85, yLim(2)-yAxisLength*0.05, textStr, 'FontSize',20);

if strcmp(scaleType, 'log')
  set(gca,'Yscale','log')
end
end


function [p, area1, area2] = ttestGroup(areaNames, data, condition)

if nargin < 3
  condition = 1;
end

areaCombos = nchoosek(1:numel(areaNames), 2);
nCombos = size(areaCombos,1);
area1 = areaNames(areaCombos(:,1));
area2 = areaNames(areaCombos(:,2));
p = zeros(1,nCombos);
for iCombo = 1:nCombos
  areaCode1 = determineArea(areaNames{areaCombos(iCombo,1)});
  areaCode2 = determineArea(areaNames{areaCombos(iCombo,2)});
  [~,p(iCombo)] = ttest2(data{condition}{areaCode1}, data{condition}{areaCode2});
end
end


function anovaOutput = runAnova(scatterGroups, areaGroups)

anovaData = [];
for iGroup = 1:numel(scatterGroups)
  anovaData = [anovaData scatterGroups{iGroup}']; %#ok<*AGROW>
  for iElement = 1:numel(scatterGroups{iGroup})
    if iGroup == 1 && iElement == 1
      anovaAreaFactor = areaGroups(iGroup);
    else
      anovaAreaFactor = [anovaAreaFactor areaGroups{iGroup}];
    end
  end
end
[~, anovaOutput] = anova1(anovaData, anovaAreaFactor);
end