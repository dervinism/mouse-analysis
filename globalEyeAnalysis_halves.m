% Run this script to assess population firing rate phase and coherence
% correlations between recording halves for area-to-pupil comparisons.
%
% The following data files are produced:
% dataDir\area2pupil\PRsFolder\phaseCorrelationsSubfolder\globalHalfPops_area2pupil.mat
%   containing compiled half phase and coherence data, as well as
%   correlation analyses data.
%
% Phase half-recording correlation figures and summary subplots are saved
%   in dataDir\area2pupil\PRsFolder\phaseCorrelationsSubfolder.

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
    rootFolder = [outputDir filesep area2pupilDir filesep PRsFolder];
  elseif strcmp(subpop, 'positive')
    rootFolder = [outputDir filesep area2pupilDir_positive filesep PRsFolder];
  elseif strcmp(subpop, 'negative')
    rootFolder = [outputDir filesep area2pupilDir_negative filesep PRsFolder];
  end
  animals = animalsOI;
elseif strcmp(repository,'uol')
  if strcmp(subpop, 'all')
    rootFolder = [outputDir filesep area2pupilDir_uol filesep PRsFolder];
  elseif strcmp(subpop, 'positive')
    rootFolder = [outputDir filesep area2pupilDir_uol_positive filesep PRsFolder];
  elseif strcmp(subpop, 'negative')
    rootFolder = [outputDir filesep area2pupilDir_uol_negative filesep PRsFolder];
  end
  animals = animalsUOLOI;
elseif strcmp(repository,'allensdk')
  if strcmp(subpop, 'all')
    rootFolder = [outputDir filesep area2pupilDir_allensdk filesep PRsFolder];
  elseif strcmp(subpop, 'positive')
    rootFolder = [outputDir filesep area2pupilDir_allensdk_positive filesep PRsFolder];
  elseif strcmp(subpop, 'negative')
    rootFolder = [outputDir filesep area2pupilDir_allensdk_negative filesep PRsFolder];
  end
  animals = animalsAllensdk;
  conditions = {'awake'};
end
mainFolderPhase = [rootFolder filesep phaseCorrelationsSubfolder];
mainFolderCoh = [rootFolder filesep coherenceCorrelationsSubfolder];

individualPhaseGraphs = false;
summaryPhaseGraphs = true;
individualCohGraphs = false;
summaryCohGraphs = true;


%% EXTRACT AND COMPILE PHASE INFO BASED ON RECORDING HALVES
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
    if animal == 1 || ~exist('areaCoh1FOIindividual', 'var')
      areaCoh1FOIindividual = {};
      areaCoh2FOIindividual = {};
      areaPhase1FOIindividual = {};
      areaPhase2FOIindividual = {};
      areaCoh1FullIndividual = {};
      areaCoh2FullIndividual = {};
      areaCohConfU1FullIndividual = {};
      areaCohConfL1FullIndividual = {};
      areaCohConfU2FullIndividual = {};
      areaCohConfL2FullIndividual = {};
      areaPhase1FullIndividual = {};
      areaPhase2FullIndividual = {};
      areaPhaseConfU1FullIndividual = {};
      areaPhaseConfL1FullIndividual = {};
      areaPhaseConfU2FullIndividual = {};
      areaPhaseConfL2FullIndividual = {};
      areaFreqFullIndividual = {};
      areaCoh1FullInterpIndividual = {};
      areaCoh2FullInterpIndividual = {};
      areaCohConfU1FullInterpIndividual = {};
      areaCohConfL1FullInterpIndividual = {};
      areaCohConfU2FullInterpIndividual = {};
      areaCohConfL2FullInterpIndividual = {};
      areaPhase1FullInterpIndividual = {};
      areaPhase2FullInterpIndividual = {};
      areaPhaseConfU1FullInterpIndividual = {};
      areaPhaseConfL1FullInterpIndividual = {};
      areaPhaseConfU2FullInterpIndividual = {};
      areaPhaseConfL2FullInterpIndividual = {};
      areaFreqFullInterpIndividual = {};
      for iCond = 1:numel(conditions)
        areaCoh1FOIindividualCond = {};
        areaCoh2FOIindividualCond = {};
        areaPhase1FOIindividualCond = {};
        areaPhase2FOIindividualCond = {};
        areaCoh1FullIndividualCond = {};
        areaCoh2FullIndividualCond = {};
        areaCohConfU1FullIndividualCond = {};
        areaCohConfL1FullIndividualCond = {};
        areaCohConfU2FullIndividualCond = {};
        areaCohConfL2FullIndividualCond = {};
        areaPhase1FullIndividualCond = {};
        areaPhase2FullIndividualCond = {};
        areaPhaseConfU1FullIndividualCond = {};
        areaPhaseConfL1FullIndividualCond = {};
        areaPhaseConfU2FullIndividualCond = {};
        areaPhaseConfL2FullIndividualCond = {};
        areaFreqFullIndividualCond = {};
        areaCoh1FullInterpIndividualCond = {};
        areaCoh2FullInterpIndividualCond = {};
        areaCohConfU1FullInterpIndividualCond = {};
        areaCohConfL1FullInterpIndividualCond = {};
        areaCohConfU2FullInterpIndividualCond = {};
        areaCohConfL2FullInterpIndividualCond = {};
        areaPhase1FullInterpIndividualCond = {};
        areaPhase2FullInterpIndividualCond = {};
        areaPhaseConfU1FullInterpIndividualCond = {};
        areaPhaseConfL1FullInterpIndividualCond = {};
        areaPhaseConfU2FullInterpIndividualCond = {};
        areaPhaseConfL2FullInterpIndividualCond = {};
        areaFreqFullInterpIndividualCond = {};
        for iArea = 1:numel(areas)
          areaCoh1FOIindividualCond{iArea} = []; %#ok<*SAGROW>
          areaCoh2FOIindividualCond{iArea} = [];
          areaPhase1FOIindividualCond{iArea} = [];
          areaPhase2FOIindividualCond{iArea} = [];
          areaCoh1FullIndividualCond{iArea} = {};
          areaCoh2FullIndividualCond{iArea} = {};
          areaCohConfU1FullIndividualCond{iArea} = {};
          areaCohConfL1FullIndividualCond{iArea} = {};
          areaCohConfU2FullIndividualCond{iArea} = {};
          areaCohConfL2FullIndividualCond{iArea} = {};
          areaPhase1FullIndividualCond{iArea} = {};
          areaPhase2FullIndividualCond{iArea} = {};
          areaPhaseConfU1FullIndividualCond{iArea} = {};
          areaPhaseConfL1FullIndividualCond{iArea} = {};
          areaPhaseConfU2FullIndividualCond{iArea} = {};
          areaPhaseConfL2FullIndividualCond{iArea} = {};
          areaFreqFullIndividualCond{iArea} = [];
          areaCoh1FullInterpIndividualCond{iArea} = [];
          areaCoh2FullInterpIndividualCond{iArea} = [];
          areaCohConfU1FullInterpIndividualCond{iArea} = [];
          areaCohConfL1FullInterpIndividualCond{iArea} = [];
          areaCohConfU2FullInterpIndividualCond{iArea} = [];
          areaCohConfL2FullInterpIndividualCond{iArea} = [];
          areaPhase1FullInterpIndividualCond{iArea} = [];
          areaPhase2FullInterpIndividualCond{iArea} = [];
          areaPhaseConfU1FullInterpIndividualCond{iArea} = [];
          areaPhaseConfL1FullInterpIndividualCond{iArea} = [];
          areaPhaseConfU2FullInterpIndividualCond{iArea} = [];
          areaPhaseConfL2FullInterpIndividualCond{iArea} = [];
          areaFreqFullInterpIndividualCond{iArea} = [];
        end
        areaCoh1FOIindividual{iCond} = areaCoh1FOIindividualCond;
        areaCoh2FOIindividual{iCond} = areaCoh2FOIindividualCond;
        areaPhase1FOIindividual{iCond} = areaPhase1FOIindividualCond;
        areaPhase2FOIindividual{iCond} = areaPhase2FOIindividualCond;
        areaCoh1FullIndividual{iCond} = areaCoh1FullIndividualCond;
        areaCoh2FullIndividual{iCond} = areaCoh2FullIndividualCond;
        areaCohConfU1FullIndividual{iCond} = areaCohConfU1FullIndividualCond;
        areaCohConfL1FullIndividual{iCond} = areaCohConfL1FullIndividualCond;
        areaCohConfU2FullIndividual{iCond} = areaCohConfU2FullIndividualCond;
        areaCohConfL2FullIndividual{iCond} = areaCohConfL2FullIndividualCond;
        areaPhase1FullIndividual{iCond} = areaPhase1FullIndividualCond;
        areaPhase2FullIndividual{iCond} = areaPhase2FullIndividualCond;
        areaPhaseConfU1FullIndividual{iCond} = areaPhaseConfU1FullIndividualCond;
        areaPhaseConfL1FullIndividual{iCond} = areaPhaseConfL1FullIndividualCond;
        areaPhaseConfU2FullIndividual{iCond} = areaPhaseConfU2FullIndividualCond;
        areaPhaseConfL2FullIndividual{iCond} = areaPhaseConfL2FullIndividualCond;
        areaFreqFullIndividual{iCond} = areaFreqFullInterpIndividualCond;
        areaCoh1FullInterpIndividual{iCond} = areaCoh1FullInterpIndividualCond;
        areaCoh2FullInterpIndividual{iCond} = areaCoh2FullInterpIndividualCond;
        areaCohConfU1FullInterpIndividual{iCond} = areaCohConfU1FullInterpIndividualCond;
        areaCohConfL1FullInterpIndividual{iCond} = areaCohConfL1FullInterpIndividualCond;
        areaCohConfU2FullInterpIndividual{iCond} = areaCohConfU2FullInterpIndividualCond;
        areaCohConfL2FullInterpIndividual{iCond} = areaCohConfL2FullInterpIndividualCond;
        areaPhase1FullInterpIndividual{iCond} = areaPhase1FullInterpIndividualCond;
        areaPhase2FullInterpIndividual{iCond} = areaPhase2FullInterpIndividualCond;
        areaPhaseConfU1FullInterpIndividual{iCond} = areaPhaseConfU1FullInterpIndividualCond;
        areaPhaseConfL1FullInterpIndividual{iCond} = areaPhaseConfL1FullInterpIndividualCond;
        areaPhaseConfU2FullInterpIndividual{iCond} = areaPhaseConfU2FullInterpIndividualCond;
        areaPhaseConfL2FullInterpIndividual{iCond} = areaPhaseConfL2FullInterpIndividualCond;
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
          
          freq = dbStruct.popData.pupil.popData.coh_halves_freq;
          freqPSD = dbStruct.popData.pupil.popData.psd_halves_freq;
          if sum(~isnan(freq)) < sum(~isnan(freqPSD))
            freq = freqPSD;
          end
          
          % Load coherence and phase values for the 1st half of the recording
          coh = dbStruct.popData.pupil.popData.coh_halves(1,:);
          coh_conf = dbStruct.popData.pupil.popData.coh_conf_halves(1,:);
          rateadjust_kappa = dbStruct.popData.pupil.popData.rateadjust_kappa_halves(1,:);
          phase1 = dbStruct.popData.pupil.popData.phase_halves(1,:);
          phase_confU = dbStruct.popData.pupil.popData.phase_conf_halves(1,:);
          phase_confL = dbStruct.popData.pupil.popData.phase_conf_halves(2,:);
          [phase1, phaseConfU1, phaseConfL1, coh1, cohConfU1, cohConfL1] = correctPhaseCoh(phase1,...
            phase_confU, phase_confL, coh, coh_conf, rateadjust_kappa);
%           [~, phaseConfU1, phaseConfL1, coh1, cohConfU1, cohConfL1] = correctPhaseCoh(phase1,...
%             phase_confU, phase_confL, coh, coh_conf, rateadjust_kappa);
          
          % Obtain and store phase and coherence values for FOI for the 1st half of the recording
          if sum(~isnan(phase1)) && sum(~isnan(coh1))
            [phase1FOI, coh1FOI] = phaseCohFOI(FOI, freq, phase1, coh1, [cohConfU1; cohConfL1]);
          else
            phase1FOI = NaN(size(FOI));
            coh1FOI = NaN(size(FOI));
          end
          
          % Load coherence and phase values for the 2nd half of the recording
          coh = dbStruct.popData.pupil.popData.coh_halves(2,:);
          coh_conf = dbStruct.popData.pupil.popData.coh_conf_halves(2,:);
          rateadjust_kappa = dbStruct.popData.pupil.popData.rateadjust_kappa_halves(2,:);
          phase2 = dbStruct.popData.pupil.popData.phase_halves(2,:);
          phase_confU = dbStruct.popData.pupil.popData.phase_conf_halves(3,:);
          phase_confL = dbStruct.popData.pupil.popData.phase_conf_halves(4,:);
          [phase2, phaseConfU2, phaseConfL2, coh2, cohConfU2, cohConfL2] = correctPhaseCoh(phase2,...
            phase_confU, phase_confL, coh, coh_conf, rateadjust_kappa);
%           [~, phaseConfU2, phaseConfL2, coh2, cohConfU2, cohConfL2] = correctPhaseCoh(phase2,...
%             phase_confU, phase_confL, coh, coh_conf, rateadjust_kappa);
          
          % Obtain and store phase and coherence values for FOI for the 2nd half of the recording
          if sum(~isnan(phase2)) && sum(~isnan(coh2))
            [phase2FOI, coh2FOI] = phaseCohFOI(FOI, freq, phase2, coh2, [cohConfU2; cohConfL2]);
          else
            phase2FOI = NaN(size(FOI));
            coh2FOI = NaN(size(FOI));
          end
          
          % Store the phase and coherence values
          areaPhase1FOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaPhase1FOIindividual{iCondPlusAll}{iAreaPlusAll}; phase1FOI];
          areaPhase2FOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaPhase1FOIindividual{iCondPlusAll}{iAreaPlusAll}; phase2FOI];
          areaCoh1FOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaCoh1FOIindividual{iCondPlusAll}{iAreaPlusAll}; coh1FOI];
          areaCoh2FOIindividual{iCondPlusAll}{iAreaPlusAll} = [areaCoh1FOIindividual{iCondPlusAll}{iAreaPlusAll}; coh2FOI];
          areaCoh1FullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCoh1FullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = coh1;
          areaCoh2FullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCoh2FullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = coh2;
          areaCohConfU1FullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCohConfU1FullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = cohConfU1;
          areaCohConfU2FullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCohConfU2FullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = cohConfU2;
          areaCohConfL1FullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCohConfL1FullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = cohConfL1;
          areaCohConfL2FullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaCohConfL2FullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = cohConfL2;
          areaPhase1FullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhase1FullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = phase1;
          areaPhase2FullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhase2FullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = phase2;
          areaPhaseConfU1FullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseConfU1FullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = phaseConfU1;
          areaPhaseConfU2FullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseConfU2FullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = phaseConfU2;
          areaPhaseConfL1FullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseConfL1FullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = phaseConfL1;
          areaPhaseConfL2FullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaPhaseConfL2FullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = phaseConfL2;
          areaFreqFullIndividual{iCondPlusAll}{iAreaPlusAll}{numel(areaFreqFullIndividual{iCondPlusAll}{iAreaPlusAll})+1} = freq;
        end
      end
    end
  end
  
  % Store interpolated phase and coherence values
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
      areaPhase1FullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaPhaseConfU1FullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaPhaseConfL1FullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaCoh1FullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaCohConfU1FullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaCohConfL1FullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaPhase2FullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaPhaseConfU2FullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaPhaseConfL2FullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaCoh2FullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaCohConfU2FullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      areaCohConfL2FullInterpIndividual{iCond}{iArea} = NaN(nRec,numel(freqCombined));
      for iRec = 1:nRec
        if ~isempty(areaPhase1FullIndividual{iCond}{iArea}) && sum(~isnan(areaPhase1FullIndividual{iCond}{iArea}{iRec}))
          areaPhase1FullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaPhase1FullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaPhaseConfU1FullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaPhaseConfU1FullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaPhaseConfL1FullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaPhaseConfL1FullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaCoh1FullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaCoh1FullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaCohConfU1FullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaCohConfU1FullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaCohConfL1FullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaCohConfL1FullIndividual{iCond}{iArea}{iRec}, freqCombined);
        end
        if ~isempty(areaPhase2FullIndividual{iCond}{iArea}) && sum(~isnan(areaPhase2FullIndividual{iCond}{iArea}{iRec}))
          areaPhase2FullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaPhase2FullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaPhaseConfU2FullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaPhaseConfU2FullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaPhaseConfL2FullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaPhaseConfL2FullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaCoh2FullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaCoh2FullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaCohConfU2FullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaCohConfU2FullIndividual{iCond}{iArea}{iRec}, freqCombined);
          areaCohConfL2FullInterpIndividual{iCond}{iArea}(iRec,:) = interp1(areaFreqFullIndividual{iCond}{iArea}{iRec},...
            areaCohConfL2FullIndividual{iCond}{iArea}{iRec}, freqCombined);
        end
      end
    end
  end
  areaFreqFullInterpIndividual = freqCombined;
end


%% SAVE OR LOAD THE DATA
filename = [rootFolder filesep 'globalHalfPRs_area2pupil.mat']; %#ok<*UNRCH>
if fullRun
  if ~exist(rootFolder, 'file')
    mkdir(rootFolder);
  end
  save(filename, 'conditions','areas','FOI','areaCoh1FOIindividual','areaCoh2FOIindividual',...
    'areaPhase1FOIindividual','areaPhase2FOIindividual','areaCoh1FullIndividual',...
    'areaCoh2FullIndividual','areaCohConfU1FullIndividual','areaCohConfL1FullIndividual',...
    'areaCohConfU2FullIndividual','areaCohConfL2FullIndividual','areaPhase1FullIndividual',...
    'areaPhase2FullIndividual','areaPhaseConfU1FullIndividual','areaPhaseConfL1FullIndividual',...
    'areaPhaseConfU2FullIndividual','areaPhaseConfL2FullIndividual','areaFreqFullInterpIndividual',...
    'areaCoh1FullInterpIndividual','areaCoh2FullInterpIndividual','areaCohConfU1FullInterpIndividual',...
    'areaCohConfL1FullInterpIndividual','areaCohConfU2FullInterpIndividual',...
    'areaCohConfL2FullInterpIndividual','areaPhase1FullInterpIndividual',...
    'areaPhase2FullInterpIndividual','areaPhaseConfU1FullInterpIndividual',...
    'areaPhaseConfL1FullInterpIndividual','areaPhaseConfU2FullInterpIndividual',...
    'areaPhaseConfL2FullInterpIndividual','areaFreqFullInterpIndividual');
else
  load(filename);
end


%% PHASE CORRELATIONS
options = struct();
options.figFolder = mainFolderPhase;
options.figSize = figSize;
options.figTitle = 'PUPIL';
options.xLabel = '1st half';
options.yLabel = '2nd half';
options.phaseLim = [phaseLim(1) 2];
options.individualGraphs = individualPhaseGraphs;
options.summaryGraphs = summaryPhaseGraphs;
options.iAreasOI = iAreasOI;
[rPhaseHalves, pvalPhaseHalves, nPhaseHalves] = halfPhaseCorr(areas, conditions(1:min([numel(conditions) 2])),...
  areaFreqFullInterpIndividual, areaPhase1FullInterpIndividual, areaPhase2FullInterpIndividual, options);


%% COHERENCE CORRELATIONS
options = struct();
options.figFolder = mainFolderCoh;
options.figSize = figSize;
options.figTitle = 'PUPIL';
options.xLabel = '1st half';
options.yLabel = '2nd half';
options.cohLim = [0 1];
options.individualGraphs = individualCohGraphs;
options.summaryGraphs = summaryCohGraphs;
options.iAreasOI = iAreasOI;
[rCohHalves, pvalCohHalves, nCohHalves] = halfCohCorr(areas, conditions(1:min([numel(conditions) 2])),...
  areaFreqFullInterpIndividual, areaCoh1FullInterpIndividual, areaCoh2FullInterpIndividual, options);


%% SAVE THE RESULTS OF CORRELATION ANALYSES
save(filename, 'rPhaseHalves','pvalPhaseHalves','nPhaseHalves','rCohHalves','pvalCohHalves','nCohHalves', '-append');