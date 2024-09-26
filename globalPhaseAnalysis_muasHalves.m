% Run this script to assess unit phase and coherence correlations between
% recording halves for cross-area comparisons.
%
% The following data files are produced:
% dataDir\caDir\MUAsFolder\phaseCorrelationsSubfolder\globalHalfMUAs_ca.mat or
%   dataDir\caDir\MUAsFolder\phaseCorrelationsSubfolder\globalHalfMUAs_ca_reverse.mat
%   containing compiled half phase and coherence data, as well as
%   correlation analyses data.
%
% Phase half-recording correlation figures and summary subplots are saved
%   in dataDir\caDir\MUAsFolder\phaseCorrelationsSubfolder or
%   dataDir\caDir\qualityMUAsFolder\phaseCorrelationsSubfolder.

clearvars -except repository subpop reverse qualityCheck allData fullRun includeRuns


%% INITIALISE PARAMETERS
params
lists

if ~exist('repository', 'var')
  repository = 'positive';
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
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep caDir_uol_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep caDir_uol_negative];
  end
  animals = animalsUOLOI;
elseif strcmp(repository,'allensdk')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep caDir_allensdk];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep caDir_allensdk_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep caDir_allensdk_negative];
  end
  animals = animalsAllensdk;
  conditions = {'awake'};
end
mainFolderPhase = [rootFolder filesep MUAsFolder filesep phaseCorrelationsSubfolder];
mainFolderCoh = [rootFolder filesep MUAsFolder filesep coherenceCorrelationsSubfolder];
areas = areas2compare;

individualPhaseGraphs = false;
summaryPhaseGraphs = false;
individualCohGraphs = false;
summaryCohGraphs = false;


%% EXTRACT AND COMPILE UNIT PHASE AND COHERENCE INFO BASED ON RECORDING HALVES
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
      if ~isfield(dbStruct_ca, 'popData') || ~isfield(dbStruct_ca.popData, 'phaseCohHalves') ||...
          isempty(dbStruct_ca.popData.phaseCohHalves)
        continue
      end
      
      % Test for exceptions
      if exceptionTest(except, seriesName1, seriesName2)
        continue
      end
      
      % Determine if population rate > 0
      breakClause = true;
      if strcmp(subpop, 'all')
        if isfield(dataStruct.seriesData, [animals{animal} '_s' seriesName1])
          [breakClause, spkDBmfr] = firingRateTest(sum(dataStruct.seriesData.([animals{animal} '_s' seriesName1]).popData.MUAsAll,1),...
            dataStruct.seriesData.([animals{animal} '_s' seriesName1]).conf.samplingParams.srData);
        end
      elseif strcmp(subpop, 'positive')
        if isfield(dataStruct.seriesData_positive, [animals{animal} '_s' seriesName1])
          [breakClause, spkDBmfr] = firingRateTest(sum(dataStruct.seriesData_positive.([animals{animal} '_s' seriesName1]).popData.MUAsAll,1),...
            dataStruct.seriesData_positive.([animals{animal} '_s' seriesName1]).conf.samplingParams.srData);
        end
      elseif strcmp(subpop, 'negative')
        if isfield(dataStruct.seriesData_negative, [animals{animal} '_s' seriesName1])
          [breakClause, spkDBmfr] = firingRateTest(sum(dataStruct.seriesData_negative.([animals{animal} '_s' seriesName1]).popData.MUAsAll,1),...
            dataStruct.seriesData_negative.([animals{animal} '_s' seriesName1]).conf.samplingParams.srData);
        end
      end
      if breakClause
        continue
      end
      breakClause = true;
      if strcmp(subpop, 'all')
        if isfield(dataStruct.seriesData, [animals{animal} '_s' seriesName2])
          [breakClause, PRmfr] = firingRateTest(sum(dataStruct.seriesData.([animals{animal} '_s' seriesName2]).popData.MUAsAll,1),...
            dataStruct.seriesData.([animals{animal} '_s' seriesName2]).conf.samplingParams.srData);
        end
      elseif strcmp(subpop, 'positive')
        if isfield(dataStruct.seriesData_positive, [animals{animal} '_s' seriesName2])
          [breakClause, PRmfr] = firingRateTest(sum(dataStruct.seriesData_positive.([animals{animal} '_s' seriesName2]).popData.MUAsAll,1),...
            dataStruct.seriesData_positive.([animals{animal} '_s' seriesName2]).conf.samplingParams.srData);
        end
      elseif strcmp(subpop, 'negative')
        if isfield(dataStruct.seriesData_negative, [animals{animal} '_s' seriesName2])
          [breakClause, PRmfr] = firingRateTest(sum(dataStruct.seriesData_negative.([animals{animal} '_s' seriesName2]).popData.MUAsAll,1),...
            dataStruct.seriesData_negative.([animals{animal} '_s' seriesName2]).conf.samplingParams.srData);
        end
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
          
          % Load and store phase and coherence data for units based on half-recordings
          for u = 1:numel(units)
            
            % Get phase and coherence values for the 1st half of the recording
            coh1 = dbStruct_ca.phaseCoh{u}.coh_halves(1,:);
            cohConf1 = dbStruct_ca.phaseCoh{u}.coh_conf_halves(1,:);
            rateadjust_kappa1 = dbStruct_ca.phaseCoh{u}.rateadjust_kappa_halves(1,:);
            phase1 = dbStruct_ca.phaseCoh{u}.phase_halves(1,:);
            phaseConfU1 = dbStruct_ca.phaseCoh{u}.phase_conf_halves(1,:);
            phaseConfL1 = dbStruct_ca.phaseCoh{u}.phase_conf_halves(2,:);
            [phase1, phaseConfU1, phaseConfL1, coh1, cohConfU1, cohConfL1] = correctPhaseCoh(phase1,...
              phaseConfU1, phaseConfL1, coh1, cohConf1, rateadjust_kappa1);
%             [~, phaseConfU1, phaseConfL1, coh1, cohConfU1, cohConfL1] = correctPhaseCoh(phase1,...
%               phaseConfU1, phaseConfL1, coh1, cohConf1, rateadjust_kappa1);
            
            % Obtain and store phase and coherence values for FOI for the 1st half of the recording
            freq = dbStruct_ca.phaseCoh{u}.coh_halves_freq;
            if sum(isnan(freq)) == numel(freq)
              freq = dbStruct_ca.phaseCoh{u}.freq(end-numel(freq)+1:end);
            end
            if sum(~isnan(phase1)) && sum(~isnan(coh1))
              [phaseFOI1, cohFOI1] = phaseCohFOI(FOI, freq, phase1, coh1, [cohConfU1; cohConfL1]);
            else
              phaseFOI1 = NaN(size(FOI));
              cohFOI1 = NaN(size(FOI));
            end
            
            % Get phase and coherence values for the 2nd half of the recording
            coh2 = dbStruct_ca.phaseCoh{u}.coh_halves(2,:);
            cohConf2 = dbStruct_ca.phaseCoh{u}.coh_conf_halves(2,:);
            rateadjust_kappa2 = dbStruct_ca.phaseCoh{u}.rateadjust_kappa_halves(2,:);
            phase2 = dbStruct_ca.phaseCoh{u}.phase_halves(2,:);
            phaseConfU2 = dbStruct_ca.phaseCoh{u}.phase_conf_halves(3,:);
            phaseConfL2 = dbStruct_ca.phaseCoh{u}.phase_conf_halves(4,:);
            [phase2, phaseConfU2, phaseConfL2, coh2, cohConfU2, cohConfL2] = correctPhaseCoh(phase2,...
              phaseConfU2, phaseConfL2, coh2, cohConf2, rateadjust_kappa2);
%             [~, phaseConfU2, phaseConfL2, coh2, cohConfU2, cohConfL2] = correctPhaseCoh(phase2,...
%               phaseConfU2, phaseConfL2, coh2, cohConf2, rateadjust_kappa2);
            
            % Obtain and store phase and coherence values for FOI for the 2nd half of the recording
            if sum(isnan(freq)) && (~sum(isnan(phase1)) || ~sum(isnan(phase2))) &&...
                (~exist('prevAnimal','var') || ~strcmpi(prevAnimal, animals{animal}))
              prevAnimal = animals{animal};
              disp(' ');
              disp([animals{animal} ' is missing frequency vector']);
              disp(' ');
            end
            if sum(~isnan(phase2)) && sum(~isnan(coh2))
              [phaseFOI2, cohFOI2] = phaseCohFOI(FOI, freq, phase2, coh2, [cohConfU2; cohConfL2]);
            else
              phaseFOI2 = NaN(size(FOI));
              cohFOI2 = NaN(size(FOI));
            end
            
            % Store the phase and coherence values
            areaPhase1FOIindividual{iCondPlusAll}{area} = [areaPhase1FOIindividual{iCondPlusAll}{area}; phaseFOI1];
            areaPhase2FOIindividual{iCondPlusAll}{area} = [areaPhase1FOIindividual{iCondPlusAll}{area}; phaseFOI2];
            areaCoh1FOIindividual{iCondPlusAll}{area} = [areaCoh1FOIindividual{iCondPlusAll}{area}; cohFOI1];
            areaCoh2FOIindividual{iCondPlusAll}{area} = [areaCoh1FOIindividual{iCondPlusAll}{area}; cohFOI2];
            areaCoh1FullIndividual{iCondPlusAll}{area}{numel(areaCoh1FullIndividual{iCondPlusAll}{area})+1} = coh1;
            areaCoh2FullIndividual{iCondPlusAll}{area}{numel(areaCoh2FullIndividual{iCondPlusAll}{area})+1} = coh2;
            areaCohConfU1FullIndividual{iCondPlusAll}{area}{numel(areaCohConfU1FullIndividual{iCondPlusAll}{area})+1} = cohConfU1;
            areaCohConfU2FullIndividual{iCondPlusAll}{area}{numel(areaCohConfU2FullIndividual{iCondPlusAll}{area})+1} = cohConfU2;
            areaCohConfL1FullIndividual{iCondPlusAll}{area}{numel(areaCohConfL1FullIndividual{iCondPlusAll}{area})+1} = cohConfL1;
            areaCohConfL2FullIndividual{iCondPlusAll}{area}{numel(areaCohConfL2FullIndividual{iCondPlusAll}{area})+1} = cohConfL2;
            areaPhase1FullIndividual{iCondPlusAll}{area}{numel(areaPhase1FullIndividual{iCondPlusAll}{area})+1} = phase1;
            areaPhase2FullIndividual{iCondPlusAll}{area}{numel(areaPhase2FullIndividual{iCondPlusAll}{area})+1} = phase2;
            areaPhaseConfU1FullIndividual{iCondPlusAll}{area}{numel(areaPhaseConfU1FullIndividual{iCondPlusAll}{area})+1} = phaseConfU1;
            areaPhaseConfU2FullIndividual{iCondPlusAll}{area}{numel(areaPhaseConfU2FullIndividual{iCondPlusAll}{area})+1} = phaseConfU2;
            areaPhaseConfL1FullIndividual{iCondPlusAll}{area}{numel(areaPhaseConfL1FullIndividual{iCondPlusAll}{area})+1} = phaseConfL1;
            areaPhaseConfL2FullIndividual{iCondPlusAll}{area}{numel(areaPhaseConfL2FullIndividual{iCondPlusAll}{area})+1} = phaseConfL2;
            areaFreqFullIndividual{iCondPlusAll}{area}{numel(areaFreqFullIndividual{iCondPlusAll}{area})+1} = freq;
            
          end
        end
      end
    end
  end
  
  areas = areasReverse;
  
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
if reverse
  filename = [fileparts(mainFolderPhase) filesep 'globalHalfMUAs_ca_reverse.mat']; %#ok<*UNRCH>
else
  filename = [fileparts(mainFolderPhase) filesep 'globalHalfMUAs_ca.mat'];
end
if fullRun
  if ~exist(fileparts(mainFolderPhase), 'file')
    mkdir(fileparts(mainFolderPhase));
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
iAreas2compareOI = find(ismember(areas,areas2compareCritical));


%% PHASE CORRELATIONS
options = struct();
options.figFolder = mainFolderPhase;
options.figSize = figSize;
options.figTitle = 'SPIKING';
options.xLabel = '1st half';
options.yLabel = '2nd half';
options.phaseLim = phaseLim;
options.individualGraphs = individualPhaseGraphs;
options.summaryGraphs = summaryPhaseGraphs;
options.iAreasOI = iAreas2compareOI;
[rPhaseHalves, pvalPhaseHalves, nPhaseHalves] = halfPhaseCorr(areas, conditions(1:min([numel(conditions) 2])),...
  areaFreqFullInterpIndividual, areaPhase1FullInterpIndividual, areaPhase2FullInterpIndividual, options);


%% COHERENCE CORRELATIONS
options = struct();
options.figFolder = mainFolderCoh;
options.figSize = figSize;
options.figTitle = 'SPIKING';
options.xLabel = '1st half';
options.yLabel = '2nd half';
options.cohLim = [0 1];
options.individualGraphs = individualCohGraphs;
options.summaryGraphs = summaryCohGraphs;
options.iAreasOI = iAreas2compareOI;
[rCohHalves, pvalCohHalves, nCohHalves] = halfCohCorr(areas, conditions(1:min([numel(conditions) 2])),...
  areaFreqFullInterpIndividual, areaCoh1FullInterpIndividual, areaCoh2FullInterpIndividual, options);


%% SAVE THE RESULTS OF CORRELATION ANALYSES
save(filename, 'rPhaseHalves','pvalPhaseHalves','nPhaseHalves','rCohHalves','pvalCohHalves','nCohHalves', '-append');