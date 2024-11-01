% Run this script to produce population rate phase and coherence frequency
% profiles and phase frequency histograms for cross-area comparisons.
%
% The following data files are produced:
% dataDir\caDir\PRsFolder\area2areaCohMats.mat and
%   dataDir\caDir\PRsFolder\area2areaCohMats_reverse.mat both contain phase and
%   coherence frequency profiles, coherence matrices, and recording counts.
%   The reverse file has area comparisons reversed.
% dataDir\caDir\PRsFolder\histosSubfolder\phaseHistos.mat and
%   dataDir\caDir\PRsFolder\histosSubfolder\phaseHistos_reverse.mat both contain phase
%   frequency profiles, phase frequency histograms, and statistical test
%   results.
%
% Phase frequency profile figures are saved in dataDir\caDir\PRsFolder\phaseFrequencyProfilesSubfolder.
% Coherence frequency profile figures are saved in dataDir\caDir\PRsFolder\coherenceFrequencyProfilesSubfolder.
% Phase sum prediction figures are saved in dataDir\caDir\PRsFolder\phaseSumPredictionsSubfolder.
% Individual phase frequency histogram figures and summary suplots are saved in dataDir\caDir\PRsFolder\histosSubfolder.
% Phase frequency maps are saved in dataDir\caDir\PRsFolder\mapsSubfolder.

clearvars -except repository subpop reverse qualityCheck allData fullRun includeRuns


%% INITIALISE PARAMETERS
params
lists

if ~exist('repository', 'var')
  repository = 'uol'; % 'uol' | 'allensdk'
end
if ~exist('subpop', 'var')
  subpop = 'all'; % 'all' | 'positive' | 'negative'
end
if ~exist('fullRun', 'var')
  fullRun = true;
end
if ~exist('reverse', 'var')
  reverse = false;
end

outputDir = [outputDir filesep includeRuns];
if strcmp(repository,'uol')
  dataDir = [dataDir_local filesep '001_uol'];
elseif strcmp(repository,'allensdk')
  dataDir = [dataDir_local filesep '002_allen'];
end
if strcmp(repository,'all')
  if strcmp(subpop, 'all')
    mainFolder = [outputDir filesep caDir filesep PRsFolder];
  elseif strcmp(subpop, 'positive')
    mainFolder = [outputDir filesep caDir_positive filesep PRsFolder];
  elseif strcmp(subpop, 'negative')
    mainFolder = [outputDir filesep caDir_negative filesep PRsFolder];
  end
  animals = animalsOI;
  xLim = freqLimUOL;
elseif strcmp(repository,'uol')
  if strcmp(subpop, 'all')
    mainFolder = [outputDir filesep caDir_uol filesep PRsFolder];
  elseif strcmp(subpop, 'positive')
    mainFolder = [outputDir filesep caDir_uol_positive filesep PRsFolder];
  elseif strcmp(subpop, 'negative')
    mainFolder = [outputDir filesep caDir_uol_negative filesep PRsFolder];
  end
  animals = animalsUOLOI;
  xLim = freqLimUOL;
elseif strcmp(repository,'allensdk')
  if strcmp(subpop, 'all')
    mainFolder = [outputDir filesep caDir_allensdk filesep PRsFolder];
  elseif strcmp(subpop, 'positive')
    mainFolder = [outputDir filesep caDir_allensdk_positive filesep PRsFolder];
  elseif strcmp(subpop, 'negative')
    mainFolder = [outputDir filesep caDir_allensdk_negative filesep PRsFolder];
  end
  animals = animalsAllensdk;
  conditions = {'awake'};
  xLim = freqLimAllen;
end
areas = areas2compare;

unitsOnlyPR = false;
drawPhaseProfiles = [true true];
drawCohProfiles = [true true];
drawPhaseSums = false;
drawPhaseHistos = [true true true];
doStats = true;


%% COMPUTE VARIABLES NEEDED FOR DISPLAYING UNIT PHASE AND COHERENCE FREQUENCY PROFILES
if fullRun
  matSize = numel(areasCohMat);
  for animal = 1:numel(animals) % Loop through animals
    disp(['processing animal ' animals{animal}]);
    if exist('allData', 'var')
      dataStruct = allData.(animals{animal}).dataStruct;
    else
      load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
    end
    animalColour = animalColours(animals(animal));
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

    % Initialise storage variables
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

    % Initialise coherence matrices
    if animal == 1 || ~exist('recMatInit', 'var')
      recMatInit = NaN(matSize,matSize,numel(FOI));
      cohMatInit = NaN(matSize,matSize,numel(FOI));
      cohConfMatInit = NaN(matSize,matSize,numel(FOI));
      mfrPR1MatInit = NaN(matSize,matSize,numel(FOI));
      mfrPR2MatInit = NaN(matSize,matSize,numel(FOI));
      for iCond = 1:numel(conditions)
        recMat{iCond} = {};
        cohMat{iCond} = {};
        cohConfUMat{iCond} = {};
        cohConfLMat{iCond} = {};
        mfrPR1Mat{iCond} = {};
        mfrPR2Mat{iCond} = {};
        seriesMat{iCond} = {};
        series1Mat{iCond} = {};
        series2Mat{iCond} = {};
        animalMat{iCond} = {};
      end
    end

    for dbCount = 1:numel(fnsData_ca) % Loop through database entries
      if strcmp(subpop, 'all')
        dbStruct_ca = dataStruct.seriesData_ca.(fnsData_ca{dbCount});
      elseif strcmp(subpop, 'positive')
        dbStruct_ca = dataStruct.seriesData_ca_positive.(fnsData_ca{dbCount});
      elseif strcmp(subpop, 'negative')
        dbStruct_ca = dataStruct.seriesData_ca_negative.(fnsData_ca{dbCount});
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
      if ~strcmpi(prevRec, recording)
        expandMatrices = true;
      end
      if ~unitsOnlyPR
        if ~isfield(dbStruct_ca.popData, 'phaseCohPop') ||...
            isempty(dbStruct_ca.popData.phaseCohPop) % this would produce an error if the population rate was empty for the first series
          continue
        end
      else
        if isempty(dbStruct_ca.popData.phaseCohPopUnitOnly)
          continue
        end
      end

      % Test for exceptions
      if exceptionTest(except, seriesName1, seriesName2)
        continue
      end

      % Determine if population rate > 0
      if strcmp(subpop, 'all')
        if isempty(dataStruct.seriesData.([animals{animal} '_s' seriesName1]))
          continue
        else
          [breakClause, spkDBmfr] = firingRateTest(sum(dataStruct.seriesData.([animals{animal} '_s' seriesName1]).popData.MUAsAll,1),...
            dataStruct.seriesData.([animals{animal} '_s' seriesName1]).conf.samplingParams.srData);
        end
      elseif strcmp(subpop, 'positive')
        if isempty(dataStruct.seriesData_positive.([animals{animal} '_s' seriesName1]))
          continue
        else
          [breakClause, spkDBmfr] = firingRateTest(sum(dataStruct.seriesData_positive.([animals{animal} '_s' seriesName1]).popData.MUAsAll,1),...
            dataStruct.seriesData_positive.([animals{animal} '_s' seriesName1]).conf.samplingParams.srData);
        end
      elseif strcmp(subpop, 'negative')
        if isempty(dataStruct.seriesData_negative.([animals{animal} '_s' seriesName1]))
          continue
        else
          [breakClause, spkDBmfr] = firingRateTest(sum(dataStruct.seriesData_negative.([animals{animal} '_s' seriesName1]).popData.MUAsAll,1),...
            dataStruct.seriesData_negative.([animals{animal} '_s' seriesName1]).conf.samplingParams.srData);
        end
      end
      if breakClause
        continue
      end
      if strcmp(subpop, 'all')
        if isempty(dataStruct.seriesData.([animals{animal} '_s' seriesName2]))
          continue
        else
          [breakClause, PRmfr] = firingRateTest(sum(dataStruct.seriesData.([animals{animal} '_s' seriesName2]).popData.MUAsAll,1),...
            dataStruct.seriesData.([animals{animal} '_s' seriesName2]).conf.samplingParams.srData);
        end
      elseif strcmp(subpop, 'positive')
        if isempty(dataStruct.seriesData_positive.([animals{animal} '_s' seriesName2]))
          continue
        else
          [breakClause, PRmfr] = firingRateTest(sum(dataStruct.seriesData_positive.([animals{animal} '_s' seriesName2]).popData.MUAsAll,1),...
            dataStruct.seriesData_positive.([animals{animal} '_s' seriesName2]).conf.samplingParams.srData);
        end
      elseif strcmp(subpop, 'negative')
        if isempty(dataStruct.seriesData_negative.([animals{animal} '_s' seriesName2]))
          continue
        else
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

      for iComp = 1:numel(comp) % Loop through non-grouped and grouped area comparisons
        area = comp(iComp);

        if numel(conditions) > 1
          condLoop = [iCond numel(conditions)];
        else
          condLoop = iCond;
        end
        for iCondPlusAll = condLoop % Loop through the main and pooled conditions

          % Assign matrices
          if expandMatrices
            recMat{iCondPlusAll}{numel(recMat{iCondPlusAll})+1} = recMatInit;
            cohMat{iCondPlusAll}{numel(cohMat{iCondPlusAll})+1} = cohMatInit;
            cohConfUMat{iCondPlusAll}{numel(cohConfUMat{iCondPlusAll})+1} = cohConfMatInit;
            cohConfLMat{iCondPlusAll}{numel(cohConfLMat{iCondPlusAll})+1} = cohConfMatInit;
            mfrPR1Mat{iCondPlusAll}{numel(mfrPR1Mat{iCondPlusAll})+1} = mfrPR1MatInit;
            mfrPR2Mat{iCondPlusAll}{numel(mfrPR2Mat{iCondPlusAll})+1} = mfrPR2MatInit;
            seriesMat{iCondPlusAll}{numel(seriesMat{iCondPlusAll})+1} = fnsData_ca{dbCount};
            series1Mat{iCondPlusAll}{numel(series1Mat{iCondPlusAll})+1} = seriesName1;
            series2Mat{iCondPlusAll}{numel(series2Mat{iCondPlusAll})+1} = seriesName2;
            animalMat{iCondPlusAll}{numel(animalMat{iCondPlusAll})+1} = animals{animal};
            if iCondPlusAll == 3
              expandMatrices = false;
            end
          end

          % Load coherence of the population rate
          if ~unitsOnlyPR
            coh = dbStruct_ca.popData.phaseCohPop.coh;
            coh_conf = dbStruct_ca.popData.phaseCohPop.coh_conf;
            rateadjust_kappa = dbStruct_ca.popData.phaseCohPop.rateadjust_kappa;
          else
            coh = dbStruct_ca.popData.phaseCohPopUnitOnly.coh;
            coh_conf = dbStruct_ca.popData.phaseCohPopUnitOnly.coh_conf;
            rateadjust_kappa = dbStruct_ca.popData.phaseCohPopUnitOnly.rateadjust_kappa;
          end
          coh_confU = (coh + coh_conf) .* rateadjust_kappa;
          coh_confL = (coh - coh_conf) .* rateadjust_kappa;
          coh_confU(coh_confL <= 0) = NaN;
          coh_confL(coh_confL <= 0) = NaN;
          coh = coh .* rateadjust_kappa;
          coh(isnan(coh_confU) | isnan(coh_confL)) = NaN;

          % Load phase of the population rate
          if ~unitsOnlyPR
            phase = bestUnwrap(dbStruct_ca.popData.phaseCohPop.phase);
            phase_confU = dbStruct_ca.popData.phaseCohPop.phase_confU;
            phase_confL = dbStruct_ca.popData.phaseCohPop.phase_confL;
            phase(isnan(phase_confU) | isnan(phase_confL) | isnan(coh) | isnan(coh_confU) | isnan(coh_confL)) = NaN;
            freq = dbStruct_ca.popData.phaseCohPop.freq;
          else
            phase = bestUnwrap(dbStruct_ca.popData.phaseCohPopUnitOnly.phase);
            phase_confU = dbStruct_ca.popData.phaseCohPopUnitOnly.phase_confU;
            phase_confL = dbStruct_ca.popData.phaseCohPopUnitOnly.phase_confL;
            phase(isnan(phase_confU) | isnan(phase_confL) | isnan(coh) | isnan(coh_confU) | isnan(coh_confL)) = NaN;
            freq = dbStruct_ca.popData.phaseCohPopUnitOnly.freq;
          end

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

          % Determine coherence matrix entry
          [areaName1, areaName2] = comparison2areas(compNames{iComp});
          [~, ~, ~, row] = determineArea(areaName1);
          [~, ~, ~, col] = determineArea(areaName2);

          % Compute coherence matrices
          if ~isempty(row) && ~isempty(col)
            for iF = 1:numel(FOI)
              recMat{iCondPlusAll}{numel(recMat{iCondPlusAll})}(col,row,iF) = 1;
              recMat{iCondPlusAll}{numel(recMat{iCondPlusAll})}(row,col,iF) = 1;
              cohMat{iCondPlusAll}{numel(cohMat{iCondPlusAll})}(col,row,iF) = cohFOI(iF);
              cohMat{iCondPlusAll}{numel(cohMat{iCondPlusAll})}(row,col,iF) = cohFOI(iF);
              cohConfUMat{iCondPlusAll}{numel(cohConfUMat{iCondPlusAll})}(col,row,iF) = cohConfUFOI(iF);
              cohConfUMat{iCondPlusAll}{numel(cohConfUMat{iCondPlusAll})}(row,col,iF) = cohConfUFOI(iF);
              cohConfLMat{iCondPlusAll}{numel(cohConfLMat{iCondPlusAll})}(col,row,iF) = cohConfLFOI(iF);
              cohConfLMat{iCondPlusAll}{numel(cohConfLMat{iCondPlusAll})}(row,col,iF) = cohConfLFOI(iF);
              for iDiagonal = 1:size(cohMat{iCondPlusAll}{numel(cohMat{iCondPlusAll})},1)
                cohMat{iCondPlusAll}{numel(cohMat{iCondPlusAll})}(iDiagonal,iDiagonal,iF) = 1;
                cohConfUMat{iCondPlusAll}{numel(cohConfUMat{iCondPlusAll})}(iDiagonal,iDiagonal,iF) = 1;
                cohConfLMat{iCondPlusAll}{numel(cohConfLMat{iCondPlusAll})}(iDiagonal,iDiagonal,iF) = 1;
              end
            end
            mfrPR1Mat{iCondPlusAll}{numel(mfrPR1Mat{iCondPlusAll})}(:,row,:) = spkDBmfr*ones(matSize,numel(FOI));
            mfrPR2Mat{iCondPlusAll}{numel(mfrPR2Mat{iCondPlusAll})}(row,:,:) = spkDBmfr*ones(matSize,numel(FOI));
            mfrPR2Mat{iCondPlusAll}{numel(mfrPR2Mat{iCondPlusAll})}(col,:,:) = PRmfr*ones(matSize,numel(FOI));
            mfrPR1Mat{iCondPlusAll}{numel(mfrPR1Mat{iCondPlusAll})}(:,col,:) = PRmfr*ones(matSize,numel(FOI));
          end
        end
      end
    end
  end

  % Calculate mean of coherence matrices
  for iCond = 1:numel(conditions)
    cohMatMean{iCond} = zeros(matSize,matSize,numel(FOI));
    tallyMat{iCond} = zeros(matSize,matSize,numel(FOI));
    cohConfUMatMean{iCond} = zeros(matSize,matSize,numel(FOI));
    cohConfLMatMean{iCond} = zeros(matSize,matSize,numel(FOI));
    for iRec = 1:numel(cohMat{iCond})
      cohMatRec = cohMat{iCond}{iRec};
      cohMatRec(isnan(cohMatRec)) = 0;
      cohMatMean{iCond} = cohMatMean{iCond} + cohMatRec;
      tallyMat{iCond} = tallyMat{iCond} + ~isnan(cohMat{iCond}{iRec});
      cohConfUMatRec = cohConfUMat{iCond}{iRec};
      cohConfUMatRec(isnan(cohConfUMatRec)) = 0;
      cohConfUfMatMean{iCond} = cohConfUMatMean{iCond} + cohConfUMatRec;
      cohConfLMatRec = cohConfLMat{iCond}{iRec};
      cohConfLMatRec(isnan(cohConfLMatRec)) = 0;
      cohConfLfMatMean{iCond} = cohConfLMatMean{iCond} + cohConfLMatRec;
    end
    cohMatMean{iCond} = cohMatMean{iCond} ./ tallyMat{iCond};
    cohConfUMatMean{iCond} = cohConfUMatMean{iCond} ./ tallyMat{iCond};
    cohConfLMatMean{iCond} = cohConfLMatMean{iCond} ./ tallyMat{iCond};
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
  filename = [mainFolder filesep 'globalPRs_ca_reverse.mat']; %#ok<*UNRCH>
else
  filename = [mainFolder filesep 'globalPRs_ca.mat'];
end
if fullRun
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  save(filename, 'cohMatMean','cohConfUMatMean','cohConfLMatMean','tallyMat','recMat','cohMat','cohConfUMat','cohConfLMat',...
    'mfrPR1Mat','mfrPR2Mat', 'seriesMat','series1Mat','series2Mat','animalMat','conditions','areas','FOI',...
    'areasCohMat','areaPhaseFOIindividual','areaCohFOIindividual','areaPhaseFullIndividual',...
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
  options.figTitle = 'Population rate phase comparisons: %s %s';
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
  options.figTitle = 'Population rate mean phase comparisons: %s';
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
  options.figTitle = 'Population rate coherence comparisons: %s %s';
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
  options.figTitle = 'Population rate mean coherence comparisons: %s';
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.cohFrequencyProfilesSubfolder = coherenceFrequencyProfilesSubfolder;
  options.freqLim = xLim;
  options.cohLim = [0 1];
  options.iAreasOI = iAreas2compareOI;
  cohFreqProfilePlotMeans(areas, conditions(1:min([numel(conditions) 2])),...
    areaFreqFullInterpIndividual, areaCohFullInterpIndividual, [], figFileName, options);
end


%% PREDICT SUMS OF PHASE FREQUENCY PROFILES
if drawPhaseSums
  options = struct();
  options.figSize = figSize;
  options.mainFolder = mainFolder;
  options.phaseSumPredictionsSubfolder = phaseSumPredictionsSubfolder;
  options.freqLim = xLim;
  options.phaseLim = [-pi pi] + [-pi/4 pi/4];
  sumPhases(areas, conditions(1:min([numel(conditions) 2])),...
    areaFreqFullInterpIndividual, areaCohFullInterpIndividual, options);
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
options.xLabelHist = '# recordings';
options.iAreasOI = iAreas2compareOI;
[phaseHistos, distributionStats] = phaseHistosPlotMaster(drawPhaseHistos, areas, conditions(1:min([numel(conditions) 2])),...
  areaFreqFullInterpIndividual, areaPhaseFullInterpIndividual, edges, options);


%% STATS ON MEAN PHASE FREQUENCY PROFILES
if doStats
  if ~strcmp(repository, 'allensdk')
    [fPEst, fWTest, strPMethod, pEst, U2, pObs, U2Obs] = phaseComparisonStats(areas, areaFreqFullInterpIndividual, FOI, areaPhaseFullInterpIndividual, [areasCritical; areas2compareCritical]);
  else
    fPEst = []; fWTest = []; strPMethod = []; pEst = []; U2 = []; pObs = []; U2Obs = [];
  end
  save(filename, 'phaseHistos','distributionStats','fPEst','fWTest','strPMethod','pEst','U2','pObs','U2Obs', '-append');
end