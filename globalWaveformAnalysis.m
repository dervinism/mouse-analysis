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
    rootFolderPositive = [dataDir filesep laDir_positive];
    rootFolderNegative = [dataDir filesep laDir_negative];
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
    rootFolderPositive = [dataDir filesep laDir_uol_positive];
    rootFolderNegative = [dataDir filesep laDir_uol_negative];
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
    rootFolderPositive = [dataDir filesep laDir_allensdk_positive];
    rootFolderNegative = [dataDir filesep laDir_allensdk_negative];
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

drawBetas = false;
drawTaus = false;
spikeWidthDistros = [true true];


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
      if ~isfield(dataStruct, 'seriesData') || isempty(dataStruct.seriesData)
        disp(['Animal ' animals{animal} ' is missing series data. Skipping to the next animal...'])
        continue
      end
    elseif strcmp(subpop, 'positive')
      if ~isfield(dataStruct, 'seriesData_positive') || isempty(dataStruct.seriesData_positive)
        disp(['Animal ' animals{animal} ' is missing postive series data. Skipping to the next animal...'])
        continue
      end
    elseif strcmp(subpop, 'negative') || isempty(dataStruct.seriesData_negative)
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
        srRecording = dataStruct.seriesData.(fnsData{1}).conf.samplingParams.srRecording;
      elseif strcmp(subpop, 'positive')
        fnsData = fieldnames(dataStruct.seriesData_positive);
        FOI = dataStruct.seriesData_positive.(fnsData{1}).conf.FOI;
        srData = dataStruct.seriesData_positive.(fnsData{1}).conf.samplingParams.srData;
        srRecording = dataStruct.seriesData_positive.(fnsData{1}).conf.samplingParams.srRecording;
      elseif strcmp(subpop, 'negative')
        fnsData = fieldnames(dataStruct.seriesData_negative);
        FOI = dataStruct.seriesData_negative.(fnsData{1}).conf.FOI;
        srData = dataStruct.seriesData_negative.(fnsData{1}).conf.samplingParams.srData;
        srRecording = dataStruct.seriesData_negative.(fnsData{1}).conf.samplingParams.srRecording;
      end
    catch
      % do nothing
    end
    if animal == 1 || ~exist('areaHalfAmpDuration', 'var')
      areaHalfAmpDuration = {};
      areaTrough2PeakTime = {};
      areaTroughVsPeakHeightRatio = {};
      for iCond = 1:numel(conditions)
        areaHalfAmpDurationCond = {};
        areaTrough2PeakTimeCond = {};
        areaTroughVsPeakHeightRatioCond = {};
        for iArea = 1:numel(areas)
          areaHalfAmpDurationCond{iArea} = []; %#ok<*SAGROW>
          areaTrough2PeakTimeCond{iArea} = [];
          areaTroughVsPeakHeightRatioCond{iArea} = [];
        end
        areaHalfAmpDuration{iCond} = areaHalfAmpDurationCond;
        areaTrough2PeakTime{iCond} = areaTrough2PeakTimeCond;
        areaTroughVsPeakHeightRatio{iCond} = areaTroughVsPeakHeightRatioCond;
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
          
          % Load and store waveform data for units
          waveforms = [];
          shankIDs = fieldnames(dbStruct.shankData);
          for sh = 1:numel(shankIDs)
            shankStruct = dbStruct.shankData.(shankIDs{sh});
            if ~isfield(shankStruct,'waveformData')
              continue
            end
            if strcmpi(repository, 'uol')
              if isempty(waveforms)
                waveforms.amplitudes = shankStruct.waveformData.amplitudes;
                waveforms.units = shankStruct.waveformData.cluIDs;
                waveforms.waveformsMaxCh = shankStruct.waveformData.maxWaveforms;
                waveforms.waveforms = shankStruct.waveformData.waveforms;
              else
                waveforms.amplitudes = [waveforms.amplitudes; shankStruct.waveformData.amplitudes];
                waveforms.units = [waveforms.units; shankStruct.waveformData.cluIDs];
                waveforms.waveformsMaxCh = [waveforms.waveformsMaxCh; shankStruct.waveformData.maxWaveforms];
                waveforms.waveforms = [waveforms.waveforms; shankStruct.waveformData.waveforms];
              end
            elseif strcmpi(repository, 'allensdk')
              if isempty(waveforms)
                waveforms.units = shankStruct.waveformData.cluIDs;
                waveforms.waveformsMaxCh = shankStruct.waveformData.maxWaveforms;
                waveforms.waveformsTimes = shankStruct.waveformData.waveformsTimes;
              else
                waveforms.units = [waveforms.units; shankStruct.waveformData.cluIDs];
                waveforms.waveformsMaxCh = [waveforms.waveformsMaxCh; shankStruct.waveformData.maxWaveforms];
                waveforms.waveformsTimes = [waveforms.waveformsTimes; shankStruct.waveformData.waveformsTimes];
              end
            end
          end
          
          % Calculate and store waveform parameters
          for u = 1:numel(qualityUnitInd)
            parameters = getSpikeParams(waveforms.waveformsMaxCh(u,:));
            parameters(parameters <= 0) = NaN;
            parameters(1:2) = parameters(1:2)*(1/srRecording);
            areaHalfAmpDuration{iCondPlusAll}{iAreaPlusAll} = [areaHalfAmpDuration{iCondPlusAll}{iAreaPlusAll}; parameters(1)];
            areaTrough2PeakTime{iCondPlusAll}{iAreaPlusAll} = [areaTrough2PeakTime{iCondPlusAll}{iAreaPlusAll}; parameters(2)];
            areaTroughVsPeakHeightRatio{iCondPlusAll}{iAreaPlusAll} = [areaTroughVsPeakHeightRatio{iCondPlusAll}{iAreaPlusAll}; parameters(3)];
          end
        end
      end
    end
  end
end


%% SAVE OR LOAD THE DATA
if qualityCheck
  filename = [mainFolder filesep 'globalUnits_quality.mat'];
else
  filename = [mainFolder filesep 'globalUnits.mat'];
end
if fullRun
  if ~exist(filename, 'file')
    save(filename, 'conditions','areas','areaHalfAmpDuration','areaTrough2PeakTime','areaTroughVsPeakHeightRatio', '-v7.3');
  else
    save(filename, 'conditions','areas','areaHalfAmpDuration','areaTrough2PeakTime','areaTroughVsPeakHeightRatio', '-append');
  end
end
load(filename); %#ok<*UNRCH>


%% SPIKE WIDTH/TROUGH-TO-PEAK TIME VS BETA EXPONENT CORRELATIONS
if ~exist([mainFolder filesep waveformsSubfolder], 'file')
  mkdir([mainFolder filesep waveformsSubfolder]);
end

if drawBetas
  if qualityCheck
    betasData = load([mainFolder filesep betaCorrelationsSubfolder filesep 'globalUnitsBeta_quality.mat']);
  else
    betasData = load([mainFolder filesep betaCorrelationsSubfolder filesep 'globalUnitsBeta.mat']);
  end
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      
      % Spike width
      spikeWidths = areaHalfAmpDuration{iCond}{iAreasOI(iArea)}.*1000;
      betas = betasData.areaBetaIndividual{iCond}{iAreasOI(iArea)};
      options.corrType = 'Spearman';
      options.figTitle = ['Unit spike half-width vs PSD exponent \beta in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
      options.figName = ['spikeWidthVsPSDbeta in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
      options.figFolder = [mainFolder filesep waveformsSubfolder];
      if strcmpi(repository, 'uol')
        options.xLim = [0 200*(1/30000)*1000];
      elseif strcmpi(repository, 'allensdk')
        options.xLim = [];
      end
      options.xLabel = 'Unit spike width (ms)';
      options.xTicks = [];
      options.yLim = [];
      options.yLabel = 'Unit PSD exponent \beta';
      options.yTicks = [];
      options.diagonal = true;
      options.figSize = figSize;
      options.saveFig = true;
      [~, rCoefSpikeWidthVsPSDbeta{iCond}{iAreasOI(iArea)}, pvalSpikeWidthVsPSDbeta{iCond}{iAreasOI(iArea)}] = corrPlot(spikeWidths, betas, options);
      
      % Spike trough to peak time
      trough2peak = areaTrough2PeakTime{iCond}{iAreasOI(iArea)}.*1000;
      options.corrType = 'Spearman';
      options.figTitle = ['Unit spike trough-to-peak time vs PSD exponent \beta in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
      options.figName = ['trough2peakVsPSDbeta in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
      options.figFolder = [mainFolder filesep waveformsSubfolder];
      if strcmpi(repository, 'uol')
        options.xLim = [0 200*(1/30000)*1000];
      elseif strcmpi(repository, 'allensdk')
        options.xLim = [];
      end
      options.xLabel = 'Unit spike trough-to-peak time (ms)';
      options.xTicks = [];
      options.yLim = [];
      options.yLabel = 'Unit PSD exponent \beta';
      options.yTicks = [];
      options.diagonal = true;
      options.figSize = figSize;
      options.saveFig = true;
      [~, rCoefTrough2peakVsPSDbeta{iCond}{iAreasOI(iArea)}, pvalTrough2peakVsPSDbeta{iCond}{iAreasOI(iArea)}] = corrPlot(trough2peak, betas, options);
    end
  end
  save(filename, 'rCoefSpikeWidthVsPSDbeta','pvalSpikeWidthVsPSDbeta','rCoefTrough2peakVsPSDbeta','pvalTrough2peakVsPSDbeta', '-append');
end


%% SPIKE WIDTH/TROUGH-TO-PEAK TIME VS AUTOCORRELATION DECAY TIME CONSTANT CORRELATIONS
if drawTaus
  if ~exist('binSize', 'var')
    binSize = [srData, 200, 100, 50, 20, 10, 8, 5];
  end
  
  % Effective taus
  for iSize = 1:numel(binSize)
    for iCond = 1
      for iArea = 1:numel(acgdUnitFitsExtra{iSize}.areas)
        areaCode = determineArea(acgdUnitFitsExtra{iSize}.areas{iArea});
        
        % Spike width
        spikeWidths = torow(areaHalfAmpDuration{iCond}{areaCode}.*1000);
        taus = torow(acgdUnitFitsExtra{iSize}.effectiveTau{iArea});
        options.corrType = 'Spearman';
        options.figTitle = ['Unit spike half-width vs AC effective decay \tau in '...
          areas{areaCode} ' during ' conditions{iCond} ' binned over ' num2str(binSize(iSize))];
        options.figName = ['spikeWidthVsACdecayTauEff in '...
          areas{areaCode} ' during ' conditions{iCond} ' binned over ' num2str(binSize(iSize))];
        options.figFolder = [mainFolder filesep waveformsSubfolder];
        if strcmpi(repository, 'uol')
          options.xLim = [0 200*(1/30000)*1000];
        elseif strcmpi(repository, 'allensdk')
          options.xLim = [];
        end
        options.xLabel = 'Unit spike width (ms)';
        options.xTicks = [];
        options.yLim = [];
        options.yLabel = '\tau_{effective}';
        options.yTicks = [];
        options.diagonal = true;
        options.figSize = figSize;
        options.saveFig = true;
        [~, rCoefSpikeWidthVsACdecayTauEff{iCond}{areaCode}{iSize}, pvalSpikeWidthVsACdecayTauEff{iCond}{areaCode}{iSize}] = corrPlot(spikeWidths, taus, options);
        
        % Spike trough-to-peak time
        trough2peak = torow(areaTrough2PeakTime{iCond}{areaCode}.*1000);
        options.corrType = 'Spearman';
        options.figTitle = ['Unit spike trough-to-peak vs AC effective decay \tau in '...
          areas{areaCode} ' during ' conditions{iCond} ' binned over ' num2str(binSize(iSize))];
        options.figName = ['trough2peakVsACdecayTauEff in '...
          areas{areaCode} ' during ' conditions{iCond} ' binned over ' num2str(binSize(iSize))];
        options.figFolder = [mainFolder filesep waveformsSubfolder];
        if strcmpi(repository, 'uol')
          options.xLim = [0 200*(1/30000)*1000];
        elseif strcmpi(repository, 'allensdk')
          options.xLim = [];
        end
        options.xLabel = 'Unit spike trough-to-peak time (ms)';
        options.xTicks = [];
        options.yLim = [];
        options.yLabel = '\tau_{effective}';
        options.yTicks = [];
        options.diagonal = true;
        options.figSize = figSize;
        options.saveFig = true;
        [~, rCoefTrough2peakVsACdecayTauEff{iCond}{areaCode}{iSize}, pvalTrough2peakVsACdecayTauEff{iCond}{areaCode}{iSize}] = corrPlot(trough2peak, taus, options);
      end
    end
  end
  
  % Fitted taus
  for iSize = 1:numel(binSize)
    for iCond = 1
      for iArea = 1:numel(acgdUnitFitsExtra{iSize}.areas)
        areaCode = determineArea(acgdUnitFitsExtra{iSize}.areas{iArea});
        
        % Spike width
        spikeWidths = torow(areaHalfAmpDuration{iCond}{areaCode}.*1000);
        taus = torow(acgdUnitFitsExtra{iSize}.fittedTau{iArea});
        options.corrType = 'Spearman';
        options.figTitle = ['Unit spike half-width vs AC fitted decay \tau in '...
          areas{areaCode} ' during ' conditions{iCond} ' binned over ' num2str(binSize(iSize))];
        options.figName = ['spikeWidthVsACdecayTauFit in '...
          areas{areaCode} ' during ' conditions{iCond} ' binned over ' num2str(binSize(iSize))];
        options.figFolder = [mainFolder filesep waveformsSubfolder];
        if strcmpi(repository, 'uol')
          options.xLim = [0 200*(1/30000)*1000];
        elseif strcmpi(repository, 'allensdk')
          options.xLim = [];
        end
        options.xLabel = 'Unit spike width (ms)';
        options.xTicks = [];
        options.yLim = [];
        options.yLabel = '\tau_{fitted}';
        options.yTicks = [];
        options.diagonal = true;
        options.figSize = figSize;
        options.saveFig = true;
        [~, rCoefSpikeWidthVsACdecayTauFit{iCond}{areaCode}{iSize}, pvalSpikeWidthVsACdecayTauFit{iCond}{areaCode}{iSize}] = corrPlot(spikeWidths, taus, options);
        
        % Spike trough-to-peak time
        trough2peak = torow(areaTrough2PeakTime{iCond}{areaCode}.*1000);
        options.corrType = 'Spearman';
        options.figTitle = ['Unit spike trough-to-peak vs AC fitted decay \tau in '...
          areas{areaCode} ' during ' conditions{iCond} ' binned over ' num2str(binSize(iSize))];
        options.figName = ['trough2peakVsACdecayTauFit in '...
          areas{areaCode} ' during ' conditions{iCond} ' binned over ' num2str(binSize(iSize))];
        options.figFolder = [mainFolder filesep waveformsSubfolder];
        if strcmpi(repository, 'uol')
          options.xLim = [0 200*(1/30000)*1000];
        elseif strcmpi(repository, 'allensdk')
          options.xLim = [];
        end
        options.xLabel = 'Unit spike trough-to-peak time (ms)';
        options.xTicks = [];
        options.yLim = [];
        options.yLabel = '\tau_{fitted}';
        options.yTicks = [];
        options.diagonal = true;
        options.figSize = figSize;
        options.saveFig = true;
        [~, rCoefTrough2peakVsACdecayTauFit{iCond}{areaCode}{iSize}, pvalTrough2peakVsACdecayTauFit{iCond}{areaCode}{iSize}] = corrPlot(trough2peak, taus, options);
      end
    end
  end
  
  save(filename, 'rCoefSpikeWidthVsACdecayTauEff','pvalSpikeWidthVsACdecayTauEff','rCoefSpikeWidthVsACdecayTauFit','pvalSpikeWidthVsACdecayTauFit',...
    'rCoefTrough2peakVsACdecayTauEff','pvalTrough2peakVsACdecayTauEff','rCoefTrough2peakVsACdecayTauFit','pvalTrough2peakVsACdecayTauFit', '-append');
end


%% SPIKE WIDTH/TROUGH-TO-PEAK TIME DISTRIBUTIONS ACROSS AREAS
if strcmp(repository, 'uol')
  areasOI = {'VB','S1','RSC','CA','DG'};
  areasOIExtra = {'Th','VB','Po','S1','RSC','CA','DG'};
elseif strcmp(repository, 'allensdk')
  areasOI = {'LGN','LP','V1','CA','DG'};
  areasOIExtra = {'Th','LGN','LP','V1','V2','CA','DG'};
end

if spikeWidthDistros(1)
  
  % Violin plots for core areas: Spike widths
  [spikeWidths, spikeWidthsMean, spikeWidthsCI95] = getDataScatter(areasOI, areaHalfAmpDuration{1}, 1000);
  [statsSpikeWidths.pval_ttest, statsSpikeWidths.area1, statsSpikeWidths.area2] = ttestGroup(areasOI, areaHalfAmpDuration);
  options.yLim = [0 0.5];
  options.yLabel = 'Unit spike width (ms)';
  fH1 = multiViolinPlots(spikeWidths, areasOI, spikeWidthsMean, spikeWidthsCI95, statsSpikeWidths, options);
  set(fH1, 'Name','Unit spike width distribution across core brain areas');
  filenameFig = [mainFolder filesep waveformsSubfolder filesep 'spikeWidthsAcrossCoreAreas' '.fig'];
  savefig(fH1, filenameFig, 'compact');
  print(fH1, [filenameFig '.png'],'-dpng','-r300');
  
  % Violin plots for core areas: Spike trough-to-peak times
  [trough2peak, trough2peakMean, trough2peakCI95] = getDataScatter(areasOI, areaTrough2PeakTime{1}, 1000);
  [statsTrough2peak.pval_ttest, statsTrough2peak.area1, statsTrough2peak.area2] = ttestGroup(areasOI, areaTrough2PeakTime);
  options.yLim = [0 1.5];
  options.yLabel = 'Unit spike trough-to-peak time (ms)';
  fH2 = multiViolinPlots(trough2peak, areasOI, trough2peakMean, trough2peakCI95, statsTrough2peak, options);
  set(fH2, 'Name','Unit spike trough-to-peak time distribution across core brain areas');
  filenameFig = [mainFolder filesep waveformsSubfolder filesep 'trough2peakAcrossCoreAreas' '.fig'];
  savefig(fH2, filenameFig, 'compact');
  print(fH2, [filenameFig '.png'],'-dpng','-r300');
  
  % Violin plots for extended areas: Spike widths
  [spikeWidths, spikeWidthsMean, spikeWidthsCI95] = getDataScatter(areasOIExtra, areaHalfAmpDuration{1}, 1000);
  [statsSpikeWidthsExtra.pval_ttest, statsSpikeWidthsExtra.area1, statsSpikeWidthsExtra.area2] = ttestGroup(areasOIExtra, areaHalfAmpDuration);
  options.yLim = [0 0.5];
  options.yLabel = 'Unit spike width (ms)';
  fH3 = multiViolinPlots(spikeWidths, areasOIExtra, spikeWidthsMean, spikeWidthsCI95, statsSpikeWidthsExtra, options);
  set(fH3, 'Name','Unit spike width distribution across extended brain areas');
  filenameFig = [mainFolder filesep waveformsSubfolder filesep 'spikeWidthsAcrossExtraAreas' '.fig'];
  savefig(fH3, filenameFig, 'compact');
  print(fH3, [filenameFig '.png'],'-dpng','-r300');
  
  % Violin plots for extended areas: Spike trough-to-peak times
  [trough2peak, trough2peakMean, trough2peakCI95] = getDataScatter(areasOIExtra, areaTrough2PeakTime{1}, 1000);
  [statsTrough2peakExtra.pval_ttest, statsTrough2peakExtra.area1, statsTrough2peakExtra.area2] = ttestGroup(areasOIExtra, areaTrough2PeakTime);
  options.yLim = [0 1.5];
  options.yLabel = 'Unit spike trough-to-peak time (ms)';
  fH4 = multiViolinPlots(trough2peak, areasOIExtra, trough2peakMean, trough2peakCI95, statsTrough2peakExtra, options);
  set(fH4, 'Name','Unit spike trough-to-peak time distribution across extended brain areas');
  filenameFig = [mainFolder filesep waveformsSubfolder filesep 'trough2peakAcrossExtraAreas' '.fig'];
  savefig(fH4, filenameFig, 'compact');
  print(fH4, [filenameFig '.png'],'-dpng','-r300');
  
  % Distribution histograms for all areas
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      
      % Spike width
      spikeWidths = areaHalfAmpDuration{iCond}{iAreasOI(iArea)}.*1000;
      nBins = 32;
      %edgesSpikeWidth = 0:ceil(max(spikeWidths))/nBins:ceil(max(spikeWidths));
      edgesSpikeWidth = 0:2/nBins:2;
      nonSignifCount = sum(isnan(spikeWidths));
      if nonSignifCount ~= numel(spikeWidths)
        spikeWidths = [nonSignifCount histcounts(spikeWidths(~isnan(spikeWidths)), edgesSpikeWidth)];
        options.figName = ['spikeWidthDistro in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = [mainFolder filesep waveformsSubfolder];
        options.xLabel = 'Unit spike width (ms)';
        options.yLabel = 'Unit count';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edgesSpikeWidth, spikeWidths, options);
      end
      
      % Spike trough to peak time
      trough2peak = areaTrough2PeakTime{iCond}{iAreasOI(iArea)}.*1000;
      nBins = 32;
      %edgesTrough2peak = 0:ceil(max(trough2peak))/nBins:ceil(max(trough2peak));
      edgesTrough2peak = 0:2/nBins:2;
      nonSignifCount = sum(isnan(trough2peak));
      if nonSignifCount ~= numel(trough2peak)
        trough2peakHist = [nonSignifCount histcounts(trough2peak(~isnan(trough2peak)), edgesTrough2peak)];
        options.figName = ['trough2peakDistro in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = [mainFolder filesep waveformsSubfolder];
        options.xLabel = 'Unit trough-to-peak time (ms)';
        options.yLabel = 'Unit count';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = true;
        options.dataScatter = trough2peak;
        histPlot(edgesTrough2peak, trough2peakHist, options);
      end
    end
  end
  
  save(filename, 'statsSpikeWidths','statsSpikeWidthsExtra','statsTrough2peak','statsTrough2peakExtra', '-append');
end


%% SPIKE WIDTH/TROUGH-TO-PEAK TIME DISTRIBUTIONS ACROSS AREAS DIVIDED BETWEEN POSITIVE AND NEGATIVE UNITS
if spikeWidthDistros(2) && strcmpi(subpop, 'all')
  if qualityCheck
    filenamePositive = [rootFolderPositive filesep qualityUnitsFolder filesep 'globalUnits_quality.mat'];
    filenameNegative = [rootFolderNegative filesep qualityUnitsFolder filesep 'globalUnits_quality.mat'];
  else
    filenamePositive = [rootFolderPositive filesep unitsFolder filesep 'globalUnits.mat'];
    filenameNegative = [rootFolderNegative filesep unitsFolder filesep filesep 'globalUnits.mat'];
  end
  positiveUnitData = load(filenamePositive);
  negativeUnitData = load(filenameNegative);
  
  % Violin plots for core areas: Spike widths
  [spikeWidthsPositive, spikeWidthsMeanPositive, spikeWidthsCI95Positive] = getDataScatter(areasOI, positiveUnitData.areaHalfAmpDuration{1}, 1000);
  [spikeWidthsNegative, spikeWidthsMeanNegative, spikeWidthsCI95Negative] = getDataScatter(areasOI, negativeUnitData.areaHalfAmpDuration{1}, 1000);
  [~, spikeWidths, spikeWidthsMean, spikeWidthsCI95, XTickLabel] = combineData(areasOI,...
    spikeWidthsPositive, spikeWidthsMeanPositive, spikeWidthsCI95Positive,...
    spikeWidthsNegative, spikeWidthsMeanNegative, spikeWidthsCI95Negative);
  [statsSpikeWidthsWithin.pval_ttest, statsSpikeWidthsWithin.area1, statsSpikeWidthsWithin.area2] = ttestGroupCombined(XTickLabel{2}, spikeWidths);
  options.yLim = [0 0.5];
  options.yLabel = 'Unit spike width (ms)';
  options.textStr = 't-test +v- p-val:   ';
  for iArea = 1:numel(statsSpikeWidthsWithin.pval_ttest)
    options.textStr = [options.textStr num2str(statsSpikeWidthsWithin.pval_ttest(iArea)) '   '];
  end
  fH5 = multiViolinPlots(spikeWidths, XTickLabel{2}, spikeWidthsMean, spikeWidthsCI95, statsSpikeWidthsWithin, options);
  set(fH5, 'Name','Unit spike width distribution within core brain areas');
  ax1 = gca;
  ax1.XTickLabel = XTickLabel{1};
  filenameFig = [mainFolder filesep waveformsSubfolder filesep 'spikeWidthsWithinCoreAreas' '.fig'];
  savefig(fH5, filenameFig, 'compact');
  print(fH5, [filenameFig '.png'],'-dpng','-r300');
  
  % Violin plots for core areas: Spike trough-to-peak times
  [trough2peakPositive, trough2peakMeanPositive, trough2peakCI95Positive] = getDataScatter(areasOI, positiveUnitData.areaTrough2PeakTime{1}, 1000);
  [trough2peakNegative, trough2peakMeanNegative, trough2peakCI95Negative] = getDataScatter(areasOI, negativeUnitData.areaTrough2PeakTime{1}, 1000);
  [~, trough2peak, trough2peakMean, trough2peakCI95, XTickLabel] = combineData(areasOI,...
    trough2peakPositive, trough2peakMeanPositive, trough2peakCI95Positive,...
    trough2peakNegative, trough2peakMeanNegative, trough2peakCI95Negative);
  [statsTrough2peakWithin.pval_ttest, statsTrough2peakWithin.area1, statsTrough2peakWithin.area2] = ttestGroupCombined(XTickLabel{2}, trough2peak);
  options.yLim = [0 1.5];
  options.yLabel = 'Unit spike trough-to-peak time (ms)';
  options.textStr = 't-test +v- p-val:   ';
  for iArea = 1:numel(statsTrough2peakWithin.pval_ttest)
    options.textStr = [options.textStr num2str(statsTrough2peakWithin.pval_ttest(iArea)) '   '];
  end
  fH6 = multiViolinPlots(trough2peak, XTickLabel{2}, trough2peakMean, trough2peakCI95, statsTrough2peakWithin, options);
  set(fH6, 'Name','Unit spike trough-to-peak time distribution within core brain areas');
  ax1 = gca;
  ax1.XTickLabel = XTickLabel{1};
  filenameFig = [mainFolder filesep waveformsSubfolder filesep 'trough2peakWithinCoreAreas' '.fig'];
  savefig(fH6, filenameFig, 'compact');
  print(fH6, [filenameFig '.png'],'-dpng','-r300');
  
  % Violin plots for extended areas: Spike widths
  [spikeWidthsPositive, spikeWidthsMeanPositive, spikeWidthsCI95Positive] = getDataScatter(areasOIExtra, positiveUnitData.areaHalfAmpDuration{1}, 1000);
  [spikeWidthsNegative, spikeWidthsMeanNegative, spikeWidthsCI95Negative] = getDataScatter(areasOIExtra, negativeUnitData.areaHalfAmpDuration{1}, 1000);
  [~, spikeWidths, spikeWidthsMean, spikeWidthsCI95, XTickLabel] = combineData(areasOIExtra,...
    spikeWidthsPositive, spikeWidthsMeanPositive, spikeWidthsCI95Positive,...
    spikeWidthsNegative, spikeWidthsMeanNegative, spikeWidthsCI95Negative);
  [statsSpikeWidthsWithinExtra.pval_ttest, statsSpikeWidthsWithinExtra.area1, statsSpikeWidthsWithinExtra.area2] = ttestGroupCombined(XTickLabel{2}, spikeWidths);
  options.yLim = [0 0.5];
  options.yLabel = 'Unit spike width (ms)';
  options.textStr = 't-test +v- p-val:   ';
  for iArea = 1:numel(statsSpikeWidthsWithinExtra.pval_ttest)
    options.textStr = [options.textStr num2str(statsSpikeWidthsWithinExtra.pval_ttest(iArea)) '   '];
  end
  fH7 = multiViolinPlots(spikeWidths, XTickLabel{2}, spikeWidthsMean, spikeWidthsCI95, statsSpikeWidthsWithinExtra, options);
  set(fH7, 'Name','Unit spike width distribution within extended brain areas');
  ax1 = gca;
  ax1.XTickLabel = XTickLabel{1};
  filenameFig = [mainFolder filesep waveformsSubfolder filesep 'spikeWidthsWithinExtraAreas' '.fig'];
  savefig(fH7, filenameFig, 'compact');
  print(fH7, [filenameFig '.png'],'-dpng','-r300');
  
  % Violin plots for extended areas: Spike trough-to-peak times
  [trough2peakPositive, trough2peakMeanPositive, trough2peakCI95Positive] = getDataScatter(areasOIExtra, positiveUnitData.areaTrough2PeakTime{1}, 1000);
  [trough2peakNegative, trough2peakMeanNegative, trough2peakCI95Negative] = getDataScatter(areasOIExtra, negativeUnitData.areaTrough2PeakTime{1}, 1000);
  [~, trough2peak, trough2peakMean, trough2peakCI95, XTickLabel] = combineData(areasOIExtra,...
    trough2peakPositive, trough2peakMeanPositive, trough2peakCI95Positive,...
    trough2peakNegative, trough2peakMeanNegative, trough2peakCI95Negative);
  [statsTrough2peakWithinExtra.pval_ttest, statsTrough2peakWithinExtra.area1, statsTrough2peakWithinExtra.area2] = ttestGroupCombined(XTickLabel{2}, trough2peak);
  options.yLim = [0 1.5];
  options.yLabel = 'Unit spike trough-to-peak time (ms)';
  options.textStr = 't-test +v- p-val:   ';
  for iArea = 1:numel(statsTrough2peakWithinExtra.pval_ttest)
    options.textStr = [options.textStr num2str(statsTrough2peakWithinExtra.pval_ttest(iArea)) '   '];
  end
  fH8 = multiViolinPlots(trough2peak, XTickLabel{2}, trough2peakMean, trough2peakCI95, statsTrough2peakWithinExtra, options);
  set(fH8, 'Name','Unit spike trough-to-peak time distribution within extended brain areas');
  ax1 = gca;
  ax1.XTickLabel = XTickLabel{1};
  filenameFig = [mainFolder filesep waveformsSubfolder filesep 'trough2peakWithinExtraAreas' '.fig'];
  savefig(fH8, filenameFig, 'compact');
  print(fH8, [filenameFig '.png'],'-dpng','-r300');
  
  save(filename, 'statsSpikeWidthsWithin','statsSpikeWidthsWithinExtra','statsTrough2peakWithin','statsTrough2peakWithinExtra', '-append');
end



%% Local functions
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
  [~,p(iCombo)] = ttest2(data{condition}{areaCode1(1)}(~isinf(data{condition}{areaCode1(1)})),...
    data{condition}{areaCode2(1)}(~isinf(data{condition}{areaCode2(1)})));
end
end


function [pSelect, area1select, area2select] = ttestGroupCombined(areaNames, data)

areaCombos = nchoosek(1:numel(areaNames), 2);
nCombos = size(areaCombos,1);
area1 = areaNames(areaCombos(:,1));
area2 = areaNames(areaCombos(:,2));
p = zeros(1,nCombos);
for iCombo = 1:nCombos
  [~,p(iCombo)] = ttest2(data{areaCombos(iCombo,1)}(~isinf(data{areaCombos(iCombo,1)})),...
    data{areaCombos(iCombo,2)}(~isinf(data{areaCombos(iCombo,2)})));
end

pSelect = [];
area1select = {};
area2select = {};
for iArea = 1:2:numel(areaNames)
  area1ind = contains(area1, areaNames{iArea}(1:end-3)) & contains(area1, 'Pos');
  area2ind = contains(area2, areaNames{iArea}(1:end-3)) & contains(area2, 'Neg');
  pSelect = [pSelect p(area1ind & area2ind)];
  area1select = [area1select area1(area1ind & area2ind)];
  area2select = [area2select area2(area1ind & area2ind)];
end
end


function [dataScatter, dataMean, dataCI95] = getDataScatter(areas, data, scaleFactor)

if nargin < 3
  scaleFactor = 1;
end

dataScatter = cell(numel(areas),1);
dataMean = NaN(1,numel(areas));
dataCI95 = NaN(2,numel(areas));
for iArea = 1:numel(areas)
  areaCode = determineArea(areas{iArea});
  dataScatter{iArea} = data{areaCode(1)}.*scaleFactor;
  [dataMean(iArea), dataCI95(:,iArea)] = datamean(data{areaCode(1)}.*scaleFactor);
end
end


function [areasCombined, spikeWidthsCombined, spikeWidthsMeanCombined, spikeWidthsCI95Combined, XTickLabel] = combineData(areas,...
  spikeWidthsPositive, spikeWidthsMeanPositive, spikeWidthsCI95Positive,...
  spikeWidthsNegative, spikeWidthsMeanNegative, spikeWidthsCI95Negative)

areasCombined = {};
spikeWidthsCombined = {};
spikeWidthsMeanCombined = [];
spikeWidthsCI95Combined = [];
XTickLabel{1} = {};
XTickLabel{2} = {};
for iArea = 1:numel(areas)
  areasCombined = [areasCombined; areas{iArea}; areas{iArea}];
  spikeWidthsCombined = [spikeWidthsCombined; spikeWidthsPositive{iArea}; spikeWidthsNegative{iArea}];
  spikeWidthsMeanCombined = [spikeWidthsMeanCombined spikeWidthsMeanPositive(iArea) spikeWidthsMeanNegative(iArea)];
  spikeWidthsCI95Combined = [spikeWidthsCI95Combined spikeWidthsCI95Positive(:,iArea) spikeWidthsCI95Negative(:,iArea)];
  XTickLabel{1} = [XTickLabel{1} [areas{iArea} '+'] [areas{iArea} '-']];
  XTickLabel{2} = [XTickLabel{2} [areas{iArea} 'Pos'] [areas{iArea} 'Neg']];
end
end