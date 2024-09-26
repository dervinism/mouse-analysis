% Run this script to produce various unit PSD exponent beta analyses figures.
%
% The following data files are produced:
% dataDir\caDir\unitsFolder\betaCorrelationsSubfolder\beta_units.mat or
%   dataDir\caDir\unitsFolder\betaCorrelationsSubfolder\beta_units_reverse.mat or
%   dataDir\caDir\unitsFolder\betaCorrelationsSubfolder\beta_units_quality.mat or
%   dataDir\caDir\unitsFolder\betaCorrelationsSubfolder\beta_units_reverse_quality.mat.
%   All of these files contain unit betas and other measures that betas are
%   being compared to.
%
% Beta analyses figures are saved in dataDir\laDir\unitsFolder\betaCorrelationsSubfolder
%   and dataDir\caDir\unitsFolder\betaCorrelationsSubfolder depending on whether within
%   area or cross area comparisons are made. The quality unit analyses are
%   save in respective quality subfolders
%   dataDir\laDir\qualityUnitsFolder\betaCorrelationsSubfolder and
%   dataDir\caDir\qualityUnitsFolder\betaCorrelationsSubfolder.

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

dataDir = [dataDir filesep includeRuns];
if strcmp(repository,'all')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep laDir];
    rootFolder_ca = [dataDir filesep caDir];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep laDir_positive];
    rootFolder_ca = [dataDir filesep caDir_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep laDir_negative];
    rootFolder_ca = [dataDir filesep caDir_negative];
  end
  animals = animalsOI;
  xLim = [FOI(end) FOI(1)];
elseif strcmp(repository,'uol')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep laDir_uol];
    rootFolder_ca = [dataDir filesep caDir_uol];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep laDir_uol_positive];
    rootFolder_ca = [dataDir filesep caDir_uol_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep laDir_uol_negative];
    rootFolder_ca = [dataDir filesep caDir_uol_negative];
  end
  animals = animalsUOLOI;
  %animals = animalsNeuropixelsUOL;
  xLim = [FOI(end) FOI(1)];
elseif strcmp(repository,'allensdk')
  if strcmp(subpop, 'all')
    rootFolder = [dataDir filesep laDir_allensdk];
    rootFolder_ca = [dataDir filesep caDir_allensdk];
  elseif strcmp(subpop, 'positive')
    rootFolder = [dataDir filesep laDir_allensdk_positive];
    rootFolder_ca = [dataDir filesep caDir_allensdk_positive];
  elseif strcmp(subpop, 'negative')
    rootFolder = [dataDir filesep laDir_allensdk_negative];
    rootFolder_ca = [dataDir filesep caDir_allensdk_negative];
  end
  animals = animalsAllensdk;
  xLim = [0.045 FOI(1)];
  conditions = {'awake'};
end
frCutoff = 0.5; % spikes per second
cutoffFreq = 0.0547; % frequency cutoff
betaWindowSize = 30; % minutes

if qualityCheck
  mainFolder = [rootFolder filesep qualityUnitsFolder filesep betaCorrelationsSubfolder]; %#ok<*UNRCH>
  mainFolder_ca = [rootFolder_ca filesep qualityUnitsFolder filesep betaCorrelationsSubfolder];
else
  mainFolder = [rootFolder filesep unitsFolder filesep betaCorrelationsSubfolder];
  mainFolder_ca = [rootFolder_ca filesep unitsFolder filesep betaCorrelationsSubfolder];
end


%% COMPUTE VARIABLES AND INITIALISE FIGURES NEEDED FOR DISPLAYING UNIT BETA PLOTS
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
    srData = dataStruct.seriesData.(fnsData{1}).conf.samplingParams.srData;
    optCoh = dataStruct.seriesData.(fnsData{1}).conf.optCoh;
    if animal == 1 || ~exist('areaBetaIndividual', 'var')
      areaBetaIndividual = {};
      areaTruncatedBetaIndividual = {};
      areaBetaWindowsIndividual = {};
      areaBetaMeanWindowsIndividual = {};
      areaBetaMatchWindowsIndividual = {};
      areaFiringRateIndividual = {};
      areaPSDIndividual = {};
      areaTruncatedPSDIndividual = {};
      areaPSDWindowsIndividual = {};
      areaPSDFreqIndividual = {};
      areaTruncatedPSDFreqIndividual = {};
      areaPSDWindowsFreqIndividual = {};
      areaSTPRIndividual = {};
      animalIDs = {};
      animalIDsWindows = {};
      animalIDsMatchWindows = {};
      unitsWindows = {};
      series = {};
      seriesWindows = {};
      seriesMatchWindows = {};
      for iCond = 1:numel(conditions)
        areaBetaIndividualCond = {};
        areaTruncatedBetaIndividualCond = {};
        areaBetaWindowsIndividualCond = {};
        areaBetaMeanWindowsIndividualCond = {};
        areaBetaMatchWindowsIndividualCond = {};
        areaFiringRateIndividualCond = {};
        areaPSDIndividualCond = {};
        areaTruncatedPSDIndividualCond = {};
        areaPSDWindowsIndividualCond = {};
        areaPSDFreqIndividualCond = {};
        areaTruncatedPSDFreqIndividualCond = {};
        areaPSDWindowsFreqIndividualCond = {};
        areaSTPRIndividualCond = {};
        animalIDsCond = {};
        animalIDsWindowsCond = {};
        animalIDsMatchWindowsCond = {};
        unitsWindowsCond = {};
        seriesCond = {};
        seriesWindowsCond = {};
        seriesMatchWindowsCond = {};
        for iArea = 1:numel(areas)
          areaBetaIndividualCond{iArea} = []; %#ok<*SAGROW>
          areaTruncatedBetaIndividualCond{iArea} = [];
          areaBetaWindowsIndividualCond{iArea} = [];
          areaBetaMeanWindowsIndividualCond{iArea} = [];
          areaBetaMatchWindowsIndividualCond{iArea} = [];
          areaFiringRateIndividualCond{iArea} = [];
          areaPSDIndividualCond{iArea} = {};
          areaTruncatedPSDIndividualCond{iArea} = {};
          areaPSDWindowsIndividualCond{iArea} = {};
          areaPSDFreqIndividualCond{iArea} = {};
          areaTruncatedPSDFreqIndividualCond{iArea} = {};
          areaPSDWindowsFreqIndividualCond{iArea} = {};
          areaSTPRIndividualCond{iArea} = {};
          animalIDsCond{iArea} = {};
          animalIDsWindowsCond{iArea} = {};
          animalIDsMatchWindowsCond{iArea} = {};
          unitsWindowsCond{iArea} = [];
          seriesCond{iArea} = {};
          seriesWindowsCond{iArea} = {};
          seriesMatchWindowsCond{iArea} = {};
        end
        areaBetaIndividual{iCond} = areaBetaIndividualCond;
        areaTruncatedBetaIndividual{iCond} = areaBetaIndividualCond;
        areaBetaWindowsIndividual{iCond} = areaBetaWindowsIndividualCond;
        areaBetaMeanWindowsIndividual{iCond} = areaBetaMeanWindowsIndividualCond;
        areaFiringRateIndividual{iCond} = areaFiringRateIndividualCond;
        areaBetaMatchWindowsIndividual{iCond} = areaBetaMatchWindowsIndividualCond;
        areaPSDIndividual{iCond} = areaPSDIndividualCond;
        areaTruncatedPSDIndividual{iCond} = areaTruncatedPSDIndividualCond;
        areaPSDWindowsIndividual{iCond} = areaPSDWindowsIndividualCond;
        areaPSDFreqIndividual{iCond} = areaPSDFreqIndividualCond;
        areaTruncatedPSDFreqIndividual{iCond} = areaTruncatedPSDFreqIndividualCond;
        areaPSDWindowsFreqIndividual{iCond} = areaPSDWindowsFreqIndividualCond;
        areaSTPRIndividual{iCond} = areaSTPRIndividualCond;
        animalIDs{iCond} = animalIDsCond;
        animalIDsWindows{iCond} = animalIDsWindowsCond;
        animalIDsMatchWindows{iCond} = animalIDsMatchWindowsCond;
        unitsWindows{iCond} = unitsWindowsCond;
        series{iCond} = seriesWindowsCond;
        seriesWindows{iCond} = seriesWindowsCond;
        seriesMatchWindows{iCond} = seriesMatchWindowsCond;
      end
    end
    
    for dbCount = 1:numel(fnsData) % Loop through db entries
      dbStruct = dataStruct.seriesData.(fnsData{dbCount});
      seriesName = seriesFromEntry(fnsData{dbCount});
      disp(['             series ' seriesName]);
      
      % Determine if series pupil data exist
      if strcmp(subpop, 'positive') || strcmp(subpop, 'negative')
        if isempty(dbStruct.popData)
          continue
        end
        if ~isfield(dbStruct.popData, 'pupil') || (isfield(dbStruct.popData, 'pupil') && isempty(dbStruct.popData.pupil))
          continue
        end
        if isempty(dbStruct.popData.pupil.popData) || isempty(dbStruct.popData.pupil.popData.phase)
          continue
        end
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
      
      betasWindows = [];
      if numel(conditions) > 1
        condLoop = [iCond numel(conditions)];
      else
        condLoop = iCond;
      end
      for iCondPlusAll = condLoop % Loop through the main and pooled conditions
        for iAreaPlusAll = area % Loop through the main and pooled areas
          close all
          
          % Load and store betas
          betas = dbStruct.popData.phaseCoh.beta(qualityUnitInd);
          areaBetaIndividual{iCondPlusAll}{iAreaPlusAll} = [areaBetaIndividual{iCondPlusAll}{iAreaPlusAll}; betas];
          
          % Load and store firing rates, truncated betas, full and truncated PSD, and associated frequency vectors
          shankIDs = fieldnames(dbStruct.shankData);
          spk = [];
          psd = {};
          psdTruncated = {};
          psdFreq = {};
          psdTruncatedFreq = {};
          betasTruncated = [];
          for sh = 1:numel(shankIDs)
            if isempty(spk)
              spk = dbStruct.shankData.(shankIDs{sh}).spk;
            else
              spk = concatenateMat(spk, dbStruct.shankData.(shankIDs{sh}).spk);
            end
            for u = 1:numel(dbStruct.shankData.(shankIDs{sh}).phaseCoh)
              psd{numel(psd)+1} = dbStruct.shankData.(shankIDs{sh}).phaseCoh{u}.psd;
              psdFreq{numel(psdFreq)+1} = dbStruct.shankData.(shankIDs{sh}).phaseCoh{u}.freq;
              psdTruncated{numel(psdTruncated)+1} = dbStruct.shankData.(shankIDs{sh}).phaseCoh{u}.psd(dbStruct.shankData.(shankIDs{sh}).phaseCoh{u}.freq >= cutoffFreq);
              psdTruncatedFreq{numel(psdTruncatedFreq)+1} = dbStruct.shankData.(shankIDs{sh}).phaseCoh{u}.freq(dbStruct.shankData.(shankIDs{sh}).phaseCoh{u}.freq >= cutoffFreq);
              betasTruncated = [betasTruncated; psdBeta(psdTruncatedFreq{end},psdTruncated{end})];
            end
          end
          firingRate = (sum(spk(qualityUnitInd,:),2).*srData)./size(spk(qualityUnitInd,:),2);
          if isempty(areaFiringRateIndividual{iCondPlusAll}{iAreaPlusAll})
            areaFiringRateIndividual{iCondPlusAll}{iAreaPlusAll} = firingRate;
            areaPSDIndividual{iCondPlusAll}{iAreaPlusAll} = psd';
            areaPSDFreqIndividual{iCondPlusAll}{iAreaPlusAll} = psdFreq';
            areaTruncatedPSDIndividual{iCondPlusAll}{iAreaPlusAll} = psdTruncated';
            areaTruncatedPSDFreqIndividual{iCondPlusAll}{iAreaPlusAll} = psdTruncatedFreq';
          else
            areaFiringRateIndividual{iCondPlusAll}{iAreaPlusAll} = [areaFiringRateIndividual{iCondPlusAll}{iAreaPlusAll}; firingRate];
            areaPSDIndividual{iCondPlusAll}{iAreaPlusAll} = [areaPSDIndividual{iCondPlusAll}{iAreaPlusAll}; psd'];
            areaPSDFreqIndividual{iCondPlusAll}{iAreaPlusAll} = [areaPSDFreqIndividual{iCondPlusAll}{iAreaPlusAll}; psdFreq'];
            areaTruncatedPSDIndividual{iCondPlusAll}{iAreaPlusAll} = [areaTruncatedPSDIndividual{iCondPlusAll}{iAreaPlusAll}; psdTruncated'];
            areaTruncatedPSDFreqIndividual{iCondPlusAll}{iAreaPlusAll} = [areaTruncatedPSDFreqIndividual{iCondPlusAll}{iAreaPlusAll}; psdTruncatedFreq'];
          end
          areaTruncatedBetaIndividual{iCondPlusAll}{iAreaPlusAll} = [areaTruncatedBetaIndividual{iCondPlusAll}{iAreaPlusAll}; betasTruncated];
          
          % Calculate betas for smaller size windows
          if isempty(betasWindows)
            [betasWindows, PSDWindows, PSDWindowsFreq, nWindows] = betaWindows(spk(qualityUnitInd,:), srData, betaWindowSize, optCoh);
            betasMeanWindows = mean(betasWindows,2);
            betasWindows = reshape(betasWindows',numel(betasWindows),1);
          end
          areaBetaWindowsIndividual{iCondPlusAll}{iAreaPlusAll} = [areaBetaWindowsIndividual{iCondPlusAll}{iAreaPlusAll}; betasWindows];
          areaBetaMeanWindowsIndividual{iCondPlusAll}{iAreaPlusAll} = [areaBetaMeanWindowsIndividual{iCondPlusAll}{iAreaPlusAll}; betasMeanWindows];
          if nWindows
            areaBetaMatchWindowsIndividual{iCondPlusAll}{iAreaPlusAll} = [areaBetaMatchWindowsIndividual{iCondPlusAll}{iAreaPlusAll}; betas];
          end
          areaPSDWindowsIndividual{iCondPlusAll}{iAreaPlusAll} = [areaPSDWindowsIndividual{iCondPlusAll}{iAreaPlusAll}; PSDWindows];
          areaPSDWindowsFreqIndividual{iCondPlusAll}{iAreaPlusAll} = [areaPSDWindowsFreqIndividual{iCondPlusAll}{iAreaPlusAll}; PSDWindowsFreq];
          
          % Load and store spike-trigerred PRs
          %           stPR = zeros(size(spk,1), size(dbStruct.shankData.shank1.phaseCoh{1}.stPRsh,2));
          %           unitCount = 0;
          %           for sh = 1:numel(shankIDs)
          %             for u = 1:numel(dbStruct.shankData.(shankIDs{sh}).phaseCoh)
          %               unitCount = unitCount + 1;
          %               stPR(unitCount,:) = dbStruct.shankData.(shankIDs{sh}).phaseCoh{1}.stPRsh;
          %             end
          %           end
          %           stPR = stPR(qualityUnitInd,:);
          %           if isempty(areaSTPRIndividual{iCondPlusAll}{iAreaPlusAll})
          %             areaSTPRIndividual{iCondPlusAll}{iAreaPlusAll} = stPR;
          %           else
          %             areaSTPRIndividual{iCondPlusAll}{iAreaPlusAll} = [areaSTPRIndividual{iCondPlusAll}{iAreaPlusAll}; stPR];
          %           end
          
          % Load and store animal IDs
          animals2add = cell(numel(betas),1);
          [animals2add{1:end}] = deal(animals{animal});
          animalIDs{iCondPlusAll}{iAreaPlusAll} = [animalIDs{iCondPlusAll}{iAreaPlusAll}; animals2add];
          
          if ~isempty(betasWindows)
            animals2add = cell(numel(betasWindows),1);
            [animals2add{1:end}] = deal(animals{animal});
            animalIDsWindows{iCondPlusAll}{iAreaPlusAll} = [animalIDsWindows{iCondPlusAll}{iAreaPlusAll}; animals2add];
          end
          
          if ~isempty(betasWindows)
            assert(numel(betas) == numel(betasMeanWindows));
            animals2add = cell(numel(betas),1);
            [animals2add{1:end}] = deal(animals{animal});
            animalIDsMatchWindows{iCondPlusAll}{iAreaPlusAll} = [animalIDsMatchWindows{iCondPlusAll}{iAreaPlusAll}; animals2add];
          end
          
          % Load and store series names
          series2add = cell(numel(betas),1);
          [series2add{1:end}] = deal(seriesName);
          series{iCondPlusAll}{iAreaPlusAll} = [series{iCondPlusAll}{iAreaPlusAll}; series2add];
          
          if ~isempty(betasWindows)
            series2add = cell(numel(betasWindows),1);
            [series2add{1:end}] = deal(seriesName);
            seriesWindows{iCondPlusAll}{iAreaPlusAll} = [seriesWindows{iCondPlusAll}{iAreaPlusAll}; series2add];
          end
          
          if ~isempty(betasWindows)
            series2add = cell(numel(betas),1);
            [series2add{1:end}] = deal(seriesName);
            seriesMatchWindows{iCondPlusAll}{iAreaPlusAll} = [seriesMatchWindows{iCondPlusAll}{iAreaPlusAll}; series2add];
          end
          
          % Load and store unit order
          if ~isempty(betasWindows)
            if isempty(unitsWindows{iCondPlusAll}{iAreaPlusAll})
              unitOrder = 1;
            else
              unitOrder = unitsWindows{iCondPlusAll}{iAreaPlusAll}(end)+1;
            end
            unitOrder = ones(size(betasWindows)).*unitOrder;
            unitsWindows{iCondPlusAll}{iAreaPlusAll} = [unitsWindows{iCondPlusAll}{iAreaPlusAll}; unitOrder];
          end
          
          % Load and store areas
          %           area2add = cell{size(spk,1),1};
          %           [~, ~, areaName] = determineArea(iArea);
          %           area2add(:) = deal(areaName);
          %           areaNames{iCondPlusAll}{iAreaPlusAll} = [areaNames{iCondPlusAll}{iAreaPlusAll}; area2add];
          
        end
      end
    end
  end
end

% Determine the file name and either save or load the data
if qualityCheck
  filename = [mainFolder filesep 'globalUnitsBeta_quality.mat'];
else
  filename = [mainFolder filesep 'globalUnitsBeta.mat'];
end
if fullRun
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  if strcmpi(repository, 'uol')
    save(filename, 'conditions','areas','FOI','animalIDs','animalIDsWindows','animalIDsMatchWindows','series','seriesWindows','seriesMatchWindows','unitsWindows',...
      'areaBetaIndividual','areaTruncatedBetaIndividual','areaBetaWindowsIndividual','areaBetaMeanWindowsIndividual','areaBetaMatchWindowsIndividual',...
      'areaPSDIndividual','areaTruncatedPSDIndividual','areaPSDWindowsIndividual',...
      'areaPSDFreqIndividual','areaTruncatedPSDFreqIndividual','areaPSDWindowsFreqIndividual','areaFiringRateIndividual');
  elseif strcmpi(repository, 'allensdk')
    save(filename, 'conditions','areas','FOI','animalIDs','animalIDsWindows','series','seriesWindows','unitsWindows',...
      'areaBetaIndividual','areaTruncatedBetaIndividual','areaBetaWindowsIndividual','areaBetaMeanWindowsIndividual',...
      'areaPSDIndividual','areaTruncatedPSDIndividual',...
      'areaPSDFreqIndividual','areaTruncatedPSDFreqIndividual','areaFiringRateIndividual');
  end
else
  load(filename);
end

%% PLOT BETAS ACROSS AREAS
if strcmp(repository, 'uol')
  areasOI = {'VB','lS1','lRSC','CA'};
  areasOIExtra = {'Th','VB','LGN','lS1','lRSC','CA','DG'};
elseif strcmp(repository, 'allensdk')
  areasOI = {'VB','LGN','V1','CA'};
  areasOIExtra = {'Th','VB','LGN','V1','CA','DG'};
end

if ~exist(mainFolder_ca, 'file')
  mkdir(mainFolder_ca);
end

% Bar plots
[fH1, pval_ttest, scatterGroups, dataMean, dataCI95] = barPlotUnits(areaBetaIndividual, areasOI, 1, [-0.5 1.5]);
set(fH1(1), 'Name','Unit PSD exponent beta across brain areas');
filenameFig = [mainFolder_ca filesep 'betaExponentAcrossAreas'];
savefig(fH1(1), filenameFig, 'compact');
print(fH1(1), [filenameFig '.png'],'-dpng','-r300');

[fH2, pval_ttestExtra, scatterGroupsExtra, dataMeanExtra, dataCI95Extra] = barPlotUnits(areaBetaIndividual, areasOIExtra, 1, [-0.5 1.5]);
set(fH2(1), 'Name','Unit PSD exponent beta across brain areas');
filenameFig = [mainFolder_ca filesep 'betaExponentAcrossAreasExtra'];
savefig(fH2(1), filenameFig, 'compact');
print(fH2(1), [filenameFig '.png'],'-dpng','-r300');

% Violin plots
options.yLim = [-0.5 1.5];
options.yLabel = 'Unit PSD exponent beta: (1/f)^{\beta}';
fH3 = multiViolinPlots(scatterGroups, areasOI, dataMean, dataCI95, pval_ttest, options);
set(fH3, 'Name','Unit PSD exponent beta across brain areas');
filenameFig = [mainFolder_ca filesep 'betaExponentAcrossAreasViolins'];
savefig(fH3, filenameFig, 'compact');
print(fH3, [filenameFig '.png'],'-dpng','-r300');

options.yLim = [-0.5 1.5];
options.yLabel = 'Unit PSD exponent beta: (1/f)^{\beta}';
fH4 = multiViolinPlots(scatterGroupsExtra, areasOIExtra, dataMeanExtra, dataCI95Extra, pval_ttestExtra, options);
set(fH4, 'Name','Unit PSD exponent beta across brain areas');
filenameFig = [mainFolder_ca filesep 'betaExponentAcrossAreasViolinsExtra'];
savefig(fH4, filenameFig, 'compact');
print(fH4, [filenameFig '.png'],'-dpng','-r300');

% ANOVA
runAnova(scatterGroups, areasOI, [filename '_betas']);
runAnova(scatterGroupsExtra, areasOIExtra, [filename '_betasExtra']);


%% PLOT WINDOWED BETAS ACROSS AREAS
% Bar plots
[fH5, pval_ttestWindows, scatterGroupsWindows, dataMeanWindows, dataCI95Windows] = barPlotUnits(areaBetaWindowsIndividual, areasOI, 1, [-0.5 1.5]);
set(fH5(1), 'Name','Windowed unit PSD exponent beta across brain areas');
filenameFig = [mainFolder_ca filesep 'betaExponentAcrossAreasWindowSize' num2str(betaWindowSize)];
savefig(fH5(1), filenameFig, 'compact');
print(fH5(1), [filenameFig '.png'],'-dpng','-r300');

[fH6, pval_ttestWindowsExtra, scatterGroupsWindowsExtra, dataMeanWindowsExtra, dataCI95WindowsExtra] = barPlotUnits(areaBetaWindowsIndividual, areasOIExtra, 1, [-0.5 1.5]);
set(fH6(1), 'Name','Windowed unit PSD exponent beta across brain areas');
filenameFig = [mainFolder_ca filesep 'betaExponentAcrossAreasWindowSizeExtra' num2str(betaWindowSize)];
savefig(fH6(1), filenameFig, 'compact');
print(fH6(1), [filenameFig '.png'],'-dpng','-r300');

% Violin plots
options.yLim = [-0.5 1.5];
options.yLabel = 'Windowed unit PSD exponent beta: (1/f)^{\beta}';
fH7 = multiViolinPlots(scatterGroupsWindows, areasOI, dataMeanWindows, dataCI95Windows, pval_ttestWindows, options);
set(fH7, 'Name','Windowed unit PSD exponent beta across brain areas');
filenameFig = [mainFolder_ca filesep 'betaExponentAcrossAreasViolinsWindowSizeExtra' num2str(betaWindowSize)];
savefig(fH7, filenameFig, 'compact');
print(fH7, [filenameFig '.png'],'-dpng','-r300');

options.yLim = [-0.5 1.5];
options.yLabel = 'Windowed unit PSD exponent beta: (1/f)^{\beta}';
fH8 = multiViolinPlots(scatterGroupsWindowsExtra, areasOIExtra, dataMeanWindowsExtra, dataCI95WindowsExtra, pval_ttestWindowsExtra, options);
set(fH8, 'Name','Windowed unit PSD exponent beta across brain areas');
filenameFig = [mainFolder_ca filesep 'betaExponentAcrossAreasViolinsWindowSizeExtra' num2str(betaWindowSize)];
savefig(fH8, filenameFig, 'compact');
print(fH8, [filenameFig '.png'],'-dpng','-r300');

% ANOVA
runAnova(scatterGroupsWindows, areasOI, [filename '_betasWindows']);
runAnova(scatterGroupsWindowsExtra, areasOIExtra, [filename '_betasWindowsExtra']);


%% PLOT MEAN WINDOWED BETAS ACROSS AREAS
% Bar plots
[fH9, pval_ttestMeanWindows, scatterGroupsMeanWindows, dataMeanMeanWindows, dataCI95MeanWindows] = barPlotUnits(areaBetaMeanWindowsIndividual, areasOI, 1, [-0.5 1.5]);
set(fH9(1), 'Name','Mean windowed unit PSD exponent beta across brain areas');
filenameFig = [mainFolder_ca filesep 'betaExponentAcrossAreasMeanWindowSize' num2str(betaWindowSize)];
savefig(fH9(1), filenameFig, 'compact');
print(fH9(1), [filenameFig '.png'],'-dpng','-r300');

[fH10, pval_ttestMeanWindowsExtra, scatterGroupsMeanWindowsExtra, dataMeanMeanWindowsExtra, dataCI95MeanWindowsExtra] = barPlotUnits(areaBetaMeanWindowsIndividual, areasOIExtra, 1, [-0.5 1.5]);
set(fH10(1), 'Name','Mean windowed unit PSD exponent beta across brain areas');
filenameFig = [mainFolder_ca filesep 'betaExponentAcrossAreasMeanWindowSizeExtra' num2str(betaWindowSize)];
savefig(fH10(1), filenameFig, 'compact');
print(fH10(1), [filenameFig '.png'],'-dpng','-r300');

% Violin plots
options.yLim = [-0.5 1.5];
options.yLabel = 'Mean windowed unit PSD exponent beta: (1/f)^{\beta}';
fH11 = multiViolinPlots(scatterGroupsMeanWindows, areasOI, dataMeanMeanWindows, dataCI95MeanWindows, pval_ttestMeanWindows, options);
set(fH11, 'Name','Mean windowed unit PSD exponent beta across brain areas');
filenameFig = [mainFolder_ca filesep 'betaExponentAcrossAreasViolinsMeanWindowSize' num2str(betaWindowSize)];
savefig(fH11, filenameFig, 'compact');
print(fH11, [filenameFig '.png'],'-dpng','-r300');

options.yLim = [-0.5 1.5];
options.yLabel = 'Mean windowed unit PSD exponent beta: (1/f)^{\beta}';
fH12 = multiViolinPlots(scatterGroupsMeanWindowsExtra, areasOIExtra, dataMeanMeanWindowsExtra, dataCI95MeanWindowsExtra, pval_ttestMeanWindowsExtra, options);
set(fH12, 'Name','Mean windowed unit PSD exponent beta across brain areas');
filenameFig = [mainFolder_ca filesep 'betaExponentAcrossAreasViolinsMeanWindowSizeExtra' num2str(betaWindowSize)];
savefig(fH12, filenameFig, 'compact');
print(fH12, [filenameFig '.png'],'-dpng','-r300');

% ANOVA
runAnova(scatterGroupsMeanWindows, areasOI, [filename '_betasMeanWindows']);
runAnova(scatterGroupsMeanWindowsExtra, areasOIExtra, [filename '_betasMeanWindowsExtra']);


%% CORRELATIONS BETWEEN REGULAR BETAS AND BETAS BASED ON 30 MIN DURATION RECORDINGS
for iCond = 1:min([2 numel(conditions)])
  for iArea = 1:numel(iAreasOI)
    betas = areaBetaMatchWindowsIndividual{iCond}{iAreasOI(iArea)};
    betasWindows = areaBetaMeanWindowsIndividual{iCond}{iAreasOI(iArea)};
    if ~isempty(betas) && ~isempty(betasWindows)
      disp(['Processing unit beta data for ' conditions{iCond} ' ' areas{iAreasOI(iArea)}...
        ' (comparison # ' num2str((iCond-1)*numel(areas) + iAreasOI(iArea)) '/' num2str(numel(conditions)*numel(areas)) ')']);
      
      % Correlations and individual graphs
      inds = ~isnan(betas) & ~isnan(betasWindows);
      [rBetaBeta{iCond}{iAreasOI(iArea)}, pvalBetaBeta{iCond}{iAreasOI(iArea)}] = corrMulti(torow(betas(inds)), torow(betasWindows(inds)), 'Spearman');
      nBetaBeta{iCond}{iAreasOI(iArea)}.significant = sum(inds);
      nBetaBeta{iCond}{iAreasOI(iArea)}.total = numel(inds);
      
      figBetaBeta = figure;
      plot(betas(inds), betasWindows(inds), '.', 'MarkerSize',5); hold on
      figTitle = ['Unit PSD beta ' num2str(betaWindowSize) ' min vs full recording correlations for ' areas{iAreasOI(iArea)} ' ' conditions{iCond}...
        ' r=' num2str(rBetaBeta{iCond}{iAreasOI(iArea)}) ' p=' num2str(pvalBetaBeta{iCond}{iAreasOI(iArea)})...
        ' n=' num2str(nBetaBeta{iCond}{iAreasOI(iArea)}.significant) '/' num2str(nBetaBeta{iCond}{iAreasOI(iArea)}.total)];
      ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
        'on', 'k', 'Unit PSD exponent (full recording)', xlim, -1:0.2:1,...
        'on', 'k', ['Unit PSD exponent (mean of ' num2str(betaWindowSize) ' min)'], ylim, -1:0.2:1);
      title(figTitle);
      figName = ['PSDbetaVsbeta' num2str(betaWindowSize) 'min_in_' areas{iAreasOI(iArea)} '_during_' conditions{iCond} '_cutoff_' num2str(frCutoff)];
      set(gcf, 'Name',figName);
      figName = [mainFolder filesep figName];
      figName = strrep(figName, ' ', '_');
      figName = strrep(figName, '.', 'p');
      
      xLim = xlim;
      xAxisLength = xLim(2)-xLim(1);
      yLim = ylim;
      yAxisLength = yLim(2)-yLim(1);
      text(xLim(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.2, ['r=' num2str(rBetaBeta{iCond}{iAreasOI(iArea)})], 'FontSize',16);
      text(xLim(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.1, ['p=' num2str(pvalBetaBeta{iCond}{iAreasOI(iArea)})], 'FontSize',16);
      text(xLim(1)+xAxisLength*0.025, yLim(2)-yAxisLength*0.05,...
        ['n=' num2str(nBetaBeta{iCond}{iAreasOI(iArea)}.significant) '/' num2str(nBetaBeta{iCond}{iAreasOI(iArea)}.total)], 'FontSize',16);
      plot([-1 2],[-1 2], ':k'); hold off
      xlim(xLim);
      ylim(yLim);
      
      label = [2 1.6];
      margin = [0.3 0.55];
      width = 1*figSize-label(1)-margin(1);
      height = 1*figSize-label(2)-margin(2);
      paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
      exportFig(figBetaBeta, [figName '.png'],'-dpng','-r300', paperSize);
      hgsave(figBetaBeta, [figName '.fig']);
      close(figBetaBeta);
    end
  end
end


%% BETA AND FIRING RATE CORRELATIONS
for iCond = 1:min([2 numel(conditions)])
  for iArea = 1:numel(iAreasOI)
    betas = areaBetaIndividual{iCond}{iAreasOI(iArea)};
    firingRate = areaFiringRateIndividual{iCond}{iAreasOI(iArea)};
    if ~isempty(betas) && ~isempty(firingRate)
      disp(['Processing unit beta data for ' conditions{iCond} ' ' areas{iAreasOI(iArea)}...
        ' (comparison # ' num2str((iCond-1)*numel(areas) + iAreasOI(iArea)) '/' num2str(numel(conditions)*numel(areas)) ')']);
      
      % Correlations and individual graphs
      inds = ~isnan(betas) & ~isnan(firingRate) & firingRate >= frCutoff;
      [rBetaFR{iCond}{iAreasOI(iArea)}, pvalBetaFR{iCond}{iAreasOI(iArea)}] = corrMulti(torow(betas(inds)), torow(firingRate(inds)), 'Spearman');
      nBetaFR{iCond}{iAreasOI(iArea)}.significant = sum(inds);
      nBetaFR{iCond}{iAreasOI(iArea)}.total = numel(inds);
      
      figBetaFR = figure;
      plot(betas(inds), firingRate(inds), '.', 'MarkerSize',5)
      set(gca,'Yscale','log')
      figTitle = ['Unit PSD beta and firing rate correlations for ' areas{iAreasOI(iArea)} ' ' conditions{iCond}...
        ' r=' num2str(rBetaFR{iCond}{iAreasOI(iArea)}) ' p=' num2str(pvalBetaFR{iCond}{iAreasOI(iArea)})...
        ' n=' num2str(nBetaFR{iCond}{iAreasOI(iArea)}.significant) '/' num2str(nBetaFR{iCond}{iAreasOI(iArea)}.total)];
      ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
        'on', 'k', 'Unit PSD exponent', xlim, [0 0.5 1],...
        'on', 'k', 'Unit firing rate', ylim, yticks);
      title(figTitle);
      figName = ['PSDbetaVsFiringRate_in_' areas{iAreasOI(iArea)} '_during_' conditions{iCond} '_cutoff_' num2str(frCutoff)];
      set(gcf, 'Name',figName);
      figName = [mainFolder filesep figName];
      figName = strrep(figName, ' ', '_');
      figName = strrep(figName, '.', 'p');
      
      xLim = xlim;
      xAxisLength = xLim(2)-xLim(1);
      yLim = ylim;
      yAxisLength = yLim(2)-yLim(1);
      text(xLim(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.2, ['r=' num2str(rBetaFR{iCond}{iAreasOI(iArea)})], 'FontSize',16);
      text(xLim(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.1, ['p=' num2str(pvalBetaFR{iCond}{iAreasOI(iArea)})], 'FontSize',16);
      text(xLim(1)+xAxisLength*0.025, yLim(2)-yAxisLength*0.05,...
        ['n=' num2str(nBetaFR{iCond}{iAreasOI(iArea)}.significant) '/' num2str(nBetaFR{iCond}{iAreasOI(iArea)}.total)], 'FontSize',16);
      
      label = [2 1.6];
      margin = [0.3 0.55];
      width = 1*figSize-label(1)-margin(1);
      height = 1*figSize-label(2)-margin(2);
      paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
      exportFig(figBetaFR, [figName '.png'],'-dpng','-r300', paperSize);
      hgsave(figBetaFR, [figName '.fig']);
      close(figBetaFR);
    end
  end
end


%% MEAN UNIT PSD FITS
fPSD = figProperties('Bar plot for units', 'normalized', [0, .005, .97, .90], 'w', 'on');
P = [];
PFit = [];
PFitTrunc = [];
PFitWindows = [];
txt = {};
txtFit = {};
txtFitTrunc = {};
txtFitWindows = {};
for iArea = 1:numel(iAreasOI)
  % Load data
  psd = areaPSDIndividual{1}{determineArea(areasOI{iAreasOI(iArea)})};
  psdFreq = areaPSDFreqIndividual{1}{determineArea(areasOI{iAreasOI(iArea)})};
  psdWindows = areaPSDWindowsIndividual{1}{determineArea(areasOI{iAreasOI(iArea)})};
  psdWindowsFreq = areaPSDWindowsFreqIndividual{1}{determineArea(areasOI{iAreasOI(iArea)})};
  
  % Interpolate PSDs and calculate fitted lines
  psdFreqCombined = [];
  beta = [];
  intercept = [];
  betaTrunc = [];
  interceptTrunc = [];
  betaWind = [];
  interceptWind = [];
  for u = 1:numel(psdFreq)
    psdFreqCombined = [psdFreqCombined psdFreq{u} FOI];
  end
  psdFreqCombined = unique(psdFreqCombined);
  psdInterp = NaN(numel(psd),numel(psdFreqCombined));
  psdWindowsInterp = NaN(numel(psdWindows),numel(psdFreqCombined));
  for u = 1:numel(psd)
    psdInterp(u,:) = interp1(psdFreq{u}, psd{u}, psdFreqCombined);
    
    [~, fitCoefs] = psdBeta(psdFreqCombined(~isnan(psdInterp(u,:))), psdInterp(u,~isnan(psdInterp(u,:))));
    beta = [beta; fitCoefs(1)];
    intercept = [intercept; fitCoefs(2)];
    
    [~, fitCoefsTrunc] = psdBeta(psdFreqCombined(psdFreqCombined >= cutoffFreq & ~isnan(psdInterp(u,:))),...
      psdInterp(u,psdFreqCombined >= cutoffFreq & ~isnan(psdInterp(u,:))));
    betaTrunc = [betaTrunc; fitCoefsTrunc(1)];
    interceptTrunc = [interceptTrunc; fitCoefsTrunc(2)];
  end
  for u = 1:numel(psdWindows)
    psdWindowsInterp(u,:) = interp1(psdWindowsFreq{u}, psdWindows{u}, psdFreqCombined);
    
    [~, fitCoefsWindows] = psdBeta(psdFreqCombined(~isnan(psdWindowsInterp(u,:))),...
      psdWindowsInterp(u,~isnan(psdWindowsInterp(u,:))));
    betaWind = [betaWind; fitCoefsWindows(1)];
    interceptWind = [interceptWind; fitCoefsWindows(2)];
  end
  
  [psdMean, psdCI95] = datamean(psdInterp);
  [~, fitCoefsPSDMean] = psdBeta(psdFreqCombined(psdFreqCombined >= 0.01268 & ~isnan(psdMean)),...
    psdMean(psdFreqCombined >= 0.01268 & ~isnan(psdMean)));
  betaPSDMean = fitCoefsPSDMean(1);
  interceptPSDMean = fitCoefsPSDMean(2);
  fitLinePSDMean = 10.^(log10(psdFreqCombined).*betaPSDMean + interceptPSDMean);
  
  [~, fitCoefsPSDMeanTrunc] = psdBeta(psdFreqCombined(psdFreqCombined >= cutoffFreq & ~isnan(psdMean)),...
    psdMean(psdFreqCombined >= cutoffFreq & ~isnan(psdMean)));
  betaPSDMeanTrunc = fitCoefsPSDMeanTrunc(1);
  interceptPSDMeanTrunc = fitCoefsPSDMeanTrunc(2);
  fitLinePSDMeanTrunc = 10.^(log10(psdFreqCombined).*betaPSDMeanTrunc + interceptPSDMeanTrunc);
  
  betaMean = mean(beta,'omitnan');
  interceptMean = mean(intercept,'omitnan');
  fitLine = 10.^(log10(psdFreqCombined).*betaMean + interceptMean);
  
  betaTruncMean = mean(betaTrunc,'omitnan');
  interceptTruncMean = mean(interceptTrunc,'omitnan');
  fitLineTrunc = 10.^(log10(psdFreqCombined).*betaTruncMean + interceptTruncMean);
  
  betaWindMean = mean(betaWind,'omitnan');
  interceptWindMean = mean(interceptWind,'omitnan');
  fitLineWindows = 10.^(log10(psdFreqCombined).*betaWindMean + interceptWindMean);
  
  if iArea == 1
    p = loglog(psdFreqCombined,psdMean, 'Color',areaColours2(determineArea(areasOI{iAreasOI(iArea)}))); hold on
  else
    p = loglog(psdFreqCombined,psdMean, 'Color',areaColours2(determineArea(areasOI{iAreasOI(iArea)})));
  end
  P = [P p];
  txt{numel(txt)+1} = [areasOI{iAreasOI(iArea)} ' mean PSD'];
  ciplot(psdMean+psdCI95(1,:), psdMean+psdCI95(2,:), psdFreqCombined, areaColours2(determineArea(areasOI{iAreasOI(iArea)})), 0.2)
  
  loglog(psdFreqCombined(psdFreqCombined >= 0.01268 & psdFreqCombined <= 1),...
    fitLinePSDMean(psdFreqCombined >= 0.01268 & psdFreqCombined <= 1),...
    'Color',areaColours2(determineArea(areasOI{iAreasOI(iArea)})));
  
  loglog(psdFreqCombined(psdFreqCombined >= cutoffFreq & psdFreqCombined <= 1),...
    fitLinePSDMeanTrunc(psdFreqCombined >= cutoffFreq & psdFreqCombined <= 1),...
    'Color',areaColours2(determineArea(areasOI{iAreasOI(iArea)})));
  
  pFit = loglog(psdFreqCombined(psdFreqCombined <= 1),fitLine(psdFreqCombined <= 1),...
    '--', 'Color',areaColours2(determineArea(areasOI{iAreasOI(iArea)})));
  PFit = [PFit pFit];
  txtFit{numel(txtFit)+1} = [areasOI{iAreasOI(iArea)} ' mean PSD Fit'];
  
  pFitTrunc = loglog(psdFreqCombined(psdFreqCombined >= cutoffFreq & psdFreqCombined <= 1),...
    fitLineTrunc(psdFreqCombined >= cutoffFreq & psdFreqCombined <= 1),...
    '-.', 'Color',areaColours2(determineArea(areasOI{iAreasOI(iArea)})));
  PFitTrunc = [PFitTrunc pFitTrunc];
  txtFitTrunc{numel(txtFitTrunc)+1} = [areasOI{iAreasOI(iArea)} ' mean PSD Fit (>=' num2str(cutoffFreq) 'Hz)'];
  
  pFitWindows = loglog(psdFreqCombined(psdFreqCombined >= cutoffFreq & psdFreqCombined <= 1),...
    fitLineWindows(psdFreqCombined >= cutoffFreq & psdFreqCombined <= 1),...
    ':', 'Color',areaColours2(determineArea(areasOI{iAreasOI(iArea)})));
  PFitWindows = [PFitWindows pFitWindows];
  txtFitWindows{numel(txtFitWindows)+1} = [areasOI{iAreasOI(iArea)} ' mean PSD Fit (' num2str(betaWindowSize) 'min)'];
end
hold off
xlabel('Frequency {Hz)');
ylabel('Power (spikes^2/Hz)');
title('Mean unit PSDs across areas');
legend([P PFit PFitTrunc PFitWindows], [txt txtFit txtFitTrunc txtFitWindows]);
legend boxoff
hgsave(fPSD, [mainFolder_ca filesep 'PSD_fits' '.fig']);



%% SAVE ADDITIONAL VARIABLES
save(filename, 'pval_ttest','scatterGroups','pval_ttestWindows','scatterGroupsWindows','pval_ttestMeanWindows','scatterGroupsMeanWindows',...
  'rBetaBeta','pvalBetaBeta','rBetaFR','pvalBetaFR', '-append');



%% Local functions
function [fH, stats, dataScatter, dataMean, dataCI95, areas] = barPlotUnits(data, areas, condition, yLim, scaleType)

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
for iArea = 1:nAreas
  areaCode = determineArea(areas{iArea});
  bars(iArea) = 1+iArea*(nBarsPerGroup+gap);
  [dataMean(iArea), dataCI95(:,iArea)] = datamean(data{condition}{areaCode});
  bar(bars(iArea), dataMean(iArea), 'FaceColor',areaColours2(areas{iArea}),...
    'EdgeColor',areaColours2(areas{iArea}), 'FaceAlpha',0.2);
  dataScatter{iArea} = data{condition}{areaCode};
end


% Scatter
for iBar = 1:numel(bars)
  scatter(bars(iBar)*ones(size(dataScatter{iBar}))', dataScatter{iBar}',...
    'MarkerEdgeColor',areaColours2(areas{iBar}), 'jitter','on');
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
  xTickPos, 'on', 'k', {'Mean unit PSD exponent beta: (1/f)^{\beta}'}, yLim, yticks);
xTickLabel = areas;
for iLabel = 1:numel(xTickLabel)
  xTickLabel{iLabel} = strrep(xTickLabel{iLabel}, 'l','');
  xTickLabel{iLabel} = strrep(xTickLabel{iLabel}, 'r','');
end
ax1.XTickLabel = xTickLabel;

xLim = xlim;
xAxisLength = xLim(2)-xLim(1);
yAxisLength = yLim(2)-yLim(1);
textStr = [];
for iArea = 1:numel(areas)
  textStr = [textStr 'v' num2str(iArea)]; %#ok<*AGROW>
end
textStr = ['t-test ' textStr(2:end) 'v1 p-val:   '];
for iArea = 1:numel(areas)-1
  textStr = [textStr num2str(stats.pval_ttest(1,ismember(stats.area1,areas{iArea}) & ismember(stats.area2,areas{iArea+1}))) '   '];
end
textStr = [textStr num2str(stats.pval_ttest(1,ismember(stats.area1,areas{1}) & ismember(stats.area2,areas{iArea+1})))];
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
  [~,p(iCombo)] = ttest2(data{condition}{areaCode1(1)}(~isinf(data{condition}{areaCode1(1)})),...
    data{condition}{areaCode2(1)}(~isinf(data{condition}{areaCode2(1)})));
end
end


function runAnova(scatterGroups, areaGroups, filename)

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
writecell(anovaOutput, [filename '_anova' '.txt']);
end


function [betas, psd, freq, nWindows] = betaWindows(spk, srData, windowSize, optPSD)

betas = [];
psd = [];
freq = [];

% Calculate the number of recording chunks
calculationWindowSize = floor(windowSize*60*srData);
nWindows = size(spk,2)/calculationWindowSize;
if nWindows > floor(nWindows)+0.95
  nWindows = ceil(nWindows);
else
  nWindows = floor(nWindows);
end

% Calculate betas
if nWindows >= 1
  betas = zeros(size(spk,1),nWindows);
  psdTemp = {};
  freqTemp = {};
  for u = 1:size(spk,1)
    for iWindow = 1:nWindows
      psdTemp{u}{iWindow} = {};
      freqTemp{u}{iWindow} = {};
    end
  end
  parfor u = 1:size(spk,1)
    disp(['               unit ' num2str(u) '/' num2str(size(spk,1))]);
    for iWindow = 1:nWindows
      [freqTemp{u}{iWindow}, psdTemp{u}{iWindow}] = freqDependentWindowCoherence(full(spk(u,(iWindow-1)*calculationWindowSize+1:min([iWindow*calculationWindowSize size(spk,2)]))), [], 1/srData, [], optPSD);
      betas(u,iWindow) = psdBeta(freqTemp{u}{iWindow}, psdTemp{u}{iWindow});
    end
  end
  psd = {};
  freq = {};
  for u = 1:size(spk,1)
    for iWindow = 1:nWindows
      psd{numel(psd)+1} = psdTemp{u}{iWindow};
      freq{numel(freq)+1} = freqTemp{u}{iWindow};
    end
  end
  psd = psd';
  freq = freq';
end
end