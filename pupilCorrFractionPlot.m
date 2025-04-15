% This script produces positive/negative unit and MUA fraction plots
% comparing different brain areas. Files and figures with bar plots, violin
% plots, and correlation plots are produced and saved in dataDir/paDir_uol
% and dataDir/paDir_allensdk folders. Tables are produces to be copied to
% Excel sheets.

clearvars -except repository subpop reverse qualityCheck allData fullRun includeRuns
params
lists

if ~exist('repository','var')
  repository = 'uol';
end
unitThr = 10; % minimum unit or MUA count per data series
alpha = 0.05; % significance level
fullRun = 1; % 1 - all, 2 - bar and violin plots and correlation between fractions in different areas onwards, 3 - Correlation tables onwards
nSampleDisplay = false; % display individual sample sizes on violin plots
violinVisibility = false; % display the violin contour
barPlots = false;
violinPlots = true;
drawCorrs = false; % correlation between fractions in different areas
drawTables = false; % correlation tables
includeRuns = 'noRun';

dataDir = [dataDir_local filesep includeRuns];
outputDir = [outputDir filesep includeRuns];
if strcmp(repository,'uol')
  mainFolder = [outputDir filesep paDir_uol filesep 'all'];
  animals = animalsUOLReduced; %animalsUOLOI;
elseif strcmp(repository,'allensdk')
  mainFolder = [outputDir filesep paDir_allensdk];
  animals = animalsAllensdk;
  conditions = {'awake'};
end


%% Load data
if ~exist('allData', 'var')
  if fullRun == 1
    if strcmp(repository,'uol')
      load([dataDir filesep 'allData_uol.mat']);
    elseif strcmp(repository,'allensdk')
      load([dataDir filesep 'allData_allensdk.mat']);
    end
  end
end



%% Calculate unit fractions for separate brain areas in different conditions
if fullRun == 1
  for animal = 1:numel(animals) % Loop through animals
    disp(['processing animal ' animals{animal}]);
    if exist('allData', 'var')
      dataStruct = allData.(animals{animal}).dataStruct;
    else
      load([dataDir filesep animals{animal} filesep animals{animal} '.mat'])
    end
    fnsData = fieldnames(dataStruct.seriesData);
    
    % Initialise storage variables
    if animal == 1
      unitTable = [];
      muaTable = [];
      
      counter = {};
      counterID = {};
      for iCond = 1:numel(conditions)
        counterCond = {};
        counterIDCond = {};
        for iArea = 1:numel(areas)
          counterCond{iArea} = [];
          counterIDCond{iArea} = {};
        end
        counter{iCond} = counterCond;
        counterID{iCond} = counterIDCond;
      end
      
      % Full recordings
      data.totalUnits = counter;
      data.positiveUnits = counter;
      %       data.positiveSignificantUnits = counter;
      data.negativeUnits = counter;
      %       data.negativeSignificantUnits = counter;
      
      data.totalMUAs = counter;
      data.positiveMUAs = counter;
      %       data.positiveSignificantMUAs = counter;
      data.negativeMUAs = counter;
      %       data.negativeSignificantMUAs = counter;
      
      data.totalUnitsThr = counter;
      data.positiveUnitsThr = counter;
      %       data.positiveSignificantUnitsThr = counter;
      data.negativeUnitsThr = counter;
      %       data.negativeSignificantUnitsThr = counter;
      
      data.totalMUAsThr = counter;
      data.positiveMUAsThr = counter;
      %       data.positiveSignificantMUAsThr = counter;
      data.negativeMUAsThr = counter;
      %       data.negativeSignificantMUAsThr = counter;
      
      %       data.neutralUnits = counter;
      %       data.neutralMUAs = counter;
      %       data.neutralUnitsThr = counter;
      %       data.neutralMUAsThr = counter;
      
      data.recordingDuration = counter;
      data.meanPupilArea = counter;
      
      data.totalUnitsRecID = counterID;
      data.positiveUnitsRecID = counterID;
      %       data.positiveSignificantUnitsRecID = counterID;
      data.negativeUnitsRecID = counterID;
      %       data.negativeSignificantUnitsRecID = counterID;
      
      data.totalMUAsRecID = counterID;
      data.positiveMUAsRecID = counterID;
      %       data.positiveSignificantMUAsRecID = counterID;
      data.negativeMUAsRecID = counterID;
      %       data.negativeSignificantMUAsRecID = counterID;
      
      data.totalUnitsThrRecID = counterID;
      data.positiveUnitsThrRecID = counterID;
      %       data.positiveSignificantUnitsThrRecID = counterID;
      data.negativeUnitsThrRecID = counterID;
      %       data.negativeSignificantUnitsThrRecID = counterID;
      
      data.totalMUAsThrRecID = counterID;
      data.positiveMUAsThrRecID = counterID;
      %       data.positiveSignificantMUAsThrRecID = counterID;
      data.negativeMUAsThrRecID = counterID;
      %       data.negativeSignificantMUAsThrRecID = counterID;
      
      %       data.neutralUnitsRecID = counterID;
      %       data.neutralMUAsRecID = counterID;
      %       data.neutralUnitsThrRecID = counterID;
      %       data.neutralMUAsThrRecID = counterID;
      
      data.recordingDurationRecID = counterID;
      data.meanPupilAreaRecID = counterID;
      
      % 10-minute recording windows
      data.totalUnits10minWindows = counter;
      data.positiveUnits10minWindows = counter;
      %       data.positiveSignificantUnits10minWindows = counter;
      data.negativeUnits10minWindows = counter;
      %       data.negativeSignificantUnits10minWindows = counter;
      
      data.totalMUAs10minWindows = counter;
      data.positiveMUAs10minWindows = counter;
      %       data.positiveSignificantMUAs10minWindows = counter;
      data.negativeMUAs10minWindows = counter;
      %       data.negativeSignificantMUAs10minWindows = counter;
      
      data.totalUnitsThr10minWindows = counter;
      data.positiveUnitsThr10minWindows = counter;
      %       data.positiveSignificantUnitsThr10minWindows = counter;
      data.negativeUnitsThr10minWindows = counter;
      %       data.negativeSignificantUnitsThr10minWindows = counter;
      
      data.totalMUAsThr10minWindows = counter;
      data.positiveMUAsThr10minWindows = counter;
      %       data.positiveSignificantMUAsThr10minWindows = counter;
      data.negativeMUAsThr10minWindows = counter;
      %       data.negativeSignificantMUAsThr10minWindows = counter;
      
      %       data.neutralUnits10minWindows = counter;
      %       data.neutralMUAs10minWindows = counter;
      %       data.neutralUnitsThr10minWindows = counter;
      %       data.neutralMUAsThr10minWindows = counter;
      
      data.meanPupilArea10minWindows = counter;
      
      data.totalUnits10minWindowsRecID = counterID;
      data.positiveUnits10minWindowsRecID = counterID;
      %       data.positiveSignificantUnits10minWindowsRecID = counterID;
      data.negativeUnits10minWindowsRecID = counterID;
      %       data.negativeSignificantUnits10minWindowsRecID = counterID;
      
      data.totalMUAs10minWindowsRecID = counterID;
      data.positiveMUAs10minWindowsRecID = counterID;
      %       data.positiveSignificantMUAs10minWindowsRecID = counterID;
      data.negativeMUAs10minWindowsRecID = counterID;
      %       data.negativeSignificantMUAs10minWindowsRecID = counterID;
      
      data.totalUnitsThr10minWindowsRecID = counterID;
      data.positiveUnitsThr10minWindowsRecID = counterID;
      %       data.positiveSignificantUnitsThr10minWindowsRecID = counterID;
      data.negativeUnitsThr10minWindowsRecID = counterID;
      %       data.negativeSignificantUnitsThr10minWindowsRecID = counterID;
      
      data.totalMUAsThr10minWindowsRecID = counterID;
      data.positiveMUAsThr10minWindowsRecID = counterID;
      %       data.positiveSignificantMUAsThr10minWindowsRecID = counterID;
      data.negativeMUAsThr10minWindowsRecID = counterID;
      %       data.negativeSignificantMUAsThr10minWindowsRecID = counterID;
      
      %       data.neutralUnits10minWindowsRecID = counterID;
      %       data.neutralMUAs10minWindowsRecID = counterID;
      %       data.neutralUnitsThr10minWindowsRecID = counterID;
      %       data.neutralMUAsThr10minWindowsRecID = counterID;
      
      data.meanPupilArea10minWindowsRecID = counterID;
      
      % 20-minute recording windows
      data.totalUnits20minWindows = counter;
      data.positiveUnits20minWindows = counter;
      %       data.positiveSignificantUnits20minWindows = counter;
      data.negativeUnits20minWindows = counter;
      %       data.negativeSignificantUnits20minWindows = counter;
      
      data.totalMUAs20minWindows = counter;
      data.positiveMUAs20minWindows = counter;
      %       data.positiveSignificantMUAs20minWindows = counter;
      data.negativeMUAs20minWindows = counter;
      %       data.negativeSignificantMUAs20minWindows = counter;
      
      data.totalUnitsThr20minWindows = counter;
      data.positiveUnitsThr20minWindows = counter;
      %       data.positiveSignificantUnitsThr20minWindows = counter;
      data.negativeUnitsThr20minWindows = counter;
      %       data.negativeSignificantUnitsThr20minWindows = counter;
      
      data.totalMUAsThr20minWindows = counter;
      data.positiveMUAsThr20minWindows = counter;
      %       data.positiveSignificantMUAsThr20minWindows = counter;
      data.negativeMUAsThr20minWindows = counter;
      %       data.negativeSignificantMUAsThr20minWindows = counter;
      
      %       data.neutralUnits20minWindows = counter;
      %       data.neutralMUAs20minWindows = counter;
      %       data.neutralUnitsThr20minWindows = counter;
      %       data.neutralMUAsThr20minWindows = counter;
      
      data.meanPupilArea20minWindows = counter;
      
      data.totalUnits20minWindowsRecID = counterID;
      data.positiveUnits20minWindowsRecID = counterID;
      %       data.positiveSignificantUnits20minWindowsRecID = counterID;
      data.negativeUnits20minWindowsRecID = counterID;
      %       data.negativeSignificantUnits20minWindowsRecID = counterID;
      
      data.totalMUAs20minWindowsRecID = counterID;
      data.positiveMUAs20minWindowsRecID = counterID;
      %       data.positiveSignificantMUAs20minWindowsRecID = counterID;
      data.negativeMUAs20minWindowsRecID = counterID;
      %       data.negativeSignificantMUAs20minWindowsRecID = counterID;
      
      data.totalUnitsThr20minWindowsRecID = counterID;
      data.positiveUnitsThr20minWindowsRecID = counterID;
      %       data.positiveSignificantUnitsThr20minWindowsRecID = counterID;
      data.negativeUnitsThr20minWindowsRecID = counterID;
      %       data.negativeSignificantUnitsThr20minWindowsRecID = counterID;
      
      data.totalMUAsThr20minWindowsRecID = counterID;
      data.positiveMUAsThr20minWindowsRecID = counterID;
      %       data.positiveSignificantMUAsThr20minWindowsRecID = counterID;
      data.negativeMUAsThr20minWindowsRecID = counterID;
      %       data.negativeSignificantMUAsThr20minWindowsRecID = counterID;
      
      %       data.neutralUnits20minWindowsRecID = counterID;
      %       data.neutralMUAs20minWindowsRecID = counterID;
      %       data.neutralUnitsThr20minWindowsRecID = counterID;
      %       data.neutralMUAsThr20minWindowsRecID = counterID;
      
      data.meanPupilArea20minWindowsRecID = counterID;
      
      % 30-minute recording windows
      data.totalUnits30minWindows = counter;
      data.positiveUnits30minWindows = counter;
      %       data.positiveSignificantUnits30minWindows = counter;
      data.negativeUnits30minWindows = counter;
      %       data.negativeSignificantUnits30minWindows = counter;
      
      data.totalMUAs30minWindows = counter;
      data.positiveMUAs30minWindows = counter;
      %       data.positiveSignificantMUAs30minWindows = counter;
      data.negativeMUAs30minWindows = counter;
      %       data.negativeSignificantMUAs30minWindows = counter;
      
      data.totalUnitsThr30minWindows = counter;
      data.positiveUnitsThr30minWindows = counter;
      %       data.positiveSignificantUnitsThr30minWindows = counter;
      data.negativeUnitsThr30minWindows = counter;
      %       data.negativeSignificantUnitsThr30minWindows = counter;
      
      data.totalMUAsThr30minWindows = counter;
      data.positiveMUAsThr30minWindows = counter;
      %       data.positiveSignificantMUAsThr30minWindows = counter;
      data.negativeMUAsThr30minWindows = counter;
      %       data.negativeSignificantMUAsThr30minWindows = counter;
      
      %       data.neutralUnits30minWindows = counter;
      %       data.neutralMUAs30minWindows = counter;
      %       data.neutralUnitsThr30minWindows = counter;
      %       data.neutralMUAsThr30minWindows = counter;
      
      data.meanPupilArea30minWindows = counter;
      
      data.totalUnits30minWindowsRecID = counterID;
      data.positiveUnits30minWindowsRecID = counterID;
      %       data.positiveSignificantUnits30minWindowsRecID = counterID;
      data.negativeUnits30minWindowsRecID = counterID;
      %       data.negativeSignificantUnits30minWindowsRecID = counterID;
      
      data.totalMUAs30minWindowsRecID = counterID;
      data.positiveMUAs30minWindowsRecID = counterID;
      %       data.positiveSignificantMUAs30minWindowsRecID = counterID;
      data.negativeMUAs30minWindowsRecID = counterID;
      %       data.negativeSignificantMUAs30minWindowsRecID = counterID;
      
      data.totalUnitsThr30minWindowsRecID = counterID;
      data.positiveUnitsThr30minWindowsRecID = counterID;
      %       data.positiveSignificantUnitsThr30minWindowsRecID = counterID;
      data.negativeUnitsThr30minWindowsRecID = counterID;
      %       data.negativeSignificantUnitsThr30minWindowsRecID = counterID;
      
      data.totalMUAsThr30minWindowsRecID = counterID;
      data.positiveMUAsThr30minWindowsRecID = counterID;
      %       data.positiveSignificantMUAsThr30minWindowsRecID = counterID;
      data.negativeMUAsThr30minWindowsRecID = counterID;
      %       data.negativeSignificantMUAsThr30minWindowsRecID = counterID;
      
      %       data.neutralUnits30minWindowsRecID = counterID;
      %       data.neutralMUAs30minWindowsRecID = counterID;
      %       data.neutralUnitsThr30minWindowsRecID = counterID;
      %       data.neutralMUAsThr30minWindowsRecID = counterID;
      
      data.meanPupilArea30minWindowsRecID = counterID;
      
      % Full recordings: IS
      %       data.positiveUnitsIS = counter;
      %       data.positiveSignificantUnitsIS = counter;
      %       data.negativeUnitsIS = counter;
      %       data.negativeSignificantUnitsIS = counter;
      %
      %       data.positiveMUAsIS = counter;
      %       data.positiveSignificantMUAsIS = counter;
      %       data.negativeMUAsIS = counter;
      %       data.negativeSignificantMUAsIS = counter;
      %
      %       data.positiveUnitsThr = counter;
      %       data.positiveSignificantUnitsThr = counter;
      %       data.negativeUnitsThr = counter;
      %       data.negativeSignificantUnitsThr = counter;
      %
      %       data.positiveMUAsThr = counter;
      %       data.positiveSignificantMUAsThr = counter;
      %       data.negativeMUAsThr = counter;
      %       data.negativeSignificantMUAsThr = counter;
      %
      %       data.neutralUnits = counter;
      %       data.neutralMUAs = counter;
      %       data.neutralUnitsThr = counter;
      %       data.neutralMUAsThr = counter;
      %
      %       data.positiveUnitsRecID = counterID;
      %       data.positiveSignificantUnitsRecID = counterID;
      %       data.negativeUnitsRecID = counterID;
      %       data.negativeSignificantUnitsRecID = counterID;
      %
      %       data.positiveMUAsRecID = counterID;
      %       data.positiveSignificantMUAsRecID = counterID;
      %       data.negativeMUAsRecID = counterID;
      %       data.negativeSignificantMUAsRecID = counterID;
      %
      %       data.positiveUnitsThrRecID = counterID;
      %       data.positiveSignificantUnitsThrRecID = counterID;
      %       data.negativeUnitsThrRecID = counterID;
      %       data.negativeSignificantUnitsThrRecID = counterID;
      %
      %       data.positiveMUAsThrRecID = counterID;
      %       data.positiveSignificantMUAsThrRecID = counterID;
      %       data.negativeMUAsThrRecID = counterID;
      %       data.negativeSignificantMUAsThrRecID = counterID;
      %
      %       data.neutralUnitsRecID = counterID;
      %       data.neutralMUAsRecID = counterID;
      %       data.neutralUnitsThrRecID = counterID;
      %       data.neutralMUAsThrRecID = counterID;
    end
    
    for dbCount = 1:numel(fnsData) % Loop through database entries
      dbStruct = dataStruct.seriesData.(fnsData{dbCount});
      seriesName = seriesFromEntry(fnsData{dbCount});
      
      % Generate unit and MUA tables
      if dbCount == 1 || ~exist('prevSeriesName', 'var')
        prevSeriesName = seriesName(1:min([numel(seriesName) 14]));
        unitTableSeries = zeros(1,15);
        muaTableSeries = zeros(1,15);
      elseif ~strcmp(prevSeriesName, seriesName(1:min([numel(seriesName) 14])))
        unitTable = [unitTable; unitTableSeries];
        unitTableSeries = zeros(1,15);
        muaTable = [muaTable; muaTableSeries];
        muaTableSeries = zeros(1,15);
        prevSeriesName = seriesName(1:min([numel(seriesName) 14]));
      end
      [~, ~, areaName] = determineArea(seriesName);
      
      if ~isempty(dbStruct.shankData.shank1.spk) && ~isempty(dbStruct.popData.spkDB)
        nUnits = size(dbStruct.shankData.shank1.spk,1);
        nMUAs = size(dbStruct.popData.spkDB,1);
      elseif isempty(dbStruct.shankData.shank1.units) && ~isempty(dbStruct.popData.spkDB)
        nUnits = 0;
        nMUAs = size(dbStruct.popData.spkDB,1);
      else
        nUnits = 0;
        nMUAs = 0;
      end
      disp([fnsData{dbCount} ' ' areaName ' ' num2str(nUnits) ' ' num2str(nMUAs)]);
      
      if strcmpi(areaName, 'mPFC') || strcmpi(areaName, 'lmPFC') || strcmpi(areaName, 'rmPFC')
        unitTableSeries(1) = unitTableSeries(1) + nUnits;
        muaTableSeries(1) = muaTableSeries(1) + nMUAs;
      elseif strcmpi(areaName, 'lV1')
        unitTableSeries(2) = unitTableSeries(2) + nUnits;
        unitTableSeries(4) = unitTableSeries(4) + nUnits;
        muaTableSeries(2) = muaTableSeries(2) + nMUAs;
        muaTableSeries(4) = muaTableSeries(4) + nMUAs;
      elseif strcmpi(areaName, 'rV1')
        unitTableSeries(3) = unitTableSeries(3) + nUnits;
        unitTableSeries(4) = unitTableSeries(4) + nUnits;
        muaTableSeries(3) = muaTableSeries(3) + nMUAs;
        muaTableSeries(4) = muaTableSeries(4) + nMUAs;
      elseif strcmpi(areaName, 'V1')
        unitTableSeries(4) = unitTableSeries(4) + nUnits;
        muaTableSeries(4) = muaTableSeries(4) + nMUAs;
      elseif strcmpi(areaName, 'lS1') || strcmpi(areaName, 'rS1') || strcmpi(areaName, 'S1')
        unitTableSeries(5) = unitTableSeries(5) + nUnits;
        muaTableSeries(5) = muaTableSeries(5) + nMUAs;
      elseif strcmpi(areaName, 'lRSC') || strcmpi(areaName, 'rRSC') || strcmpi(areaName, 'RSC')
        unitTableSeries(6) = unitTableSeries(6) + nUnits;
        muaTableSeries(6) = muaTableSeries(6) + nMUAs;
      elseif strcmpi(areaName, 'lVB') || strcmpi(areaName, 'rVB') || strcmpi(areaName, 'VB') ||...
          strcmpi(areaName, 'lVB1') || strcmpi(areaName, 'rVB1') || strcmpi(areaName, 'VB1') ||...
          strcmpi(areaName, 'lVB2') || strcmpi(areaName, 'rVB2') || strcmpi(areaName, 'VB2')
        unitTableSeries(7) = unitTableSeries(7) + nUnits;
        unitTableSeries(11) = unitTableSeries(11) + nUnits;
        muaTableSeries(7) = muaTableSeries(7) + nMUAs;
        muaTableSeries(11) = muaTableSeries(11) + nMUAs;
      elseif strcmpi(areaName, 'lPo') || strcmpi(areaName, 'rPo') || strcmpi(areaName, 'Po') ||...
          strcmpi(areaName, 'lPo1') || strcmpi(areaName, 'rPo1') || strcmpi(areaName, 'Po1') ||...
          strcmpi(areaName, 'lPo2') || strcmpi(areaName, 'rPo2') || strcmpi(areaName, 'Po2')
        unitTableSeries(8) = unitTableSeries(8) + nUnits;
        unitTableSeries(11) = unitTableSeries(11) + nUnits;
        muaTableSeries(8) = muaTableSeries(8) + nMUAs;
        muaTableSeries(11) = muaTableSeries(11) + nMUAs;
      elseif strcmpi(areaName, 'lLP') || strcmpi(areaName, 'rLP') || strcmpi(areaName, 'LP') ||...
          strcmpi(areaName, 'lLP1') || strcmpi(areaName, 'rLP1') || strcmpi(areaName, 'LP1') ||...
          strcmpi(areaName, 'lLP2') || strcmpi(areaName, 'rLP2') || strcmpi(areaName, 'LP2')
        unitTableSeries(9) = unitTableSeries(9) + nUnits;
        unitTableSeries(11) = unitTableSeries(11) + nUnits;
        muaTableSeries(9) = muaTableSeries(9) + nMUAs;
        muaTableSeries(11) = muaTableSeries(11) + nMUAs;
      elseif strcmpi(areaName, 'lLGN') || strcmpi(areaName, 'rLGN') || strcmpi(areaName, 'LGN') ||...
          strcmpi(areaName, 'lLGN1') || strcmpi(areaName, 'rLGN1') || strcmpi(areaName, 'LGN1') ||...
          strcmpi(areaName, 'lLGN2') || strcmpi(areaName, 'rLGN2') || strcmpi(areaName, 'LGN2')
        unitTableSeries(10) = unitTableSeries(10) + nUnits;
        unitTableSeries(11) = unitTableSeries(11) + nUnits;
        muaTableSeries(10) = muaTableSeries(10) + nMUAs;
        muaTableSeries(11) = muaTableSeries(11) + nMUAs;
      elseif strcmpi(areaName, 'lCA1') || strcmpi(areaName, 'rCA1') || strcmpi(areaName, 'CA1')
        unitTableSeries(12) = unitTableSeries(12) + nUnits;
        unitTableSeries(14) = unitTableSeries(14) + nUnits;
        muaTableSeries(12) = muaTableSeries(12) + nMUAs;
        muaTableSeries(14) = muaTableSeries(14) + nMUAs;
      elseif strcmpi(areaName, 'lCA3') || strcmpi(areaName, 'rCA3') || strcmpi(areaName, 'CA3')
        unitTableSeries(13) = unitTableSeries(13) + nUnits;
        unitTableSeries(14) = unitTableSeries(14) + nUnits;
        muaTableSeries(13) = muaTableSeries(13) + nMUAs;
        muaTableSeries(14) = muaTableSeries(14) + nMUAs;
      elseif strcmpi(areaName, 'lDG') || strcmpi(areaName, 'rDG') || strcmpi(areaName, 'DG')
        unitTableSeries(15) = unitTableSeries(15) + nUnits;
        muaTableSeries(15) = muaTableSeries(15) + nMUAs;
      end
      
      if dbCount == numel(fnsData)
        unitTable = [unitTable; unitTableSeries];
        muaTable = [muaTable; muaTableSeries];
      end
      
      % Determine if series pupil data exist
      if isempty(dbStruct.popData)
        continue
      end
      if ~isfield(dbStruct.popData, 'pupil') || (isfield(dbStruct.popData, 'pupil') && isempty(dbStruct.popData.pupil))
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
      
      % Get positively and negatively correlated unit and MUA data
      if isfield(dataStruct, 'seriesData_positive') && isfield(dataStruct.seriesData_positive, fnsData{dbCount})
        dbStructPositive = dataStruct.seriesData_positive.(fnsData{dbCount});
      else
        dbStructPositive = [];
      end
      if isfield(dataStruct, 'seriesData_negative') && isfield(dataStruct.seriesData_negative, fnsData{dbCount})
        dbStructNegative = dataStruct.seriesData_negative.(fnsData{dbCount});
      else
        dbStructNegative = [];
      end
      if isempty(dbStructPositive) && isempty(dbStructNegative)
        continue
      end
      dbStruct = dataStruct.seriesData.(fnsData{dbCount});
      entryNameEye = [animals{animal} '_s' seriesName(1:14)];
      eyeDataDB = dataStruct.eyeData.(entryNameEye);
      
      for iAreaPlusAll = area % Loop through the main and pooled areas. This is the core of data extraction
        if numel(conditions) > 1
          condLoop = [iCond numel(conditions)];
        else
          condLoop = iCond;
        end
        for iCondPlusAll = condLoop % Loop through the main and pooled conditions
          
          if ~isempty(dbStructPositive) && ~isempty(dbStructPositive.shankData) && ~isempty(dbStructPositive.shankData.shank1) % We are only interested in single shank recordings
            unitMetadataPositive = dbStructPositive.shankData.shank1.unitMetadata;
          else
            unitMetadataPositive = [];
          end
          if ~isempty(dbStructNegative) && ~isempty(dbStructNegative.shankData) && ~isempty(dbStructNegative.shankData.shank1)
            unitMetadataNegative = dbStructNegative.shankData.shank1.unitMetadata;
          else
            unitMetadataNegative = [];
          end
          srData = dbStruct.conf.samplingParams.srData;
          commonPeriod = combinePeriods(dbStruct.dbSeries.period, eyeDataDB.period, srData);
          [~, times] = pupilFilt(eyeDataDB, srData, dbStruct.popData.MUAsAll, 2*srData, commonPeriod, srData);
          data.recordingDuration{iCondPlusAll}{iAreaPlusAll} = [data.recordingDuration{iCondPlusAll}{iAreaPlusAll};...
            numel(times)/srData];
          data.recordingDurationRecID{iCondPlusAll}{iAreaPlusAll} = [data.recordingDurationRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          data.meanPupilArea{iCondPlusAll}{iAreaPlusAll} = [data.meanPupilArea{iCondPlusAll}{iAreaPlusAll}; dbStruct.popData.meanPupilArea];
          data.meanPupilAreaRecID{iCondPlusAll}{iAreaPlusAll} = [data.meanPupilAreaRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          for iWindow = 1:numel(dbStruct.popData.meanPupilArea10minWindows)
            data.meanPupilArea10minWindows{iCondPlusAll}{iAreaPlusAll} = [data.meanPupilArea10minWindows{iCondPlusAll}{iAreaPlusAll};...
              dbStruct.popData.meanPupilArea10minWindows(iWindow)];
            data.meanPupilArea10minWindowsRecID{iCondPlusAll}{iAreaPlusAll} = [data.meanPupilArea10minWindowsRecID{iCondPlusAll}{iAreaPlusAll};...
              [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
          end
          for iWindow = 1:numel(dbStruct.popData.meanPupilArea20minWindows)
            data.meanPupilArea20minWindows{iCondPlusAll}{iAreaPlusAll} = [data.meanPupilArea20minWindows{iCondPlusAll}{iAreaPlusAll};...
              dbStruct.popData.meanPupilArea20minWindows(iWindow)];
            data.meanPupilArea20minWindowsRecID{iCondPlusAll}{iAreaPlusAll} = [data.meanPupilArea20minWindowsRecID{iCondPlusAll}{iAreaPlusAll};...
              [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
          end
          for iWindow = 1:numel(dbStruct.popData.meanPupilArea30minWindows)
            data.meanPupilArea30minWindows{iCondPlusAll}{iAreaPlusAll} = [data.meanPupilArea30minWindows{iCondPlusAll}{iAreaPlusAll};...
              dbStruct.popData.meanPupilArea30minWindows(iWindow)];
            data.meanPupilArea30minWindowsRecID{iCondPlusAll}{iAreaPlusAll} = [data.meanPupilArea30minWindowsRecID{iCondPlusAll}{iAreaPlusAll};...
              [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
          end
          
          allUnits = size(unitMetadataPositive,1) + size(unitMetadataNegative,1);
          data.totalUnits{iCondPlusAll}{iAreaPlusAll} = [data.totalUnits{iCondPlusAll}{iAreaPlusAll};...
            allUnits]; %#ok<*SAGROW>
          data.positiveUnits{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnits{iCondPlusAll}{iAreaPlusAll};...
            size(unitMetadataPositive,1)];
          %           if ~isempty(unitMetadataPositive)
          %             plusUnits = sum(unitMetadataPositive(:,11) < alpha);
          %           else
          %             plusUnits = 0;
          %           end
          %           data.positiveSignificantUnits{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnits{iCondPlusAll}{iAreaPlusAll};...
          %             plusUnits];
          data.negativeUnits{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnits{iCondPlusAll}{iAreaPlusAll};...
            size(unitMetadataNegative,1)];
          %           if ~isempty(unitMetadataNegative)
          %             minusUnits = sum(unitMetadataNegative(:,11) < alpha);
          %           else
          %             minusUnits = 0;
          %           end
          %           data.negativeSignificantUnits{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnits{iCondPlusAll}{iAreaPlusAll};...
          %             minusUnits];
          %           data.neutralUnits{iCondPlusAll}{iAreaPlusAll} = [data.neutralUnits{iCondPlusAll}{iAreaPlusAll};...
          %             allUnits - plusUnits - minusUnits];
          data.totalUnitsRecID{iCondPlusAll}{iAreaPlusAll} = [data.totalUnitsRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          data.positiveUnitsRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          %           data.positiveSignificantUnitsRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsRecID{iCondPlusAll}{iAreaPlusAll};...
          %             seriesName(1:min([numel(seriesName) 14]))];
          data.negativeUnitsRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          %           data.negativeSignificantUnitsRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsRecID{iCondPlusAll}{iAreaPlusAll};...
          %             seriesName(1:min([numel(seriesName) 14]))];
          %           data.neutralUnitsRecID{iCondPlusAll}{iAreaPlusAll} = [data.neutralUnitsRecID{iCondPlusAll}{iAreaPlusAll};...
          %             seriesName(1:min([numel(seriesName) 14]))];
          
          if ~isempty(dbStructPositive) && ~isempty(dbStructNegative)
            data.totalMUAs{iCondPlusAll}{iAreaPlusAll} = [data.totalMUAs{iCondPlusAll}{iAreaPlusAll};...
              numel(dbStructPositive.popData.spkDB_units) + numel(dbStructNegative.popData.spkDB_units)];
            data.positiveMUAs{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAs{iCondPlusAll}{iAreaPlusAll};...
              numel(dbStructPositive.popData.spkDB_units)];
            %             data.positiveSignificantMUAs{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAs{iCondPlusAll}{iAreaPlusAll};...
            %               sum(dbStructPositive.popData.pvalSpearman < alpha)];
            data.negativeMUAs{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAs{iCondPlusAll}{iAreaPlusAll};...
              numel(dbStructNegative.popData.spkDB_units)];
            %             data.negativeSignificantMUAs{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAs{iCondPlusAll}{iAreaPlusAll};...
            %               sum(dbStructNegative.popData.pvalSpearman < alpha)];
            %             data.neutralMUAs{iCondPlusAll}{iAreaPlusAll} = [data.neutralMUAs{iCondPlusAll}{iAreaPlusAll};...
            %               numel(dbStructPositive.popData.spkDB_units) + numel(dbStructNegative.popData.spkDB_units) - ...
            %               sum(dbStructPositive.popData.pvalSpearman < alpha) - sum(dbStructNegative.popData.pvalSpearman < alpha)];
            data.totalMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.totalMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.positiveMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            %             data.positiveSignificantMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
            %               seriesName(1:min([numel(seriesName) 14]))];
            data.negativeMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            %             data.negativeSignificantMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
            %               seriesName(1:min([numel(seriesName) 14]))];
            %             data.neutralMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.neutralMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
            %               seriesName(1:min([numel(seriesName) 14]))];
          elseif ~isempty(dbStructPositive) && isempty(dbStructNegative)
            data.totalMUAs{iCondPlusAll}{iAreaPlusAll} = [data.totalMUAs{iCondPlusAll}{iAreaPlusAll};...
              numel(dbStructPositive.popData.spkDB_units)];
            data.positiveMUAs{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAs{iCondPlusAll}{iAreaPlusAll};...
              numel(dbStructPositive.popData.spkDB_units)];
            %             data.positiveSignificantMUAs{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAs{iCondPlusAll}{iAreaPlusAll};...
            %               sum(dbStructPositive.popData.pvalSpearman < alpha)];
            %             data.neutralMUAs{iCondPlusAll}{iAreaPlusAll} = [data.neutralMUAs{iCondPlusAll}{iAreaPlusAll};...
            %               numel(dbStructPositive.popData.spkDB_units) - sum(dbStructPositive.popData.pvalSpearman < alpha)];
            data.totalMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.totalMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.positiveMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            %             data.positiveSignificantMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
            %               seriesName(1:min([numel(seriesName) 14]))];
            %             data.neutralMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.neutralMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
            %               seriesName(1:min([numel(seriesName) 14]))];
          elseif isempty(dbStructPositive) && ~isempty(dbStructNegative)
            data.totalMUAs{iCondPlusAll}{iAreaPlusAll} = [data.totalMUAs{iCondPlusAll}{iAreaPlusAll};...
              numel(dbStructNegative.popData.spkDB_units)];
            data.negativeMUAs{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAs{iCondPlusAll}{iAreaPlusAll};...
              numel(dbStructNegative.popData.spkDB_units)];
            %             data.negativeSignificantMUAs{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAs{iCondPlusAll}{iAreaPlusAll};...
            %               sum(dbStructNegative.popData.pvalSpearman < alpha)];
            %             data.neutralMUAs{iCondPlusAll}{iAreaPlusAll} = [data.neutralMUAs{iCondPlusAll}{iAreaPlusAll};...
            %               numel(dbStructNegative.popData.spkDB_units) - sum(dbStructNegative.popData.pvalSpearman < alpha)];
            data.totalMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.totalMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.negativeMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            %             data.negativeSignificantMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
            %               seriesName(1:min([numel(seriesName) 14]))];
            %             data.neutralMUAsRecID{iCondPlusAll}{iAreaPlusAll} = [data.neutralMUAsRecID{iCondPlusAll}{iAreaPlusAll};...
            %               seriesName(1:min([numel(seriesName) 14]))];
          end
          
          if size(unitMetadataPositive,1) + size(unitMetadataNegative,1) >= unitThr
            data.totalUnitsThr{iCondPlusAll}{iAreaPlusAll} = [data.totalUnitsThr{iCondPlusAll}{iAreaPlusAll};...
              allUnits]; %#ok<*SAGROW>
            data.positiveUnitsThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsThr{iCondPlusAll}{iAreaPlusAll};...
              size(unitMetadataPositive,1)];
            %             data.positiveSignificantUnitsThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsThr{iCondPlusAll}{iAreaPlusAll};...
            %               plusUnits];
            data.negativeUnitsThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsThr{iCondPlusAll}{iAreaPlusAll};...
              size(unitMetadataNegative,1)];
            %             data.negativeSignificantUnitsThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsThr{iCondPlusAll}{iAreaPlusAll};...
            %               minusUnits];
            %             data.neutralUnitsThr{iCondPlusAll}{iAreaPlusAll} = [data.neutralUnitsThr{iCondPlusAll}{iAreaPlusAll};...
            %               allUnits - plusUnits - minusUnits];
            data.totalUnitsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.totalUnitsThrRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.positiveUnitsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsThrRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            %             data.positiveSignificantUnitsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsThrRecID{iCondPlusAll}{iAreaPlusAll};...
            %               seriesName(1:min([numel(seriesName) 14]))];
            data.negativeUnitsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsThrRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            %             data.negativeSignificantUnitsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsThrRecID{iCondPlusAll}{iAreaPlusAll};...
            %               seriesName(1:min([numel(seriesName) 14]))];
            %             data.neutralUnitsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.neutralUnitsThrRecID{iCondPlusAll}{iAreaPlusAll};...
            %               seriesName(1:min([numel(seriesName) 14]))];
          end
          
          if (~isempty(dbStructPositive) || ~isempty(dbStructNegative)) &&...
              numel(dbStructPositive.popData.spkDB_units) + numel(dbStructNegative.popData.spkDB_units) >= unitThr
            if ~isempty(dbStructPositive) && ~isempty(dbStructNegative)
              data.totalMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.totalMUAsThr{iCondPlusAll}{iAreaPlusAll};...
                numel(dbStructPositive.popData.spkDB_units) + numel(dbStructNegative.popData.spkDB_units)];
              data.positiveMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsThr{iCondPlusAll}{iAreaPlusAll};...
                numel(dbStructPositive.popData.spkDB_units)];
              %               data.positiveSignificantMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsThr{iCondPlusAll}{iAreaPlusAll};...
              %                 sum(dbStructPositive.popData.pvalSpearman < alpha)];
              data.negativeMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsThr{iCondPlusAll}{iAreaPlusAll};...
                numel(dbStructNegative.popData.spkDB_units)];
              %               data.negativeSignificantMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsThr{iCondPlusAll}{iAreaPlusAll};...
              %                 sum(dbStructNegative.popData.pvalSpearman < alpha)];
              %               data.neutralMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.neutralMUAsThr{iCondPlusAll}{iAreaPlusAll};...
              %                 numel(dbStructPositive.popData.spkDB_units) + numel(dbStructNegative.popData.spkDB_units) - ...
              %                 sum(dbStructPositive.popData.pvalSpearman < alpha) - sum(dbStructNegative.popData.pvalSpearman < alpha)];
              data.totalMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.totalMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.positiveMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              %               data.positiveSignificantMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
              %                 seriesName(1:min([numel(seriesName) 14]))];
              data.negativeMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              %               data.negativeSignificantMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
              %                 seriesName(1:min([numel(seriesName) 14]))];
              %               data.neutralMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.neutralMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
              %                 seriesName(1:min([numel(seriesName) 14]))];
            elseif ~isempty(dbStructPositive) && isempty(dbStructNegative)
              data.totalMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.totalMUAsThr{iCondPlusAll}{iAreaPlusAll};...
                numel(dbStructPositive.popData.spkDB_units)];
              data.positiveMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsThr{iCondPlusAll}{iAreaPlusAll};...
                numel(dbStructPositive.popData.spkDB_units)];
              %               data.positiveSignificantMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsThr{iCondPlusAll}{iAreaPlusAll};...
              %                 sum(dbStructPositive.popData.pvalSpearman < alpha)];
              %               data.neutralMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.neutralMUAsThr{iCondPlusAll}{iAreaPlusAll};...
              %                 numel(dbStructPositive.popData.spkDB_units) - sum(dbStructPositive.popData.pvalSpearman < alpha)];
              data.totalMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.totalMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.positiveMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              %               data.positiveSignificantMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
              %                 seriesName(1:min([numel(seriesName) 14]))];
              %               data.neutralMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.neutralMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
              %                 seriesName(1:min([numel(seriesName) 14]))];
            elseif isempty(dbStructPositive) && ~isempty(dbStructNegative)
              data.totalMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.totalMUAsThr{iCondPlusAll}{iAreaPlusAll};...
                numel(dbStructNegative.popData.spkDB_units)];
              data.negativeMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsThr{iCondPlusAll}{iAreaPlusAll};...
                numel(dbStructNegative.popData.spkDB_units)];
              %               data.negativeSignificantMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsThr{iCondPlusAll}{iAreaPlusAll};...
              %                 sum(dbStructNegative.popData.pvalSpearman < alpha)];
              %               data.neutralMUAsThr{iCondPlusAll}{iAreaPlusAll} = [data.neutralMUAsThr{iCondPlusAll}{iAreaPlusAll};...
              %                 numel(dbStructNegative.popData.spkDB_units) - sum(dbStructNegative.popData.pvalSpearman < alpha)];
              data.totalMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.totalMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.negativeMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              %               data.negativeSignificantMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
              %                 seriesName(1:min([numel(seriesName) 14]))];
              %               data.neutralMUAsThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.neutralMUAsThrRecID{iCondPlusAll}{iAreaPlusAll};...
              %                 seriesName(1:min([numel(seriesName) 14]))];
            end
          end
          
          % Smaller windows. Similar script as for full recordings (lines 540-767 but for smaller recording windows).
          data = assignValues2Windows(dbStruct, 10, data, iCondPlusAll, iAreaPlusAll, seriesName, unitThr, alpha);
          data = assignValues2Windows(dbStruct, 20, data, iCondPlusAll, iAreaPlusAll, seriesName, unitThr, alpha);
          data = assignValues2Windows(dbStruct, 30, data, iCondPlusAll, iAreaPlusAll, seriesName, unitThr, alpha);
        end
      end
    end
  end
  
  % Calculate mean fractions and 95% confidence limits on mean fractions
  for iCond = 1:numel(conditions)
    for iArea = 1:numel(areas)
      [data.positiveUnitsFrCI95{iCond}{iArea}, ~, ~, data.positiveUnitsFr{iCond}{iArea}] = calc95CI((data.positiveUnits{iCond}{iArea}./data.totalUnits{iCond}{iArea}));
      %       [data.positiveSignificantUnitsFrCI95{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFr{iCond}{iArea}] = calc95CI((data.positiveSignificantUnits{iCond}{iArea}./data.totalUnits{iCond}{iArea}));
      [data.negativeUnitsFrCI95{iCond}{iArea}, ~, ~, data.negativeUnitsFr{iCond}{iArea}] = calc95CI((data.negativeUnits{iCond}{iArea}./data.totalUnits{iCond}{iArea}));
      %       [data.negativeSignificantUnitsFrCI95{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFr{iCond}{iArea}] = calc95CI((data.negativeSignificantUnits{iCond}{iArea}./data.totalUnits{iCond}{iArea}));
      [data.positiveMUAsFrCI95{iCond}{iArea}, ~, ~, data.positiveMUAsFr{iCond}{iArea}] = calc95CI((data.positiveMUAs{iCond}{iArea}./data.totalMUAs{iCond}{iArea}));
      %       [data.positiveSignificantMUAsFrCI95{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFr{iCond}{iArea}] = calc95CI((data.positiveSignificantMUAs{iCond}{iArea}./data.totalMUAs{iCond}{iArea}));
      [data.negativeMUAsFrCI95{iCond}{iArea}, ~, ~, data.negativeMUAsFr{iCond}{iArea}] = calc95CI((data.negativeMUAs{iCond}{iArea}./data.totalMUAs{iCond}{iArea}));
      %       [data.negativeSignificantMUAsFrCI95{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFr{iCond}{iArea}] = calc95CI((data.negativeSignificantMUAs{iCond}{iArea}./data.totalMUAs{iCond}{iArea}));
      [data.positiveUnitsFrCI95Thr{iCond}{iArea}, ~, ~, data.positiveUnitsFrThr{iCond}{iArea}] = calc95CI((data.positiveUnitsThr{iCond}{iArea}./data.totalUnitsThr{iCond}{iArea}));
      %       [data.positiveSignificantUnitsFrCI95Thr{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrThr{iCond}{iArea}] = calc95CI((data.positiveSignificantUnitsThr{iCond}{iArea}./data.totalUnitsThr{iCond}{iArea}));
      [data.negativeUnitsFrCI95Thr{iCond}{iArea}, ~, ~, data.negativeUnitsFrThr{iCond}{iArea}] = calc95CI((data.negativeUnitsThr{iCond}{iArea}./data.totalUnitsThr{iCond}{iArea}));
      %       [data.negativeSignificantUnitsFrCI95Thr{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrThr{iCond}{iArea}] = calc95CI((data.negativeSignificantUnitsThr{iCond}{iArea}./data.totalUnitsThr{iCond}{iArea}));
      [data.positiveMUAsFrCI95Thr{iCond}{iArea}, ~, ~, data.positiveMUAsFrThr{iCond}{iArea}] = calc95CI((data.positiveMUAsThr{iCond}{iArea}./data.totalMUAsThr{iCond}{iArea}));
      %       [data.positiveSignificantMUAsFrCI95Thr{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrThr{iCond}{iArea}] = calc95CI((data.positiveSignificantMUAsThr{iCond}{iArea}./data.totalMUAsThr{iCond}{iArea}));
      [data.negativeMUAsFrCI95Thr{iCond}{iArea}, ~, ~, data.negativeMUAsFrThr{iCond}{iArea}] = calc95CI((data.negativeMUAsThr{iCond}{iArea}./data.totalMUAsThr{iCond}{iArea}));
      %       [data.negativeSignificantMUAsFrCI95Thr{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrThr{iCond}{iArea}] = calc95CI((data.negativeSignificantMUAsThr{iCond}{iArea}./data.totalMUAsThr{iCond}{iArea}));
      
      [data.positiveUnitsFrCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFr10minWindows{iCond}{iArea}] = calc95CI((data.positiveUnits10minWindows{iCond}{iArea}./data.totalUnits10minWindows{iCond}{iArea}));
      %       [data.positiveSignificantUnitsFrCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFr10minWindows{iCond}{iArea}] = calc95CI((data.positiveSignificantUnits10minWindows{iCond}{iArea}./data.totalUnits10minWindows{iCond}{iArea}));
      [data.negativeUnitsFrCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFr10minWindows{iCond}{iArea}] = calc95CI((data.negativeUnits10minWindows{iCond}{iArea}./data.totalUnits10minWindows{iCond}{iArea}));
      %       [data.negativeSignificantUnitsFrCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFr10minWindows{iCond}{iArea}] = calc95CI((data.negativeSignificantUnits10minWindows{iCond}{iArea}./data.totalUnits10minWindows{iCond}{iArea}));
      [data.positiveMUAsFrCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFr10minWindows{iCond}{iArea}] = calc95CI((data.positiveMUAs10minWindows{iCond}{iArea}./data.totalMUAs10minWindows{iCond}{iArea}));
      %       [data.positiveSignificantMUAsFrCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFr10minWindows{iCond}{iArea}] = calc95CI((data.positiveSignificantMUAs10minWindows{iCond}{iArea}./data.totalMUAs10minWindows{iCond}{iArea}));
      [data.negativeMUAsFrCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFr10minWindows{iCond}{iArea}] = calc95CI((data.negativeMUAs10minWindows{iCond}{iArea}./data.totalMUAs10minWindows{iCond}{iArea}));
      %       [data.negativeSignificantMUAsFrCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFr10minWindows{iCond}{iArea}] = calc95CI((data.negativeSignificantMUAs10minWindows{iCond}{iArea}./data.totalMUAs10minWindows{iCond}{iArea}));
      [data.positiveUnitsFrCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrThr10minWindows{iCond}{iArea}] = calc95CI((data.positiveUnitsThr{iCond}{iArea}./data.totalUnitsThr{iCond}{iArea}));
      %       [data.positiveSignificantUnitsFrCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrThr10minWindows{iCond}{iArea}] = calc95CI((data.positiveSignificantUnitsThr10minWindows{iCond}{iArea}./data.totalUnitsThr10minWindows{iCond}{iArea}));
      [data.negativeUnitsFrCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrThr10minWindows{iCond}{iArea}] = calc95CI((data.negativeUnitsThr10minWindows{iCond}{iArea}./data.totalUnitsThr10minWindows{iCond}{iArea}));
      %       [data.negativeSignificantUnitsFrCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrThr10minWindows{iCond}{iArea}] = calc95CI((data.negativeSignificantUnitsThr10minWindows{iCond}{iArea}./data.totalUnitsThr10minWindows{iCond}{iArea}));
      [data.positiveMUAsFrCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrThr10minWindows{iCond}{iArea}] = calc95CI((data.positiveMUAsThr10minWindows{iCond}{iArea}./data.totalMUAsThr10minWindows{iCond}{iArea}));
      %       [data.positiveSignificantMUAsFrCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrThr10minWindows{iCond}{iArea}] = calc95CI((data.positiveSignificantMUAsThr10minWindows{iCond}{iArea}./data.totalMUAsThr10minWindows{iCond}{iArea}));
      [data.negativeMUAsFrCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrThr10minWindows{iCond}{iArea}] = calc95CI((data.negativeMUAsThr10minWindows{iCond}{iArea}./data.totalMUAsThr10minWindows{iCond}{iArea}));
      %       [data.negativeSignificantMUAsFrCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrThr10minWindows{iCond}{iArea}] = calc95CI((data.negativeSignificantMUAsThr10minWindows{iCond}{iArea}./data.totalMUAsThr10minWindows{iCond}{iArea}));
      
      [data.positiveUnitsFrCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFr20minWindows{iCond}{iArea}] = calc95CI((data.positiveUnits20minWindows{iCond}{iArea}./data.totalUnits20minWindows{iCond}{iArea}));
      %       [data.positiveSignificantUnitsFrCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFr20minWindows{iCond}{iArea}] = calc95CI((data.positiveSignificantUnits20minWindows{iCond}{iArea}./data.totalUnits20minWindows{iCond}{iArea}));
      [data.negativeUnitsFrCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFr20minWindows{iCond}{iArea}] = calc95CI((data.negativeUnits20minWindows{iCond}{iArea}./data.totalUnits20minWindows{iCond}{iArea}));
      %       [data.negativeSignificantUnitsFrCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFr20minWindows{iCond}{iArea}] = calc95CI((data.negativeSignificantUnits20minWindows{iCond}{iArea}./data.totalUnits20minWindows{iCond}{iArea}));
      [data.positiveMUAsFrCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFr20minWindows{iCond}{iArea}] = calc95CI((data.positiveMUAs20minWindows{iCond}{iArea}./data.totalMUAs20minWindows{iCond}{iArea}));
      %       [data.positiveSignificantMUAsFrCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFr20minWindows{iCond}{iArea}] = calc95CI((data.positiveSignificantMUAs20minWindows{iCond}{iArea}./data.totalMUAs20minWindows{iCond}{iArea}));
      [data.negativeMUAsFrCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFr20minWindows{iCond}{iArea}] = calc95CI((data.negativeMUAs20minWindows{iCond}{iArea}./data.totalMUAs20minWindows{iCond}{iArea}));
      %       [data.negativeSignificantMUAsFrCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFr20minWindows{iCond}{iArea}] = calc95CI((data.negativeSignificantMUAs20minWindows{iCond}{iArea}./data.totalMUAs20minWindows{iCond}{iArea}));
      [data.positiveUnitsFrCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrThr20minWindows{iCond}{iArea}] = calc95CI((data.positiveUnitsThr{iCond}{iArea}./data.totalUnitsThr{iCond}{iArea}));
      %       [data.positiveSignificantUnitsFrCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrThr20minWindows{iCond}{iArea}] = calc95CI((data.positiveSignificantUnitsThr20minWindows{iCond}{iArea}./data.totalUnitsThr20minWindows{iCond}{iArea}));
      [data.negativeUnitsFrCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrThr20minWindows{iCond}{iArea}] = calc95CI((data.negativeUnitsThr20minWindows{iCond}{iArea}./data.totalUnitsThr20minWindows{iCond}{iArea}));
      %       [data.negativeSignificantUnitsFrCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrThr20minWindows{iCond}{iArea}] = calc95CI((data.negativeSignificantUnitsThr20minWindows{iCond}{iArea}./data.totalUnitsThr20minWindows{iCond}{iArea}));
      [data.positiveMUAsFrCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrThr20minWindows{iCond}{iArea}] = calc95CI((data.positiveMUAsThr20minWindows{iCond}{iArea}./data.totalMUAsThr20minWindows{iCond}{iArea}));
      %       [data.positiveSignificantMUAsFrCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrThr20minWindows{iCond}{iArea}] = calc95CI((data.positiveSignificantMUAsThr20minWindows{iCond}{iArea}./data.totalMUAsThr20minWindows{iCond}{iArea}));
      [data.negativeMUAsFrCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrThr20minWindows{iCond}{iArea}] = calc95CI((data.negativeMUAsThr20minWindows{iCond}{iArea}./data.totalMUAsThr20minWindows{iCond}{iArea}));
      %       [data.negativeSignificantMUAsFrCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrThr20minWindows{iCond}{iArea}] = calc95CI((data.negativeSignificantMUAsThr20minWindows{iCond}{iArea}./data.totalMUAsThr20minWindows{iCond}{iArea}));
      
      [data.positiveUnitsFrCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFr30minWindows{iCond}{iArea}] = calc95CI((data.positiveUnits30minWindows{iCond}{iArea}./data.totalUnits30minWindows{iCond}{iArea}));
      %       [data.positiveSignificantUnitsFrCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFr30minWindows{iCond}{iArea}] = calc95CI((data.positiveSignificantUnits30minWindows{iCond}{iArea}./data.totalUnits30minWindows{iCond}{iArea}));
      [data.negativeUnitsFrCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFr30minWindows{iCond}{iArea}] = calc95CI((data.negativeUnits30minWindows{iCond}{iArea}./data.totalUnits30minWindows{iCond}{iArea}));
      %       [data.negativeSignificantUnitsFrCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFr30minWindows{iCond}{iArea}] = calc95CI((data.negativeSignificantUnits30minWindows{iCond}{iArea}./data.totalUnits30minWindows{iCond}{iArea}));
      [data.positiveMUAsFrCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFr30minWindows{iCond}{iArea}] = calc95CI((data.positiveMUAs30minWindows{iCond}{iArea}./data.totalMUAs30minWindows{iCond}{iArea}));
      %       [data.positiveSignificantMUAsFrCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFr30minWindows{iCond}{iArea}] = calc95CI((data.positiveSignificantMUAs30minWindows{iCond}{iArea}./data.totalMUAs30minWindows{iCond}{iArea}));
      [data.negativeMUAsFrCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFr30minWindows{iCond}{iArea}] = calc95CI((data.negativeMUAs30minWindows{iCond}{iArea}./data.totalMUAs30minWindows{iCond}{iArea}));
      %       [data.negativeSignificantMUAsFrCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFr30minWindows{iCond}{iArea}] = calc95CI((data.negativeSignificantMUAs30minWindows{iCond}{iArea}./data.totalMUAs30minWindows{iCond}{iArea}));
      [data.positiveUnitsFrCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrThr30minWindows{iCond}{iArea}] = calc95CI((data.positiveUnitsThr{iCond}{iArea}./data.totalUnitsThr{iCond}{iArea}));
      %       [data.positiveSignificantUnitsFrCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrThr30minWindows{iCond}{iArea}] = calc95CI((data.positiveSignificantUnitsThr30minWindows{iCond}{iArea}./data.totalUnitsThr30minWindows{iCond}{iArea}));
      [data.negativeUnitsFrCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrThr30minWindows{iCond}{iArea}] = calc95CI((data.negativeUnitsThr30minWindows{iCond}{iArea}./data.totalUnitsThr30minWindows{iCond}{iArea}));
      %       [data.negativeSignificantUnitsFrCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrThr30minWindows{iCond}{iArea}] = calc95CI((data.negativeSignificantUnitsThr30minWindows{iCond}{iArea}./data.totalUnitsThr30minWindows{iCond}{iArea}));
      [data.positiveMUAsFrCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrThr30minWindows{iCond}{iArea}] = calc95CI((data.positiveMUAsThr30minWindows{iCond}{iArea}./data.totalMUAsThr30minWindows{iCond}{iArea}));
      %       [data.positiveSignificantMUAsFrCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrThr30minWindows{iCond}{iArea}] = calc95CI((data.positiveSignificantMUAsThr30minWindows{iCond}{iArea}./data.totalMUAsThr30minWindows{iCond}{iArea}));
      [data.negativeMUAsFrCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrThr30minWindows{iCond}{iArea}] = calc95CI((data.negativeMUAsThr30minWindows{iCond}{iArea}./data.totalMUAsThr30minWindows{iCond}{iArea}));
      %       [data.negativeSignificantMUAsFrCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrThr30minWindows{iCond}{iArea}] = calc95CI((data.negativeSignificantMUAsThr30minWindows{iCond}{iArea}./data.totalMUAsThr30minWindows{iCond}{iArea}));
    end
  end
  
  % Save processed data
  if ~exist(mainFolder, 'dir')
    mkdir(mainFolder)
  end
  save([mainFolder filesep 'pupilCorrFractions.mat'],...
    'repository','unitThr','alpha','animals','conditions','areas','data','unitTable','muaTable', '-v7.3');
else
  load([mainFolder filesep 'pupilCorrFractions.mat']);
end


if fullRun <= 2
  
  %% Areas of interest
  if strcmp(repository, 'uol')
    area1 = 'VB'; area2 = 'S1'; area3 = 'RSC'; area4 = 'CA';
  elseif strcmp(repository, 'allensdk')
    area1 = 'VB'; area2 = 'LGN'; area3 = 'V1'; area4 = 'CA';
  end
  areasOI = {area1; area2; area3; area4};
  
  
  %% Bar plots for units
  if barPlots
    [fH, scatterGroups, areaGroups, ~, ~, corrGroups] = barPlotUnits(data, area1, area2, area3, area4); %#ok<*UNRCH>
    set(fH(1), 'Name','Proportion of units positively and negatively correlated to pupil size in different brain areas');
    filename = [mainFolder filesep 'pupilCorrFractionsUnits'];
    savefig(fH(1), filename, 'compact');
    print(fH(1), [filename '.png'],'-dpng','-r300');
    %close all
    
    figure(fH(2));
    %   yLim = ylim;
    %   if strcmp(repository, 'allensdk')
    %     ylim([yLim(1) 13]);
    %   else
    %     ylim([yLim(1) 10]);
    %   end
    %xLabel = get(gca, 'xlabel');
    yLabel = get(gca, 'ylabel');
    ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 40, 4/3, 2, [0.005 0], 'out',...
      'on', 'k', {}, xlim, xticks,...
      'on', 'k', {yLabel.String{1}}, ylim, yticks);
    set(fH(2), 'Name','Proportion of units positively and negatively correlated to pupil size in different brain areas');
    
    figFileName = [mainFolder filesep 'pupilCorrFractionsUnitsRedux'];
    figSize = 15;
    label = [4.5 1.8];
    margin = [0 0.8];
    width = ((100-23)/(147-119))*figSize-label(1)-margin(1);
    height = figSize-label(2)-margin(2);
    paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
    savefig(gcf, [figFileName '.fig'], 'compact');
    exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);
    exportFig(gcf, [figFileName '.eps'],'-depsc','-r1200', paperSize);
    %close all
    
    
    %% ANOVA for units:
    runAnova(scatterGroups, corrGroups, areaGroups, filename)
    
    
    %% Bar plots for MUAs:
    [fH, scatterGroups, areaGroups, ~, ~, corrGroups] = barPlotMUAs(data, area1, area2, area3, area4);
    set(fH(1), 'Name','Proportion of units and MUAs positively and negatively correlated to pupil size in different brain areas');
    filename = [mainFolder filesep 'pupilCorrFractionsMUAs'];
    savefig(fH(1), filename, 'compact');
    print(fH(1), [filename '.png'],'-dpng','-r300');
    %close all
    
    figure(fH(2));
    %   yLim = ylim;
    %   if strcmp(repository, 'allensdk')
    %     ylim([yLim(1) 13]);
    %   else
    %     ylim([yLim(1) 10]);
    %   end
    %xLabel = get(gca, 'xlabel');
    yLabel = get(gca, 'ylabel');
    ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 40, 4/3, 2, [0.005 0], 'out',...
      'on', 'k', {}, xlim, xticks,...
      'on', 'k', {yLabel.String{1}}, ylim, yticks);
    set(fH(2), 'Name','Proportion of units positively and negatively correlated to pupil size in different brain areas');
    
    figFileName = [mainFolder filesep 'pupilCorrFractionsMUAsRedux'];
    figSize = 15;
    label = [4.5 1.8];
    margin = [0 0.8];
    width = ((100-23)/(147-119))*figSize-label(1)-margin(1);
    height = figSize-label(2)-margin(2);
    paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
    savefig(gcf, [figFileName '.fig'], 'compact');
    exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);
    exportFig(gcf, [figFileName '.eps'],'-depsc','-r1200', paperSize);
    %close all
    
    
    %% ANOVA for MUAs:
    runAnova(scatterGroups, corrGroups, areaGroups, filename)
  end
  
  
  %% Fraction violin plots for units:
  if violinPlots
    if strcmp(repository, 'uol')
      area1 = 'Th'; area2 = 'VB'; area3 = 'Po'; area4 = 'S1'; area5 = 'RSC'; area6 = 'CA'; area7 = 'DG'; area8 = 'Cx'; area9 = 'Hp';
      areasOIFull = {area1; area4; area5; area6; area7};
      areasOIMerged = {area1; area8; area9};
    elseif strcmp(repository, 'allensdk')
      area1 = 'Th'; area2 = 'LGN'; area3 = 'LP'; area4 = 'V1'; area5 = 'V2'; area6 = 'CA'; area7 = 'DG'; area8 = 'VIS'; area9 = 'Hp';
      areasOIFull = {area1; area4; area5; area6; area7};
      areasOIMerged = {area1; area8; area9};
    end
    
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    violinAreas = {};
    violinData = {};
    options.nSample = {};
    iCond = 1;
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Pos'];
      violinData{(iArea*2)-1} = data.positiveUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(violinData{(iArea*2)-1});
      options.nSample{(iArea*2)-1} = data.totalUnitsThr{iCond}{area};
      violinAreas{(iArea*2)} = [areasOIFull{iArea} 'Neg'];
      violinData{(iArea*2)} = data.negativeUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(violinData{(iArea*2)});
      options.nSample{(iArea*2)} = data.totalUnitsThr{iCond}{area};
      statsFractionsUnits{iArea} = meanTest([violinData{(iArea*2)-1} violinData{(iArea*2)}], 'ANOVARM', 'off');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    if ~nSampleDisplay
      options.nSample = [];
    end
    options.violinVisibility = violinVisibility;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 12.75]);
    end
    set(fH, 'Name','Proportion of units positively and negatively correlated to pupil size in different brain areas');
    filename = [mainFolder filesep 'pupilCorrFractionsUnitsViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
    
    dataMean = zeros(1, numel(areasOIMerged)*2);
    dataCI95 = zeros(2, numel(areasOIMerged)*2);
    violinAreas = {};
    violinData = {};
    options.nSample = {};
    iCond = 1;
    for iArea = 1:numel(areasOIMerged)
      area = determineArea(areasOIMerged{iArea});
      violinAreas{(iArea*2)-1} = [areasOIMerged{iArea} 'Pos'];
      violinData{(iArea*2)-1} = data.positiveUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(violinData{(iArea*2)-1});
      options.nSample{(iArea*2)-1} = data.totalUnitsThr{iCond}{area};
      violinAreas{(iArea*2)} = [areasOIMerged{iArea} 'Neg'];
      violinData{(iArea*2)} = data.negativeUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(violinData{(iArea*2)});
      options.nSample{(iArea*2)} = data.totalUnitsThr{iCond}{area};
      statsFractionsUnitsMerged{iArea} = meanTest([violinData{(iArea*2)-1} violinData{(iArea*2)}], 'ANOVARM', 'off');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    if ~nSampleDisplay
      options.nSample = [];
    end
    options.violinVisibility = violinVisibility;
    fH2 = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 6.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 6.75]);
    end
    set(fH2, 'Name','Proportion of units positively and negatively correlated to pupil size in different brain areas');
    filename = [mainFolder filesep 'pupilCorrFractionsUnitsMergedViolins'];
    savefig(fH2, filename, 'compact');
    print(fH2, [filename '.png'],'-dpng','-r300');
    
    
    %% Fraction violin plots for MUAs:
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    violinAreas = {};
    violinData = {};
    options.nSample = {};
    iCond = 1;
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Pos'];
      violinData{(iArea*2)-1} = data.positiveMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(violinData{(iArea*2)-1});
      options.nSample{(iArea*2)-1} = data.totalMUAsThr{iCond}{area};
      violinAreas{(iArea*2)} = [areasOIFull{iArea} 'Neg'];
      violinData{(iArea*2)} = data.negativeMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(violinData{(iArea*2)});
      options.nSample{(iArea*2)} = data.totalMUAsThr{iCond}{area};
      statsFractionsMUAs{iArea} = meanTest([violinData{(iArea*2)-1} violinData{(iArea*2)}], 'ANOVARM', 'off');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    if ~nSampleDisplay
      options.nSample = [];
    end
    options.violinVisibility = violinVisibility;
    fH3 = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 12.75]);
    end
    set(fH3, 'Name','Proportion of units and MUAs positively and negatively correlated to pupil size in different brain areas');
    filename = [mainFolder filesep 'pupilCorrFractionsMUAsViolins'];
    savefig(fH3, filename, 'compact');
    print(fH3, [filename '.png'],'-dpng','-r300');
    
    dataMean = zeros(1, numel(areasOIMerged)*2);
    dataCI95 = zeros(2, numel(areasOIMerged)*2);
    violinAreas = {};
    violinData = {};
    options.nSample = {};
    iCond = 1;
    for iArea = 1:numel(areasOIMerged)
      area = determineArea(areasOIMerged{iArea});
      violinAreas{(iArea*2)-1} = [areasOIMerged{iArea} 'Pos'];
      violinData{(iArea*2)-1} = data.positiveMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(violinData{(iArea*2)-1});
      options.nSample{(iArea*2)-1} = data.totalMUAsThr{iCond}{area};
      violinAreas{(iArea*2)} = [areasOIMerged{iArea} 'Neg'];
      violinData{(iArea*2)} = data.negativeMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(violinData{(iArea*2)});
      options.nSample{(iArea*2)} = data.totalMUAsThr{iCond}{area};
      statsFractionsMUAsMerged{iArea} = meanTest([violinData{(iArea*2)-1} violinData{(iArea*2)}], 'ANOVARM', 'off');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    if ~nSampleDisplay
      options.nSample = [];
    end
    options.violinVisibility = violinVisibility;
    fH4 = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 6.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 6.75]);
    end
    set(fH4, 'Name','Proportion of units and MUAs positively and negatively correlated to pupil size in different brain areas');
    filename = [mainFolder filesep 'pupilCorrFractionsMUAsMergedViolins'];
    savefig(fH4, filename, 'compact');
    print(fH4, [filename '.png'],'-dpng','-r300');
    
    
    %% Partial (positive only) fraction violin plots for units:
    dataMean = zeros(1, numel(areasOIFull));
    dataCI95 = zeros(2, numel(areasOIFull));
    violinAreas = {};
    violinData = {};
    options.nSample = {};
    iCond = 1;
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      violinAreas{iArea} = [areasOIFull{iArea} 'Pos'];
      violinData{iArea} = data.positiveUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
      [dataMean(iArea), dataCI95(:,iArea)] = datamean(violinData{iArea});
      options.nSample{iArea} = data.totalUnitsThr{iCond}{area};
    end
    statsFractionsUnitsPositive = meanTest2(violinData, 'ttest');
    violinDataSubtracted = violinData;
    for iArea = 1:numel(violinDataSubtracted)
      violinDataSubtracted{iArea} = violinDataSubtracted{iArea} - 0.5;
    end
    statsFractionsUnitsPositiveZero = meanTestZero(violinDataSubtracted, 'ttest');
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    if ~nSampleDisplay
      options.nSample = [];
    end
    options.violinVisibility = violinVisibility;
    fH5 = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 6.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 6.75]);
    end
    set(fH5, 'Name','Proportion of units positively correlated to pupil size in different brain areas');
    filename = [mainFolder filesep 'pupilCorrFractionsUnitsViolinsPositive'];
    savefig(fH5, filename, 'compact');
    print(fH5, [filename '.png'],'-dpng','-r300');
    
    dataMean = zeros(1, numel(areasOIMerged));
    dataCI95 = zeros(2, numel(areasOIMerged));
    violinAreas = {};
    violinData = {};
    options.nSample = {};
    iCond = 1;
    for iArea = 1:numel(areasOIMerged)
      area = determineArea(areasOIMerged{iArea});
      violinAreas{iArea} = [areasOIMerged{iArea} 'Pos'];
      violinData{iArea} = data.positiveUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
      [dataMean(iArea), dataCI95(:,iArea)] = datamean(violinData{iArea});
      options.nSample{iArea} = data.totalUnitsThr{iCond}{area};
    end
    statsFractionsUnitsPositiveMerged = meanTest2(violinData, 'ttest');
    violinDataSubtracted = violinData;
    for iArea = 1:numel(violinDataSubtracted)
      violinDataSubtracted{iArea} = violinDataSubtracted{iArea} - 0.5;
    end
    statsFractionsUnitsPositiveZeroMerged = meanTestZero(violinDataSubtracted, 'ttest');
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    if ~nSampleDisplay
      options.nSample = [];
    end
    options.violinVisibility = violinVisibility;
    fH6 = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 3.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 3.75]);
    end
    set(fH6, 'Name','Proportion of units positively correlated to pupil size in different brain areas');
    filename = [mainFolder filesep 'pupilCorrFractionsUnitsMergedViolinsPositive'];
    savefig(fH6, filename, 'compact');
    print(fH6, [filename '.png'],'-dpng','-r300');
    
    if strcmpi(repository, 'uol')
      dataMean = zeros(1, numel(areasOIFull));
      dataCI95 = zeros(2, numel(areasOIFull));
      violinAreas = {};
      violinData = {};
      options.nSample = {};
      iCond = 2;
      for iArea = 1:numel(areasOIFull)
        area = determineArea(areasOIFull{iArea});
        violinAreas{iArea} = [areasOIFull{iArea} 'Pos'];
        violinData{iArea} = data.positiveUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
        [dataMean(iArea), dataCI95(:,iArea)] = datamean(violinData{iArea});
        options.nSample{iArea} = data.totalUnitsThr{iCond}{area};
      end
      statsFractionsUnitsPositiveAnaest = meanTest2(violinData, 'ttest');
      violinDataSubtracted = violinData;
      for iArea = 1:numel(violinDataSubtracted)
        violinDataSubtracted{iArea} = violinDataSubtracted{iArea} - 0.5;
      end
      statsFractionsUnitsPositiveZeroAnaest = meanTestZero(violinDataSubtracted, 'ttest');
      options.yLim = [0 1];
      options.yLabel = 'Fraction';
      options.showNotches = false;
      options.medianPlot = false;
      if ~nSampleDisplay
        options.nSample = [];
      end
      options.violinVisibility = violinVisibility;
      fH7 = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
      if strcmp(repository, 'uol')
        xlim([0.25 6.75]);
      elseif strcmp(repository, 'allensdk')
        xlim([0.25 6.75]);
      end
      set(fH7, 'Name','Proportion of units positively correlated to pupil size in different brain areas');
      filename = [mainFolder filesep 'pupilCorrFractionsUnitsViolinsPositiveAnaesthesia'];
      savefig(fH7, filename, 'compact');
      print(fH7, [filename '.png'],'-dpng','-r300');
      
      dataMean = zeros(1, numel(areasOIMerged));
      dataCI95 = zeros(2, numel(areasOIMerged));
      violinAreas = {};
      violinData = {};
      options.nSample = {};
      iCond = 2;
      for iArea = 1:numel(areasOIMerged)
        area = determineArea(areasOIMerged{iArea});
        violinAreas{iArea} = [areasOIMerged{iArea} 'Pos'];
        violinData{iArea} = data.positiveUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
        [dataMean(iArea), dataCI95(:,iArea)] = datamean(violinData{iArea});
        options.nSample{iArea} = data.totalUnitsThr{iCond}{area};
      end
      statsFractionsUnitsPositiveAnaestMerged = meanTest2(violinData, 'ttest');
      violinDataSubtracted = violinData;
      for iArea = 1:numel(violinDataSubtracted)
        violinDataSubtracted{iArea} = violinDataSubtracted{iArea} - 0.5;
      end
      statsFractionsUnitsPositiveZeroAnaestMerged = meanTestZero(violinDataSubtracted, 'ttest');
      options.yLim = [0 1];
      options.yLabel = 'Fraction';
      options.showNotches = false;
      options.medianPlot = false;
      if ~nSampleDisplay
        options.nSample = [];
      end
      options.violinVisibility = violinVisibility;
      fH8 = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
      if strcmp(repository, 'uol')
        xlim([0.25 3.75]);
      elseif strcmp(repository, 'allensdk')
        xlim([0.25 3.75]);
      end
      set(fH8, 'Name','Proportion of units positively correlated to pupil size in different brain areas');
      filename = [mainFolder filesep 'pupilCorrFractionsUnitsMergedViolinsPositiveAnaesthesia'];
      savefig(fH8, filename, 'compact');
      print(fH8, [filename '.png'],'-dpng','-r300');
    end
    
    
    %% Partial (positive only) fraction violin plots for MUAs:
    dataMean = zeros(1, numel(areasOIFull));
    dataCI95 = zeros(2, numel(areasOIFull));
    violinAreas = {};
    violinData = {};
    options.nSample{iArea} = {};
    iCond = 1;
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      violinAreas{iArea} = [areasOIFull{iArea} 'Pos'];
      violinData{iArea} = data.positiveMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
      [dataMean(iArea), dataCI95(:,iArea)] = datamean(violinData{iArea});
      options.nSample{iArea} = data.totalMUAsThr{iCond}{area};
    end
    statsFractionsMUAsPositive = meanTest2(violinData, 'ttest'); % we care about this
    violinDataSubtracted = violinData;
    for iArea = 1:numel(violinDataSubtracted)
      violinDataSubtracted{iArea} = violinDataSubtracted{iArea} - 0.5;
    end
    statsFractionsMUAsPositiveZero = meanTestZero(violinDataSubtracted, 'ttest'); % we care about this
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    options.markerFaceAlpha = 1;
    if ~nSampleDisplay
      options.nSample = [];
    end
    options.violinVisibility = violinVisibility;
    fH9 = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 6.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 6.75]);
    end
    yticks([0 0.2 0.4 0.6 0.8 1]);
    set(fH9, 'Name','Proportion of units and MUAs positively correlated to pupil size in different brain areas');
    filename = [mainFolder filesep 'pupilCorrFractionsMUAsViolinsPositive'];
    savefig(fH9, filename, 'compact');
    print(fH9, [filename '.png'],'-dpng','-r300');
    
    dataMean = zeros(1, numel(areasOIMerged));
    dataCI95 = zeros(2, numel(areasOIMerged));
    violinAreas = {};
    violinData = {};
    options.nSample = {};
    iCond = 1;
    for iArea = 1:numel(areasOIMerged)
      area = determineArea(areasOIMerged{iArea});
      violinAreas{iArea} = [areasOIMerged{iArea} 'Pos'];
      violinData{iArea} = data.positiveMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
      [dataMean(iArea), dataCI95(:,iArea)] = datamean(violinData{iArea});
      options.nSample{iArea} = data.totalMUAsThr{iCond}{area};
    end
    statsFractionsMUAsPositiveMerged = meanTest2(violinData, 'ttest'); % we care about this
    violinDataSubtracted = violinData;
    for iArea = 1:numel(violinDataSubtracted)
      violinDataSubtracted{iArea} = violinDataSubtracted{iArea} - 0.5;
    end
    statsFractionsMUAsPositiveZeroMerged = meanTestZero(violinDataSubtracted, 'ttest'); % we care about this
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    options.markerFaceAlpha = 1;
    if ~nSampleDisplay
      options.nSample = [];
    end
    options.violinVisibility = violinVisibility;
    fH10 = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 3.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 3.75]);
    end
    yticks([0 0.2 0.4 0.6 0.8 1]);
    set(fH10, 'Name','Proportion of units and MUAs positively correlated to pupil size in different brain areas');
    filename = [mainFolder filesep 'pupilCorrFractionsMUAsMergedViolinsPositive'];
    savefig(fH10, filename, 'compact');
    print(fH10, [filename '.png'],'-dpng','-r300');
    
    if strcmpi(repository, 'uol')
      dataMean = zeros(1, numel(areasOIFull));
      dataCI95 = zeros(2, numel(areasOIFull));
      violinAreas = {};
      violinData = {};
      options.nSample = {};
      iCond = 2;
      for iArea = 1:numel(areasOIFull)
        area = determineArea(areasOIFull{iArea});
        violinAreas{iArea} = [areasOIFull{iArea} 'Pos'];
        violinData{iArea} = data.positiveMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
        [dataMean(iArea), dataCI95(:,iArea)] = datamean(violinData{iArea});
        options.nSample{iArea} = data.totalMUAsThr{iCond}{area};
      end
      statsFractionsMUAsPositiveAnaest = meanTest2(violinData, 'ttest');
      violinDataSubtracted = violinData;
      for iArea = 1:numel(violinDataSubtracted)
        violinDataSubtracted{iArea} = violinDataSubtracted{iArea} - 0.5;
      end
      statsFractionsMUAsPositiveZeroAnaest = meanTestZero(violinDataSubtracted, 'ttest');
      options.yLim = [0 1];
      options.yLabel = 'Fraction';
      options.showNotches = false;
      options.medianPlot = false;
      options.markerFaceAlpha = 1;
      if ~nSampleDisplay
        options.nSample = [];
      end
      options.violinVisibility = violinVisibility;
      fH11 = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
      if strcmp(repository, 'uol')
        xlim([0.25 6.75]);
      elseif strcmp(repository, 'allensdk')
        xlim([0.25 6.75]);
      end
      yticks([0 0.2 0.4 0.6 0.8 1]);
      set(fH11, 'Name','Proportion of units and MUAs positively correlated to pupil size in different brain areas');
      filename = [mainFolder filesep 'pupilCorrFractionsMUAsViolinsPositiveAnaesthesia'];
      savefig(fH11, filename, 'compact');
      print(fH11, [filename '.png'],'-dpng','-r300');
      
      dataMean = zeros(1, numel(areasOIMerged));
      dataCI95 = zeros(2, numel(areasOIMerged));
      violinAreas = {};
      violinData = {};
      options.nSample = {};
      iCond = 2;
      for iArea = 1:numel(areasOIMerged)
        area = determineArea(areasOIMerged{iArea});
        violinAreas{iArea} = [areasOIMerged{iArea} 'Pos'];
        violinData{iArea} = data.positiveMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
        [dataMean(iArea), dataCI95(:,iArea)] = datamean(violinData{iArea});
        options.nSample{iArea} = data.totalMUAsThr{iCond}{area};
      end
      statsFractionsMUAsPositiveAnaestMerged = meanTest2(violinData, 'ttest');
      violinDataSubtracted = violinData;
      for iArea = 1:numel(violinDataSubtracted)
        violinDataSubtracted{iArea} = violinDataSubtracted{iArea} - 0.5;
      end
      statsFractionsMUAsPositiveZeroAnaestMerged = meanTestZero(violinDataSubtracted, 'ttest');
      options.yLim = [0 1];
      options.yLabel = 'Fraction';
      options.showNotches = false;
      options.medianPlot = false;
      options.markerFaceAlpha = 1;
      if ~nSampleDisplay
        options.nSample = [];
      end
      options.violinVisibility = violinVisibility;
      fH12 = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
      if strcmp(repository, 'uol')
        xlim([0.25 3.75]);
      elseif strcmp(repository, 'allensdk')
        xlim([0.25 3.75]);
      end
      yticks([0 0.2 0.4 0.6 0.8 1]);
      set(fH12, 'Name','Proportion of units and MUAs positively correlated to pupil size in different brain areas');
      filename = [mainFolder filesep 'pupilCorrFractionsMUAsMergedViolinsPositiveAnaesthesia'];
      savefig(fH12, filename, 'compact');
      print(fH12, [filename '.png'],'-dpng','-r300');
    end
    
    
%     %% Significant fraction violin plots for positive units comparing wakefulness to anaesthesia:
%     if strcmp(repository, 'uol')
%       area1 = 'Th'; area2 = 'VB'; area3 = 'Po'; area4 = 'S1'; area5 = 'RSC'; area6 = 'DG'; area7 = 'CA';
%       areasOIFull = {area1; area2; area3; area4; area5; area6; area7};
%       dataMean = zeros(1, numel(areasOIFull)*2);
%       dataCI95 = zeros(2, numel(areasOIFull)*2);
%       for iArea = 1:numel(areasOIFull)
%         area = determineArea(areasOIFull{iArea});
%         iCond = 1;
%         violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
%         violinData{(iArea*2)-1} = data.positiveSignificantUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
%         [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(violinData{(iArea*2)-1});
%         iCond = 2;
%         violinAreas{(iArea*2)} = [areasOIFull{iArea} 'Anaest'];
%         violinData{(iArea*2)} = data.positiveSignificantUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
%         [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(violinData{(iArea*2)});
%         [~, statsFractionsUnitsPositiveSignificant{iArea}.p] = ttest2(violinData{(iArea*2)-1}, violinData{(iArea*2)});
%       end
%       options.yLim = [0 1];
%       options.yLabel = 'Fraction';
%       options.showNotches = false;
%       options.medianPlot = false;
%       fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
%       if strcmp(repository, 'uol')
%         xlim([0.25 12.75]);
%       elseif strcmp(repository, 'allensdk')
%         xlim([0.25 10.75]);
%       end
%       set(fH, 'Name','Proportion of significant positive units in different brain areas across conditions');
%       filename = [mainFolder filesep 'pupilCorrFractionsSignificantPositiveUnitsViolins'];
%       savefig(fH, filename, 'compact');
%       print(fH, [filename '.png'],'-dpng','-r300');
%     else
%       statsFractionsUnitsPositiveSignificant = [];
%     end
%     
%     
%     %% Significant fraction violin plots for positive units+MUAs comparing wakefulness to anaesthesia:
%     if strcmp(repository, 'uol')
%       area1 = 'Th'; area2 = 'VB'; area3 = 'Po'; area4 = 'S1'; area5 = 'RSC'; area6 = 'DG'; area7 = 'CA';
%       areasOIFull = {area1; area2; area3; area4; area5; area6; area7};
%       dataMean = zeros(1, numel(areasOIFull)*2);
%       dataCI95 = zeros(2, numel(areasOIFull)*2);
%       for iArea = 1:numel(areasOIFull)
%         area = determineArea(areasOIFull{iArea});
%         iCond = 1;
%         violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
%         violinData{(iArea*2)-1} = data.positiveSignificantMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
%         [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(violinData{(iArea*2)-1});
%         iCond = 2;
%         violinAreas{(iArea*2)} = [areasOIFull{iArea} 'Anaest'];
%         violinData{(iArea*2)} = data.positiveSignificantMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
%         [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(violinData{(iArea*2)});
%         [~, statsFractionsMUAsPositiveSignificant{iArea}.p] = ttest2(violinData{(iArea*2)-1}, violinData{(iArea*2)});
%       end
%       options.yLim = [0 1];
%       options.yLabel = 'Fraction';
%       options.showNotches = false;
%       options.medianPlot = false;
%       fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
%       if strcmp(repository, 'uol')
%         xlim([0.25 12.75]);
%       elseif strcmp(repository, 'allensdk')
%         xlim([0.25 10.75]);
%       end
%       set(fH, 'Name','Proportion of significant positive units+MUAs in different brain areas across conditions');
%       filename = [mainFolder filesep 'pupilCorrFractionsSignificantPositiveMUAsViolins'];
%       savefig(fH, filename, 'compact');
%       print(fH, [filename '.png'],'-dpng','-r300');
%     else
%       statsFractionsMUAsPositiveSignificant = [];
%     end
%     
%     
%     %% Significant fraction violin plots for negative units comparing wakefulness to anaesthesia:
%     if strcmp(repository, 'uol')
%       area1 = 'Th'; area2 = 'VB'; area3 = 'Po'; area4 = 'S1'; area5 = 'RSC'; area6 = 'DG'; area7 = 'CA';
%       areasOIFull = {area1; area2; area3; area4; area5; area6; area7};
%       dataMean = zeros(1, numel(areasOIFull)*2);
%       dataCI95 = zeros(2, numel(areasOIFull)*2);
%       for iArea = 1:numel(areasOIFull)
%         area = determineArea(areasOIFull{iArea});
%         iCond = 1;
%         violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
%         violinData{(iArea*2)-1} = data.negativeSignificantUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
%         [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(violinData{(iArea*2)-1});
%         iCond = 2;
%         violinAreas{(iArea*2)} = [areasOIFull{iArea} 'Anaest'];
%         violinData{(iArea*2)} = data.negativeSignificantUnitsThr{iCond}{area}./data.totalUnitsThr{iCond}{area};
%         [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(violinData{(iArea*2)});
%         [~, statsFractionsUnitsNegativeSignificant{iArea}.p] = ttest2(violinData{(iArea*2)-1}, violinData{(iArea*2)});
%       end
%       options.yLim = [0 1];
%       options.yLabel = 'Fraction';
%       options.showNotches = false;
%       options.medianPlot = false;
%       fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
%       if strcmp(repository, 'uol')
%         xlim([0.25 12.75]);
%       elseif strcmp(repository, 'allensdk')
%         xlim([0.25 10.75]);
%       end
%       set(fH, 'Name','Proportion of significant negative units in different brain areas across conditions');
%       filename = [mainFolder filesep 'pupilCorrFractionsSignificantNegativeUnitsViolins'];
%       savefig(fH, filename, 'compact');
%       print(fH, [filename '.png'],'-dpng','-r300');
%     else
%       statsFractionsUnitsNegativeSignificant = [];
%     end
%     
%     
%     %% Significant fraction violin plots for negative units+MUAs comparing wakefulness to anaesthesia:
%     if strcmp(repository, 'uol')
%       area1 = 'Th'; area2 = 'VB'; area3 = 'Po'; area4 = 'S1'; area5 = 'RSC'; area6 = 'DG'; area7 = 'CA';
%       areasOIFull = {area1; area2; area3; area4; area5; area6; area7};
%       dataMean = zeros(1, numel(areasOIFull)*2);
%       dataCI95 = zeros(2, numel(areasOIFull)*2);
%       for iArea = 1:numel(areasOIFull)
%         area = determineArea(areasOIFull{iArea});
%         iCond = 1;
%         violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
%         violinData{(iArea*2)-1} = data.negativeSignificantMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
%         [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(violinData{(iArea*2)-1});
%         iCond = 2;
%         violinAreas{(iArea*2)} = [areasOIFull{iArea} 'Anaest'];
%         violinData{(iArea*2)} = data.negativeSignificantMUAsThr{iCond}{area}./data.totalMUAsThr{iCond}{area};
%         [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(violinData{(iArea*2)});
%         [~, statsFractionsMUAsNegativeSignificant{iArea}.p] = ttest2(violinData{(iArea*2)-1}, violinData{(iArea*2)});
%       end
%       options.yLim = [0 1];
%       options.yLabel = 'Fraction';
%       options.showNotches = false;
%       options.medianPlot = false;
%       fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
%       if strcmp(repository, 'uol')
%         xlim([0.25 12.75]);
%       elseif strcmp(repository, 'allensdk')
%         xlim([0.25 10.75]);
%       end
%       set(fH, 'Name','Proportion of significant negative units+MUAs in different brain areas across conditions');
%       filename = [mainFolder filesep 'pupilCorrFractionsSignificantNegativeMUAsViolins'];
%       savefig(fH, filename, 'compact');
%       print(fH, [filename '.png'],'-dpng','-r300');
%     else
%       statsFractionsMUAsNegativeSignificant = [];
%     end
  end
  
  
  %% Correlation between fractions in different brain areas
  if drawCorrs
    windowSize = [];
    [~, rPositiveUnits, pvalPositiveUnits] = corrPlot(areasOI, data, 'positiveUnits', windowSize, mainFolder);
    %   [~, rPositiveSignificantUnits, pvalPositiveSignificantUnits] = corrPlot(areasOI, data, 'positiveSignificantUnits', windowSize, mainFolder);
    [~, rNegativeUnits, pvalNegativeUnits] = corrPlot(areasOI, data, 'negativeUnits', windowSize, mainFolder);
    %   [~, rNegativeSignificantUnits, pvalNegativeSignificantUnits] = corrPlot(areasOI, data, 'negativeSignificantUnits', windowSize, mainFolder);
    %   [~, rNeutralUnits, pvalNeutralUnits] = corrPlot(areasOI, data, 'neutralUnits', windowSize, mainFolder);
    [~, rPositiveMUAs, pvalPositiveMUAs] = corrPlot(areasOI, data, 'positiveMUAs', windowSize, mainFolder);
    %   [~, rPositiveSignificantMUAs, pvalPositiveSignificantMUAs] = corrPlot(areasOI, data, 'positiveSignificantMUAs', windowSize, mainFolder);
    [~, rNegativeMUAs, pvalNegativeMUAs] = corrPlot(areasOI, data, 'negativeMUAs', windowSize, mainFolder);
    %   [~, rNegativeSignificantMUAs, pvalNegativeSignificantMUAs] = corrPlot(areasOI, data, 'negativeSignificantMUAs', windowSize, mainFolder);
    %   [~, rNeutralMUAs, pvalNeutralMUAs] = corrPlot(areasOI, data, 'neutralMUAs', windowSize, mainFolder);
    
    windowSize = 10;
    [~, rPositiveUnits10minWindows, pvalPositiveUnits10minWindows] = corrPlot(areasOI, data, 'positiveUnits', windowSize, mainFolder);
    %   [~, rPositiveSignificantUnits10minWindows, pvalPositiveSignificantUnits10minWindows] = corrPlot(areasOI, data, 'positiveSignificantUnits', windowSize, mainFolder);
    [~, rNegativeUnits10minWindows, pvalNegativeUnits10minWindows] = corrPlot(areasOI, data, 'negativeUnits', windowSize, mainFolder);
    %   [~, rNegativeSignificantUnits10minWindows, pvalNegativeSignificantUnits10minWindows] = corrPlot(areasOI, data, 'negativeSignificantUnits', windowSize, mainFolder);
    %   [~, rNeutralUnits10minWindows, pvalNeutralUnits10minWindows] = corrPlot(areasOI, data, 'neutralUnits', windowSize, mainFolder);
    [~, rPositiveMUAs10minWindows, pvalPositiveMUAs10minWindows] = corrPlot(areasOI, data, 'positiveMUAs', windowSize, mainFolder);
    %   [~, rPositiveSignificantMUAs10minWindows, pvalPositiveSignificantMUAs10minWindows] = corrPlot(areasOI, data, 'positiveSignificantMUAs', windowSize, mainFolder);
    [~, rNegativeMUAs10minWindows, pvalNegativeMUAs10minWindows] = corrPlot(areasOI, data, 'negativeMUAs', windowSize, mainFolder);
    %   [~, rNegativeSignificantMUAs10minWindows, pvalNegativeSignificantMUAs10minWindows] = corrPlot(areasOI, data, 'negativeSignificantMUAs', windowSize, mainFolder);
    %   [~, rNeutralMUAs10minWindows, pvalNeutralMUAs10minWindows] = corrPlot(areasOI, data, 'neutralMUAs', windowSize, mainFolder);
    
    windowSize = 20;
    [~, rPositiveUnits20minWindows, pvalPositiveUnits20minWindows] = corrPlot(areasOI, data, 'positiveUnits', windowSize, mainFolder);
    %   [~, rPositiveSignificantUnits20minWindows, pvalPositiveSignificantUnits20minWindows] = corrPlot(areasOI, data, 'positiveSignificantUnits', windowSize, mainFolder);
    [~, rNegativeUnits20minWindows, pvalNegativeUnits20minWindows] = corrPlot(areasOI, data, 'negativeUnits', windowSize, mainFolder);
    %   [~, rNegativeSignificantUnits20minWindows, pvalNegativeSignificantUnits20minWindows] = corrPlot(areasOI, data, 'negativeSignificantUnits', windowSize, mainFolder);
    %   [~, rNeutralUnits20minWindows, pvalNeutralUnits20minWindows] = corrPlot(areasOI, data, 'neutralUnits', windowSize, mainFolder);
    [~, rPositiveMUAs20minWindows, pvalPositiveMUAs20minWindows] = corrPlot(areasOI, data, 'positiveMUAs', windowSize, mainFolder);
    %   [~, rPositiveSignificantMUAs20minWindows, pvalPositiveSignificantMUAs20minWindows] = corrPlot(areasOI, data, 'positiveSignificantMUAs', windowSize, mainFolder);
    [~, rNegativeMUAs20minWindows, pvalNegativeMUAs20minWindows] = corrPlot(areasOI, data, 'negativeMUAs', windowSize, mainFolder);
    %   [~, rNegativeSignificantMUAs20minWindows, pvalNegativeSignificantMUAs20minWindows] = corrPlot(areasOI, data, 'negativeSignificantMUAs', windowSize, mainFolder);
    %   [~, rNeutralMUAs20minWindows, pvalNeutralMUAs20minWindows] = corrPlot(areasOI, data, 'neutralMUAs', windowSize, mainFolder);
    
    windowSize = 30;
    [~, rPositiveUnits30minWindows, pvalPositiveUnits30minWindows] = corrPlot(areasOI, data, 'positiveUnits', windowSize, mainFolder);
    %   [~, rPositiveSignificantUnits30minWindows, pvalPositiveSignificantUnits30minWindows] = corrPlot(areasOI, data, 'positiveSignificantUnits', windowSize, mainFolder);
    [~, rNegativeUnits30minWindows, pvalNegativeUnits30minWindows] = corrPlot(areasOI, data, 'negativeUnits', windowSize, mainFolder);
    %   [~, rNegativeSignificantUnits30minWindows, pvalNegativeSignificantUnits30minWindows] = corrPlot(areasOI, data, 'negativeSignificantUnits', windowSize, mainFolder);
    %   [~, rNeutralUnits30minWindows, pvalNeutralUnits30minWindows] = corrPlot(areasOI, data, 'neutralUnits', windowSize, mainFolder);
    [~, rPositiveMUAs30minWindows, pvalPositiveMUAs30minWindows] = corrPlot(areasOI, data, 'positiveMUAs', windowSize, mainFolder);
    %   [~, rPositiveSignificantMUAs30minWindows, pvalPositiveSignificantMUAs30minWindows] = corrPlot(areasOI, data, 'positiveSignificantMUAs', windowSize, mainFolder);
    [~, rNegativeMUAs30minWindows, pvalNegativeMUAs30minWindows] = corrPlot(areasOI, data, 'negativeMUAs', windowSize, mainFolder);
    %   [~, rNegativeSignificantMUAs30minWindows, pvalNegativeSignificantMUAs30minWindows] = corrPlot(areasOI, data, 'negativeSignificantMUAs', windowSize, mainFolder);
    %   [~, rNeutralMUAs30minWindows, pvalNeutralMUAs30minWindows] = corrPlot(areasOI, data, 'neutralMUAs', windowSize, mainFolder);
    
    
    %% Correlation between fractions and pupil size
    windowSize = [];
    [~, rPositiveUnitsPupil, pvalPositiveUnitsPupil] = corrPlotPupil(areasOI, data, 'positiveUnits', windowSize, mainFolder);
    %   [~, rPositiveSignificantUnitsPupil, pvalPositiveSignificantUnitsPupil] = corrPlotPupil(areasOI, data, 'positiveSignificantUnits', windowSize, mainFolder);
    [~, rNegativeUnitsPupil, pvalNegativeUnitsPupil] = corrPlotPupil(areasOI, data, 'negativeUnits', windowSize, mainFolder);
    %   [~, rNegativeSignificantUnitsPupil, pvalNegativeSignificantUnitsPupil] = corrPlotPupil(areasOI, data, 'negativeSignificantUnits', windowSize, mainFolder);
    %   [~, rNeutralUnitsPupil, pvalNeutralUnitsPupil] = corrPlotPupil(areasOI, data, 'neutralUnits', windowSize, mainFolder);
    [~, rPositiveMUAsPupil, pvalPositiveMUAsPupil] = corrPlotPupil(areasOI, data, 'positiveMUAs', windowSize, mainFolder);
    %   [~, rPositiveSignificantMUAsPupil, pvalPositiveSignificantMUAsPupil] = corrPlotPupil(areasOI, data, 'positiveSignificantMUAs', windowSize, mainFolder);
    [~, rNegativeMUAsPupil, pvalNegativeMUAsPupil] = corrPlotPupil(areasOI, data, 'negativeMUAs', windowSize, mainFolder);
    %   [~, rNegativeSignificantMUAsPupil, pvalNegativeSignificantMUAsPupil] = corrPlotPupil(areasOI, data, 'negativeSignificantMUAs', windowSize, mainFolder);
    %   [~, rNeutralMUAsPupil, pvalNeutralMUAsPupil] = corrPlotPupil(areasOI, data, 'neutralMUAs', windowSize, mainFolder);
    
    windowSize = 10;
    [~, rPositiveUnitsPupil10minWindows, pvalPositiveUnitsPupil10minWindows] = corrPlotPupil(areasOI, data, 'positiveUnits', windowSize, mainFolder);
    %   [~, rPositiveSignificantUnitsPupil10minWindows, pvalPositiveSignificantUnitsPupil10minWindows] = corrPlotPupil(areasOI, data, 'positiveSignificantUnits', windowSize, mainFolder);
    [~, rNegativeUnitsPupil10minWindows, pvalNegativeUnitsPupil10minWindows] = corrPlotPupil(areasOI, data, 'negativeUnits', windowSize, mainFolder);
    %   [~, rNegativeSignificantUnitsPupil10minWindows, pvalNegativeSignificantUnitsPupil10minWindows] = corrPlotPupil(areasOI, data, 'negativeSignificantUnits', windowSize, mainFolder);
    %   [~, rNeutralUnitsPupil10minWindows, pvalNeutralUnitsPupil10minWindows] = corrPlotPupil(areasOI, data, 'neutralUnits', windowSize, mainFolder);
    [~, rPositiveMUAsPupil10minWindows, pvalPositiveMUAsPupil10minWindows] = corrPlotPupil(areasOI, data, 'positiveMUAs', windowSize, mainFolder);
    %   [~, rPositiveSignificantMUAsPupil10minWindows, pvalPositiveSignificantMUAsPupil10minWindows] = corrPlotPupil(areasOI, data, 'positiveSignificantMUAs', windowSize, mainFolder);
    [~, rNegativeMUAsPupil10minWindows, pvalNegativeMUAsPupil10minWindows] = corrPlotPupil(areasOI, data, 'negativeMUAs', windowSize, mainFolder);
    %   [~, rNegativeSignificantMUAsPupil10minWindows, pvalNegativeSignificantMUAsPupil10minWindows] = corrPlotPupil(areasOI, data, 'negativeSignificantMUAs', windowSize, mainFolder);
    %   [~, rNeutralMUAsPupil10minWindows, pvalNeutralMUAsPupil10minWindows] = corrPlotPupil(areasOI, data, 'neutralMUAs', windowSize, mainFolder);
    
    windowSize = 20;
    [~, rPositiveUnitsPupil20minWindows, pvalPositiveUnitsPupil20minWindows] = corrPlotPupil(areasOI, data, 'positiveUnits', windowSize, mainFolder);
    %   [~, rPositiveSignificantUnitsPupil20minWindows, pvalPositiveSignificantUnitsPupil20minWindows] = corrPlotPupil(areasOI, data, 'positiveSignificantUnits', windowSize, mainFolder);
    [~, rNegativeUnitsPupil20minWindows, pvalNegativeUnitsPupil20minWindows] = corrPlotPupil(areasOI, data, 'negativeUnits', windowSize, mainFolder);
    %   [~, rNegativeSignificantUnitsPupil20minWindows, pvalNegativeSignificantUnitsPupil20minWindows] = corrPlotPupil(areasOI, data, 'negativeSignificantUnits', windowSize, mainFolder);
    %   [~, rNeutralUnitsPupil20minWindows, pvalNeutralUnitsPupil20minWindows] = corrPlotPupil(areasOI, data, 'neutralUnits', windowSize, mainFolder);
    [~, rPositiveMUAsPupil20minWindows, pvalPositiveMUAsPupil20minWindows] = corrPlotPupil(areasOI, data, 'positiveMUAs', windowSize, mainFolder);
    %   [~, rPositiveSignificantMUAsPupil20minWindows, pvalPositiveSignificantMUAsPupil20minWindows] = corrPlotPupil(areasOI, data, 'positiveSignificantMUAs', windowSize, mainFolder);
    [~, rNegativeMUAsPupil20minWindows, pvalNegativeMUAsPupil20minWindows] = corrPlotPupil(areasOI, data, 'negativeMUAs', windowSize, mainFolder);
    %   [~, rNegativeSignificantMUAsPupil20minWindows, pvalNegativeSignificantMUAsPupil20minWindows] = corrPlotPupil(areasOI, data, 'negativeSignificantMUAs', windowSize, mainFolder);
    %   [~, rNeutralMUAsPupil20minWindows, pvalNeutralMUAsPupil20minWindows] = corrPlotPupil(areasOI, data, 'neutralMUAs', windowSize, mainFolder);
    
    windowSize = 30;
    [~, rPositiveUnitsPupil30minWindows, pvalPositiveUnitsPupil30minWindows] = corrPlotPupil(areasOI, data, 'positiveUnits', windowSize, mainFolder);
    %   [~, rPositiveSignificantUnitsPupil30minWindows, pvalPositiveSignificantUnitsPupil30minWindows] = corrPlotPupil(areasOI, data, 'positiveSignificantUnits', windowSize, mainFolder);
    [~, rNegativeUnitsPupil30minWindows, pvalNegativeUnitsPupil30minWindows] = corrPlotPupil(areasOI, data, 'negativeUnits', windowSize, mainFolder);
    %   [~, rNegativeSignificantUnitsPupil30minWindows, pvalNegativeSignificantUnitsPupil30minWindows] = corrPlotPupil(areasOI, data, 'negativeSignificantUnits', windowSize, mainFolder);
    %   [~, rNeutralUnitsPupil30minWindows, pvalNeutralUnitsPupil30minWindows] = corrPlotPupil(areasOI, data, 'neutralUnits', windowSize, mainFolder);
    [~, rPositiveMUAsPupil30minWindows, pvalPositiveMUAsPupil30minWindows] = corrPlotPupil(areasOI, data, 'positiveMUAs', windowSize, mainFolder);
    %   [~, rPositiveSignificantMUAsPupil30minWindows, pvalPositiveSignificantMUAsPupil30minWindows] = corrPlotPupil(areasOI, data, 'positiveSignificantMUAs', windowSize, mainFolder);
    [~, rNegativeMUAsPupil30minWindows, pvalNegativeMUAsPupil30minWindows] = corrPlotPupil(areasOI, data, 'negativeMUAs', windowSize, mainFolder);
    %   [~, rNegativeSignificantMUAsPupil30minWindows, pvalNegativeSignificantMUAsPupil30minWindows] = corrPlotPupil(areasOI, data, 'negativeSignificantMUAs', windowSize, mainFolder);
    %   [~, rNeutralMUAsPupil30minWindows, pvalNeutralMUAsPupil30minWindows] = corrPlotPupil(areasOI, data, 'neutralMUAs', windowSize, mainFolder);
  end
  
%   try
%     save([mainFolder filesep 'pupilCorrFractions.mat'], 'statsFractionsUnits', 'statsFractionsMUAs',...
%       'rPositiveUnits', 'pvalPositiveUnits', 'rPositiveSignificantUnits', 'pvalPositiveSignificantUnits', 'rNegativeUnits', 'pvalNegativeUnits','rNegativeSignificantUnits',...
%       'pvalNegativeSignificantUnits', 'rNeutralUnits', 'pvalNeutralUnits', 'rPositiveMUAs', 'pvalPositiveMUAs', 'rPositiveSignificantMUAs', 'pvalPositiveSignificantMUAs',...
%       'rNegativeMUAs', 'pvalNegativeMUAs', 'rNegativeSignificantMUAs', 'pvalNegativeSignificantMUAs', 'rNeutralMUAs', 'pvalNeutralMUAs',...
%       'rPositiveUnits10minWindows', 'pvalPositiveUnits10minWindows', 'rPositiveSignificantUnits10minWindows', 'pvalPositiveSignificantUnits10minWindows',...
%       'rNegativeUnits10minWindows', 'pvalNegativeUnits10minWindows','rNegativeSignificantUnits10minWindows', 'pvalNegativeSignificantUnits10minWindows',...
%       'rNeutralUnits10minWindows', 'pvalNeutralUnits10minWindows', 'rPositiveMUAs10minWindows', 'pvalPositiveMUAs10minWindows',...
%       'rPositiveSignificantMUAs10minWindows', 'pvalPositiveSignificantMUAs10minWindows', 'rNegativeMUAs10minWindows', 'pvalNegativeMUAs10minWindows',...
%       'rNegativeSignificantMUAs10minWindows', 'pvalNegativeSignificantMUAs10minWindows', 'rNeutralMUAs10minWindows', 'pvalNeutralMUAs10minWindows',...
%       'rPositiveUnits20minWindows', 'pvalPositiveUnits20minWindows', 'rPositiveSignificantUnits20minWindows', 'pvalPositiveSignificantUnits20minWindows',...
%       'rNegativeUnits20minWindows', 'pvalNegativeUnits20minWindows','rNegativeSignificantUnits20minWindows', 'pvalNegativeSignificantUnits20minWindows',...
%       'rNeutralUnits20minWindows', 'pvalNeutralUnits20minWindows', 'rPositiveMUAs20minWindows', 'pvalPositiveMUAs20minWindows',...
%       'rPositiveSignificantMUAs20minWindows', 'pvalPositiveSignificantMUAs20minWindows', 'rNegativeMUAs20minWindows', 'pvalNegativeMUAs20minWindows',...
%       'rNegativeSignificantMUAs20minWindows', 'pvalNegativeSignificantMUAs20minWindows', 'rNeutralMUAs20minWindows', 'pvalNeutralMUAs20minWindows',...
%       'rPositiveUnits30minWindows', 'pvalPositiveUnits30minWindows', 'rPositiveSignificantUnits30minWindows', 'pvalPositiveSignificantUnits30minWindows',...
%       'rNegativeUnits30minWindows', 'pvalNegativeUnits30minWindows','rNegativeSignificantUnits30minWindows', 'pvalNegativeSignificantUnits30minWindows',...
%       'rNeutralUnits30minWindows', 'pvalNeutralUnits30minWindows', 'rPositiveMUAs30minWindows', 'pvalPositiveMUAs30minWindows',...
%       'rPositiveSignificantMUAs30minWindows', 'pvalPositiveSignificantMUAs30minWindows', 'rNegativeMUAs30minWindows', 'pvalNegativeMUAs30minWindows',...
%       'rNegativeSignificantMUAs30minWindows', 'pvalNegativeSignificantMUAs30minWindows', 'rNeutralMUAs30minWindows', 'pvalNeutralMUAs30minWindows',...
%       'rPositiveUnitsPupil', 'pvalPositiveUnitsPupil', 'rPositiveSignificantUnitsPupil', 'pvalPositiveSignificantUnitsPupil', 'rNegativeUnitsPupil', 'pvalNegativeUnitsPupil','rNegativeSignificantUnitsPupil',...
%       'pvalNegativeSignificantUnitsPupil', 'rNeutralUnitsPupil', 'pvalNeutralUnitsPupil', 'rPositiveMUAsPupil', 'pvalPositiveMUAsPupil', 'rPositiveSignificantMUAsPupil', 'pvalPositiveSignificantMUAsPupil',...
%       'rNegativeMUAsPupil', 'pvalNegativeMUAsPupil', 'rNegativeSignificantMUAsPupil', 'pvalNegativeSignificantMUAsPupil', 'rNeutralMUAsPupil', 'pvalNeutralMUAsPupil',...
%       'rPositiveUnitsPupil10minWindows', 'pvalPositiveUnitsPupil10minWindows', 'rPositiveSignificantUnitsPupil10minWindows', 'pvalPositiveSignificantUnitsPupil10minWindows',...
%       'rNegativeUnitsPupil10minWindows', 'pvalNegativeUnitsPupil10minWindows','rNegativeSignificantUnitsPupil10minWindows', 'pvalNegativeSignificantUnitsPupil10minWindows',...
%       'rNeutralUnitsPupil10minWindows', 'pvalNeutralUnitsPupil10minWindows', 'rPositiveMUAsPupil10minWindows', 'pvalPositiveMUAsPupil10minWindows',...
%       'rPositiveSignificantMUAsPupil10minWindows', 'pvalPositiveSignificantMUAsPupil10minWindows', 'rNegativeMUAsPupil10minWindows', 'pvalNegativeMUAsPupil10minWindows',...
%       'rNegativeSignificantMUAsPupil10minWindows', 'pvalNegativeSignificantMUAsPupil10minWindows', 'rNeutralMUAsPupil10minWindows', 'pvalNeutralMUAsPupil10minWindows',...
%       'rPositiveUnitsPupil20minWindows', 'pvalPositiveUnitsPupil20minWindows', 'rPositiveSignificantUnitsPupil20minWindows', 'pvalPositiveSignificantUnitsPupil20minWindows',...
%       'rNegativeUnitsPupil20minWindows', 'pvalNegativeUnitsPupil20minWindows','rNegativeSignificantUnitsPupil20minWindows', 'pvalNegativeSignificantUnitsPupil20minWindows',...
%       'rNeutralUnitsPupil20minWindows', 'pvalNeutralUnitsPupil20minWindows', 'rPositiveMUAsPupil20minWindows', 'pvalPositiveMUAsPupil20minWindows',...
%       'rPositiveSignificantMUAsPupil20minWindows', 'pvalPositiveSignificantMUAsPupil20minWindows', 'rNegativeMUAsPupil20minWindows', 'pvalNegativeMUAsPupil20minWindows',...
%       'rNegativeSignificantMUAsPupil20minWindows', 'pvalNegativeSignificantMUAsPupil20minWindows', 'rNeutralMUAsPupil20minWindows', 'pvalNeutralMUAsPupil20minWindows',...
%       'rPositiveUnitsPupil30minWindows', 'pvalPositiveUnitsPupil30minWindows', 'rPositiveSignificantUnitsPupil30minWindows', 'pvalPositiveSignificantUnitsPupil30minWindows',...
%       'rNegativeUnitsPupil30minWindows', 'pvalNegativeUnitsPupil30minWindows','rNegativeSignificantUnitsPupil30minWindows', 'pvalNegativeSignificantUnitsPupil30minWindows',...
%       'rNeutralUnitsPupil30minWindows', 'pvalNeutralUnitsPupil30minWindows', 'rPositiveMUAsPupil30minWindows', 'pvalPositiveMUAsPupil30minWindows',...
%       'rPositiveSignificantMUAsPupil30minWindows', 'pvalPositiveSignificantMUAsPupil30minWindows', 'rNegativeMUAsPupil30minWindows', 'pvalNegativeMUAsPupil30minWindows',...
%       'rNegativeSignificantMUAsPupil30minWindows', 'pvalNegativeSignificantMUAsPupil30minWindows', 'rNeutralMUAsPupil30minWindows', 'pvalNeutralMUAsPupil30minWindows',...
%       'statsFractionsUnits','statsFractionsMUAs','statsFractionsUnitsPositive','statsFractionsUnitsPositiveZero','statsFractionsMUAsPositive','statsFractionsMUAsPositiveZero',...
%       'statsFractionsUnitsPositiveAnaest','statsFractionsUnitsPositiveZeroAnaest','statsFractionsMUAsPositiveAnaest','statsFractionsMUAsPositiveZeroAnaest',...
%       'statsFractionsUnitsPositiveSignificant','statsFractionsUnitsNegativeSignificant','statsFractionsMUAsPositiveSignificant','statsFractionsMUAsNegativeSignificant',...
%       'areasOI', '-append');
%   catch
%     save([mainFolder filesep 'all' filesep 'pupilCorrFractions.mat'], 'statsFractionsUnits', 'statsFractionsMUAs',...
%       'rPositiveUnits', 'pvalPositiveUnits', 'rPositiveSignificantUnits', 'pvalPositiveSignificantUnits', 'rNegativeUnits', 'pvalNegativeUnits','rNegativeSignificantUnits',...
%       'pvalNegativeSignificantUnits', 'rNeutralUnits', 'pvalNeutralUnits', 'rPositiveMUAs', 'pvalPositiveMUAs', 'rPositiveSignificantMUAs', 'pvalPositiveSignificantMUAs',...
%       'rNegativeMUAs', 'pvalNegativeMUAs', 'rNegativeSignificantMUAs', 'pvalNegativeSignificantMUAs', 'rNeutralMUAs', 'pvalNeutralMUAs',...
%       'rPositiveUnits10minWindows', 'pvalPositiveUnits10minWindows', 'rPositiveSignificantUnits10minWindows', 'pvalPositiveSignificantUnits10minWindows',...
%       'rNegativeUnits10minWindows', 'pvalNegativeUnits10minWindows','rNegativeSignificantUnits10minWindows', 'pvalNegativeSignificantUnits10minWindows',...
%       'rNeutralUnits10minWindows', 'pvalNeutralUnits10minWindows', 'rPositiveMUAs10minWindows', 'pvalPositiveMUAs10minWindows',...
%       'rPositiveSignificantMUAs10minWindows', 'pvalPositiveSignificantMUAs10minWindows', 'rNegativeMUAs10minWindows', 'pvalNegativeMUAs10minWindows',...
%       'rNegativeSignificantMUAs10minWindows', 'pvalNegativeSignificantMUAs10minWindows', 'rNeutralMUAs10minWindows', 'pvalNeutralMUAs10minWindows',...
%       'rPositiveUnits20minWindows', 'pvalPositiveUnits20minWindows', 'rPositiveSignificantUnits20minWindows', 'pvalPositiveSignificantUnits20minWindows',...
%       'rNegativeUnits20minWindows', 'pvalNegativeUnits20minWindows','rNegativeSignificantUnits20minWindows', 'pvalNegativeSignificantUnits20minWindows',...
%       'rNeutralUnits20minWindows', 'pvalNeutralUnits20minWindows', 'rPositiveMUAs20minWindows', 'pvalPositiveMUAs20minWindows',...
%       'rPositiveSignificantMUAs20minWindows', 'pvalPositiveSignificantMUAs20minWindows', 'rNegativeMUAs20minWindows', 'pvalNegativeMUAs20minWindows',...
%       'rNegativeSignificantMUAs20minWindows', 'pvalNegativeSignificantMUAs20minWindows', 'rNeutralMUAs20minWindows', 'pvalNeutralMUAs20minWindows',...
%       'rPositiveUnits30minWindows', 'pvalPositiveUnits30minWindows', 'rPositiveSignificantUnits30minWindows', 'pvalPositiveSignificantUnits30minWindows',...
%       'rNegativeUnits30minWindows', 'pvalNegativeUnits30minWindows','rNegativeSignificantUnits30minWindows', 'pvalNegativeSignificantUnits30minWindows',...
%       'rNeutralUnits30minWindows', 'pvalNeutralUnits30minWindows', 'rPositiveMUAs30minWindows', 'pvalPositiveMUAs30minWindows',...
%       'rPositiveSignificantMUAs30minWindows', 'pvalPositiveSignificantMUAs30minWindows', 'rNegativeMUAs30minWindows', 'pvalNegativeMUAs30minWindows',...
%       'rNegativeSignificantMUAs30minWindows', 'pvalNegativeSignificantMUAs30minWindows', 'rNeutralMUAs30minWindows', 'pvalNeutralMUAs30minWindows',...
%       'rPositiveUnitsPupil', 'pvalPositiveUnitsPupil', 'rPositiveSignificantUnitsPupil', 'pvalPositiveSignificantUnitsPupil', 'rNegativeUnitsPupil', 'pvalNegativeUnitsPupil','rNegativeSignificantUnitsPupil',...
%       'pvalNegativeSignificantUnitsPupil', 'rNeutralUnitsPupil', 'pvalNeutralUnitsPupil', 'rPositiveMUAsPupil', 'pvalPositiveMUAsPupil', 'rPositiveSignificantMUAsPupil', 'pvalPositiveSignificantMUAsPupil',...
%       'rNegativeMUAsPupil', 'pvalNegativeMUAsPupil', 'rNegativeSignificantMUAsPupil', 'pvalNegativeSignificantMUAsPupil', 'rNeutralMUAsPupil', 'pvalNeutralMUAsPupil',...
%       'rPositiveUnitsPupil10minWindows', 'pvalPositiveUnitsPupil10minWindows', 'rPositiveSignificantUnitsPupil10minWindows', 'pvalPositiveSignificantUnitsPupil10minWindows',...
%       'rNegativeUnitsPupil10minWindows', 'pvalNegativeUnitsPupil10minWindows','rNegativeSignificantUnitsPupil10minWindows', 'pvalNegativeSignificantUnitsPupil10minWindows',...
%       'rNeutralUnitsPupil10minWindows', 'pvalNeutralUnitsPupil10minWindows', 'rPositiveMUAsPupil10minWindows', 'pvalPositiveMUAsPupil10minWindows',...
%       'rPositiveSignificantMUAsPupil10minWindows', 'pvalPositiveSignificantMUAsPupil10minWindows', 'rNegativeMUAsPupil10minWindows', 'pvalNegativeMUAsPupil10minWindows',...
%       'rNegativeSignificantMUAsPupil10minWindows', 'pvalNegativeSignificantMUAsPupil10minWindows', 'rNeutralMUAsPupil10minWindows', 'pvalNeutralMUAsPupil10minWindows',...
%       'rPositiveUnitsPupil20minWindows', 'pvalPositiveUnitsPupil20minWindows', 'rPositiveSignificantUnitsPupil20minWindows', 'pvalPositiveSignificantUnitsPupil20minWindows',...
%       'rNegativeUnitsPupil20minWindows', 'pvalNegativeUnitsPupil20minWindows','rNegativeSignificantUnitsPupil20minWindows', 'pvalNegativeSignificantUnitsPupil20minWindows',...
%       'rNeutralUnitsPupil20minWindows', 'pvalNeutralUnitsPupil20minWindows', 'rPositiveMUAsPupil20minWindows', 'pvalPositiveMUAsPupil20minWindows',...
%       'rPositiveSignificantMUAsPupil20minWindows', 'pvalPositiveSignificantMUAsPupil20minWindows', 'rNegativeMUAsPupil20minWindows', 'pvalNegativeMUAsPupil20minWindows',...
%       'rNegativeSignificantMUAsPupil20minWindows', 'pvalNegativeSignificantMUAsPupil20minWindows', 'rNeutralMUAsPupil20minWindows', 'pvalNeutralMUAsPupil20minWindows',...
%       'rPositiveUnitsPupil30minWindows', 'pvalPositiveUnitsPupil30minWindows', 'rPositiveSignificantUnitsPupil30minWindows', 'pvalPositiveSignificantUnitsPupil30minWindows',...
%       'rNegativeUnitsPupil30minWindows', 'pvalNegativeUnitsPupil30minWindows','rNegativeSignificantUnitsPupil30minWindows', 'pvalNegativeSignificantUnitsPupil30minWindows',...
%       'rNeutralUnitsPupil30minWindows', 'pvalNeutralUnitsPupil30minWindows', 'rPositiveMUAsPupil30minWindows', 'pvalPositiveMUAsPupil30minWindows',...
%       'rPositiveSignificantMUAsPupil30minWindows', 'pvalPositiveSignificantMUAsPupil30minWindows', 'rNegativeMUAsPupil30minWindows', 'pvalNegativeMUAsPupil30minWindows',...
%       'rNegativeSignificantMUAsPupil30minWindows', 'pvalNegativeSignificantMUAsPupil30minWindows', 'rNeutralMUAsPupil30minWindows', 'pvalNeutralMUAsPupil30minWindows',...
%       'statsFractionsUnits','statsFractionsMUAs','statsFractionsUnitsPositive','statsFractionsUnitsPositiveZero','statsFractionsMUAsPositive','statsFractionsMUAsPositiveZero',...
%       'statsFractionsUnitsPositiveAnaest','statsFractionsUnitsPositiveZeroAnaest','statsFractionsMUAsPositiveAnaest','statsFractionsMUAsPositiveZeroAnaest',...
%       'statsFractionsUnitsPositiveSignificant','statsFractionsUnitsNegativeSignificant','statsFractionsMUAsPositiveSignificant','statsFractionsMUAsNegativeSignificant',...
%       'areasOI', '-append');
%   end
  if drawCorrs
    if strcmpi(repository,'uol')
      try
        save([mainFolder filesep 'pupilCorrFractions.mat'], 'statsFractionsUnits', 'statsFractionsMUAs',...
          'rPositiveUnits', 'pvalPositiveUnits', 'rNegativeUnits', 'pvalNegativeUnits','rPositiveMUAs', 'pvalPositiveMUAs', 'rNegativeMUAs', 'pvalNegativeMUAs',...
          'rPositiveUnits10minWindows', 'pvalPositiveUnits10minWindows', 'rNegativeUnits10minWindows', 'pvalNegativeUnits10minWindows',...
          'rPositiveMUAs10minWindows', 'pvalPositiveMUAs10minWindows','rNegativeMUAs10minWindows', 'pvalNegativeMUAs10minWindows',...
          'rPositiveUnits20minWindows', 'pvalPositiveUnits20minWindows', 'rNegativeUnits20minWindows', 'pvalNegativeUnits20minWindows',...
          'rPositiveMUAs20minWindows', 'pvalPositiveMUAs20minWindows', 'rNegativeMUAs20minWindows', 'pvalNegativeMUAs20minWindows',...
          'rPositiveUnits30minWindows', 'pvalPositiveUnits30minWindows', 'rNegativeUnits30minWindows', 'pvalNegativeUnits30minWindows',...
          'rPositiveMUAs30minWindows', 'pvalPositiveMUAs30minWindows','rNegativeMUAs30minWindows', 'pvalNegativeMUAs30minWindows',...
          'rPositiveUnitsPupil', 'pvalPositiveUnitsPupil', 'rNegativeUnitsPupil', 'pvalNegativeUnitsPupil',...
          'rPositiveMUAsPupil', 'pvalPositiveMUAsPupil', 'rNegativeMUAsPupil', 'pvalNegativeMUAsPupil',...
          'rPositiveUnitsPupil10minWindows', 'pvalPositiveUnitsPupil10minWindows','rNegativeUnitsPupil10minWindows', 'pvalNegativeUnitsPupil10minWindows',...
          'rPositiveMUAsPupil10minWindows', 'pvalPositiveMUAsPupil10minWindows', 'rNegativeMUAsPupil10minWindows', 'pvalNegativeMUAsPupil10minWindows',...
          'rPositiveUnitsPupil20minWindows', 'pvalPositiveUnitsPupil20minWindows', 'rNegativeUnitsPupil20minWindows', 'pvalNegativeUnitsPupil20minWindows',...
          'rPositiveMUAsPupil20minWindows', 'pvalPositiveMUAsPupil20minWindows', 'rNegativeMUAsPupil20minWindows', 'pvalNegativeMUAsPupil20minWindows',...
          'rPositiveUnitsPupil30minWindows', 'pvalPositiveUnitsPupil30minWindows', 'rNegativeUnitsPupil30minWindows', 'pvalNegativeUnitsPupil30minWindows',...
          'rPositiveMUAsPupil30minWindows', 'pvalPositiveMUAsPupil30minWindows', 'rNegativeMUAsPupil30minWindows', 'pvalNegativeMUAsPupil30minWindows',...
          'statsFractionsUnits','statsFractionsMUAs','statsFractionsUnitsPositive','statsFractionsUnitsPositiveZero','statsFractionsMUAsPositive','statsFractionsMUAsPositiveZero',...
          'statsFractionsUnitsPositiveAnaest','statsFractionsUnitsPositiveZeroAnaest','statsFractionsMUAsPositiveAnaest','statsFractionsMUAsPositiveZeroAnaest',...
          'statsFractionsUnitsMerged','statsFractionsMUAsMerged','statsFractionsUnitsPositiveMerged','statsFractionsUnitsPositiveZeroMerged','statsFractionsMUAsPositiveMerged','statsFractionsMUAsPositiveZeroMerged',...
          'statsFractionsUnitsPositiveAnaestMerged','statsFractionsUnitsPositiveZeroAnaestMerged','statsFractionsMUAsPositiveAnaestMerged','statsFractionsMUAsPositiveZeroAnaestMerged',...
          'areasOI', '-append');
      catch
        save([mainFolder filesep 'all' filesep 'pupilCorrFractions.mat'], 'statsFractionsUnits', 'statsFractionsMUAs',...
          'rPositiveUnits', 'pvalPositiveUnits', 'rNegativeUnits', 'pvalNegativeUnits','rPositiveMUAs', 'pvalPositiveMUAs', 'rNegativeMUAs', 'pvalNegativeMUAs',...
          'rPositiveUnits10minWindows', 'pvalPositiveUnits10minWindows', 'rNegativeUnits10minWindows', 'pvalNegativeUnits10minWindows',...
          'rPositiveMUAs10minWindows', 'pvalPositiveMUAs10minWindows','rNegativeMUAs10minWindows', 'pvalNegativeMUAs10minWindows',...
          'rPositiveUnits20minWindows', 'pvalPositiveUnits20minWindows', 'rNegativeUnits20minWindows', 'pvalNegativeUnits20minWindows',...
          'rPositiveMUAs20minWindows', 'pvalPositiveMUAs20minWindows', 'rNegativeMUAs20minWindows', 'pvalNegativeMUAs20minWindows',...
          'rPositiveUnits30minWindows', 'pvalPositiveUnits30minWindows', 'rNegativeUnits30minWindows', 'pvalNegativeUnits30minWindows',...
          'rPositiveMUAs30minWindows', 'pvalPositiveMUAs30minWindows','rNegativeMUAs30minWindows', 'pvalNegativeMUAs30minWindows',...
          'rPositiveUnitsPupil', 'pvalPositiveUnitsPupil', 'rNegativeUnitsPupil', 'pvalNegativeUnitsPupil',...
          'rPositiveMUAsPupil', 'pvalPositiveMUAsPupil', 'rNegativeMUAsPupil', 'pvalNegativeMUAsPupil',...
          'rPositiveUnitsPupil10minWindows', 'pvalPositiveUnitsPupil10minWindows','rNegativeUnitsPupil10minWindows', 'pvalNegativeUnitsPupil10minWindows',...
          'rPositiveMUAsPupil10minWindows', 'pvalPositiveMUAsPupil10minWindows', 'rNegativeMUAsPupil10minWindows', 'pvalNegativeMUAsPupil10minWindows',...
          'rPositiveUnitsPupil20minWindows', 'pvalPositiveUnitsPupil20minWindows', 'rNegativeUnitsPupil20minWindows', 'pvalNegativeUnitsPupil20minWindows',...
          'rPositiveMUAsPupil20minWindows', 'pvalPositiveMUAsPupil20minWindows', 'rNegativeMUAsPupil20minWindows', 'pvalNegativeMUAsPupil20minWindows',...
          'rPositiveUnitsPupil30minWindows', 'pvalPositiveUnitsPupil30minWindows', 'rNegativeUnitsPupil30minWindows', 'pvalNegativeUnitsPupil30minWindows',...
          'rPositiveMUAsPupil30minWindows', 'pvalPositiveMUAsPupil30minWindows', 'rNegativeMUAsPupil30minWindows', 'pvalNegativeMUAsPupil30minWindows',...
          'statsFractionsUnits','statsFractionsMUAs','statsFractionsUnitsPositive','statsFractionsUnitsPositiveZero','statsFractionsMUAsPositive','statsFractionsMUAsPositiveZero',...
          'statsFractionsUnitsPositiveAnaest','statsFractionsUnitsPositiveZeroAnaest','statsFractionsMUAsPositiveAnaest','statsFractionsMUAsPositiveZeroAnaest',...
          'statsFractionsUnitsMerged','statsFractionsMUAsMerged','statsFractionsUnitsPositiveMerged','statsFractionsUnitsPositiveZeroMerged','statsFractionsMUAsPositiveMerged','statsFractionsMUAsPositiveZeroMerged',...
          'statsFractionsUnitsPositiveAnaestMerged','statsFractionsUnitsPositiveZeroAnaestMerged','statsFractionsMUAsPositiveAnaestMerged','statsFractionsMUAsPositiveZeroAnaestMerged',...
          'areasOI', '-append');
      end
    else
      try
        save([mainFolder filesep 'pupilCorrFractions.mat'], 'statsFractionsUnits', 'statsFractionsMUAs',...
          'rPositiveUnits', 'pvalPositiveUnits', 'rNegativeUnits', 'pvalNegativeUnits','rPositiveMUAs', 'pvalPositiveMUAs', 'rNegativeMUAs', 'pvalNegativeMUAs',...
          'rPositiveUnits10minWindows', 'pvalPositiveUnits10minWindows', 'rNegativeUnits10minWindows', 'pvalNegativeUnits10minWindows',...
          'rPositiveMUAs10minWindows', 'pvalPositiveMUAs10minWindows','rNegativeMUAs10minWindows', 'pvalNegativeMUAs10minWindows',...
          'rPositiveUnits20minWindows', 'pvalPositiveUnits20minWindows', 'rNegativeUnits20minWindows', 'pvalNegativeUnits20minWindows',...
          'rPositiveMUAs20minWindows', 'pvalPositiveMUAs20minWindows', 'rNegativeMUAs20minWindows', 'pvalNegativeMUAs20minWindows',...
          'rPositiveUnits30minWindows', 'pvalPositiveUnits30minWindows', 'rNegativeUnits30minWindows', 'pvalNegativeUnits30minWindows',...
          'rPositiveMUAs30minWindows', 'pvalPositiveMUAs30minWindows','rNegativeMUAs30minWindows', 'pvalNegativeMUAs30minWindows',...
          'rPositiveUnitsPupil', 'pvalPositiveUnitsPupil', 'rNegativeUnitsPupil', 'pvalNegativeUnitsPupil',...
          'rPositiveMUAsPupil', 'pvalPositiveMUAsPupil', 'rNegativeMUAsPupil', 'pvalNegativeMUAsPupil',...
          'rPositiveUnitsPupil10minWindows', 'pvalPositiveUnitsPupil10minWindows','rNegativeUnitsPupil10minWindows', 'pvalNegativeUnitsPupil10minWindows',...
          'rPositiveMUAsPupil10minWindows', 'pvalPositiveMUAsPupil10minWindows', 'rNegativeMUAsPupil10minWindows', 'pvalNegativeMUAsPupil10minWindows',...
          'rPositiveUnitsPupil20minWindows', 'pvalPositiveUnitsPupil20minWindows', 'rNegativeUnitsPupil20minWindows', 'pvalNegativeUnitsPupil20minWindows',...
          'rPositiveMUAsPupil20minWindows', 'pvalPositiveMUAsPupil20minWindows', 'rNegativeMUAsPupil20minWindows', 'pvalNegativeMUAsPupil20minWindows',...
          'rPositiveUnitsPupil30minWindows', 'pvalPositiveUnitsPupil30minWindows', 'rNegativeUnitsPupil30minWindows', 'pvalNegativeUnitsPupil30minWindows',...
          'rPositiveMUAsPupil30minWindows', 'pvalPositiveMUAsPupil30minWindows', 'rNegativeMUAsPupil30minWindows', 'pvalNegativeMUAsPupil30minWindows',...
          'statsFractionsUnits','statsFractionsMUAs','statsFractionsUnitsPositive','statsFractionsUnitsPositiveZero','statsFractionsMUAsPositive','statsFractionsMUAsPositiveZero',...
          'statsFractionsUnitsMerged','statsFractionsMUAsMerged','statsFractionsUnitsPositiveMerged','statsFractionsUnitsPositiveZeroMerged','statsFractionsMUAsPositiveMerged','statsFractionsMUAsPositiveZeroMerged',...
          'areasOI', '-append');
      catch
        save([mainFolder filesep 'all' filesep 'pupilCorrFractions.mat'], 'statsFractionsUnits', 'statsFractionsMUAs',...
          'rPositiveUnits', 'pvalPositiveUnits', 'rNegativeUnits', 'pvalNegativeUnits','rPositiveMUAs', 'pvalPositiveMUAs', 'rNegativeMUAs', 'pvalNegativeMUAs',...
          'rPositiveUnits10minWindows', 'pvalPositiveUnits10minWindows', 'rNegativeUnits10minWindows', 'pvalNegativeUnits10minWindows',...
          'rPositiveMUAs10minWindows', 'pvalPositiveMUAs10minWindows','rNegativeMUAs10minWindows', 'pvalNegativeMUAs10minWindows',...
          'rPositiveUnits20minWindows', 'pvalPositiveUnits20minWindows', 'rNegativeUnits20minWindows', 'pvalNegativeUnits20minWindows',...
          'rPositiveMUAs20minWindows', 'pvalPositiveMUAs20minWindows', 'rNegativeMUAs20minWindows', 'pvalNegativeMUAs20minWindows',...
          'rPositiveUnits30minWindows', 'pvalPositiveUnits30minWindows', 'rNegativeUnits30minWindows', 'pvalNegativeUnits30minWindows',...
          'rPositiveMUAs30minWindows', 'pvalPositiveMUAs30minWindows','rNegativeMUAs30minWindows', 'pvalNegativeMUAs30minWindows',...
          'rPositiveUnitsPupil', 'pvalPositiveUnitsPupil', 'rNegativeUnitsPupil', 'pvalNegativeUnitsPupil',...
          'rPositiveMUAsPupil', 'pvalPositiveMUAsPupil', 'rNegativeMUAsPupil', 'pvalNegativeMUAsPupil',...
          'rPositiveUnitsPupil10minWindows', 'pvalPositiveUnitsPupil10minWindows','rNegativeUnitsPupil10minWindows', 'pvalNegativeUnitsPupil10minWindows',...
          'rPositiveMUAsPupil10minWindows', 'pvalPositiveMUAsPupil10minWindows', 'rNegativeMUAsPupil10minWindows', 'pvalNegativeMUAsPupil10minWindows',...
          'rPositiveUnitsPupil20minWindows', 'pvalPositiveUnitsPupil20minWindows', 'rNegativeUnitsPupil20minWindows', 'pvalNegativeUnitsPupil20minWindows',...
          'rPositiveMUAsPupil20minWindows', 'pvalPositiveMUAsPupil20minWindows', 'rNegativeMUAsPupil20minWindows', 'pvalNegativeMUAsPupil20minWindows',...
          'rPositiveUnitsPupil30minWindows', 'pvalPositiveUnitsPupil30minWindows', 'rNegativeUnitsPupil30minWindows', 'pvalNegativeUnitsPupil30minWindows',...
          'rPositiveMUAsPupil30minWindows', 'pvalPositiveMUAsPupil30minWindows', 'rNegativeMUAsPupil30minWindows', 'pvalNegativeMUAsPupil30minWindows',...
          'statsFractionsUnits','statsFractionsMUAs','statsFractionsUnitsPositive','statsFractionsUnitsPositiveZero','statsFractionsMUAsPositive','statsFractionsMUAsPositiveZero',...
          'statsFractionsUnitsMerged','statsFractionsMUAsMerged','statsFractionsUnitsPositiveMerged','statsFractionsUnitsPositiveZeroMerged','statsFractionsMUAsPositiveMerged','statsFractionsMUAsPositiveZeroMerged',...
          'areasOI', '-append');
      end
    end
  else
    if strcmpi(repository,'uol')
      try
        save([mainFolder filesep 'pupilCorrFractions.mat'], 'statsFractionsUnits', 'statsFractionsMUAs',...
          'statsFractionsUnits','statsFractionsMUAs','statsFractionsUnitsPositive','statsFractionsUnitsPositiveZero','statsFractionsMUAsPositive','statsFractionsMUAsPositiveZero',...
          'statsFractionsUnitsPositiveAnaest','statsFractionsUnitsPositiveZeroAnaest','statsFractionsMUAsPositiveAnaest','statsFractionsMUAsPositiveZeroAnaest',...
          'statsFractionsUnitsMerged','statsFractionsMUAsMerged','statsFractionsUnitsPositiveMerged','statsFractionsUnitsPositiveZeroMerged','statsFractionsMUAsPositiveMerged','statsFractionsMUAsPositiveZeroMerged',...
          'statsFractionsUnitsPositiveAnaestMerged','statsFractionsUnitsPositiveZeroAnaestMerged','statsFractionsMUAsPositiveAnaestMerged','statsFractionsMUAsPositiveZeroAnaestMerged',...
          'areasOI', '-append');
      catch
        save([mainFolder filesep 'all' filesep 'pupilCorrFractions.mat'], 'statsFractionsUnits', 'statsFractionsMUAs',...
          'statsFractionsUnits','statsFractionsMUAs','statsFractionsUnitsPositive','statsFractionsUnitsPositiveZero','statsFractionsMUAsPositive','statsFractionsMUAsPositiveZero',...
          'statsFractionsUnitsPositiveAnaest','statsFractionsUnitsPositiveZeroAnaest','statsFractionsMUAsPositiveAnaest','statsFractionsMUAsPositiveZeroAnaest',...
          'statsFractionsUnitsMerged','statsFractionsMUAsMerged','statsFractionsUnitsPositiveMerged','statsFractionsUnitsPositiveZeroMerged','statsFractionsMUAsPositiveMerged','statsFractionsMUAsPositiveZeroMerged',...
          'statsFractionsUnitsPositiveAnaestMerged','statsFractionsUnitsPositiveZeroAnaestMerged','statsFractionsMUAsPositiveAnaestMerged','statsFractionsMUAsPositiveZeroAnaestMerged',...
          'areasOI', '-append');
      end
    else
      try
        save([mainFolder filesep 'pupilCorrFractions.mat'], 'statsFractionsUnits', 'statsFractionsMUAs',...
          'statsFractionsUnits','statsFractionsMUAs','statsFractionsUnitsPositive','statsFractionsUnitsPositiveZero','statsFractionsMUAsPositive','statsFractionsMUAsPositiveZero',...
          'statsFractionsUnitsMerged','statsFractionsMUAsMerged','statsFractionsUnitsPositiveMerged','statsFractionsUnitsPositiveZeroMerged','statsFractionsMUAsPositiveMerged','statsFractionsMUAsPositiveZeroMerged',...
          'areasOI', '-append');
      catch
        save([mainFolder filesep 'all' filesep 'pupilCorrFractions.mat'], 'statsFractionsUnits', 'statsFractionsMUAs',...
          'statsFractionsUnits','statsFractionsMUAs','statsFractionsUnitsPositive','statsFractionsUnitsPositiveZero','statsFractionsMUAsPositive','statsFractionsMUAsPositiveZero',...
          'statsFractionsUnitsMerged','statsFractionsMUAsMerged','statsFractionsUnitsPositiveMerged','statsFractionsUnitsPositiveZeroMerged','statsFractionsMUAsPositiveMerged','statsFractionsMUAsPositiveZeroMerged',...
          'areasOI', '-append');
      end
    end
  end
else
  if strcmp(repository, 'allensdk')
    load([mainFolder filesep 'pupilCorrFractions.mat']);
  elseif strcmp(repository, 'uol')
    load([mainFolder filesep 'all' filesep 'pupilCorrFractions.mat']);
  end
end


if fullRun <= 3 && drawTables
  
  %% Produce tables
  
  % MUAs 5sec: S1 v RSC
  if strcmp(repository, 'uol')
    inds = findComboEntries('lS1', 'lRSC', areasOI);
    positiveMUAsTable = [rPositiveMUAs10minWindows(inds) rPositiveMUAs20minWindows(inds) rPositiveMUAs30minWindows(inds) rPositiveMUAs(inds)];
    positiveMUAsTable = [positiveMUAsTable; positiveMUAsTable.^2];
    positiveMUAsTable = [positiveMUAsTable; [pvalPositiveMUAs10minWindows(inds) pvalPositiveMUAs20minWindows(inds) pvalPositiveMUAs30minWindows(inds) pvalPositiveMUAs(inds)]];
    negativeMUAsTable = [rPositiveMUAs10minWindows(inds) rPositiveMUAs20minWindows(inds) rPositiveMUAs30minWindows(inds) rPositiveMUAs(inds)];
    negativeMUAsTable = [negativeMUAsTable; negativeMUAsTable.^2];
    negativeMUAsTable = [negativeMUAsTable; [pvalPositiveMUAs10minWindows(inds) pvalPositiveMUAs20minWindows(inds) pvalPositiveMUAs30minWindows(inds) pvalPositiveMUAs(inds)]];
    MUAsTable_S1vRSC = [positiveMUAsTable; negativeMUAsTable];
  end
  
  % MUAs 5sec: S1 v CA
  if strcmp(repository, 'uol')
    inds = findComboEntries('lS1', 'CA', areasOI);
    positiveMUAsTable = [rPositiveMUAs10minWindows(inds) rPositiveMUAs20minWindows(inds) rPositiveMUAs30minWindows(inds) rPositiveMUAs(inds)];
    positiveMUAsTable = [positiveMUAsTable; positiveMUAsTable.^2];
    positiveMUAsTable = [positiveMUAsTable; [pvalPositiveMUAs10minWindows(inds) pvalPositiveMUAs20minWindows(inds) pvalPositiveMUAs30minWindows(inds) pvalPositiveMUAs(inds)]];
    negativeMUAsTable = [rPositiveMUAs10minWindows(inds) rPositiveMUAs20minWindows(inds) rPositiveMUAs30minWindows(inds) rPositiveMUAs(inds)];
    negativeMUAsTable = [negativeMUAsTable; negativeMUAsTable.^2];
    negativeMUAsTable = [negativeMUAsTable; [pvalPositiveMUAs10minWindows(inds) pvalPositiveMUAs20minWindows(inds) pvalPositiveMUAs30minWindows(inds) pvalPositiveMUAs(inds)]];
    MUAsTable_S1vCA = [positiveMUAsTable; negativeMUAsTable];
  end
  
  % MUAs 5sec: RSC v CA
  if strcmp(repository, 'uol')
    inds = findComboEntries('lRSC', 'CA', areasOI);
    positiveMUAsTable = [rPositiveMUAs10minWindows(inds) rPositiveMUAs20minWindows(inds) rPositiveMUAs30minWindows(inds) rPositiveMUAs(inds)];
    positiveMUAsTable = [positiveMUAsTable; positiveMUAsTable.^2];
    positiveMUAsTable = [positiveMUAsTable; [pvalPositiveMUAs10minWindows(inds) pvalPositiveMUAs20minWindows(inds) pvalPositiveMUAs30minWindows(inds) pvalPositiveMUAs(inds)]];
    negativeMUAsTable = [rPositiveMUAs10minWindows(inds) rPositiveMUAs20minWindows(inds) rPositiveMUAs30minWindows(inds) rPositiveMUAs(inds)];
    negativeMUAsTable = [negativeMUAsTable; negativeMUAsTable.^2];
    negativeMUAsTable = [negativeMUAsTable; [pvalPositiveMUAs10minWindows(inds) pvalPositiveMUAs20minWindows(inds) pvalPositiveMUAs30minWindows(inds) pvalPositiveMUAs(inds)]];
    MUAsTable_RSCvCA = [positiveMUAsTable; negativeMUAsTable];
  end
  
  % MUAs 5sec: V1 v CA
  if strcmp(repository, 'allensdk')
    inds = findComboEntries('V1', 'CA', areasOI);
    positiveMUAsTable = [rPositiveMUAs10minWindows(inds) rPositiveMUAs20minWindows(inds) rPositiveMUAs30minWindows(inds) rPositiveMUAs(inds)];
    positiveMUAsTable = [positiveMUAsTable; positiveMUAsTable.^2];
    positiveMUAsTable = [positiveMUAsTable; [pvalPositiveMUAs10minWindows(inds) pvalPositiveMUAs20minWindows(inds) pvalPositiveMUAs30minWindows(inds) pvalPositiveMUAs(inds)]];
    negativeMUAsTable = [rPositiveMUAs10minWindows(inds) rPositiveMUAs20minWindows(inds) rPositiveMUAs30minWindows(inds) rPositiveMUAs(inds)];
    negativeMUAsTable = [negativeMUAsTable; negativeMUAsTable.^2];
    negativeMUAsTable = [negativeMUAsTable; [pvalPositiveMUAs10minWindows(inds) pvalPositiveMUAs20minWindows(inds) pvalPositiveMUAs30minWindows(inds) pvalPositiveMUAs(inds)]];
    MUAsTable_V1vCA = [positiveMUAsTable; negativeMUAsTable];
  end
  
  % Units 5sec: S1 v RSC
  if strcmp(repository, 'uol')
    inds = findComboEntries('lS1', 'lRSC', areasOI);
    positiveUnitsTable = [rPositiveUnits10minWindows(inds) rPositiveUnits20minWindows(inds) rPositiveUnits30minWindows(inds) rPositiveUnits(inds)];
    positiveUnitsTable = [positiveUnitsTable; positiveUnitsTable.^2];
    positiveUnitsTable = [positiveUnitsTable; [pvalPositiveUnits10minWindows(inds) pvalPositiveUnits20minWindows(inds) pvalPositiveUnits30minWindows(inds) pvalPositiveUnits(inds)]];
    negativeUnitsTable = [rPositiveUnits10minWindows(inds) rPositiveUnits20minWindows(inds) rPositiveUnits30minWindows(inds) rPositiveUnits(inds)];
    negativeUnitsTable = [negativeUnitsTable; negativeUnitsTable.^2];
    negativeUnitsTable = [negativeUnitsTable; [pvalPositiveUnits10minWindows(inds) pvalPositiveUnits20minWindows(inds) pvalPositiveUnits30minWindows(inds) pvalPositiveUnits(inds)]];
    unitsTable_S1vRSC = [positiveUnitsTable; negativeUnitsTable];
  end
  
  % Units 5sec: S1 v CA
  if strcmp(repository, 'uol')
    inds = findComboEntries('lS1', 'CA', areasOI);
    positiveUnitsTable = [rPositiveUnits10minWindows(inds) rPositiveUnits20minWindows(inds) rPositiveUnits30minWindows(inds) rPositiveUnits(inds)];
    positiveUnitsTable = [positiveUnitsTable; positiveUnitsTable.^2];
    positiveUnitsTable = [positiveUnitsTable; [pvalPositiveUnits10minWindows(inds) pvalPositiveUnits20minWindows(inds) pvalPositiveUnits30minWindows(inds) pvalPositiveUnits(inds)]];
    negativeUnitsTable = [rPositiveUnits10minWindows(inds) rPositiveUnits20minWindows(inds) rPositiveUnits30minWindows(inds) rPositiveUnits(inds)];
    negativeUnitsTable = [negativeUnitsTable; negativeUnitsTable.^2];
    negativeUnitsTable = [negativeUnitsTable; [pvalPositiveUnits10minWindows(inds) pvalPositiveUnits20minWindows(inds) pvalPositiveUnits30minWindows(inds) pvalPositiveUnits(inds)]];
    unitsTable_S1vCA = [positiveUnitsTable; negativeUnitsTable];
  end
  
  % Units 5sec: RSC v CA
  if strcmp(repository, 'uol')
    inds = findComboEntries('lRSC', 'CA', areasOI);
    positiveUnitsTable = [rPositiveUnits10minWindows(inds) rPositiveUnits20minWindows(inds) rPositiveUnits30minWindows(inds) rPositiveUnits(inds)];
    positiveUnitsTable = [positiveUnitsTable; positiveUnitsTable.^2];
    positiveUnitsTable = [positiveUnitsTable; [pvalPositiveUnits10minWindows(inds) pvalPositiveUnits20minWindows(inds) pvalPositiveUnits30minWindows(inds) pvalPositiveUnits(inds)]];
    negativeUnitsTable = [rPositiveUnits10minWindows(inds) rPositiveUnits20minWindows(inds) rPositiveUnits30minWindows(inds) rPositiveUnits(inds)];
    negativeUnitsTable = [negativeUnitsTable; negativeUnitsTable.^2];
    negativeUnitsTable = [negativeUnitsTable; [pvalPositiveUnits10minWindows(inds) pvalPositiveUnits20minWindows(inds) pvalPositiveUnits30minWindows(inds) pvalPositiveUnits(inds)]];
    unitsTable_RSCvCA = [positiveUnitsTable; negativeUnitsTable];
  end
  
  % Units 5sec: V1 v CA
  if strcmp(repository, 'allensdk')
    inds = findComboEntries('V1', 'CA', areasOI);
    positiveUnitsTable = [rPositiveUnits10minWindows(inds) rPositiveUnits20minWindows(inds) rPositiveUnits30minWindows(inds) rPositiveUnits(inds)];
    positiveUnitsTable = [positiveUnitsTable; positiveUnitsTable.^2];
    positiveUnitsTable = [positiveUnitsTable; [pvalPositiveUnits10minWindows(inds) pvalPositiveUnits20minWindows(inds) pvalPositiveUnits30minWindows(inds) pvalPositiveUnits(inds)]];
    negativeUnitsTable = [rPositiveUnits10minWindows(inds) rPositiveUnits20minWindows(inds) rPositiveUnits30minWindows(inds) rPositiveUnits(inds)];
    negativeUnitsTable = [negativeUnitsTable; negativeUnitsTable.^2];
    negativeUnitsTable = [negativeUnitsTable; [pvalPositiveUnits10minWindows(inds) pvalPositiveUnits20minWindows(inds) pvalPositiveUnits30minWindows(inds) pvalPositiveUnits(inds)]];
    unitsTable_V1vCA = [positiveUnitsTable; negativeUnitsTable];
  end
  
  % MUAs vs pupil: S1
  if strcmp(repository, 'uol')
    inds = find(ismember(areasOI, 'lS1'));
    positiveMUAsPupilTable = [rPositiveMUAsPupil10minWindows(inds) rPositiveMUAsPupil20minWindows(inds) rPositiveMUAsPupil30minWindows(inds) rPositiveMUAsPupil(inds)];
    positiveMUAsPupilTable = [positiveMUAsPupilTable; positiveMUAsPupilTable.^2];
    positiveMUAsPupilTable = [positiveMUAsPupilTable; [pvalPositiveMUAsPupil10minWindows(inds) pvalPositiveMUAsPupil20minWindows(inds) pvalPositiveMUAsPupil30minWindows(inds) pvalPositiveMUAsPupil(inds)]];
    negativeMUAsPupilTable = [rPositiveMUAsPupil10minWindows(inds) rPositiveMUAsPupil20minWindows(inds) rPositiveMUAsPupil30minWindows(inds) rPositiveMUAsPupil(inds)];
    negativeMUAsPupilTable = [negativeMUAsPupilTable; negativeMUAsPupilTable.^2];
    negativeMUAsPupilTable = [negativeMUAsPupilTable; [pvalPositiveMUAsPupil10minWindows(inds) pvalPositiveMUAsPupil20minWindows(inds) pvalPositiveMUAsPupil30minWindows(inds) pvalPositiveMUAsPupil(inds)]];
    MUAsPupilTable_S1 = [positiveMUAsPupilTable; negativeMUAsPupilTable];
  end
  
  % MUAs vs pupil: RSC
  if strcmp(repository, 'uol')
    inds = find(ismember(areasOI, 'lRSC'));
    positiveMUAsPupilTable = [rPositiveMUAsPupil10minWindows(inds) rPositiveMUAsPupil20minWindows(inds) rPositiveMUAsPupil30minWindows(inds) rPositiveMUAsPupil(inds)];
    positiveMUAsPupilTable = [positiveMUAsPupilTable; positiveMUAsPupilTable.^2];
    positiveMUAsPupilTable = [positiveMUAsPupilTable; [pvalPositiveMUAsPupil10minWindows(inds) pvalPositiveMUAsPupil20minWindows(inds) pvalPositiveMUAsPupil30minWindows(inds) pvalPositiveMUAsPupil(inds)]];
    negativeMUAsPupilTable = [rPositiveMUAsPupil10minWindows(inds) rPositiveMUAsPupil20minWindows(inds) rPositiveMUAsPupil30minWindows(inds) rPositiveMUAsPupil(inds)];
    negativeMUAsPupilTable = [negativeMUAsPupilTable; negativeMUAsPupilTable.^2];
    negativeMUAsPupilTable = [negativeMUAsPupilTable; [pvalPositiveMUAsPupil10minWindows(inds) pvalPositiveMUAsPupil20minWindows(inds) pvalPositiveMUAsPupil30minWindows(inds) pvalPositiveMUAsPupil(inds)]];
    MUAsPupilTable_RSC = [positiveMUAsPupilTable; negativeMUAsPupilTable];
  end
  
  % MUAs vs pupil: CA
  inds = find(ismember(areasOI, 'CA'));
  positiveMUAsPupilTable = [rPositiveMUAsPupil10minWindows(inds) rPositiveMUAsPupil20minWindows(inds) rPositiveMUAsPupil30minWindows(inds) rPositiveMUAsPupil(inds)];
  positiveMUAsPupilTable = [positiveMUAsPupilTable; positiveMUAsPupilTable.^2];
  positiveMUAsPupilTable = [positiveMUAsPupilTable; [pvalPositiveMUAsPupil10minWindows(inds) pvalPositiveMUAsPupil20minWindows(inds) pvalPositiveMUAsPupil30minWindows(inds) pvalPositiveMUAsPupil(inds)]];
  negativeMUAsPupilTable = [rPositiveMUAsPupil10minWindows(inds) rPositiveMUAsPupil20minWindows(inds) rPositiveMUAsPupil30minWindows(inds) rPositiveMUAsPupil(inds)];
  negativeMUAsPupilTable = [negativeMUAsPupilTable; negativeMUAsPupilTable.^2];
  negativeMUAsPupilTable = [negativeMUAsPupilTable; [pvalPositiveMUAsPupil10minWindows(inds) pvalPositiveMUAsPupil20minWindows(inds) pvalPositiveMUAsPupil30minWindows(inds) pvalPositiveMUAsPupil(inds)]];
  MUAsPupilTable_CA = [positiveMUAsPupilTable; negativeMUAsPupilTable];
  
  % MUAs vs pupil: V1
  if strcmp(repository, 'allensdk')
    inds = find(ismember(areasOI, 'V1'));
    positiveMUAsPupilTable = [rPositiveMUAsPupil10minWindows(inds) rPositiveMUAsPupil20minWindows(inds) rPositiveMUAsPupil30minWindows(inds) rPositiveMUAsPupil(inds)];
    positiveMUAsPupilTable = [positiveMUAsPupilTable; positiveMUAsPupilTable.^2];
    positiveMUAsPupilTable = [positiveMUAsPupilTable; [pvalPositiveMUAsPupil10minWindows(inds) pvalPositiveMUAsPupil20minWindows(inds) pvalPositiveMUAsPupil30minWindows(inds) pvalPositiveMUAsPupil(inds)]];
    negativeMUAsPupilTable = [rPositiveMUAsPupil10minWindows(inds) rPositiveMUAsPupil20minWindows(inds) rPositiveMUAsPupil30minWindows(inds) rPositiveMUAsPupil(inds)];
    negativeMUAsPupilTable = [negativeMUAsPupilTable; negativeMUAsPupilTable.^2];
    negativeMUAsPupilTable = [negativeMUAsPupilTable; [pvalPositiveMUAsPupil10minWindows(inds) pvalPositiveMUAsPupil20minWindows(inds) pvalPositiveMUAsPupil30minWindows(inds) pvalPositiveMUAsPupil(inds)]];
    MUAsPupilTable_V1 = [positiveMUAsPupilTable; negativeMUAsPupilTable];
  end
  
  % Units vs pupil: S1
  if strcmp(repository, 'uol')
    inds = find(ismember(areasOI, 'lS1'));
    positiveUnitsPupilTable = [rPositiveUnitsPupil10minWindows(inds) rPositiveUnitsPupil20minWindows(inds) rPositiveUnitsPupil30minWindows(inds) rPositiveUnitsPupil(inds)];
    positiveUnitsPupilTable = [positiveUnitsPupilTable; positiveUnitsPupilTable.^2];
    positiveUnitsPupilTable = [positiveUnitsPupilTable; [pvalPositiveUnitsPupil10minWindows(inds) pvalPositiveUnitsPupil20minWindows(inds) pvalPositiveUnitsPupil30minWindows(inds) pvalPositiveUnitsPupil(inds)]];
    negativeUnitsPupilTable = [rPositiveUnitsPupil10minWindows(inds) rPositiveUnitsPupil20minWindows(inds) rPositiveUnitsPupil30minWindows(inds) rPositiveUnitsPupil(inds)];
    negativeUnitsPupilTable = [negativeUnitsPupilTable; negativeUnitsPupilTable.^2];
    negativeUnitsPupilTable = [negativeUnitsPupilTable; [pvalPositiveUnitsPupil10minWindows(inds) pvalPositiveUnitsPupil20minWindows(inds) pvalPositiveUnitsPupil30minWindows(inds) pvalPositiveUnitsPupil(inds)]];
    unitsPupilTable_S1 = [positiveUnitsPupilTable; negativeUnitsPupilTable];
  end
  
  % Units vs pupil: RSC
  if strcmp(repository, 'uol')
    inds = find(ismember(areasOI, 'lRSC'));
    positiveUnitsPupilTable = [rPositiveUnitsPupil10minWindows(inds) rPositiveUnitsPupil20minWindows(inds) rPositiveUnitsPupil30minWindows(inds) rPositiveUnitsPupil(inds)];
    positiveUnitsPupilTable = [positiveUnitsPupilTable; positiveUnitsPupilTable.^2];
    positiveUnitsPupilTable = [positiveUnitsPupilTable; [pvalPositiveUnitsPupil10minWindows(inds) pvalPositiveUnitsPupil20minWindows(inds) pvalPositiveUnitsPupil30minWindows(inds) pvalPositiveUnitsPupil(inds)]];
    negativeUnitsPupilTable = [rPositiveUnitsPupil10minWindows(inds) rPositiveUnitsPupil20minWindows(inds) rPositiveUnitsPupil30minWindows(inds) rPositiveUnitsPupil(inds)];
    negativeUnitsPupilTable = [negativeUnitsPupilTable; negativeUnitsPupilTable.^2];
    negativeUnitsPupilTable = [negativeUnitsPupilTable; [pvalPositiveUnitsPupil10minWindows(inds) pvalPositiveUnitsPupil20minWindows(inds) pvalPositiveUnitsPupil30minWindows(inds) pvalPositiveUnitsPupil(inds)]];
    unitsPupilTable_RSC = [positiveUnitsPupilTable; negativeUnitsPupilTable];
  end
  
  % Units vs pupil: CA
  inds = find(ismember(areasOI, 'CA'));
  positiveUnitsPupilTable = [rPositiveUnitsPupil10minWindows(inds) rPositiveUnitsPupil20minWindows(inds) rPositiveUnitsPupil30minWindows(inds) rPositiveUnitsPupil(inds)];
  positiveUnitsPupilTable = [positiveUnitsPupilTable; positiveUnitsPupilTable.^2];
  positiveUnitsPupilTable = [positiveUnitsPupilTable; [pvalPositiveUnitsPupil10minWindows(inds) pvalPositiveUnitsPupil20minWindows(inds) pvalPositiveUnitsPupil30minWindows(inds) pvalPositiveUnitsPupil(inds)]];
  negativeUnitsPupilTable = [rPositiveUnitsPupil10minWindows(inds) rPositiveUnitsPupil20minWindows(inds) rPositiveUnitsPupil30minWindows(inds) rPositiveUnitsPupil(inds)];
  negativeUnitsPupilTable = [negativeUnitsPupilTable; negativeUnitsPupilTable.^2];
  negativeUnitsPupilTable = [negativeUnitsPupilTable; [pvalPositiveUnitsPupil10minWindows(inds) pvalPositiveUnitsPupil20minWindows(inds) pvalPositiveUnitsPupil30minWindows(inds) pvalPositiveUnitsPupil(inds)]];
  unitsPupilTable_CA = [positiveUnitsPupilTable; negativeUnitsPupilTable];
  
  % Units vs pupil: V1
  if strcmp(repository, 'allensdk')
    inds = find(ismember(areasOI, 'V1'));
    positiveUnitsPupilTable = [rPositiveUnitsPupil10minWindows(inds) rPositiveUnitsPupil20minWindows(inds) rPositiveUnitsPupil30minWindows(inds) rPositiveUnitsPupil(inds)];
    positiveUnitsPupilTable = [positiveUnitsPupilTable; positiveUnitsPupilTable.^2];
    positiveUnitsPupilTable = [positiveUnitsPupilTable; [pvalPositiveUnitsPupil10minWindows(inds) pvalPositiveUnitsPupil20minWindows(inds) pvalPositiveUnitsPupil30minWindows(inds) pvalPositiveUnitsPupil(inds)]];
    negativeUnitsPupilTable = [rPositiveUnitsPupil10minWindows(inds) rPositiveUnitsPupil20minWindows(inds) rPositiveUnitsPupil30minWindows(inds) rPositiveUnitsPupil(inds)];
    negativeUnitsPupilTable = [negativeUnitsPupilTable; negativeUnitsPupilTable.^2];
    negativeUnitsPupilTable = [negativeUnitsPupilTable; [pvalPositiveUnitsPupil10minWindows(inds) pvalPositiveUnitsPupil20minWindows(inds) pvalPositiveUnitsPupil30minWindows(inds) pvalPositiveUnitsPupil(inds)]];
    unitsPupilTable_V1 = [positiveUnitsPupilTable; negativeUnitsPupilTable];
  end
  
  if strcmp(repository, 'uol')
    MUAsTable = [MUAsTable_S1vRSC MUAsTable_S1vCA MUAsTable_RSCvCA];
    unitsTable = [unitsTable_S1vRSC unitsTable_S1vCA unitsTable_RSCvCA];
    MUAsPupilTable = [MUAsPupilTable_S1 MUAsPupilTable_RSC MUAsPupilTable_CA];
    unitsPupilTable = [unitsPupilTable_S1 unitsPupilTable_RSC unitsPupilTable_CA];
    pupilTable = [MUAsPupilTable unitsPupilTable];
  end
  
  
  %% Save the tables
  if strcmp(repository, 'uol')
    save([mainFolder filesep 'pupilCorrFractions.mat'],...
      'MUAsTable_S1vRSC', 'MUAsTable_S1vCA', 'MUAsTable_RSCvCA',...
      'unitsTable_S1vRSC', 'unitsTable_S1vCA', 'unitsTable_RSCvCA',...
      'MUAsPupilTable_S1', 'MUAsPupilTable_RSC', 'MUAsPupilTable_CA',...
      'unitsPupilTable_S1', 'unitsPupilTable_RSC', 'unitsPupilTable_CA',...
      'MUAsTable', 'unitsTable', 'MUAsPupilTable', 'unitsPupilTable',...
      'pupilTable', '-append');
  else
    save([mainFolder filesep 'pupilCorrFractions.mat'],...
      'MUAsTable_V1vCA', 'unitsTable_V1vCA', 'MUAsPupilTable_V1', 'unitsPupilTable_CA',...
      '-append');
  end
  
% else
%   if strcmp(repository, 'allensdk')
%     load([mainFolder filesep 'pupilCorrFractions.mat']);
%   elseif strcmp(repository, 'uol')
%     load([mainFolder filesep 'all' filesep 'pupilCorrFractions.mat']);
%   end
end



%% Local functions
function [dataMean, dataScatter, areaGroup, colourGroup, corrGroup] = barGroup(areaName, colourNumber, barPosition, data)

areaCode = determineArea(areaName);
colour = matlabColours(colourNumber);
bar(1+barPosition, data.positiveUnitsFrThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour, 'FaceAlpha',0.2);
bar(2+barPosition, data.positiveSignificantUnitsFrThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
bar(3+barPosition, data.negativeUnitsFrThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3), 'FaceAlpha',0.2);
bar(4+barPosition, data.negativeSignificantUnitsFrThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
dataMean = [data.positiveUnitsFrThr{1}{areaCode} data.positiveSignificantUnitsFrThr{1}{areaCode}...
  data.negativeUnitsFrThr{1}{areaCode} data.negativeSignificantUnitsFrThr{1}{areaCode}];
dataScatter = {data.positiveUnitsThr{1}{areaCode}./data.totalUnitsThr{1}{areaCode},...
  data.positiveSignificantUnitsThr{1}{areaCode}./data.totalUnitsThr{1}{areaCode},...
  data.negativeUnitsThr{1}{areaCode}./data.totalUnitsThr{1}{areaCode},...
  data.negativeSignificantUnitsThr{1}{areaCode}./data.totalUnitsThr{1}{areaCode}};
areaGroup = {areaName, areaName, areaName, areaName};
colourGroup = [colour; colour; colour; colour];
corrGroup = {'positive','positive','negative','negative'};
end


function [dataMean, dataScatter, areaGroup, colourGroup, corrGroup] = barGroupRedux(areaName, colourNumber, barPosition, data)

areaCode = determineArea(areaName);
colour = matlabColours(colourNumber);
bar(1+barPosition, data.positiveUnitsFrThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
bar(2+barPosition, data.negativeUnitsFrThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
dataMean = [data.positiveUnitsFrThr{1}{areaCode} data.negativeUnitsFrThr{1}{areaCode}];
dataScatter = {data.positiveUnitsThr{1}{areaCode}./data.totalUnitsThr{1}{areaCode},...
  data.negativeUnitsThr{1}{areaCode}./data.totalUnitsThr{1}{areaCode}};
areaGroup = {areaName, areaName};
colourGroup = [colour; colour];
corrGroup = {'positive','negative'};
end


function [dataMean, dataScatter, areaGroup, colourGroup, corrGroup] = barGroupMUAs(areaName, colourNumber, barPosition, data)

areaCode = determineArea(areaName);
colour = matlabColours(colourNumber);
bar(1+barPosition, data.positiveMUAsFrThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour, 'FaceAlpha',0.2);
bar(2+barPosition, data.positiveSignificantMUAsFrThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
bar(3+barPosition, data.negativeMUAsFrThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3), 'FaceAlpha',0.2);
bar(4+barPosition, data.negativeSignificantMUAsFrThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
dataMean = [data.positiveMUAsFrThr{1}{areaCode} data.positiveSignificantMUAsFrThr{1}{areaCode}...
  data.negativeMUAsFrThr{1}{areaCode} data.negativeSignificantMUAsFrThr{1}{areaCode}];
dataScatter = {data.positiveMUAsThr{1}{areaCode}./data.totalMUAsThr{1}{areaCode},...
  data.positiveSignificantMUAsThr{1}{areaCode}./data.totalMUAsThr{1}{areaCode},...
  data.negativeMUAsThr{1}{areaCode}./data.totalMUAsThr{1}{areaCode},...
  data.negativeSignificantMUAsThr{1}{areaCode}./data.totalMUAsThr{1}{areaCode}};
areaGroup = {areaName, areaName, areaName, areaName};
colourGroup = [colour; colour; colour; colour];
corrGroup = {'positive','positive','negative','negative'};
end


function [dataMean, dataScatter, areaGroup, colourGroup, corrGroup] = barGroupMUAsRedux(areaName, colourNumber, barPosition, data)

areaCode = determineArea(areaName);
colour = matlabColours(colourNumber);
bar(1+barPosition, data.positiveMUAsFrThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
bar(2+barPosition, data.negativeMUAsFrThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
dataMean = [data.positiveMUAsFrThr{1}{areaCode} data.negativeMUAsFrThr{1}{areaCode}];
dataScatter = {data.positiveMUAsThr{1}{areaCode}./data.totalMUAsThr{1}{areaCode},...
  data.negativeMUAsThr{1}{areaCode}./data.totalMUAsThr{1}{areaCode}};
areaGroup = {areaName, areaName};
colourGroup = [colour; colour];
corrGroup = {'positive','negative'};
end


function [p, pSig] = ttestGroup(areaName, data)

areaCode = determineArea(areaName);
[~,p] = ttest2(data.positiveUnitsThr{1}{areaCode}./data.totalUnitsThr{1}{areaCode},...
  data.negativeUnitsThr{1}{areaCode}./data.totalUnitsThr{1}{areaCode});
[~,pSig] = ttest2(data.positiveSignificantUnitsThr{1}{areaCode}./data.totalUnitsThr{1}{areaCode},...
  data.negativeSignificantUnitsThr{1}{areaCode}./data.totalUnitsThr{1}{areaCode});
end


function [p, pSig] = ttestGroupMUAs(areaName, data)

areaCode = determineArea(areaName);
[~,p] = ttest2(data.positiveMUAsThr{1}{areaCode}./data.totalMUAsThr{1}{areaCode},...
  data.negativeMUAsThr{1}{areaCode}./data.totalMUAsThr{1}{areaCode});
[~,pSig] = ttest2(data.positiveSignificantMUAsThr{1}{areaCode}./data.totalMUAsThr{1}{areaCode},...
  data.negativeSignificantMUAsThr{1}{areaCode}./data.totalMUAsThr{1}{areaCode});
end


function [fH, scatterGroups, areaGroups, colourGroups, colorCode, corrGroups] = barPlotUnits(data, area1, area2, area3, area4, yLim)

[pVB, pVBSig] = ttestGroup(area1, data);
[plS1, plS1Sig] = ttestGroup(area2, data);
[plRSC, plRSCSig] = ttestGroup(area3, data);
[pCA, pCASig] = ttestGroup(area4, data);


%% Draw combined bar graphs
fH1 = figProperties('Bar plot for units', 'normalized', [0, .005, .97, .90], 'w', 'on');
hold on


% Data bars
gap = 1;
nBarsPerGroup = 4;
[barGroup1, scatterGroup1, areaGroup1, colourGroup1, corrGroup1] = barGroup(area1, 1, 0*(nBarsPerGroup+gap), data);
if strcmp(area2,'LGN') || strcmp(area2,'lLGN') || strcmp(area2,'rLGN')
  [barGroup2, scatterGroup2, areaGroup2, colourGroup2, corrGroup2] = barGroup(area2, 2, 1*(nBarsPerGroup+gap), data);
elseif strcmp(area2,'S1') || strcmp(area2,'lS1') || strcmp(area2,'rS1')
  [barGroup2, scatterGroup2, areaGroup2, colourGroup2, corrGroup2] = barGroup(area2, 13, 1*(nBarsPerGroup+gap), data);
end
if strcmp(area3,'V1') || strcmp(area3,'lV1') || strcmp(area3,'rV1')
  [barGroup3, scatterGroup3, areaGroup3, colourGroup3, corrGroup3] = barGroup(area3, 11, 2*(nBarsPerGroup+gap), data);
elseif strcmp(area3,'RSC') || strcmp(area3,'lRSC') || strcmp(area3,'rRSC')
  [barGroup3, scatterGroup3, areaGroup3, colourGroup3, corrGroup3] = barGroup(area3, 3, 2*(nBarsPerGroup+gap), data);
end
[barGroup4, scatterGroup4, areaGroup4, colourGroup4, corrGroup4] = barGroup(area4, 4, 3*(nBarsPerGroup+gap), data);


% Error bars
gaps = 5:5:15;
bars = sort([1 1+gaps 2 2+gaps 3 3+gaps 4 4+gaps]);
dataMean = [barGroup1 barGroup2 barGroup3 barGroup4];
err = [data.positiveUnitsFrCI95Thr{1}{determineArea(area1)} data.positiveSignificantUnitsFrCI95Thr{1}{determineArea(area1)}...
  data.negativeUnitsFrCI95Thr{1}{determineArea(area1)} data.negativeSignificantUnitsFrCI95Thr{1}{determineArea(area1)}...
  data.positiveUnitsFrCI95Thr{1}{determineArea(area2)} data.positiveSignificantUnitsFrCI95Thr{1}{determineArea(area2)}...
  data.negativeUnitsFrCI95Thr{1}{determineArea(area2)} data.negativeSignificantUnitsFrCI95Thr{1}{determineArea(area2)}...
  data.positiveUnitsFrCI95Thr{1}{determineArea(area3)} data.positiveSignificantUnitsFrCI95Thr{1}{determineArea(area3)}...
  data.negativeUnitsFrCI95Thr{1}{determineArea(area3)} data.negativeSignificantUnitsFrCI95Thr{1}{determineArea(area3)}...
  data.positiveUnitsFrCI95Thr{1}{determineArea(area4)} data.positiveSignificantUnitsFrCI95Thr{1}{determineArea(area4)}...
  data.negativeUnitsFrCI95Thr{1}{determineArea(area4)} data.negativeSignificantUnitsFrCI95Thr{1}{determineArea(area4)}];
err(:, ~dataMean) = 0;

if numel(bars) == numel(dataMean)
  er = errorbar(bars,dataMean,err(2,:),err(1,:));
  er.Color = [0 0 0];
  er.LineStyle = 'none';
end


% Scatter
scatterGroups = [scatterGroup1 scatterGroup2 scatterGroup3 scatterGroup4];
areaGroups = [areaGroup1 areaGroup2 areaGroup3 areaGroup4];
colourGroups = [colourGroup1; colourGroup2; colourGroup3; colourGroup4;];
colorCode = [30 30 30]./255;
corrGroups = [corrGroup1 corrGroup2 corrGroup3 corrGroup4];
for iBar = 1:numel(bars)
  scatter(bars(iBar)*ones(size(scatterGroups{iBar}))', scatterGroups{iBar}',...
    'MarkerEdgeColor',colorCode, 'jitter','on'); %colourGroups(iBar,:));
end


% Graph adjustments
xTickPos = [gaps 20]-2.5;
if nargin < 6 || isempty(yLim)
  yLim = ylim;
end
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Arial', 25, 4/3, 2, [0.003 0.003], 'out', 'on', 'k', {}, [],...
  xTickPos, 'on', 'k', {'Proportion'}, yLim, yticks);
xTickLabel = {area1,area2,area3,area4};
for iLabel = 1:numel(xTickLabel)
  xTickLabel{iLabel} = strrep(xTickLabel{iLabel}, 'l','');
  xTickLabel{iLabel} = strrep(xTickLabel{iLabel}, 'r','');
end
ax1.XTickLabel = xTickLabel;

xLim = xlim;
xAxisLength = xLim(2)-xLim(1);
yAxisLength = yLim(2)-yLim(1);
textStr = ['t-test +v- p-val for all units: ' num2str(pVB) '(' area1 ') ' num2str(plS1) '(' area2 ') '...
  num2str(plRSC) '(' area3 ') ' num2str(pCA) '(' area4 ') '];
text(xLim(2)-xAxisLength*0.85, yLim(2)-yAxisLength*0.05, textStr, 'FontSize',20);
textStr = ['t-test +v- p-val for significant units only: ' num2str(pVBSig) '(' area1 ') ' num2str(plS1Sig) '(' area2 ') '...
  num2str(plRSCSig) '(' area3 ') ' num2str(pCASig) '(' area4 ') '];
text(xLim(2)-xAxisLength*0.85, yLim(2)-yAxisLength*0.01, textStr, 'FontSize',20);
text(1-0.2,0.03*yLim(2), 'All', 'FontSize',20)
text(2-0.25,0.03*yLim(2), 'Sig', 'FontSize',20)
for iBar = bars([1 2 5 6 9 10 13 14])
  text(iBar-0.075,-0.01*yLim(2), '+', 'FontSize',12)
  text(iBar+2-0.075,-0.01*yLim(2), '-', 'FontSize',12)
end

ylim([0 1.05]);


%% Draw bar graphs without bars for significant units only
fH2 = figProperties('Bar plot for units', 'normalized', [0, .005, .97, .90], 'w', 'on');
hold on


% Data bars
gap = 1;
nBarsPerGroup = 2;
[barGroupRedux1, scatterGroupRedux1] = barGroupMUAsRedux(area1, 1, 0*(nBarsPerGroup+gap), data);
if strcmp(area2,'LGN') || strcmp(area2,'lLGN') || strcmp(area2,'rLGN')
  [barGroupRedux2, scatterGroupRedux2] = barGroupMUAsRedux(area2, 2, 1*(nBarsPerGroup+gap), data);
elseif strcmp(area2,'S1') || strcmp(area2,'lS1') || strcmp(area2,'rS1')
  [barGroupRedux2, scatterGroupRedux2] = barGroupMUAsRedux(area2, 13, 1*(nBarsPerGroup+gap), data);
end
if strcmp(area3,'V1') || strcmp(area3,'lV1') || strcmp(area3,'rV1')
  [barGroupRedux3, scatterGroupRedux3] = barGroupMUAsRedux(area3, 11, 2*(nBarsPerGroup+gap), data);
elseif strcmp(area3,'RSC') || strcmp(area3,'lRSC') || strcmp(area3,'rRSC')
  [barGroupRedux3, scatterGroupRedux3] = barGroupMUAsRedux(area3, 3, 2*(nBarsPerGroup+gap), data);
end
[barGroupRedux4, scatterGroupRedux4] = barGroupMUAsRedux(area4, 4, 3*(nBarsPerGroup+gap), data);


% Error bars
gaps = 3:3:9;
bars = sort([1 1+gaps 2 2+gaps]);
dataMean = [barGroupRedux1 barGroupRedux2 barGroupRedux3 barGroupRedux4];
err = [data.positiveUnitsFrCI95Thr{1}{determineArea(area1)}...
  data.negativeUnitsFrCI95Thr{1}{determineArea(area1)}...
  data.positiveUnitsFrCI95Thr{1}{determineArea(area2)}...
  data.negativeUnitsFrCI95Thr{1}{determineArea(area2)}...
  data.positiveUnitsFrCI95Thr{1}{determineArea(area3)}...
  data.negativeUnitsFrCI95Thr{1}{determineArea(area3)}...
  data.positiveUnitsFrCI95Thr{1}{determineArea(area4)}...
  data.negativeUnitsFrCI95Thr{1}{determineArea(area4)}];

if numel(bars) == numel(dataMean)
  er = errorbar(bars,dataMean,err(2,:),err(1,:));
  er.Color = [0 0 0];
  er.LineStyle = 'none';
end


% Scatter
scatterGroupsRedux = [scatterGroupRedux1 scatterGroupRedux2 scatterGroupRedux3 scatterGroupRedux4];
for iBar = 1:numel(bars)
  scatter(bars(iBar)*ones(size(scatterGroupsRedux{iBar}))', scatterGroupsRedux{iBar}',...
    'MarkerEdgeColor',colorCode, 'jitter','on'); %colourGroups(iBar,:));
end


% Graph adjustments
xTickPos = [gaps 12]-1.5;
if nargin < 6 || isempty(yLim)
  yLim = ylim;
end
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Arial', 25, 4/3, 2, [0.003 0.003], 'out', 'on', 'k', {}, [],...
  xTickPos, 'on', 'k', {'Proportion'}, yLim, yticks);
ax1.XTickLabel = xTickLabel;

for iBar = bars([1 3 5 7])
  text(iBar-0.075,-0.01*yLim(2), '+', 'FontSize',20)
  text(iBar+1-0.075,-0.01*yLim(2), '-', 'FontSize',20)
end

ylim([0 1.05]);

fH = [fH1 fH2];
end


function [fH, scatterGroups, areaGroups, colourGroups, colorCode, corrGroups] = barPlotMUAs(data, area1, area2, area3, area4, yLim)

[pVB, pVBSig] = ttestGroupMUAs(area1, data);
[plS1, plS1Sig] = ttestGroupMUAs(area2, data);
[plRSC, plRSCSig] = ttestGroupMUAs(area3, data);
[pCA, pCASig] = ttestGroupMUAs(area4, data);


%% Draw the bar graphs
fH1 = figProperties('Bar plot for MUAs', 'normalized', [0, .005, .97, .90], 'w', 'on');
hold on


% Data bars
gap = 1;
nBarsPerGroup = 4;
[barGroup1, scatterGroup1, areaGroup1, colourGroup1, corrGroup1] = barGroupMUAs(area1, 1, 0*(nBarsPerGroup+gap), data);
if strcmp(area2,'LGN') || strcmp(area2,'lLGN') || strcmp(area2,'rLGN')
  [barGroup2, scatterGroup2, areaGroup2, colourGroup2, corrGroup2] = barGroupMUAs(area2, 2, 1*(nBarsPerGroup+gap), data);
elseif strcmp(area2,'S1') || strcmp(area2,'lS1') || strcmp(area2,'rS1')
  [barGroup2, scatterGroup2, areaGroup2, colourGroup2, corrGroup2] = barGroupMUAs(area2, 13, 1*(nBarsPerGroup+gap), data);
end
if strcmp(area3,'V1') || strcmp(area3,'lV1') || strcmp(area3,'rV1')
  [barGroup3, scatterGroup3, areaGroup3, colourGroup3, corrGroup3] = barGroupMUAs(area3, 11, 2*(nBarsPerGroup+gap), data);
elseif strcmp(area3,'RSC') || strcmp(area3,'lRSC') || strcmp(area3,'rRSC')
  [barGroup3, scatterGroup3, areaGroup3, colourGroup3, corrGroup3] = barGroupMUAs(area3, 3, 2*(nBarsPerGroup+gap), data);
end
[barGroup4, scatterGroup4, areaGroup4, colourGroup4, corrGroup4] = barGroupMUAs(area4, 4, 3*(nBarsPerGroup+gap), data);


% Error bars
gaps = 5:5:15;
bars = sort([1 1+gaps 2 2+gaps 3 3+gaps 4 4+gaps]);
dataMean = [barGroup1 barGroup2 barGroup3 barGroup4];
err = [data.positiveMUAsFrCI95Thr{1}{determineArea(area1)} data.positiveSignificantMUAsFrCI95Thr{1}{determineArea(area1)}...
  data.negativeMUAsFrCI95Thr{1}{determineArea(area1)} data.negativeSignificantMUAsFrCI95Thr{1}{determineArea(area1)}...
  data.positiveMUAsFrCI95Thr{1}{determineArea(area2)} data.positiveSignificantMUAsFrCI95Thr{1}{determineArea(area2)}...
  data.negativeMUAsFrCI95Thr{1}{determineArea(area2)} data.negativeSignificantMUAsFrCI95Thr{1}{determineArea(area2)}...
  data.positiveMUAsFrCI95Thr{1}{determineArea(area3)} data.positiveSignificantMUAsFrCI95Thr{1}{determineArea(area3)}...
  data.negativeMUAsFrCI95Thr{1}{determineArea(area3)} data.negativeSignificantMUAsFrCI95Thr{1}{determineArea(area3)}...
  data.positiveMUAsFrCI95Thr{1}{determineArea(area4)} data.positiveSignificantMUAsFrCI95Thr{1}{determineArea(area4)}...
  data.negativeMUAsFrCI95Thr{1}{determineArea(area4)} data.negativeSignificantMUAsFrCI95Thr{1}{determineArea(area4)}];
err(:, ~dataMean) = 0;

if numel(bars) == numel(dataMean)
  er = errorbar(bars,dataMean,err(2,:),err(1,:));
  er.Color = [0 0 0];
  er.LineStyle = 'none';
end


% Scatter
scatterGroups = [scatterGroup1 scatterGroup2 scatterGroup3 scatterGroup4];
areaGroups = [areaGroup1 areaGroup2 areaGroup3 areaGroup4];
colourGroups = [colourGroup1; colourGroup2; colourGroup3; colourGroup4;];
colorCode = [30 30 30]./255;
corrGroups = [corrGroup1 corrGroup2 corrGroup3 corrGroup4];
for iBar = 1:numel(bars)
  scatter(bars(iBar)*ones(size(scatterGroups{iBar}))', scatterGroups{iBar}',...
    'MarkerEdgeColor',colorCode, 'jitter','on'); %colourGroups(iBar,:));
end


% Graph adjustments
xTickPos = [gaps 20]-2.5;
if nargin < 6 || isempty(yLim)
  yLim = ylim;
end
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Arial', 25, 4/3, 2, [0.003 0.003], 'out', 'on', 'k', {}, [],...
  xTickPos, 'on', 'k', {'Proportion'}, yLim, yticks);
xTickLabel = {area1,area2,area3,area4};
for iLabel = 1:numel(xTickLabel)
  xTickLabel{iLabel} = strrep(xTickLabel{iLabel}, 'l','');
  xTickLabel{iLabel} = strrep(xTickLabel{iLabel}, 'r','');
end
ax1.XTickLabel = xTickLabel;

xLim = xlim;
xAxisLength = xLim(2)-xLim(1);
yAxisLength = yLim(2)-yLim(1);
textStr = ['t-test +v- p-val for all units: ' num2str(pVB) '(' area1 ') ' num2str(plS1) '(' area2 ') '...
  num2str(plRSC) '(' area3 ') ' num2str(pCA) '(' area4 ') '];
text(xLim(2)-xAxisLength*0.85, yLim(2)-yAxisLength*0.05, textStr, 'FontSize',20);
textStr = ['t-test +v- p-val for significant units only: ' num2str(pVBSig) '(' area1 ') ' num2str(plS1Sig) '(' area2 ') '...
  num2str(plRSCSig) '(' area3 ') ' num2str(pCASig) '(' area4 ') '];
text(xLim(2)-xAxisLength*0.85, yLim(2)-yAxisLength*0.01, textStr, 'FontSize',20);
text(1-0.2,0.03*yLim(2), 'All', 'FontSize',20)
text(2-0.25,0.03*yLim(2), 'Sig', 'FontSize',20)
for iBar = bars([1 2 5 6 9 10 13 14])
  text(iBar-0.075,-0.01*yLim(2), '+', 'FontSize',12)
  text(iBar+2-0.075,-0.01*yLim(2), '-', 'FontSize',12)
end

ylim([0 1.05]);


%% Draw bar graphs without bars for significant units+MUAs only
fH2 = figProperties('Bar plot for MUAs', 'normalized', [0, .005, .97, .90], 'w', 'on');
hold on


% Data bars
gap = 1;
nBarsPerGroup = 2;
[barGroupRedux1, scatterGroupRedux1] = barGroupMUAsRedux(area1, 1, 0*(nBarsPerGroup+gap), data);
if strcmp(area2,'LGN') || strcmp(area2,'lLGN') || strcmp(area2,'rLGN')
  [barGroupRedux2, scatterGroupRedux2] = barGroupMUAsRedux(area2, 2, 1*(nBarsPerGroup+gap), data);
elseif strcmp(area2,'S1') || strcmp(area2,'lS1') || strcmp(area2,'rS1')
  [barGroupRedux2, scatterGroupRedux2] = barGroupMUAsRedux(area2, 13, 1*(nBarsPerGroup+gap), data);
end
if strcmp(area3,'V1') || strcmp(area3,'lV1') || strcmp(area3,'rV1')
  [barGroupRedux3, scatterGroupRedux3] = barGroupMUAsRedux(area3, 11, 2*(nBarsPerGroup+gap), data);
elseif strcmp(area3,'RSC') || strcmp(area3,'lRSC') || strcmp(area3,'rRSC')
  [barGroupRedux3, scatterGroupRedux3] = barGroupMUAsRedux(area3, 3, 2*(nBarsPerGroup+gap), data);
end
[barGroupRedux4, scatterGroupRedux4] = barGroupMUAsRedux(area4, 4, 3*(nBarsPerGroup+gap), data);


% Error bars
gaps = 3:3:9;
bars = sort([1 1+gaps 2 2+gaps]);
dataMean = [barGroupRedux1 barGroupRedux2 barGroupRedux3 barGroupRedux4];
err = [data.positiveMUAsFrCI95Thr{1}{determineArea(area1)}...
  data.negativeMUAsFrCI95Thr{1}{determineArea(area1)}...
  data.positiveMUAsFrCI95Thr{1}{determineArea(area2)}...
  data.negativeMUAsFrCI95Thr{1}{determineArea(area2)}...
  data.positiveMUAsFrCI95Thr{1}{determineArea(area3)}...
  data.negativeMUAsFrCI95Thr{1}{determineArea(area3)}...
  data.positiveMUAsFrCI95Thr{1}{determineArea(area4)}...
  data.negativeMUAsFrCI95Thr{1}{determineArea(area4)}];
err(:, ~dataMean) = 0;

if numel(bars) == numel(dataMean)
  er = errorbar(bars,dataMean,err(2,:),err(1,:));
  er.Color = [0 0 0];
  er.LineStyle = 'none';
end


% Scatter
scatterGroupsRedux = [scatterGroupRedux1 scatterGroupRedux2 scatterGroupRedux3 scatterGroupRedux4];
for iBar = 1:numel(bars)
  scatter(bars(iBar)*ones(size(scatterGroupsRedux{iBar}))', scatterGroupsRedux{iBar}',...
    'MarkerEdgeColor',colorCode, 'jitter','on'); %colourGroups(iBar,:));
end


% Graph adjustments
xTickPos = [gaps 12]-1.5;
if nargin < 6 || isempty(yLim)
  yLim = ylim;
end
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Arial', 25, 4/3, 2, [0.003 0.003], 'out', 'on', 'k', {}, [],...
  xTickPos, 'on', 'k', {'Proportion'}, yLim, yticks);
ax1.XTickLabel = xTickLabel;

for iBar = bars([1 3 5 7])
  text(iBar-0.075,-0.01*yLim(2), '+', 'FontSize',20)
  text(iBar+1-0.075,-0.01*yLim(2), '-', 'FontSize',20)
end

ylim([0 1.05]);

fH = [fH1 fH2];
end


function runAnova(scatterGroups, corrGroups, areaGroups, filename)

anovaData = [];
for iGroup = 1:numel(scatterGroups)
  anovaData = [anovaData scatterGroups{iGroup}']; %#ok<*AGROW>
  for iElement = 1:numel(scatterGroups{iGroup})
    if iGroup == 1 && iElement == 1
      anovaCorrDirectionFactor = corrGroups(iGroup);
      anovaAreaFactor = areaGroups(iGroup);
    else
      anovaCorrDirectionFactor = [anovaCorrDirectionFactor corrGroups{iGroup}];
      anovaAreaFactor = [anovaAreaFactor areaGroups{iGroup}];
    end
  end
end
[~, anovaOutput] = anovan(anovaData(1:2:end), {anovaCorrDirectionFactor(1:2:end),anovaAreaFactor(1:2:end)},...
  'model','interaction', 'varnames',{'corrSign','area'});
writecell(anovaOutput, [filename '_anova_all' '.txt'])
[~, anovaOutput] = anovan(anovaData(2:2:end), {anovaCorrDirectionFactor(2:2:end),anovaAreaFactor(2:2:end)},...
  'model','interaction', 'varnames',{'corrSign','area'});
writecell(anovaOutput, [filename '_anova_sig' '.txt'])
end


function [fH, r, pval] = corrPlot(areas, data, type, windowSize, outputDir)

if ~isempty(windowSize)
  windowSize = [num2str(windowSize) 'minWindows'];
end
areaCombos = nchoosek(1:numel(areas),2);
fH = zeros(1,size(areaCombos,1));
r = zeros(1,size(areaCombos,1));
pval = zeros(1,size(areaCombos,2));
for iCombo = 1:size(areaCombos,1)
  areaName1 = areas{areaCombos(iCombo,1)};
  areaName2 = areas{areaCombos(iCombo,2)};
  areaCode1 = determineArea(areaName1);
  areaCode2 = determineArea(areaName2);
  if strcmpi(type, 'positiveUnits')
    areaInd1 = ismember(data.(['positiveUnitsThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveUnitsThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveUnitsThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveUnitsThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveUnitsThr' windowSize]){1}{areaCode1}(areaInd1)./data.(['totalUnitsThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveUnitsThr' windowSize]){1}{areaCode2}(areaInd2)./data.(['totalUnitsThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveUnitsThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveUnitsThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Positive unit fraction in ' areaName1], ['Positive unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Positive unit fraction in ' areaName1 ' v ' areaName2 ': Full recordings'];
    else
      figName = ['Positive unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveSignificantUnits')
    areaInd1 = ismember(data.(['positiveSignificantUnitsThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveSignificantUnitsThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveSignificantUnitsThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveSignificantUnitsThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveSignificantUnitsThr' windowSize]){1}{areaCode1}(areaInd1)./data.(['totalUnitsThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveSignificantUnitsThr' windowSize]){1}{areaCode2}(areaInd2)./data.(['totalUnitsThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveSignificantUnitsThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveSignificantUnitsThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Positive significant unit fraction in ' areaName1], ['Positive significant unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Positive significant unit fraction in ' areaName1 ' v ' areaName2 ': Full recordings'];
    else
      figName = ['Positive significant unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeUnits')
    areaInd1 = ismember(data.(['negativeUnitsThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeUnitsThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeUnitsThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeUnitsThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeUnitsThr' windowSize]){1}{areaCode1}(areaInd1)./data.(['totalUnitsThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeUnitsThr' windowSize]){1}{areaCode2}(areaInd2)./data.(['totalUnitsThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeUnitsThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeUnitsThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Negative unit fraction in ' areaName1], ['Negative unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Negative unit fraction in ' areaName1 ' v ' areaName2 ': Full recordings'];
    else
      figName = ['Negative unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeSignificantUnits')
    areaInd1 = ismember(data.(['negativeSignificantUnitsThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeSignificantUnitsThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeSignificantUnitsThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeSignificantUnitsThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeSignificantUnitsThr' windowSize]){1}{areaCode1}(areaInd1)./data.(['totalUnitsThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeSignificantUnitsThr' windowSize]){1}{areaCode2}(areaInd2)./data.(['totalUnitsThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeSignificantUnitsThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeSignificantUnitsThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Negative significant unit fraction in ' areaName1], ['Negative significant unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Negative significant unit fraction in ' areaName1 ' v ' areaName2 ': Full recordings'];
    else
      figName = ['Negative significant unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'neutralUnits')
    areaInd1 = ismember(data.(['neutralUnitsThr' windowSize 'RecID']){1}{areaCode1}, data.(['neutralUnitsThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['neutralUnitsThr' windowSize 'RecID']){1}{areaCode2}, data.(['neutralUnitsThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['neutralUnitsThr' windowSize]){1}{areaCode1}(areaInd1)./data.(['totalUnitsThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['neutralUnitsThr' windowSize]){1}{areaCode2}(areaInd2)./data.(['totalUnitsThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['neutralUnitsThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['neutralUnitsThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Neutral unit fraction in ' areaName1], ['Neutral unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Neutral unit fraction in ' areaName1 ' v ' areaName2 ': Full recordings'];
    else
      figName = ['Neutral unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveMUAs')
    areaInd1 = ismember(data.(['positiveMUAsThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveMUAsThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveMUAsThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveMUAsThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveMUAsThr' windowSize]){1}{areaCode1}(areaInd1)./data.(['totalMUAsThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveMUAsThr' windowSize]){1}{areaCode2}(areaInd2)./data.(['totalMUAsThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveMUAsThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveMUAsThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Positive unit+MUA fraction in ' areaName1], ['Positive unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Positive unit+MUA fraction in ' areaName1 ' v ' areaName2 ': Full recordings'];
    else
      figName = ['Positive unit+MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveSignificantMUAs')
    areaInd1 = ismember(data.(['positiveSignificantMUAsThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveSignificantMUAsThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveSignificantMUAsThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveSignificantMUAsThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveSignificantMUAsThr' windowSize]){1}{areaCode1}(areaInd1)./data.(['totalMUAsThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveSignificantMUAsThr' windowSize]){1}{areaCode2}(areaInd2)./data.(['totalMUAsThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveSignificantMUAsThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveSignificantMUAsThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Positive significant unit+MUA fraction in ' areaName1], ['Positive significant unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Positive significant unit+MUA fraction in ' areaName1 ' v ' areaName2 ': Full recordings'];
    else
      figName = ['Positive significant unit+MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeMUAs')
    areaInd1 = ismember(data.(['negativeMUAsThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeMUAsThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeMUAsThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeMUAsThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeMUAsThr' windowSize]){1}{areaCode1}(areaInd1)./data.(['totalMUAsThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeMUAsThr' windowSize]){1}{areaCode2}(areaInd2)./data.(['totalMUAsThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeMUAsThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeMUAsThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Negative unit+MUA fraction in ' areaName1], ['Negative unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Negative unit+MUA fraction in ' areaName1 ' v ' areaName2 ': Full recordings'];
    else
      figName = ['Negative unit+MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeSignificantMUAs')
    areaInd1 = ismember(data.(['negativeSignificantMUAsThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeSignificantMUAsThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeSignificantMUAsThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeSignificantMUAsThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeSignificantMUAsThr' windowSize]){1}{areaCode1}(areaInd1)./data.(['totalMUAsThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeSignificantMUAsThr' windowSize]){1}{areaCode2}(areaInd2)./data.(['totalMUAsThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeSignificantMUAsThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeSignificantMUAsThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Negative significant unit+MUA fraction in ' areaName1], ['Negative significant unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Negative significant unit+MUA fraction in ' areaName1 ' v ' areaName2 ': Full recordings'];
    else
      figName = ['Negative significant unit+MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'neutralMUAs')
    areaInd1 = ismember(data.(['neutralMUAsThr' windowSize 'RecID']){1}{areaCode1}, data.(['neutralMUAsThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['neutralMUAsThr' windowSize 'RecID']){1}{areaCode2}, data.(['neutralMUAsThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['neutralMUAsThr' windowSize]){1}{areaCode1}(areaInd1)./data.(['totalMUAsThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['neutralMUAsThr' windowSize]){1}{areaCode2}(areaInd2)./data.(['totalMUAsThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['neutralMUAsThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['neutralMUAsThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Neutral unit+MUA fraction in ' areaName1], ['Neutral unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Neutral unit+MUA fraction in ' areaName1 ' v ' areaName2 ': Full recordings'];
    else
      figName = ['Neutral unit+MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  end
  
  xLim = xlim;
  xAxisLength = xLim(2)-xLim(1);
  yLim = ylim;
  yAxisLength = yLim(2)-yLim(1);
  textStr = ['r=' num2str(r(iCombo))];
  text(xLim(2)-xAxisLength*0.15, yLim(1)+yAxisLength*0.10, textStr, 'FontSize',20);
  textStr = ['p=' num2str(pval(iCombo))];
  text(xLim(2)-xAxisLength*0.15, yLim(1)+yAxisLength*0.05, textStr, 'FontSize',20);
  
  xLabel = get(gca, 'xlabel');
  yLabel = get(gca, 'ylabel');
  ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 20, 4/3, 2, [0.005 0], 'out',...
    'on', 'k', {xLabel.String}, xlim, xticks,...
    'on', 'k', {yLabel.String}, ylim, yticks);
  set(gcf,'color','w');
  
  %title(figName);
  set(fH(iCombo), 'Name',figName);
  if isempty(areaInd1) && isempty(areaInd2)
    continue
  end
  figName = strrep(figName, ':', '_');
  filename = [outputDir filesep figName];
  filename = strrep(filename, ' ', '_');
  figSize = 15;
  label = [2.75 2.75];
  margin = [0.5 0.5];
  width = figSize-label(1)-margin(1);
  height = figSize-label(2)-margin(2);
  paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
  savefig(fH(iCombo), filename, 'compact');
  exportFig(fH(iCombo), [filename '.png'],'-dpng','-r300', paperSize);
end
end


function [fH, r, pval] = corrPlotPupil(areas, data, type, windowSize, outputDir)

if ~isempty(windowSize)
  windowSize = [num2str(windowSize) 'minWindows'];
end
areaCombos = (1:numel(areas))';
fH = zeros(1,size(areaCombos,1));
r = zeros(1,size(areaCombos,1));
pval = zeros(1,size(areaCombos,2));
for iCombo = 1:size(areaCombos,1)
  areaName = areas{areaCombos(iCombo,1)};
  areaCode = determineArea(areaName);
  if strcmpi(type, 'positiveUnits')
    pupilInd = ismember(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}, data.(['positiveUnitsThr' windowSize 'RecID']){1}{areaCode(1)});
    areaVec = data.(['positiveUnitsThr' windowSize]){1}{areaCode(1)}./data.(['totalUnitsThr' windowSize]){1}{areaCode(1)};
    pupilVec = data.(['meanPupilArea' windowSize]){1}{areaCode(1)}(pupilInd);
    areaVec = reduce2unique(data.(['positiveUnitsThr' windowSize 'RecID']){1}{areaCode(1)}, areaVec);
    pupilVec = reduce2unique(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}(pupilInd), pupilVec);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec', pupilVec', 'Spearman');
    fH(iCombo) = xyPlot(pupilVec, areaVec, 1, 'Mean pupil area size (a.u.)', ['Positive unit fraction in ' areaName]);
    if isempty(windowSize)
      figName = ['Positive unit fraction in ' areaName ' v mean pupil area size: Full recordings'];
    else
      figName = ['Positive unit fraction in ' areaName ' v mean pupil area size: ' windowSize];
    end
  elseif strcmpi(type, 'positiveSignificantUnits')
    pupilInd = ismember(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}, data.(['positiveSignificantUnitsThr' windowSize 'RecID']){1}{areaCode(1)});
    areaVec = data.(['positiveSignificantUnitsThr' windowSize]){1}{areaCode(1)}./data.(['totalUnitsThr' windowSize]){1}{areaCode(1)};
    pupilVec = data.(['meanPupilArea' windowSize]){1}{areaCode(1)}(pupilInd);
    areaVec = reduce2unique(data.(['positiveSignificantUnitsThr' windowSize 'RecID']){1}{areaCode(1)}, areaVec);
    pupilVec = reduce2unique(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}(pupilInd), pupilVec);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec', pupilVec', 'Spearman');
    fH(iCombo) = xyPlot(pupilVec, areaVec, 1, 'Mean pupil area size (a.u.)', ['Positive Significant unit fraction in ' areaName]);
    if isempty(windowSize)
      figName = ['Positive Significant unit fraction in ' areaName ' v mean pupil area size: Full recordings'];
    else
      figName = ['Positive Significant unit fraction in ' areaName ' v mean pupil area size: ' windowSize];
    end
  elseif strcmpi(type, 'negativeUnits')
    pupilInd = ismember(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}, data.(['negativeUnitsThr' windowSize 'RecID']){1}{areaCode(1)});
    areaVec = data.(['negativeUnitsThr' windowSize]){1}{areaCode(1)}./data.(['totalUnitsThr' windowSize]){1}{areaCode(1)};
    pupilVec = data.(['meanPupilArea' windowSize]){1}{areaCode(1)}(pupilInd);
    areaVec = reduce2unique(data.(['negativeUnitsThr' windowSize 'RecID']){1}{areaCode(1)}, areaVec);
    pupilVec = reduce2unique(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}(pupilInd), pupilVec);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec', pupilVec', 'Spearman');
    fH(iCombo) = xyPlot(pupilVec, areaVec, 1, 'Mean pupil area size (a.u.)', ['Negative unit fraction in ' areaName]);
    if isempty(windowSize)
      figName = ['Negative unit fraction in ' areaName ' v mean pupil area size: Full recordings'];
    else
      figName = ['Negative unit fraction in ' areaName ' v mean pupil area size: ' windowSize];
    end
  elseif strcmpi(type, 'negativeSignificantUnits')
    pupilInd = ismember(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}, data.(['negativeSignificantUnitsThr' windowSize 'RecID']){1}{areaCode(1)});
    areaVec = data.(['negativeSignificantUnitsThr' windowSize]){1}{areaCode(1)}./data.(['totalUnitsThr' windowSize]){1}{areaCode(1)};
    pupilVec = data.(['meanPupilArea' windowSize]){1}{areaCode(1)}(pupilInd);
    areaVec = reduce2unique(data.(['negativeSignificantUnitsThr' windowSize 'RecID']){1}{areaCode(1)}, areaVec);
    pupilVec = reduce2unique(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}(pupilInd), pupilVec);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec', pupilVec', 'Spearman');
    fH(iCombo) = xyPlot(pupilVec, areaVec, 1, 'Mean pupil area size (a.u.)', ['Negative Significant unit fraction in ' areaName]);
    if isempty(windowSize)
      figName = ['Negative Significant unit fraction in ' areaName ' v mean pupil area size: Full recordings'];
    else
      figName = ['Negative Significant unit fraction in ' areaName ' v mean pupil area size: ' windowSize];
    end
  elseif strcmpi(type, 'neutralUnits')
    pupilInd = ismember(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}, data.(['neutralUnitsThr' windowSize 'RecID']){1}{areaCode(1)});
    areaVec = data.(['neutralUnitsThr' windowSize]){1}{areaCode(1)}./data.(['totalUnitsThr' windowSize]){1}{areaCode(1)};
    pupilVec = data.(['meanPupilArea' windowSize]){1}{areaCode(1)}(pupilInd);
    areaVec = reduce2unique(data.(['neutralUnitsThr' windowSize 'RecID']){1}{areaCode(1)}, areaVec);
    pupilVec = reduce2unique(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}(pupilInd), pupilVec);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec', pupilVec', 'Spearman');
    fH(iCombo) = xyPlot(pupilVec, areaVec, 1, 'Mean pupil area size (a.u.)', ['Neutral unit fraction in ' areaName]);
    if isempty(windowSize)
      figName = ['Neutral unit fraction in ' areaName ' v mean pupil area size: Full recordings'];
    else
      figName = ['Neutral unit fraction in ' areaName ' v mean pupil area size: ' windowSize];
    end
  elseif strcmpi(type, 'positiveMUAs')
    pupilInd = ismember(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}, data.(['positiveMUAsThr' windowSize 'RecID']){1}{areaCode(1)});
    areaVec = data.(['positiveMUAsThr' windowSize]){1}{areaCode(1)}./data.(['totalMUAsThr' windowSize]){1}{areaCode(1)};
    pupilVec = data.(['meanPupilArea' windowSize]){1}{areaCode(1)}(pupilInd);
    areaVec = reduce2unique(data.(['positiveMUAsThr' windowSize 'RecID']){1}{areaCode(1)}, areaVec);
    pupilVec = reduce2unique(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}(pupilInd), pupilVec);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec', pupilVec', 'Spearman');
    fH(iCombo) = xyPlot(pupilVec, areaVec, 1, 'Mean pupil area size (a.u.)', ['Positive unit+MUA fraction in ' areaName]);
    if isempty(windowSize)
      figName = ['Positive unit+MUA fraction in ' areaName ' v mean pupil area size: Full recordings'];
    else
      figName = ['Positive unit+MUA fraction in ' areaName ' v mean pupil area size: ' windowSize];
    end
  elseif strcmpi(type, 'positiveSignificantMUAs')
    pupilInd = ismember(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}, data.(['positiveSignificantMUAsThr' windowSize 'RecID']){1}{areaCode(1)});
    areaVec = data.(['positiveSignificantMUAsThr' windowSize]){1}{areaCode(1)}./data.(['totalMUAsThr' windowSize]){1}{areaCode(1)};
    pupilVec = data.(['meanPupilArea' windowSize]){1}{areaCode(1)}(pupilInd);
    areaVec = reduce2unique(data.(['positiveSignificantMUAsThr' windowSize 'RecID']){1}{areaCode(1)}, areaVec);
    pupilVec = reduce2unique(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}(pupilInd), pupilVec);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec', pupilVec', 'Spearman');
    fH(iCombo) = xyPlot(pupilVec, areaVec, 1, 'Mean pupil area size (a.u.)', ['Positive Significant unit+MUA fraction in ' areaName]);
    if isempty(windowSize)
      figName = ['Positive Significant unit+MUA fraction in ' areaName ' v mean pupil area size: Full recordings'];
    else
      figName = ['Positive Significant unit+MUA fraction in ' areaName ' v mean pupil area size: ' windowSize];
    end
  elseif strcmpi(type, 'negativeMUAs')
    pupilInd = ismember(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}, data.(['negativeMUAsThr' windowSize 'RecID']){1}{areaCode(1)});
    areaVec = data.(['negativeMUAsThr' windowSize]){1}{areaCode(1)}./data.(['totalMUAsThr' windowSize]){1}{areaCode(1)};
    pupilVec = data.(['meanPupilArea' windowSize]){1}{areaCode(1)}(pupilInd);
    areaVec = reduce2unique(data.(['negativeMUAsThr' windowSize 'RecID']){1}{areaCode(1)}, areaVec);
    pupilVec = reduce2unique(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}(pupilInd), pupilVec);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec', pupilVec', 'Spearman');
    fH(iCombo) = xyPlot(pupilVec, areaVec, 1, 'Mean pupil area size (a.u.)', ['Negative unit+MUA fraction in ' areaName]);
    if isempty(windowSize)
      figName = ['Negative unit+MUA fraction in ' areaName ' v mean pupil area size: Full recordings'];
    else
      figName = ['Negative unit+MUA fraction in ' areaName ' v mean pupil area size: ' windowSize];
    end
  elseif strcmpi(type, 'negativeSignificantMUAs')
    pupilInd = ismember(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}, data.(['negativeSignificantMUAsThr' windowSize 'RecID']){1}{areaCode(1)});
    areaVec = data.(['negativeSignificantMUAsThr' windowSize]){1}{areaCode(1)}./data.(['totalMUAsThr' windowSize]){1}{areaCode(1)};
    pupilVec = data.(['meanPupilArea' windowSize]){1}{areaCode(1)}(pupilInd);
    areaVec = reduce2unique(data.(['negativeSignificantMUAsThr' windowSize 'RecID']){1}{areaCode(1)}, areaVec);
    pupilVec = reduce2unique(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}(pupilInd), pupilVec);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec', pupilVec', 'Spearman');
    fH(iCombo) = xyPlot(pupilVec, areaVec, 1, 'Mean pupil area size (a.u.)', ['Negative Significant unit+MUA fraction in ' areaName]);
    if isempty(windowSize)
      figName = ['Negative Significant unit+MUA fraction in ' areaName ' v mean pupil area size: Full recordings'];
    else
      figName = ['Negative Significant unit+MUA fraction in ' areaName ' v mean pupil area size: ' windowSize];
    end
  elseif strcmpi(type, 'neutralMUAs')
    pupilInd = ismember(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}, data.(['neutralMUAsThr' windowSize 'RecID']){1}{areaCode(1)});
    areaVec = data.(['neutralMUAsThr' windowSize]){1}{areaCode(1)}./data.(['totalMUAsThr' windowSize]){1}{areaCode(1)};
    pupilVec = data.(['meanPupilArea' windowSize]){1}{areaCode(1)}(pupilInd);
    areaVec = reduce2unique(data.(['neutralMUAsThr' windowSize 'RecID']){1}{areaCode(1)}, areaVec);
    pupilVec = reduce2unique(data.(['meanPupilArea' windowSize 'RecID']){1}{areaCode(1)}(pupilInd), pupilVec);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec', pupilVec', 'Spearman');
    fH(iCombo) = xyPlot(pupilVec, areaVec, 1, 'Mean pupil area size (a.u.)', ['Neutral unit+MUA fraction in ' areaName]);
    if isempty(windowSize)
      figName = ['Neutral unit+MUA fraction in ' areaName ' v mean pupil area size: Full recordings'];
    else
      figName = ['Neutral unit+MUA fraction in ' areaName ' v mean pupil area size: ' windowSize];
    end
  end
  
  xLim = xlim;
  xAxisLength = xLim(2)-xLim(1);
  yLim = ylim;
  yAxisLength = yLim(2)-yLim(1);
  textStr = ['r=' num2str(r(iCombo))];
  text(xLim(2)-xAxisLength*0.15, yLim(1)+yAxisLength*0.10, textStr, 'FontSize',20);
  textStr = ['p=' num2str(pval(iCombo))];
  text(xLim(2)-xAxisLength*0.15, yLim(1)+yAxisLength*0.05, textStr, 'FontSize',20);
  
  xLabel = get(gca, 'xlabel');
  yLabel = get(gca, 'ylabel');
  ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 20, 4/3, 2, [0.005 0], 'out',...
    'on', 'k', {xLabel.String}, xlim, xticks,...
    'on', 'k', {yLabel.String}, ylim, yticks);
  set(gcf,'color','w');
  
  %title(figName);
  set(fH(iCombo), 'Name',figName);
  if isempty(pupilInd)
    continue
  end
  figName = strrep(figName, ':', '_');
  filename = [outputDir filesep figName];
  filename = strrep(filename, ' ', '_');
  figSize = 15;
  label = [2.75 2.75];
  margin = [0.5 0.5];
  width = figSize-label(1)-margin(1);
  height = figSize-label(2)-margin(2);
  paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
  savefig(fH(iCombo), filename, 'compact');
  exportFig(fH(iCombo), [filename '.png'],'-dpng','-r300', paperSize);
end
end


function uniqueVec = reduce2unique(nameVec, dataVec)

[~, uniqueElements] = unique(nameVec, 'stable');
uniqueVec = zeros(numel(uniqueElements),1);
for i = 1:numel(uniqueElements)
  if uniqueElements(i)+1 <= size(nameVec,1) && strcmp(nameVec(uniqueElements(i)), nameVec(uniqueElements(i)+1))
    uniqueVec(i) = mean(dataVec(uniqueElements(i):uniqueElements(i)+1));
  else
    uniqueVec(i) = dataVec(uniqueElements(i));
  end
end
end


function data = assignValues2Windows(dbStruct, windowDuration, data, iCondPlusAll, iAreaPlusAll, seriesName, unitThr, alpha) %#ok<*INUSD,*DEFNU>

% units
type = [num2str(windowDuration) 'minWindows'];
shankIDs = fieldnames(dbStruct.shankData);
nUnitsWindows = zeros(numel(dbStruct.popData.(['meanPupilArea' type])),1);
for iWindow = 1:numel(dbStruct.popData.(['meanPupilArea' type]))
  rSpearmanWindow = [];
  pvalSpearmanWindow = [];
  for sh = 1:numel(shankIDs)
    if ~isempty(dbStruct.shankData.(shankIDs{sh}).(['rSpearman' type]))
      rSpearmanWindow = [rSpearmanWindow; dbStruct.shankData.(shankIDs{sh}).(['rSpearman' type]){iWindow}'];
      pvalSpearmanWindow = [pvalSpearmanWindow; dbStruct.shankData.(shankIDs{sh}).(['pvalSpearman' type]){iWindow}'];
    end
  end
  nUnitsWindows(iWindow) = numel(rSpearmanWindow);
  data.(['totalUnits' type]){iCondPlusAll}{iAreaPlusAll} =...
    [data.(['totalUnits' type]){iCondPlusAll}{iAreaPlusAll}; numel(rSpearmanWindow)];
  data.(['totalUnits' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
    [data.(['totalUnits' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
    [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
  data.(['positiveUnits' type]){iCondPlusAll}{iAreaPlusAll} =...
    [data.(['positiveUnits' type]){iCondPlusAll}{iAreaPlusAll}; sum(rSpearmanWindow >= 0)];
  data.(['positiveUnits' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
    [data.(['positiveUnits' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
    [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
%   data.(['positiveSignificantUnits' type]){iCondPlusAll}{iAreaPlusAll} =...
%     [data.(['positiveSignificantUnits' type]){iCondPlusAll}{iAreaPlusAll};...
%     sum(rSpearmanWindow >= 0 & pvalSpearmanWindow < alpha)];
%   data.(['positiveSignificantUnits' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
%     [data.(['positiveSignificantUnits' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
%     [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
  data.(['negativeUnits' type]){iCondPlusAll}{iAreaPlusAll} =...
    [data.(['negativeUnits' type]){iCondPlusAll}{iAreaPlusAll}; sum(rSpearmanWindow < 0)];
  data.(['negativeUnits' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
    [data.(['negativeUnits' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
    [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
%   data.(['negativeSignificantUnits' type]){iCondPlusAll}{iAreaPlusAll} =...
%     [data.(['negativeSignificantUnits' type]){iCondPlusAll}{iAreaPlusAll};...
%     sum(rSpearmanWindow < 0 & pvalSpearmanWindow < alpha)];
%   data.(['negativeSignificantUnits' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
%     [data.(['negativeSignificantUnits' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
%     [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
%   data.(['neutralUnits' type]){iCondPlusAll}{iAreaPlusAll} =...
%     [data.(['neutralUnits' type]){iCondPlusAll}{iAreaPlusAll};
%     sum(rSpearmanWindow < 0 & pvalSpearmanWindow >= alpha)];
%   data.(['neutralUnits' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
%     [data.(['neutralUnits' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
%     [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
  if nUnitsWindows(iWindow) >= unitThr
    data.(['totalUnitsThr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['totalUnitsThr' type]){iCondPlusAll}{iAreaPlusAll}; numel(rSpearmanWindow)];
    data.(['totalUnitsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['totalUnitsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    data.(['positiveUnitsThr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['positiveUnitsThr' type]){iCondPlusAll}{iAreaPlusAll}; sum(rSpearmanWindow >= 0)];
    data.(['positiveUnitsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['positiveUnitsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
%     data.(['positiveSignificantUnitsThr' type]){iCondPlusAll}{iAreaPlusAll} =...
%       [data.(['positiveSignificantUnitsThr' type]){iCondPlusAll}{iAreaPlusAll};...
%       sum(rSpearmanWindow >= 0 & pvalSpearmanWindow < alpha)];
%     data.(['positiveSignificantUnitsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
%       [data.(['positiveSignificantUnitsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
%       [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    data.(['negativeUnitsThr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['negativeUnitsThr' type]){iCondPlusAll}{iAreaPlusAll}; sum(rSpearmanWindow < 0)];
    data.(['negativeUnitsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['negativeUnitsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
%     data.(['negativeSignificantUnitsThr' type]){iCondPlusAll}{iAreaPlusAll} =...
%       [data.(['negativeSignificantUnitsThr' type]){iCondPlusAll}{iAreaPlusAll};...
%       sum(rSpearmanWindow < 0 & pvalSpearmanWindow < alpha)];
%     data.(['negativeSignificantUnitsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
%       [data.(['negativeSignificantUnitsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
%       [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
%     data.(['neutralUnitsThr' type]){iCondPlusAll}{iAreaPlusAll} =...
%       [data.(['neutralUnitsThr' type]){iCondPlusAll}{iAreaPlusAll};
%       sum(rSpearmanWindow < 0 & pvalSpearmanWindow >= alpha)];
%     data.(['neutralUnitsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
%       [data.(['neutralUnitsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
%       [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
  end
end

% MUAs
for iWindow = 1:numel(dbStruct.popData.(['meanPupilArea' type]))
  rSpearmanWindow = dbStruct.popData.(['rSpearman' type]){iWindow}';
  pvalSpearmanWindow = dbStruct.popData.(['pvalSpearman' type]){iWindow}'; %#ok<*NASGU>
  data.(['totalMUAs' type]){iCondPlusAll}{iAreaPlusAll} =...
    [data.(['totalMUAs' type]){iCondPlusAll}{iAreaPlusAll}; numel(rSpearmanWindow)];
  data.(['totalMUAs' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
    [data.(['totalMUAs' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
    [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
  data.(['positiveMUAs' type]){iCondPlusAll}{iAreaPlusAll} =...
    [data.(['positiveMUAs' type]){iCondPlusAll}{iAreaPlusAll}; sum(rSpearmanWindow >= 0)];
  data.(['positiveMUAs' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
    [data.(['positiveMUAs' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
    [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
%   data.(['positiveSignificantMUAs' type]){iCondPlusAll}{iAreaPlusAll} =...
%     [data.(['positiveSignificantMUAs' type]){iCondPlusAll}{iAreaPlusAll};...
%     sum(rSpearmanWindow >= 0 & pvalSpearmanWindow < alpha)];
%   data.(['positiveSignificantMUAs' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
%     [data.(['positiveSignificantMUAs' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
%     [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
  data.(['negativeMUAs' type]){iCondPlusAll}{iAreaPlusAll} =...
    [data.(['negativeMUAs' type]){iCondPlusAll}{iAreaPlusAll}; sum(rSpearmanWindow < 0)];
  data.(['negativeMUAs' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
    [data.(['negativeMUAs' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
    [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
%   data.(['negativeSignificantMUAs' type]){iCondPlusAll}{iAreaPlusAll} =...
%     [data.(['negativeSignificantMUAs' type]){iCondPlusAll}{iAreaPlusAll};...
%     sum(rSpearmanWindow < 0 & pvalSpearmanWindow < alpha)];
%   data.(['negativeSignificantMUAs' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
%     [data.(['negativeSignificantMUAs' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
%     [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
%   data.(['neutralMUAs' type]){iCondPlusAll}{iAreaPlusAll} =...
%     [data.(['neutralMUAs' type]){iCondPlusAll}{iAreaPlusAll};
%     sum(rSpearmanWindow < 0 & pvalSpearmanWindow >= alpha)];
%   data.(['neutralMUAs' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
%     [data.(['neutralMUAs' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
%     [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
  if nUnitsWindows(iWindow) >= unitThr
    data.(['totalMUAsThr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['totalMUAsThr' type]){iCondPlusAll}{iAreaPlusAll}; numel(rSpearmanWindow)];
    data.(['totalMUAsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['totalMUAsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    data.(['positiveMUAsThr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['positiveMUAsThr' type]){iCondPlusAll}{iAreaPlusAll}; sum(rSpearmanWindow >= 0)];
    data.(['positiveMUAsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['positiveMUAsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
%     data.(['positiveSignificantMUAsThr' type]){iCondPlusAll}{iAreaPlusAll} =...
%       [data.(['positiveSignificantMUAsThr' type]){iCondPlusAll}{iAreaPlusAll};...
%       sum(rSpearmanWindow >= 0 & pvalSpearmanWindow < alpha)];
%     data.(['positiveSignificantMUAsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
%       [data.(['positiveSignificantMUAsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
%       [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    data.(['negativeMUAsThr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['negativeMUAsThr' type]){iCondPlusAll}{iAreaPlusAll}; sum(rSpearmanWindow < 0)];
    data.(['negativeMUAsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['negativeMUAsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
%     data.(['negativeSignificantMUAsThr' type]){iCondPlusAll}{iAreaPlusAll} =...
%       [data.(['negativeSignificantMUAsThr' type]){iCondPlusAll}{iAreaPlusAll};...
%       sum(rSpearmanWindow < 0 & pvalSpearmanWindow < alpha)];
%     data.(['negativeSignificantMUAsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
%       [data.(['negativeSignificantMUAsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
%       [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
%     data.(['neutralMUAsThr' type]){iCondPlusAll}{iAreaPlusAll} =...
%       [data.(['neutralMUAsThr' type]){iCondPlusAll}{iAreaPlusAll};
%       sum(rSpearmanWindow < 0 & pvalSpearmanWindow >= alpha)];
%     data.(['neutralMUAsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
%       [data.(['neutralMUAsThr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
%       [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
  end
end
end


function comboEntries = findComboEntries(areaName1, areaName2, areas)

areaCombos = nchoosek(1:numel(areas),2);

area1rank = find(ismember(areas, areaName1));
area2rank = find(ismember(areas, areaName2));
area1col1 = areaCombos(:,1) == area1rank;
area2col1 = areaCombos(:,1) == area2rank;
areasCol1 = area1col1 | area2col1;
area1col2 = areaCombos(:,2) == area1rank;
area2col2 = areaCombos(:,2) == area2rank;
areasCol2 = area1col2 | area2col2;
comboEntries = areasCol1 & areasCol2;
end