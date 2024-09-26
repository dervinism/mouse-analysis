% This script produces positive/negative unit and MUA fraction plots
% comparing different brain areas on ultra-slow (US), infra-slow (IS), and
% slow (S) timescales. Nested configurations, where sub-fractions of
% positive IS cells within positive US (USIS) cells or vice-versa  are
% considered, are also compared. Files and figures with bar plots, violin
% plots, and correlation plots are produced and saved in dataDir/paDir_uol
% and dataDir/paDir_allensdk folders. Tables are produces to be copied to
% Excel sheets.

clearvars -except allData
params
lists

repository = 'uol';
unitThr = 10; % minimum unit or MUA count per data series
alpha = 0.05; % significance level
fullRun = 1; % 1 - all, 2 - bar and violin plots and correlation between fractions in different areas onwards, 3 - Correlation tables onwards
includeRuns = 'run';

dataDir = [dataDir filesep includeRuns];
if strcmp(repository,'uol')
  mainFolder = [dataDir filesep paDir_uol];
  animals = animalsUOLOI;
elseif strcmp(repository,'allensdk')
  mainFolder = [dataDir filesep paDir_allensdk];
  animals = animalsAllensdk;
  conditions = {'awake'};
end


%% Load data
if ~exist('allData', 'var')
  params
  if fullRun == 1
    if strcmp(repository,'all')
      load([dataDir filesep 'allData.mat']);
    elseif strcmp(repository,'uol')
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
      
      % Correlations
      data.rSpearmanUnitsS = counter;
      data.pvalSpearmanUnitsS = counter;
      data.rSpearmanUnitsIS = counter;
      data.pvalSpearmanUnitsIS = counter;
      data.rSpearmanUnitsUS = counter;
      data.pvalSpearmanUnitsUS = counter;
      
      data.rSpearmanMUAsS = counter;
      data.pvalSpearmanMUAsS = counter;
      data.rSpearmanMUAsIS = counter;
      data.pvalSpearmanMUAsIS = counter;
      data.rSpearmanMUAsUS = counter;
      data.pvalSpearmanMUAsUS = counter;
      
      data.rSpearmanUnitsSThr = counter;
      data.pvalSpearmanUnitsSThr = counter;
      data.rSpearmanUnitsISThr = counter;
      data.pvalSpearmanUnitsISThr = counter;
      data.rSpearmanUnitsUSThr = counter;
      data.pvalSpearmanUnitsUSThr = counter;
      
      data.rSpearmanMUAsSThr = counter;
      data.pvalSpearmanMUAsSThr = counter;
      data.rSpearmanMUAsISThr = counter;
      data.pvalSpearmanMUAsISThr = counter;
      data.rSpearmanMUAsUSThr = counter;
      data.pvalSpearmanMUAsUSThr = counter;
      
      data.rSpearmanUnitsSRecID = counterID;
      data.pvalSpearmanUnitsSRecID = counterID;
      data.rSpearmanUnitsISRecID = counterID;
      data.pvalSpearmanUnitsISRecID = counterID;
      data.rSpearmanUnitsUSRecID = counterID;
      data.pvalSpearmanUnitsUSRecID = counterID;
      
      data.rSpearmanMUAsSRecID = counterID;
      data.pvalSpearmanMUAsSRecID = counterID;
      data.rSpearmanMUAsISRecID = counterID;
      data.pvalSpearmanMUAsISRecID = counterID;
      data.rSpearmanMUAsUSRecID = counterID;
      data.pvalSpearmanMUAsUSRecID = counterID;
      
      data.rSpearmanUnitsSThrRecID = counterID;
      data.pvalSpearmanUnitsSThrRecID = counterID;
      data.rSpearmanUnitsISThrRecID = counterID;
      data.pvalSpearmanUnitsISThrRecID = counterID;
      data.rSpearmanUnitsUSThrRecID = counterID;
      data.pvalSpearmanUnitsUSThrRecID = counterID;
      
      data.rSpearmanMUAsSThrRecID = counterID;
      data.pvalSpearmanMUAsSThrRecID = counterID;
      data.rSpearmanMUAsISThrRecID = counterID;
      data.pvalSpearmanMUAsISThrRecID = counterID;
      data.rSpearmanMUAsUSThrRecID = counterID;
      data.pvalSpearmanMUAsUSThrRecID = counterID;
      
      % Full recordings
      data.positiveUnitsFrIS = counter;
      data.positiveUnitsFrUS = counter;
      data.negativeUnitsFrIS = counter;
      data.negativeUnitsFrUS = counter;
      data.positiveUnitsFrUSIS = counter;
      data.negativeUnitsFrUSIS = counter;
      
      data.positiveMUAsFrIS = counter;
      data.positiveMUAsFrUS = counter;
      data.negativeMUAsFrIS = counter;
      data.negativeMUAsFrUS = counter;
      data.positiveMUAsFrUSIS = counter;
      data.negativeMUAsFrUSIS = counter;
      
      data.positiveUnitsFrISThr = counter;
      data.positiveUnitsFrUSThr = counter;
      data.negativeUnitsFrISThr = counter;
      data.negativeUnitsFrUSThr = counter;
      data.positiveUnitsFrUSISThr = counter;
      data.negativeUnitsFrUSISThr = counter;
      
      data.positiveMUAsFrISThr = counter;
      data.positiveMUAsFrUSThr = counter;
      data.negativeMUAsFrISThr = counter;
      data.negativeMUAsFrUSThr = counter;
      data.positiveMUAsFrUSISThr = counter;
      data.negativeMUAsFrUSISThr = counter;
      
      data.positiveUnitsFrISRecID = counterID;
      data.positiveUnitsFrUSRecID = counterID;
      data.negativeUnitsFrISRecID = counterID;
      data.negativeUnitsFrUSRecID = counterID;
      data.positiveUnitsFrUSISRecID = counterID;
      data.negativeUnitsFrUSISRecID = counterID;
      
      data.positiveMUAsFrISRecID = counterID;
      data.positiveMUAsFrUSRecID = counterID;
      data.negativeMUAsFrISRecID = counterID;
      data.negativeMUAsFrUSRecID = counterID;
      data.positiveMUAsFrUSISRecID = counterID;
      data.negativeMUAsFrUSISRecID = counterID;
      
      data.positiveUnitsFrISThrRecID = counterID;
      data.positiveUnitsFrUSThrRecID = counterID;
      data.negativeUnitsFrISThrRecID = counterID;
      data.negativeUnitsFrUSThrRecID = counterID;
      data.positiveUnitsFrUSISThrRecID = counterID;
      data.negativeUnitsFrUSISThrRecID = counterID;
      
      data.positiveMUAsFrISThrRecID = counterID;
      data.positiveMUAsFrUSThrRecID = counterID;
      data.negativeMUAsFrISThrRecID = counterID;
      data.negativeMUAsFrUSThrRecID = counterID;
      data.positiveMUAsFrUSISThrRecID = counterID;
      data.negativeMUAsFrUSISThrRecID = counterID;
      
      data.positiveSignificantUnitsFrIS = counter;
      data.positiveSignificantUnitsFrUS = counter;
      data.negativeSignificantUnitsFrIS = counter;
      data.negativeSignificantUnitsFrUS = counter;
      data.positiveSignificantUnitsFrUSIS = counter;
      data.negativeSignificantUnitsFrUSIS = counter;
      
      data.positiveSignificantMUAsFrIS = counter;
      data.positiveSignificantMUAsFrUS = counter;
      data.negativeSignificantMUAsFrIS = counter;
      data.negativeSignificantMUAsFrUS = counter;
      data.positiveSignificantMUAsFrUSIS = counter;
      data.negativeSignificantMUAsFrUSIS = counter;
      
      data.positiveSignificantUnitsFrISThr = counter;
      data.positiveSignificantUnitsFrUSThr = counter;
      data.negativeSignificantUnitsFrISThr = counter;
      data.negativeSignificantUnitsFrUSThr = counter;
      data.positiveSignificantUnitsFrUSISThr = counter;
      data.negativeSignificantUnitsFrUSISThr = counter;
      
      data.positiveSignificantMUAsFrISThr = counter;
      data.positiveSignificantMUAsFrUSThr = counter;
      data.negativeSignificantMUAsFrISThr = counter;
      data.negativeSignificantMUAsFrUSThr = counter;
      data.positiveSignificantMUAsFrUSISThr = counter;
      data.negativeSignificantMUAsFrUSISThr = counter;
      
      data.positiveSignificantUnitsFrISRecID = counterID;
      data.positiveSignificantUnitsFrUSRecID = counterID;
      data.negativeSignificantUnitsFrISRecID = counterID;
      data.negativeSignificantUnitsFrUSRecID = counterID;
      data.positiveSignificantUnitsFrUSISRecID = counterID;
      data.negativeSignificantUnitsFrUSISRecID = counterID;
      
      data.positiveSignificantMUAsFrISRecID = counterID;
      data.positiveSignificantMUAsFrUSRecID = counterID;
      data.negativeSignificantMUAsFrISRecID = counterID;
      data.negativeSignificantMUAsFrUSRecID = counterID;
      data.positiveSignificantMUAsFrUSISRecID = counterID;
      data.negativeSignificantMUAsFrUSISRecID = counterID;
      
      data.positiveSignificantUnitsFrISThrRecID = counterID;
      data.positiveSignificantUnitsFrUSThrRecID = counterID;
      data.negativeSignificantUnitsFrISThrRecID = counterID;
      data.negativeSignificantUnitsFrUSThrRecID = counterID;
      data.positiveSignificantUnitsFrUSISThrRecID = counterID;
      data.negativeSignificantUnitsFrUSISThrRecID = counterID;
      
      data.positiveSignificantMUAsFrISThrRecID = counterID;
      data.positiveSignificantMUAsFrUSThrRecID = counterID;
      data.negativeSignificantMUAsFrISThrRecID = counterID;
      data.negativeSignificantMUAsFrUSThrRecID = counterID;
      data.positiveSignificantMUAsFrUSISThrRecID = counterID;
      data.negativeSignificantMUAsFrUSISThrRecID = counterID;
      
      % 10-minute windows
      data.positiveUnitsFrIS10minWindows = counter;
      data.positiveUnitsFrUS10minWindows = counter;
      data.negativeUnitsFrIS10minWindows = counter;
      data.negativeUnitsFrUS10minWindows = counter;
      data.positiveUnitsFrUSIS10minWindows = counter;
      data.negativeUnitsFrUSIS10minWindows = counter;
      
      data.positiveMUAsFrIS10minWindows = counter;
      data.positiveMUAsFrUS10minWindows = counter;
      data.negativeMUAsFrIS10minWindows = counter;
      data.negativeMUAsFrUS10minWindows = counter;
      data.positiveMUAsFrUSIS10minWindows = counter;
      data.negativeMUAsFrUSIS10minWindows = counter;
      
      data.positiveUnitsFrISThr10minWindows = counter;
      data.positiveUnitsFrUSThr10minWindows = counter;
      data.negativeUnitsFrISThr10minWindows = counter;
      data.negativeUnitsFrUSThr10minWindows = counter;
      data.positiveUnitsFrUSISThr10minWindows = counter;
      data.negativeUnitsFrUSISThr10minWindows = counter;
      
      data.positiveMUAsFrISThr10minWindows = counter;
      data.positiveMUAsFrUSThr10minWindows = counter;
      data.negativeMUAsFrISThr10minWindows = counter;
      data.negativeMUAsFrUSThr10minWindows = counter;
      data.positiveMUAsFrUSISThr10minWindows = counter;
      data.negativeMUAsFrUSISThr10minWindows = counter;
      
      data.positiveUnitsFrIS10minWindowsRecID = counterID;
      data.positiveUnitsFrUS10minWindowsRecID = counterID;
      data.negativeUnitsFrIS10minWindowsRecID = counterID;
      data.negativeUnitsFrUS10minWindowsRecID = counterID;
      data.positiveUnitsFrUSIS10minWindowsRecID = counterID;
      data.negativeUnitsFrUSIS10minWindowsRecID = counterID;
      
      data.positiveMUAsFrIS10minWindowsRecID = counterID;
      data.positiveMUAsFrUS10minWindowsRecID = counterID;
      data.negativeMUAsFrIS10minWindowsRecID = counterID;
      data.negativeMUAsFrUS10minWindowsRecID = counterID;
      data.positiveMUAsFrUSIS10minWindowsRecID = counterID;
      data.negativeMUAsFrUSIS10minWindowsRecID = counterID;
      
      data.positiveUnitsFrISThr10minWindowsRecID = counterID;
      data.positiveUnitsFrUSThr10minWindowsRecID = counterID;
      data.negativeUnitsFrISThr10minWindowsRecID = counterID;
      data.negativeUnitsFrUSThr10minWindowsRecID = counterID;
      data.positiveUnitsFrUSISThr10minWindowsRecID = counterID;
      data.negativeUnitsFrUSISThr10minWindowsRecID = counterID;
      
      data.positiveMUAsFrISThr10minWindowsRecID = counterID;
      data.positiveMUAsFrUSThr10minWindowsRecID = counterID;
      data.negativeMUAsFrISThr10minWindowsRecID = counterID;
      data.negativeMUAsFrUSThr10minWindowsRecID = counterID;
      data.positiveMUAsFrUSISThr10minWindowsRecID = counterID;
      data.negativeMUAsFrUSISThr10minWindowsRecID = counterID;
      
      data.positiveSignificantUnitsFrIS10minWindows = counter;
      data.positiveSignificantUnitsFrUS10minWindows = counter;
      data.negativeSignificantUnitsFrIS10minWindows = counter;
      data.negativeSignificantUnitsFrUS10minWindows = counter;
      data.positiveSignificantUnitsFrUSIS10minWindows = counter;
      data.negativeSignificantUnitsFrUSIS10minWindows = counter;
      
      data.positiveSignificantMUAsFrIS10minWindows = counter;
      data.positiveSignificantMUAsFrUS10minWindows = counter;
      data.negativeSignificantMUAsFrIS10minWindows = counter;
      data.negativeSignificantMUAsFrUS10minWindows = counter;
      data.positiveSignificantMUAsFrUSIS10minWindows = counter;
      data.negativeSignificantMUAsFrUSIS10minWindows = counter;
      
      data.positiveSignificantUnitsFrISThr10minWindows = counter;
      data.positiveSignificantUnitsFrUSThr10minWindows = counter;
      data.negativeSignificantUnitsFrISThr10minWindows = counter;
      data.negativeSignificantUnitsFrUSThr10minWindows = counter;
      data.positiveSignificantUnitsFrUSISThr10minWindows = counter;
      data.negativeSignificantUnitsFrUSISThr10minWindows = counter;
      
      data.positiveSignificantMUAsFrISThr10minWindows = counter;
      data.positiveSignificantMUAsFrUSThr10minWindows = counter;
      data.negativeSignificantMUAsFrISThr10minWindows = counter;
      data.negativeSignificantMUAsFrUSThr10minWindows = counter;
      data.positiveSignificantMUAsFrUSISThr10minWindows = counter;
      data.negativeSignificantMUAsFrUSISThr10minWindows = counter;
      
      data.positiveSignificantUnitsFrIS10minWindowsRecID = counterID;
      data.positiveSignificantUnitsFrUS10minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrIS10minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrUS10minWindowsRecID = counterID;
      data.positiveSignificantUnitsFrUSIS10minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrUSIS10minWindowsRecID = counterID;
      
      data.positiveSignificantMUAsFrIS10minWindowsRecID = counterID;
      data.positiveSignificantMUAsFrUS10minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrIS10minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrUS10minWindowsRecID = counterID;
      data.positiveSignificantMUAsFrUSIS10minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrUSIS10minWindowsRecID = counterID;
      
      data.positiveSignificantUnitsFrISThr10minWindowsRecID = counterID;
      data.positiveSignificantUnitsFrUSThr10minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrISThr10minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrUSThr10minWindowsRecID = counterID;
      data.positiveSignificantUnitsFrUSISThr10minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrUSISThr10minWindowsRecID = counterID;
      
      data.positiveSignificantMUAsFrISThr10minWindowsRecID = counterID;
      data.positiveSignificantMUAsFrUSThr10minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrISThr10minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrUSThr10minWindowsRecID = counterID;
      data.positiveSignificantMUAsFrUSISThr10minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrUSISThr10minWindowsRecID = counterID;
      
      % 20-minute windows
      data.positiveUnitsFrIS20minWindows = counter;
      data.positiveUnitsFrUS20minWindows = counter;
      data.negativeUnitsFrIS20minWindows = counter;
      data.negativeUnitsFrUS20minWindows = counter;
      data.positiveUnitsFrUSIS20minWindows = counter;
      data.negativeUnitsFrUSIS20minWindows = counter;
      
      data.positiveMUAsFrIS20minWindows = counter;
      data.positiveMUAsFrUS20minWindows = counter;
      data.negativeMUAsFrIS20minWindows = counter;
      data.negativeMUAsFrUS20minWindows = counter;
      data.positiveMUAsFrUSIS20minWindows = counter;
      data.negativeMUAsFrUSIS20minWindows = counter;
      
      data.positiveUnitsFrISThr20minWindows = counter;
      data.positiveUnitsFrUSThr20minWindows = counter;
      data.negativeUnitsFrISThr20minWindows = counter;
      data.negativeUnitsFrUSThr20minWindows = counter;
      data.positiveUnitsFrUSISThr20minWindows = counter;
      data.negativeUnitsFrUSISThr20minWindows = counter;
      
      data.positiveMUAsFrISThr20minWindows = counter;
      data.positiveMUAsFrUSThr20minWindows = counter;
      data.negativeMUAsFrISThr20minWindows = counter;
      data.negativeMUAsFrUSThr20minWindows = counter;
      data.positiveMUAsFrUSISThr20minWindows = counter;
      data.negativeMUAsFrUSISThr20minWindows = counter;
      
      data.positiveUnitsFrIS20minWindowsRecID = counterID;
      data.positiveUnitsFrUS20minWindowsRecID = counterID;
      data.negativeUnitsFrIS20minWindowsRecID = counterID;
      data.negativeUnitsFrUS20minWindowsRecID = counterID;
      data.positiveUnitsFrUSIS20minWindowsRecID = counterID;
      data.negativeUnitsFrUSIS20minWindowsRecID = counterID;
      
      data.positiveMUAsFrIS20minWindowsRecID = counterID;
      data.positiveMUAsFrUS20minWindowsRecID = counterID;
      data.negativeMUAsFrIS20minWindowsRecID = counterID;
      data.negativeMUAsFrUS20minWindowsRecID = counterID;
      data.positiveMUAsFrUSIS20minWindowsRecID = counterID;
      data.negativeMUAsFrUSIS20minWindowsRecID = counterID;
      
      data.positiveUnitsFrISThr20minWindowsRecID = counterID;
      data.positiveUnitsFrUSThr20minWindowsRecID = counterID;
      data.negativeUnitsFrISThr20minWindowsRecID = counterID;
      data.negativeUnitsFrUSThr20minWindowsRecID = counterID;
      data.positiveUnitsFrUSISThr20minWindowsRecID = counterID;
      data.negativeUnitsFrUSISThr20minWindowsRecID = counterID;
      
      data.positiveMUAsFrISThr20minWindowsRecID = counterID;
      data.positiveMUAsFrUSThr20minWindowsRecID = counterID;
      data.negativeMUAsFrISThr20minWindowsRecID = counterID;
      data.negativeMUAsFrUSThr20minWindowsRecID = counterID;
      data.positiveMUAsFrUSISThr20minWindowsRecID = counterID;
      data.negativeMUAsFrUSISThr20minWindowsRecID = counterID;
      
      data.positiveSignificantUnitsFrIS20minWindows = counter;
      data.positiveSignificantUnitsFrUS20minWindows = counter;
      data.negativeSignificantUnitsFrIS20minWindows = counter;
      data.negativeSignificantUnitsFrUS20minWindows = counter;
      data.positiveSignificantUnitsFrUSIS20minWindows = counter;
      data.negativeSignificantUnitsFrUSIS20minWindows = counter;
      
      data.positiveSignificantMUAsFrIS20minWindows = counter;
      data.positiveSignificantMUAsFrUS20minWindows = counter;
      data.negativeSignificantMUAsFrIS20minWindows = counter;
      data.negativeSignificantMUAsFrUS20minWindows = counter;
      data.positiveSignificantMUAsFrUSIS20minWindows = counter;
      data.negativeSignificantMUAsFrUSIS20minWindows = counter;
      
      data.positiveSignificantUnitsFrISThr20minWindows = counter;
      data.positiveSignificantUnitsFrUSThr20minWindows = counter;
      data.negativeSignificantUnitsFrISThr20minWindows = counter;
      data.negativeSignificantUnitsFrUSThr20minWindows = counter;
      data.positiveSignificantUnitsFrUSISThr20minWindows = counter;
      data.negativeSignificantUnitsFrUSISThr20minWindows = counter;
      
      data.positiveSignificantMUAsFrISThr20minWindows = counter;
      data.positiveSignificantMUAsFrUSThr20minWindows = counter;
      data.negativeSignificantMUAsFrISThr20minWindows = counter;
      data.negativeSignificantMUAsFrUSThr20minWindows = counter;
      data.positiveSignificantMUAsFrUSISThr20minWindows = counter;
      data.negativeSignificantMUAsFrUSISThr20minWindows = counter;
      
      data.positiveSignificantUnitsFrIS20minWindowsRecID = counterID;
      data.positiveSignificantUnitsFrUS20minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrIS20minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrUS20minWindowsRecID = counterID;
      data.positiveSignificantUnitsFrUSIS20minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrUSIS20minWindowsRecID = counterID;
      
      data.positiveSignificantMUAsFrIS20minWindowsRecID = counterID;
      data.positiveSignificantMUAsFrUS20minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrIS20minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrUS20minWindowsRecID = counterID;
      data.positiveSignificantMUAsFrUSIS20minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrUSIS20minWindowsRecID = counterID;
      
      data.positiveSignificantUnitsFrISThr20minWindowsRecID = counterID;
      data.positiveSignificantUnitsFrUSThr20minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrISThr20minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrUSThr20minWindowsRecID = counterID;
      data.positiveSignificantUnitsFrUSISThr20minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrUSISThr20minWindowsRecID = counterID;
      
      data.positiveSignificantMUAsFrISThr20minWindowsRecID = counterID;
      data.positiveSignificantMUAsFrUSThr20minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrISThr20minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrUSThr20minWindowsRecID = counterID;
      data.positiveSignificantMUAsFrUSISThr20minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrUSISThr20minWindowsRecID = counterID;
      
      % 30-minute windows
      data.positiveUnitsFrIS30minWindows = counter;
      data.positiveUnitsFrUS30minWindows = counter;
      data.negativeUnitsFrIS30minWindows = counter;
      data.negativeUnitsFrUS30minWindows = counter;
      data.positiveUnitsFrUSIS30minWindows = counter;
      data.negativeUnitsFrUSIS30minWindows = counter;
      
      data.positiveMUAsFrIS30minWindows = counter;
      data.positiveMUAsFrUS30minWindows = counter;
      data.negativeMUAsFrIS30minWindows = counter;
      data.negativeMUAsFrUS30minWindows = counter;
      data.positiveMUAsFrUSIS30minWindows = counter;
      data.negativeMUAsFrUSIS30minWindows = counter;
      
      data.positiveUnitsFrISThr30minWindows = counter;
      data.positiveUnitsFrUSThr30minWindows = counter;
      data.negativeUnitsFrISThr30minWindows = counter;
      data.negativeUnitsFrUSThr30minWindows = counter;
      data.positiveUnitsFrUSISThr30minWindows = counter;
      data.negativeUnitsFrUSISThr30minWindows = counter;
      
      data.positiveMUAsFrISThr30minWindows = counter;
      data.positiveMUAsFrUSThr30minWindows = counter;
      data.negativeMUAsFrISThr30minWindows = counter;
      data.negativeMUAsFrUSThr30minWindows = counter;
      data.positiveMUAsFrUSISThr30minWindows = counter;
      data.negativeMUAsFrUSISThr30minWindows = counter;
      
      data.positiveUnitsFrIS30minWindowsRecID = counterID;
      data.positiveUnitsFrUS30minWindowsRecID = counterID;
      data.negativeUnitsFrIS30minWindowsRecID = counterID;
      data.negativeUnitsFrUS30minWindowsRecID = counterID;
      data.positiveUnitsFrUSIS30minWindowsRecID = counterID;
      data.negativeUnitsFrUSIS30minWindowsRecID = counterID;
      
      data.positiveMUAsFrIS30minWindowsRecID = counterID;
      data.positiveMUAsFrUS30minWindowsRecID = counterID;
      data.negativeMUAsFrIS30minWindowsRecID = counterID;
      data.negativeMUAsFrUS30minWindowsRecID = counterID;
      data.positiveMUAsFrUSIS30minWindowsRecID = counterID;
      data.negativeMUAsFrUSIS30minWindowsRecID = counterID;
      
      data.positiveUnitsFrISThr30minWindowsRecID = counterID;
      data.positiveUnitsFrUSThr30minWindowsRecID = counterID;
      data.negativeUnitsFrISThr30minWindowsRecID = counterID;
      data.negativeUnitsFrUSThr30minWindowsRecID = counterID;
      data.positiveUnitsFrUSISThr30minWindowsRecID = counterID;
      data.negativeUnitsFrUSISThr30minWindowsRecID = counterID;
      
      data.positiveMUAsFrISThr30minWindowsRecID = counterID;
      data.positiveMUAsFrUSThr30minWindowsRecID = counterID;
      data.negativeMUAsFrISThr30minWindowsRecID = counterID;
      data.negativeMUAsFrUSThr30minWindowsRecID = counterID;
      data.positiveMUAsFrUSISThr30minWindowsRecID = counterID;
      data.negativeMUAsFrUSISThr30minWindowsRecID = counterID;
      
      data.positiveSignificantUnitsFrIS30minWindows = counter;
      data.positiveSignificantUnitsFrUS30minWindows = counter;
      data.negativeSignificantUnitsFrIS30minWindows = counter;
      data.negativeSignificantUnitsFrUS30minWindows = counter;
      data.positiveSignificantUnitsFrUSIS30minWindows = counter;
      data.negativeSignificantUnitsFrUSIS30minWindows = counter;
      
      data.positiveSignificantMUAsFrIS30minWindows = counter;
      data.positiveSignificantMUAsFrUS30minWindows = counter;
      data.negativeSignificantMUAsFrIS30minWindows = counter;
      data.negativeSignificantMUAsFrUS30minWindows = counter;
      data.positiveSignificantMUAsFrUSIS30minWindows = counter;
      data.negativeSignificantMUAsFrUSIS30minWindows = counter;
      
      data.positiveSignificantUnitsFrISThr30minWindows = counter;
      data.positiveSignificantUnitsFrUSThr30minWindows = counter;
      data.negativeSignificantUnitsFrISThr30minWindows = counter;
      data.negativeSignificantUnitsFrUSThr30minWindows = counter;
      data.positiveSignificantUnitsFrUSISThr30minWindows = counter;
      data.negativeSignificantUnitsFrUSISThr30minWindows = counter;
      
      data.positiveSignificantMUAsFrISThr30minWindows = counter;
      data.positiveSignificantMUAsFrUSThr30minWindows = counter;
      data.negativeSignificantMUAsFrISThr30minWindows = counter;
      data.negativeSignificantMUAsFrUSThr30minWindows = counter;
      data.positiveSignificantMUAsFrUSISThr30minWindows = counter;
      data.negativeSignificantMUAsFrUSISThr30minWindows = counter;
      
      data.positiveSignificantUnitsFrIS30minWindowsRecID = counterID;
      data.positiveSignificantUnitsFrUS30minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrIS30minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrUS30minWindowsRecID = counterID;
      data.positiveSignificantUnitsFrUSIS30minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrUSIS30minWindowsRecID = counterID;
      
      data.positiveSignificantMUAsFrIS30minWindowsRecID = counterID;
      data.positiveSignificantMUAsFrUS30minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrIS30minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrUS30minWindowsRecID = counterID;
      data.positiveSignificantMUAsFrUSIS30minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrUSIS30minWindowsRecID = counterID;
      
      data.positiveSignificantUnitsFrISThr30minWindowsRecID = counterID;
      data.positiveSignificantUnitsFrUSThr30minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrISThr30minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrUSThr30minWindowsRecID = counterID;
      data.positiveSignificantUnitsFrUSISThr30minWindowsRecID = counterID;
      data.negativeSignificantUnitsFrUSISThr30minWindowsRecID = counterID;
      
      data.positiveSignificantMUAsFrISThr30minWindowsRecID = counterID;
      data.positiveSignificantMUAsFrUSThr30minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrISThr30minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrUSThr30minWindowsRecID = counterID;
      data.positiveSignificantMUAsFrUSISThr30minWindowsRecID = counterID;
      data.negativeSignificantMUAsFrUSISThr30minWindowsRecID = counterID;
    end
    
    for dbCount = 1:numel(fnsData) % Loop through database entries
      dbStruct = dataStruct.seriesData.(fnsData{dbCount});
      seriesName = seriesFromEntry(fnsData{dbCount});
      
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
      
      for iAreaPlusAll = area % Loop through the main and pooled areas
        if numel(conditions) > 1
          condLoop = [iCond numel(conditions)];
        else
          condLoop = iCond;
        end
        for iCondPlusAll = condLoop % Loop through the main and pooled conditions
          
          % Units
          rSpearmanS = [];
          rSpearmanIS = [];
          rSpearmanUS = [];
          for sh = 1:numel(fieldnames(dbStruct.shankData))
            if ~isempty(dbStruct.shankData.(['shank' num2str(sh)]).spk)
              rSpearmanS = [rSpearmanS; torow(dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanSS)'];
              rSpearmanIS = [rSpearmanIS; torow(dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanISIS)'];
              rSpearmanUS = [rSpearmanUS; torow(dbStruct.shankData.(['shank' num2str(sh)]).rSpearmanUSUS)'];
            end
          end
          data.rSpearmanUnitsS{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanUnitsS{iCondPlusAll}{iAreaPlusAll}; torow(rSpearmanS)'];
          data.rSpearmanUnitsIS{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanUnitsIS{iCondPlusAll}{iAreaPlusAll}; torow(rSpearmanIS)'];
          data.rSpearmanUnitsUS{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanUnitsUS{iCondPlusAll}{iAreaPlusAll}; torow(rSpearmanUS)'];
          data.rSpearmanUnitsSRecID{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanUnitsSRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          data.rSpearmanUnitsISRecID{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanUnitsISRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          data.rSpearmanUnitsUSRecID{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanUnitsUSRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          if numel(rSpearmanS) >= unitThr
            data.rSpearmanUnitsSThr{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanUnitsSThr{iCondPlusAll}{iAreaPlusAll}; torow(rSpearmanS)'];
            data.rSpearmanUnitsISThr{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanUnitsISThr{iCondPlusAll}{iAreaPlusAll}; torow(rSpearmanIS)'];
            data.rSpearmanUnitsUSThr{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanUnitsUSThr{iCondPlusAll}{iAreaPlusAll}; torow(rSpearmanUS)'];
            for i = 1:numel(rSpearmanIS)
              data.rSpearmanUnitsSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanUnitsSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.rSpearmanUnitsISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanUnitsISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.rSpearmanUnitsUSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanUnitsUSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
            end
          end
          
          if ~isempty(rSpearmanIS) && ~isempty(rSpearmanUS)
            data.positiveUnitsFrIS{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsFrIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanIS >= 0)/numel(rSpearmanIS)];
            data.positiveUnitsFrISRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsFrISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.negativeUnitsFrIS{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsFrIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanIS < 0)/numel(rSpearmanIS)];
            data.negativeUnitsFrISRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsFrISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            if numel(rSpearmanIS) >= unitThr
              data.positiveUnitsFrISThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsFrISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanIS >= 0)/numel(rSpearmanIS)];
              data.positiveUnitsFrISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsFrISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.negativeUnitsFrISThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsFrISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanIS < 0)/numel(rSpearmanIS)];
              data.negativeUnitsFrISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsFrISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
            end
            data.positiveUnitsFrUS{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsFrUS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS >= 0)/numel(rSpearmanUS)];
            data.positiveUnitsFrUSRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsFrUSRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.negativeUnitsFrUS{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsFrUS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS < 0)/numel(rSpearmanUS)];
            data.negativeUnitsFrUSRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsFrUSRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.positiveUnitsFrUSIS{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsFrUSIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS >= 0 & rSpearmanIS >= 0)/sum(rSpearmanUS >= 0)];
            data.positiveUnitsFrUSISRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsFrUSISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.negativeUnitsFrUSIS{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsFrUSIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS < 0 & rSpearmanIS < 0)/sum(rSpearmanUS < 0)];
            data.negativeUnitsFrUSISRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsFrUSISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            if numel(rSpearmanUS) >= unitThr
              data.positiveUnitsFrUSThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsFrUSThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS >= 0)/numel(rSpearmanUS)];
              data.positiveUnitsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.negativeUnitsFrUSThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsFrUSThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS < 0)/numel(rSpearmanUS)];
              data.negativeUnitsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.positiveUnitsFrUSISThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsFrUSISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS >= 0 & rSpearmanIS >= 0)/sum(rSpearmanUS >= 0)];
              data.positiveUnitsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveUnitsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.negativeUnitsFrUSISThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsFrUSISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS < 0 & rSpearmanIS < 0)/sum(rSpearmanUS < 0)];
              data.negativeUnitsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeUnitsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
            end
          end
          
          % Significant units
          pvalSpearmanS = [];
          pvalSpearmanIS = [];
          pvalSpearmanUS = [];
          for sh = 1:numel(fieldnames(dbStruct.shankData))
            if ~isempty(dbStruct.shankData.(['shank' num2str(sh)]).spk)
              pvalSpearmanS = [pvalSpearmanS; torow(dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanSS)'];
              pvalSpearmanIS = [pvalSpearmanIS; torow(dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanISIS)'];
              pvalSpearmanUS = [pvalSpearmanUS; torow(dbStruct.shankData.(['shank' num2str(sh)]).pvalSpearmanUSUS)'];
            end
          end
          data.pvalSpearmanUnitsS{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanUnitsS{iCondPlusAll}{iAreaPlusAll}; torow(pvalSpearmanS)'];
          data.pvalSpearmanUnitsIS{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanUnitsIS{iCondPlusAll}{iAreaPlusAll}; torow(pvalSpearmanIS)'];
          data.pvalSpearmanUnitsUS{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanUnitsUS{iCondPlusAll}{iAreaPlusAll}; torow(pvalSpearmanUS)'];
          data.pvalSpearmanUnitsSRecID{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanUnitsSRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          data.pvalSpearmanUnitsISRecID{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanUnitsISRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          data.pvalSpearmanUnitsUSRecID{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanUnitsUSRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          if numel(pvalSpearmanS) >= unitThr
            data.pvalSpearmanUnitsSThr{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanUnitsSThr{iCondPlusAll}{iAreaPlusAll}; torow(pvalSpearmanS)'];
            data.pvalSpearmanUnitsISThr{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanUnitsISThr{iCondPlusAll}{iAreaPlusAll}; torow(pvalSpearmanIS)'];
            data.pvalSpearmanUnitsUSThr{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanUnitsUSThr{iCondPlusAll}{iAreaPlusAll}; torow(pvalSpearmanUS)'];
            for i = 1:numel(pvalSpearmanS)
              data.pvalSpearmanUnitsSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanUnitsSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.pvalSpearmanUnitsISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanUnitsISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.pvalSpearmanUnitsUSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanUnitsUSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
            end
          end
          
          if ~isempty(rSpearmanIS) && ~isempty(rSpearmanUS)
            data.positiveSignificantUnitsFrIS{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsFrIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanIS >= 0 & pvalSpearmanIS < alpha)/numel(rSpearmanIS)];
            data.positiveSignificantUnitsFrISRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsFrISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.negativeSignificantUnitsFrIS{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsFrIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanIS < 0 & pvalSpearmanIS < alpha)/numel(rSpearmanIS)];
            data.negativeSignificantUnitsFrISRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsFrISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            if numel(rSpearmanIS) >= unitThr
              data.positiveSignificantUnitsFrISThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsFrISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanIS >= 0 & pvalSpearmanIS < alpha)/numel(rSpearmanIS)];
              data.positiveSignificantUnitsFrISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsFrISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.negativeSignificantUnitsFrISThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsFrISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanIS < 0 & pvalSpearmanIS < alpha)/numel(rSpearmanIS)];
              data.negativeSignificantUnitsFrISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsFrISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
            end
            data.positiveSignificantUnitsFrUS{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsFrUS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS >= 0 & pvalSpearmanUS < alpha)/numel(rSpearmanUS)];
            data.positiveSignificantUnitsFrUSRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsFrUSRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.negativeSignificantUnitsFrUS{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsFrUS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS < 0 & pvalSpearmanUS < alpha)/numel(rSpearmanUS)];
            data.negativeSignificantUnitsFrUSRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsFrUSRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.positiveSignificantUnitsFrUSIS{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsFrUSIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS >= 0 & rSpearmanIS >= 0 & pvalSpearmanIS < alpha & pvalSpearmanUS < alpha)/sum(rSpearmanUS >= 0)];
            data.positiveSignificantUnitsFrUSISRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsFrUSISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.negativeSignificantUnitsFrUSIS{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsFrUSIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS < 0 & rSpearmanIS < 0 & pvalSpearmanIS < alpha & pvalSpearmanUS < alpha)/sum(rSpearmanUS < 0)];
            data.negativeSignificantUnitsFrUSISRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsFrUSISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            if numel(rSpearmanUS) >= unitThr
              data.positiveSignificantUnitsFrUSThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsFrUSThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS >= 0 & pvalSpearmanUS < alpha)/numel(rSpearmanUS)];
              data.positiveSignificantUnitsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.negativeSignificantUnitsFrUSThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsFrUSThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS < 0 & pvalSpearmanUS < alpha)/numel(rSpearmanUS)];
              data.negativeSignificantUnitsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.positiveSignificantUnitsFrUSISThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsFrUSISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS >= 0 & rSpearmanIS >= 0 & pvalSpearmanIS < alpha & pvalSpearmanUS < alpha)/sum(rSpearmanUS >= 0)];
              data.positiveSignificantUnitsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantUnitsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.negativeSignificantUnitsFrUSISThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsFrUSISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS < 0 & rSpearmanIS < 0 & pvalSpearmanIS < alpha & pvalSpearmanUS < alpha)/sum(rSpearmanUS < 0)];
              data.negativeSignificantUnitsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantUnitsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
            end
          end
          
          % MUAs
          rSpearmanS = [];
          rSpearmanIS = [];
          rSpearmanUS = [];
          if ~isempty(dbStruct.popData.spkDB)
            rSpearmanS = torow(dbStruct.popData.rSpearmanSS)';
            rSpearmanIS = torow(dbStruct.popData.rSpearmanISIS)';
            rSpearmanUS = torow(dbStruct.popData.rSpearmanUSUS)';
          end
          data.rSpearmanMUAsS{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanMUAsS{iCondPlusAll}{iAreaPlusAll}; torow(rSpearmanS)'];
          data.rSpearmanMUAsIS{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanMUAsIS{iCondPlusAll}{iAreaPlusAll}; torow(rSpearmanIS)'];
          data.rSpearmanMUAsUS{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanMUAsUS{iCondPlusAll}{iAreaPlusAll}; torow(rSpearmanUS)'];
          data.rSpearmanMUAsSRecID{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanMUAsSRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          data.rSpearmanMUAsISRecID{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanMUAsISRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          data.rSpearmanMUAsUSRecID{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanMUAsUSRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          if numel(rSpearmanS) >= unitThr
            data.rSpearmanMUAsSThr{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanMUAsSThr{iCondPlusAll}{iAreaPlusAll}; torow(rSpearmanS)'];
            data.rSpearmanMUAsISThr{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanMUAsISThr{iCondPlusAll}{iAreaPlusAll}; torow(rSpearmanIS)'];
            data.rSpearmanMUAsUSThr{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanMUAsUSThr{iCondPlusAll}{iAreaPlusAll}; torow(rSpearmanUS)'];
            for i = 1:numel(rSpearmanS)
              data.rSpearmanMUAsSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanMUAsSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.rSpearmanMUAsISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanMUAsISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.rSpearmanMUAsUSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.rSpearmanMUAsUSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
            end
          end
          
          if ~isempty(rSpearmanIS) && ~isempty(rSpearmanUS)
            data.positiveMUAsFrIS{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsFrIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanIS >= 0)/numel(rSpearmanIS)];
            data.positiveMUAsFrISRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsFrISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.negativeMUAsFrIS{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsFrIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanIS < 0)/numel(rSpearmanIS)];
            data.negativeMUAsFrISRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsFrISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            if numel(rSpearmanIS) >= unitThr
              data.positiveMUAsFrISThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsFrISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanIS >= 0)/numel(rSpearmanIS)];
              data.positiveMUAsFrISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsFrISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.negativeMUAsFrISThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsFrISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanIS < 0)/numel(rSpearmanIS)];
              data.negativeMUAsFrISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsFrISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
            end
            data.positiveMUAsFrUS{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsFrUS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS >= 0)/numel(rSpearmanUS)];
            data.positiveMUAsFrUSRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsFrUSRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.negativeMUAsFrUS{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsFrUS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS < 0)/numel(rSpearmanUS)];
            data.negativeMUAsFrUSRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsFrUSRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.positiveMUAsFrUSIS{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsFrUSIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS >= 0 & rSpearmanIS >= 0)/sum(rSpearmanUS >= 0)];
            data.positiveMUAsFrUSISRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsFrUSISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.negativeMUAsFrUSIS{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsFrUSIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS < 0 & rSpearmanIS < 0)/sum(rSpearmanUS < 0)];
            data.negativeMUAsFrUSISRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsFrUSISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            if numel(rSpearmanUS) >= unitThr
              data.positiveMUAsFrUSThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsFrUSThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS >= 0)/numel(rSpearmanUS)];
              data.positiveMUAsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.negativeMUAsFrUSThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsFrUSThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS < 0)/numel(rSpearmanUS)];
              data.negativeMUAsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.positiveMUAsFrUSISThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsFrUSISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS >= 0 & rSpearmanIS >= 0)/sum(rSpearmanUS >= 0)];
              data.positiveMUAsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveMUAsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.negativeMUAsFrUSISThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsFrUSISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS < 0 & rSpearmanIS < 0)/sum(rSpearmanUS < 0)];
              data.negativeMUAsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeMUAsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
            end
          end
          
          % Significant MUAs
          pvalSpearmanS = [];
          pvalSpearmanIS = [];
          pvalSpearmanUS = [];
          if ~isempty(dbStruct.popData.spkDB)
            pvalSpearmanS = torow(dbStruct.popData.pvalSpearmanSS)';
            pvalSpearmanIS = torow(dbStruct.popData.pvalSpearmanISIS)';
            pvalSpearmanUS = torow(dbStruct.popData.pvalSpearmanUSUS)';
          end
          data.pvalSpearmanMUAsS{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanMUAsS{iCondPlusAll}{iAreaPlusAll}; torow(pvalSpearmanS)'];
          data.pvalSpearmanMUAsIS{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanMUAsIS{iCondPlusAll}{iAreaPlusAll}; torow(pvalSpearmanIS)'];
          data.pvalSpearmanMUAsUS{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanMUAsUS{iCondPlusAll}{iAreaPlusAll}; torow(pvalSpearmanUS)'];
          data.pvalSpearmanMUAsSRecID{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanMUAsSRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          data.pvalSpearmanMUAsISRecID{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanMUAsISRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          data.pvalSpearmanMUAsUSRecID{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanMUAsUSRecID{iCondPlusAll}{iAreaPlusAll};...
            seriesName(1:min([numel(seriesName) 14]))];
          if numel(pvalSpearmanS) >= unitThr
            data.pvalSpearmanMUAsSThr{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanMUAsSThr{iCondPlusAll}{iAreaPlusAll}; torow(pvalSpearmanS)'];
            data.pvalSpearmanMUAsISThr{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanMUAsISThr{iCondPlusAll}{iAreaPlusAll}; torow(pvalSpearmanIS)'];
            data.pvalSpearmanMUAsUSThr{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanMUAsUSThr{iCondPlusAll}{iAreaPlusAll}; torow(pvalSpearmanUS)'];
            for i = 1:numel(pvalSpearmanIS)
              data.pvalSpearmanMUAsSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanMUAsSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.pvalSpearmanMUAsISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanMUAsISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.pvalSpearmanMUAsUSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.pvalSpearmanMUAsUSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
            end
          end
          
          if ~isempty(rSpearmanIS) && ~isempty(rSpearmanUS)
            data.positiveSignificantMUAsFrIS{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsFrIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanIS >= 0 & pvalSpearmanIS < alpha)/numel(rSpearmanIS)];
            data.positiveSignificantMUAsFrISRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsFrISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.negativeSignificantMUAsFrIS{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsFrIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanIS < 0 & pvalSpearmanIS < alpha)/numel(rSpearmanIS)];
            data.negativeSignificantMUAsFrISRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsFrISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            if numel(rSpearmanIS) >= unitThr
              data.positiveSignificantMUAsFrISThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsFrISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanIS >= 0 & pvalSpearmanIS < alpha)/numel(rSpearmanIS)];
              data.positiveSignificantMUAsFrISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsFrISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.negativeSignificantMUAsFrISThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsFrISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanIS < 0 & pvalSpearmanIS < alpha)/numel(rSpearmanIS)];
              data.negativeSignificantMUAsFrISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsFrISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
            end
            data.positiveSignificantMUAsFrUS{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsFrUS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS >= 0 & pvalSpearmanUS < alpha)/numel(rSpearmanUS)];
            data.positiveSignificantMUAsFrUSRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsFrUSRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.negativeSignificantMUAsFrUS{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsFrUS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS < 0 & pvalSpearmanUS < alpha)/numel(rSpearmanUS)];
            data.negativeSignificantMUAsFrUSRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsFrUSRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.positiveSignificantMUAsFrUSIS{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsFrUSIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS >= 0 & rSpearmanIS >= 0 & pvalSpearmanIS < alpha & pvalSpearmanUS < alpha)/sum(rSpearmanUS >= 0)];
            data.positiveSignificantMUAsFrUSISRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsFrUSISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            data.negativeSignificantMUAsFrUSIS{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsFrUSIS{iCondPlusAll}{iAreaPlusAll};...
              sum(rSpearmanUS < 0 & rSpearmanIS < 0 & pvalSpearmanIS < alpha & pvalSpearmanUS < alpha)/sum(rSpearmanUS < 0)];
            data.negativeSignificantMUAsFrUSISRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsFrUSISRecID{iCondPlusAll}{iAreaPlusAll};...
              seriesName(1:min([numel(seriesName) 14]))];
            if numel(rSpearmanUS) >= unitThr
              data.positiveSignificantMUAsFrUSThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsFrUSThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS >= 0 & pvalSpearmanUS < alpha)/numel(rSpearmanUS)];
              data.positiveSignificantMUAsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.negativeSignificantMUAsFrUSThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsFrUSThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS < 0 & pvalSpearmanUS < alpha)/numel(rSpearmanUS)];
              data.negativeSignificantMUAsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsFrUSThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.positiveSignificantMUAsFrUSISThr{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsFrUSISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS >= 0 & rSpearmanIS >= 0 & pvalSpearmanIS < alpha & pvalSpearmanUS < alpha)/sum(rSpearmanUS >= 0)];
              data.positiveSignificantMUAsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.positiveSignificantMUAsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
              data.negativeSignificantMUAsFrUSISThr{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsFrUSISThr{iCondPlusAll}{iAreaPlusAll};...
                sum(rSpearmanUS < 0 & rSpearmanIS < 0 & pvalSpearmanIS < alpha & pvalSpearmanUS < alpha)/sum(rSpearmanUS < 0)];
              data.negativeSignificantMUAsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll} = [data.negativeSignificantMUAsFrUSISThrRecID{iCondPlusAll}{iAreaPlusAll};...
                seriesName(1:min([numel(seriesName) 14]))];
            end
          end
          
          % Smaller windows
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
      [data.positiveUnitsFrISCI95{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanIS{iCond}{iArea}] = calc95CI(data.positiveUnitsFrIS{iCond}{iArea}); %#ok<*SAGROW>
      [data.negativeUnitsFrISCI95{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanIS{iCond}{iArea}] = calc95CI(data.negativeUnitsFrIS{iCond}{iArea});
      [data.positiveUnitsFrUSCI95{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUS{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUS{iCond}{iArea});
      [data.negativeUnitsFrUSCI95{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUS{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUS{iCond}{iArea});
      [data.positiveUnitsFrUSISCI95{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUSIS{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUSIS{iCond}{iArea});
      [data.negativeUnitsFrUSISCI95{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUSIS{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUSIS{iCond}{iArea});
      [data.positiveMUAsFrISCI95{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanIS{iCond}{iArea}] = calc95CI(data.positiveMUAsFrIS{iCond}{iArea});
      [data.negativeMUAsFrISCI95{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanIS{iCond}{iArea}] = calc95CI(data.negativeMUAsFrIS{iCond}{iArea});
      [data.positiveMUAsFrUSCI95{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUS{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUS{iCond}{iArea});
      [data.negativeMUAsFrUSCI95{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUS{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUS{iCond}{iArea});
      [data.positiveMUAsFrUSISCI95{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUSIS{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUSIS{iCond}{iArea});
      [data.negativeMUAsFrUSISCI95{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUSIS{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUSIS{iCond}{iArea});
      [data.positiveUnitsFrISCI95Thr{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanISThr{iCond}{iArea}] = calc95CI(data.positiveUnitsFrISThr{iCond}{iArea});
      [data.negativeUnitsFrISCI95Thr{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanISThr{iCond}{iArea}] = calc95CI(data.negativeUnitsFrISThr{iCond}{iArea});
      [data.positiveUnitsFrUSCI95Thr{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUSThr{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUSThr{iCond}{iArea});
      [data.negativeUnitsFrUSCI95Thr{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUSThr{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUSThr{iCond}{iArea});
      [data.positiveUnitsFrUSISCI95Thr{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUSISThr{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUSISThr{iCond}{iArea});
      [data.negativeUnitsFrUSISCI95Thr{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUSISThr{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUSISThr{iCond}{iArea});
      [data.positiveMUAsFrISCI95Thr{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanISThr{iCond}{iArea}] = calc95CI(data.positiveMUAsFrISThr{iCond}{iArea});
      [data.negativeMUAsFrISCI95Thr{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanISThr{iCond}{iArea}] = calc95CI(data.negativeMUAsFrISThr{iCond}{iArea});
      [data.positiveMUAsFrUSCI95Thr{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUSThr{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUSThr{iCond}{iArea});
      [data.negativeMUAsFrUSCI95Thr{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUSThr{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUSThr{iCond}{iArea});
      [data.positiveMUAsFrUSISCI95Thr{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUSISThr{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUSISThr{iCond}{iArea});
      [data.negativeMUAsFrUSISCI95Thr{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUSISThr{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUSISThr{iCond}{iArea});
      [data.positiveSignificantUnitsFrISCI95{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanIS{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrIS{iCond}{iArea});
      [data.negativeSignificantUnitsFrISCI95{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanIS{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrIS{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSCI95{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUS{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUS{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSCI95{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUS{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUS{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSISCI95{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUSIS{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUSIS{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSISCI95{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUSIS{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUSIS{iCond}{iArea});
      [data.positiveSignificantMUAsFrISCI95{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanIS{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrIS{iCond}{iArea});
      [data.negativeSignificantMUAsFrISCI95{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanIS{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrIS{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSCI95{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUS{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUS{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSCI95{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUS{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUS{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSISCI95{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUSIS{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUSIS{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSISCI95{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUSIS{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUSIS{iCond}{iArea});
      [data.positiveSignificantUnitsFrISCI95Thr{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanISThr{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrISThr{iCond}{iArea});
      [data.negativeSignificantUnitsFrISCI95Thr{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanISThr{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrISThr{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSCI95Thr{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUSThr{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUSThr{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSCI95Thr{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUSThr{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUSThr{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSISCI95Thr{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUSISThr{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUSISThr{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSISCI95Thr{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUSISThr{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUSISThr{iCond}{iArea});
      [data.positiveSignificantMUAsFrISCI95Thr{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanISThr{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrISThr{iCond}{iArea});
      [data.negativeSignificantMUAsFrISCI95Thr{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanISThr{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrISThr{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSCI95Thr{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUSThr{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUSThr{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSCI95Thr{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUSThr{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUSThr{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSISCI95Thr{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUSISThr{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUSISThr{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSISCI95Thr{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUSISThr{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUSISThr{iCond}{iArea});
      
      [data.positiveUnitsFrISCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanIS10minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrIS10minWindows{iCond}{iArea});
      [data.negativeUnitsFrISCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanIS10minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrIS10minWindows{iCond}{iArea});
      [data.positiveUnitsFrUSCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUS10minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUS10minWindows{iCond}{iArea});
      [data.negativeUnitsFrUSCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUS10minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUS10minWindows{iCond}{iArea});
      [data.positiveUnitsFrUSISCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUSIS10minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUSIS10minWindows{iCond}{iArea});
      [data.negativeUnitsFrUSISCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUSIS10minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUSIS10minWindows{iCond}{iArea});
      [data.positiveMUAsFrISCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanIS10minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrIS10minWindows{iCond}{iArea});
      [data.negativeMUAsFrISCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanIS10minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrIS10minWindows{iCond}{iArea});
      [data.positiveMUAsFrUSCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUS10minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUS10minWindows{iCond}{iArea});
      [data.negativeMUAsFrUSCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUS10minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUS10minWindows{iCond}{iArea});
      [data.positiveMUAsFrUSISCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUSIS10minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUSIS10minWindows{iCond}{iArea});
      [data.negativeMUAsFrUSISCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUSIS10minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUSIS10minWindows{iCond}{iArea});
      [data.positiveUnitsFrISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanISThr10minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrISThr10minWindows{iCond}{iArea});
      [data.negativeUnitsFrISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanISThr10minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrISThr10minWindows{iCond}{iArea});
      [data.positiveUnitsFrUSCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUSThr10minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUSThr10minWindows{iCond}{iArea});
      [data.negativeUnitsFrUSCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUSThr10minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUSThr10minWindows{iCond}{iArea});
      [data.positiveUnitsFrUSISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUSISThr10minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUSISThr10minWindows{iCond}{iArea});
      [data.negativeUnitsFrUSISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUSISThr10minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUSISThr10minWindows{iCond}{iArea});
      [data.positiveMUAsFrISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanISThr10minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrISThr10minWindows{iCond}{iArea});
      [data.negativeMUAsFrISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanISThr10minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrISThr10minWindows{iCond}{iArea});
      [data.positiveMUAsFrUSCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUSThr10minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUSThr10minWindows{iCond}{iArea});
      [data.negativeMUAsFrUSCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUSThr10minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUSThr10minWindows{iCond}{iArea});
      [data.positiveMUAsFrUSISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUSISThr10minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUSISThr10minWindows{iCond}{iArea});
      [data.negativeMUAsFrUSISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUSISThr10minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUSISThr10minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrISCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanIS10minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrIS10minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrISCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanIS10minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrIS10minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUS10minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUS10minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUS10minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUS10minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSISCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUSIS10minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUSIS10minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSISCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUSIS10minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUSIS10minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrISCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanIS10minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrIS10minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrISCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanIS10minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrIS10minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUS10minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUS10minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUS10minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUS10minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSISCI9510minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUSIS10minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUSIS10minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSISCI9510minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUSIS10minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUSIS10minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanISThr10minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrISThr10minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanISThr10minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrISThr10minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUSThr10minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUSThr10minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUSThr10minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUSThr10minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUSISThr10minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUSISThr10minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUSISThr10minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUSISThr10minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanISThr10minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrISThr10minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanISThr10minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrISThr10minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUSThr10minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUSThr10minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUSThr10minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUSThr10minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUSISThr10minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUSISThr10minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSISCI95Thr10minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUSISThr10minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUSISThr10minWindows{iCond}{iArea});
      
      [data.positiveUnitsFrISCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanIS20minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrIS20minWindows{iCond}{iArea});
      [data.negativeUnitsFrISCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanIS20minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrIS20minWindows{iCond}{iArea});
      [data.positiveUnitsFrUSCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUS20minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUS20minWindows{iCond}{iArea});
      [data.negativeUnitsFrUSCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUS20minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUS20minWindows{iCond}{iArea});
      [data.positiveUnitsFrUSISCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUSIS20minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUSIS20minWindows{iCond}{iArea});
      [data.negativeUnitsFrUSISCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUSIS20minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUSIS20minWindows{iCond}{iArea});
      [data.positiveMUAsFrISCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanIS20minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrIS20minWindows{iCond}{iArea});
      [data.negativeMUAsFrISCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanIS20minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrIS20minWindows{iCond}{iArea});
      [data.positiveMUAsFrUSCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUS20minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUS20minWindows{iCond}{iArea});
      [data.negativeMUAsFrUSCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUS20minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUS20minWindows{iCond}{iArea});
      [data.positiveMUAsFrUSISCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUSIS20minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUSIS20minWindows{iCond}{iArea});
      [data.negativeMUAsFrUSISCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUSIS20minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUSIS20minWindows{iCond}{iArea});
      [data.positiveUnitsFrISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanISThr20minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrISThr20minWindows{iCond}{iArea});
      [data.negativeUnitsFrISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanISThr20minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrISThr20minWindows{iCond}{iArea});
      [data.positiveUnitsFrUSCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUSThr20minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUSThr20minWindows{iCond}{iArea});
      [data.negativeUnitsFrUSCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUSThr20minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUSThr20minWindows{iCond}{iArea});
      [data.positiveUnitsFrUSISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUSISThr20minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUSISThr20minWindows{iCond}{iArea});
      [data.negativeUnitsFrUSISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUSISThr20minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUSISThr20minWindows{iCond}{iArea});
      [data.positiveMUAsFrISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanISThr20minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrISThr20minWindows{iCond}{iArea});
      [data.negativeMUAsFrISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanISThr20minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrISThr20minWindows{iCond}{iArea});
      [data.positiveMUAsFrUSCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUSThr20minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUSThr20minWindows{iCond}{iArea});
      [data.negativeMUAsFrUSCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUSThr20minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUSThr20minWindows{iCond}{iArea});
      [data.positiveMUAsFrUSISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUSISThr20minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUSISThr20minWindows{iCond}{iArea});
      [data.negativeMUAsFrUSISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUSISThr20minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUSISThr20minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrISCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanIS20minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrIS20minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrISCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanIS20minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrIS20minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUS20minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUS20minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUS20minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUS20minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSISCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUSIS20minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUSIS20minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSISCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUSIS20minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUSIS20minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrISCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanIS20minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrIS20minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrISCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanIS20minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrIS20minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUS20minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUS20minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUS20minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUS20minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSISCI9520minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUSIS20minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUSIS20minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSISCI9520minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUSIS20minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUSIS20minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanISThr20minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrISThr20minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanISThr20minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrISThr20minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUSThr20minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUSThr20minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUSThr20minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUSThr20minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUSISThr20minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUSISThr20minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUSISThr20minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUSISThr20minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanISThr20minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrISThr20minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanISThr20minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrISThr20minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUSThr20minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUSThr20minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUSThr20minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUSThr20minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUSISThr20minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUSISThr20minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSISCI95Thr20minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUSISThr20minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUSISThr20minWindows{iCond}{iArea});
      
      [data.positiveUnitsFrISCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanIS30minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrIS30minWindows{iCond}{iArea});
      [data.negativeUnitsFrISCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanIS30minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrIS30minWindows{iCond}{iArea});
      [data.positiveUnitsFrUSCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUS30minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUS30minWindows{iCond}{iArea});
      [data.negativeUnitsFrUSCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUS30minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUS30minWindows{iCond}{iArea});
      [data.positiveUnitsFrUSISCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUSIS30minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUSIS30minWindows{iCond}{iArea});
      [data.negativeUnitsFrUSISCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUSIS30minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUSIS30minWindows{iCond}{iArea});
      [data.positiveMUAsFrISCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanIS30minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrIS30minWindows{iCond}{iArea});
      [data.negativeMUAsFrISCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanIS30minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrIS30minWindows{iCond}{iArea});
      [data.positiveMUAsFrUSCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUS30minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUS30minWindows{iCond}{iArea});
      [data.negativeMUAsFrUSCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUS30minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUS30minWindows{iCond}{iArea});
      [data.positiveMUAsFrUSISCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUSIS30minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUSIS30minWindows{iCond}{iArea});
      [data.negativeMUAsFrUSISCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUSIS30minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUSIS30minWindows{iCond}{iArea});
      [data.positiveUnitsFrISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanISThr30minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrISThr30minWindows{iCond}{iArea});
      [data.negativeUnitsFrISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanISThr30minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrISThr30minWindows{iCond}{iArea});
      [data.positiveUnitsFrUSCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUSThr30minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUSThr30minWindows{iCond}{iArea});
      [data.negativeUnitsFrUSCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUSThr30minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUSThr30minWindows{iCond}{iArea});
      [data.positiveUnitsFrUSISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveUnitsFrMeanUSISThr30minWindows{iCond}{iArea}] = calc95CI(data.positiveUnitsFrUSISThr30minWindows{iCond}{iArea});
      [data.negativeUnitsFrUSISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeUnitsFrMeanUSISThr30minWindows{iCond}{iArea}] = calc95CI(data.negativeUnitsFrUSISThr30minWindows{iCond}{iArea});
      [data.positiveMUAsFrISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanISThr30minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrISThr30minWindows{iCond}{iArea});
      [data.negativeMUAsFrISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanISThr30minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrISThr30minWindows{iCond}{iArea});
      [data.positiveMUAsFrUSCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUSThr30minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUSThr30minWindows{iCond}{iArea});
      [data.negativeMUAsFrUSCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUSThr30minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUSThr30minWindows{iCond}{iArea});
      [data.positiveMUAsFrUSISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveMUAsFrMeanUSISThr30minWindows{iCond}{iArea}] = calc95CI(data.positiveMUAsFrUSISThr30minWindows{iCond}{iArea});
      [data.negativeMUAsFrUSISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeMUAsFrMeanUSISThr30minWindows{iCond}{iArea}] = calc95CI(data.negativeMUAsFrUSISThr30minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrISCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanIS30minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrIS30minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrISCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanIS30minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrIS30minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUS30minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUS30minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUS30minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUS30minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSISCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUSIS30minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUSIS30minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSISCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUSIS30minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUSIS30minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrISCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanIS30minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrIS30minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrISCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanIS30minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrIS30minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUS30minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUS30minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUS30minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUS30minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSISCI9530minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUSIS30minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUSIS30minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSISCI9530minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUSIS30minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUSIS30minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanISThr30minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrISThr30minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanISThr30minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrISThr30minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUSThr30minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUSThr30minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUSThr30minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUSThr30minWindows{iCond}{iArea});
      [data.positiveSignificantUnitsFrUSISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantUnitsFrMeanUSISThr30minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantUnitsFrUSISThr30minWindows{iCond}{iArea});
      [data.negativeSignificantUnitsFrUSISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantUnitsFrMeanUSISThr30minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantUnitsFrUSISThr30minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanISThr30minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrISThr30minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanISThr30minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrISThr30minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUSThr30minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUSThr30minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUSThr30minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUSThr30minWindows{iCond}{iArea});
      [data.positiveSignificantMUAsFrUSISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.positiveSignificantMUAsFrMeanUSISThr30minWindows{iCond}{iArea}] = calc95CI(data.positiveSignificantMUAsFrUSISThr30minWindows{iCond}{iArea});
      [data.negativeSignificantMUAsFrUSISCI95Thr30minWindows{iCond}{iArea}, ~, ~, data.negativeSignificantMUAsFrMeanUSISThr30minWindows{iCond}{iArea}] = calc95CI(data.negativeSignificantMUAsFrUSISThr30minWindows{iCond}{iArea});
    end
  end
  
  % Save processed data
  try
    save([mainFolder filesep 'pupilCorrFractionsFilter.mat'],...
      'unitThr','alpha','repository','animals','conditions','areas','data', '-v7.3');
  catch
    save([mainFolder filesep 'all' filesep 'pupilCorrFractionsFilter.mat'],...
      'unitThr','alpha','repository','animals','conditions','areas','data', '-v7.3');
  end
else
  try
    load([mainFolder filesep 'pupilCorrFractionsFilter.mat']); %#ok<*UNRCH>
  catch
    load([mainFolder filesep 'all' filesep 'pupilCorrFractionsFilter.mat']);
  end
end


if fullRun <= 2
  
  %% Areas of interest
  if strcmp(repository, 'uol')
    area1 = 'VB'; area2 = 'lS1'; area3 = 'lRSC'; area4 = 'CA';
  elseif strcmp(repository, 'allensdk')
    area1 = 'VB'; area2 = 'LGN'; area3 = 'V1'; area4 = 'CA';
  end
  areasOI = {area1; area2; area3; area4};
  
  
  %% Bar plots for units IS
  [fH, scatterGroups, areaGroups, ~, ~, corrGroups] = barPlotUnits(data, area1, area2, area3, area4, [], 'IS');
  set(fH, 'Name','Proportion of units positively and negatively correlated to pupil size on the infra-slow timescale in different brain areas');
  filename = [mainFolder filesep 'pupilCorrFractionsUnitsIS'];
  savefig(fH, filename, 'compact');
  print(fH, [filename '.png'],'-dpng','-r300');
  %close all
  
  
  %% ANOVA for units IS:
  runAnova(scatterGroups, corrGroups, areaGroups, filename)
  
  
  %% Bar plots for MUAs IS:
  [fH, scatterGroups, areaGroups, ~, ~, corrGroups] = barPlotMUAs(data, area1, area2, area3, area4, [], 'IS');
  set(fH, 'Name','Proportion of units and MUAs positively and negatively correlated to pupil size on the infra-slow timescale in different brain areas');
  filename = [mainFolder filesep 'pupilCorrFractionsMUAsIS'];
  savefig(fH, filename, 'compact');
  print(fH, [filename '.png'],'-dpng','-r300');
  %close all
  
  
  %% ANOVA for MUAs IS:
  runAnova(scatterGroups, corrGroups, areaGroups, filename)
  
  
  %% Bar plots for units US
  [fH, scatterGroups, areaGroups, ~, ~, corrGroups] = barPlotUnits(data, area1, area2, area3, area4, [], 'US');
  set(fH, 'Name','Proportion of units positively and negatively correlated to pupil size on the ultra-slow timescale in different brain areas');
  filename = [mainFolder filesep 'pupilCorrFractionsUnitsUS'];
  savefig(fH, filename, 'compact');
  print(fH, [filename '.png'],'-dpng','-r300');
  %close all
  
  
  %% ANOVA for units US:
  runAnova(scatterGroups, corrGroups, areaGroups, filename)
  
  
  %% Bar plots for MUAs US:
  [fH, scatterGroups, areaGroups, ~, ~, corrGroups] = barPlotMUAs(data, area1, area2, area3, area4, [], 'US');
  set(fH, 'Name','Proportion of units and MUAs positively and negatively correlated to pupil size on the ultra-slow timescale in different brain areas');
  filename = [mainFolder filesep 'pupilCorrFractionsMUAsUS'];
  savefig(fH, filename, 'compact');
  print(fH, [filename '.png'],'-dpng','-r300');
  %close all
  
  
  %% ANOVA for MUAs US:
  runAnova(scatterGroups, corrGroups, areaGroups, filename)
  
  
  %% Bar plots for units USIS
  [fH, scatterGroups, areaGroups, ~, ~, corrGroups] = barPlotUnits(data, area1, area2, area3, area4, [], 'USIS');
  set(fH, 'Name','Proportion of positive US_I_S/positive US and negative US_I_S/negative US units in different brain areas');
  filename = [mainFolder filesep 'pupilCorrFractionsUnitsUSIS'];
  savefig(fH, filename, 'compact');
  print(fH, [filename '.png'],'-dpng','-r300');
  %close all
  
  
  %% ANOVA for units USIS:
  runAnova(scatterGroups, corrGroups, areaGroups, filename)
  
  
  %% Bar plots for MUAs USIS:
  [fH, scatterGroups, areaGroups, ~, ~, corrGroups] = barPlotMUAs(data, area1, area2, area3, area4, [], 'USIS');
  set(fH, 'Name','Proportion of positive US_I_S/positive US and negative US_I_S/negative US units and MUAs in different brain areas');
  filename = [mainFolder filesep 'pupilCorrFractionsMUAsUSIS'];
  savefig(fH, filename, 'compact');
  print(fH, [filename '.png'],'-dpng','-r300');
  %close all
  
  
  %% ANOVA for MUAs USIS:
  runAnova(scatterGroups, corrGroups, areaGroups, filename)
  
  
  %% Significant fraction violin plots for negative units comparing wakefulness to anaesthesia: S
  if strcmp(repository, 'uol')
    area1 = 'Th'; area2 = 'VB'; area3 = 'lS1'; area4 = 'lRSC'; area5 = 'CA'; area6 = 'DG';
    areasOIFull = {area1; area2; area3; area4; area5; area6};
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      iCond = 1;
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
      [~, ~, ~, violinData{(iArea*2)-1}] = reduce2unique2(data.rSpearmanUnitsSThrRecID{iCond}{area},...
        data.rSpearmanUnitsSThr{iCond}{area}, data.pvalSpearmanUnitsSThr{iCond}{area});
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(torow(violinData{(iArea*2)-1})');
      iCond = 2;
      violinAreas{iArea*2} = [areasOIFull{iArea} 'Anaest'];
      [~, ~, ~, violinData{iArea*2}] = reduce2unique2(data.rSpearmanUnitsSThrRecID{iCond}{area},...
        data.rSpearmanUnitsSThr{iCond}{area}, data.pvalSpearmanUnitsSThr{iCond}{area});
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(torow(violinData{(iArea*2)})');
      [~, statsFractionsUnitsNegativeSSignificant{iArea}.p] = ttest2(torow(violinData{(iArea*2)-1})',...
        torow(violinData{(iArea*2)})');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 10.75]);
    end
    set(fH, 'Name','Proportion of slow significant negative units in different brain areas across conditions');
    filename = [mainFolder filesep 'pupilCorrFractionsSignificantNegativeUnitsSViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
  else
    statsFractionsUnitsNegativeSSignificant = [];
  end
  
  
  %% Significant fraction violin plots for negative MUAs comparing wakefulness to anaesthesia: S
  if strcmp(repository, 'uol')
    area1 = 'Th'; area2 = 'VB'; area3 = 'lS1'; area4 = 'lRSC'; area5 = 'CA'; area6 = 'DG';
    areasOIFull = {area1; area2; area3; area4; area5; area6};
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      iCond = 1;
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
      [~, ~, ~, violinData{(iArea*2)-1}] = reduce2unique2(data.rSpearmanMUAsSThrRecID{iCond}{area},...
        data.rSpearmanMUAsSThr{iCond}{area}, data.pvalSpearmanMUAsSThr{iCond}{area});
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(torow(violinData{(iArea*2)-1})');
      iCond = 2;
      violinAreas{iArea*2} = [areasOIFull{iArea} 'Anaest'];
      [~, ~, ~, violinData{iArea*2}] = reduce2unique2(data.rSpearmanMUAsSThrRecID{iCond}{area},...
        data.rSpearmanMUAsSThr{iCond}{area}, data.pvalSpearmanMUAsSThr{iCond}{area});
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(torow(violinData{(iArea*2)})');
      [~, statsFractionsMUAsNegativeSSignificant{iArea}.p] = ttest2(torow(violinData{(iArea*2)-1})',...
        torow(violinData{(iArea*2)})');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 10.75]);
    end
    set(fH, 'Name','Proportion of slow significant negative MUAs in different brain areas across conditions');
    filename = [mainFolder filesep 'pupilCorrFractionsSignificantNegativeMUAsSViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
  else
    statsFractionsMUAsNegativeSSignificant = [];
  end
  
  
  %% Significant fraction violin plots for negative units comparing wakefulness to anaesthesia: IS
  if strcmp(repository, 'uol')
    area1 = 'Th'; area2 = 'VB'; area3 = 'lS1'; area4 = 'lRSC'; area5 = 'CA'; area6 = 'DG';
    areasOIFull = {area1; area2; area3; area4; area5; area6};
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      iCond = 1;
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
      [~, ~, ~, violinData{(iArea*2)-1}] = reduce2unique2(data.rSpearmanUnitsISThrRecID{iCond}{area},...
        data.rSpearmanUnitsISThr{iCond}{area}, data.pvalSpearmanUnitsISThr{iCond}{area});
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(torow(violinData{(iArea*2)-1})');
      iCond = 2;
      violinAreas{iArea*2} = [areasOIFull{iArea} 'Anaest'];
      [~, ~, ~, violinData{iArea*2}] = reduce2unique2(data.rSpearmanUnitsISThrRecID{iCond}{area},...
        data.rSpearmanUnitsISThr{iCond}{area}, data.pvalSpearmanUnitsISThr{iCond}{area});
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(torow(violinData{(iArea*2)})');
      [~, statsFractionsUnitsNegativeISSignificant{iArea}.p] = ttest2(torow(violinData{(iArea*2)-1})',...
        torow(violinData{(iArea*2)})');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 10.75]);
    end
    set(fH, 'Name','Proportion of IS significant negative units in different brain areas across conditions');
    filename = [mainFolder filesep 'pupilCorrFractionsSignificantNegativeUnitsISViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
  else
    statsFractionsUnitsNegativeISSignificant = [];
  end
  
  
  %% Significant fraction violin plots for negative MUAs comparing wakefulness to anaesthesia: IS
  if strcmp(repository, 'uol')
    area1 = 'Th'; area2 = 'VB'; area3 = 'lS1'; area4 = 'lRSC'; area5 = 'CA'; area6 = 'DG';
    areasOIFull = {area1; area2; area3; area4; area5; area6};
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      iCond = 1;
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
      [~, ~, ~, violinData{(iArea*2)-1}] = reduce2unique2(data.rSpearmanMUAsISThrRecID{iCond}{area},...
        data.rSpearmanMUAsISThr{iCond}{area}, data.pvalSpearmanMUAsISThr{iCond}{area});
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(torow(violinData{(iArea*2)-1})');
      iCond = 2;
      violinAreas{iArea*2} = [areasOIFull{iArea} 'Anaest'];
      [~, ~, ~, violinData{iArea*2}] = reduce2unique2(data.rSpearmanMUAsISThrRecID{iCond}{area},...
        data.rSpearmanMUAsISThr{iCond}{area}, data.pvalSpearmanMUAsISThr{iCond}{area});
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(torow(violinData{(iArea*2)})');
      [~, statsFractionsMUAsNegativeISSignificant{iArea}.p] = ttest2(torow(violinData{(iArea*2)-1})',...
        torow(violinData{(iArea*2)})');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 10.75]);
    end
    set(fH, 'Name','Proportion of IS significant negative MUAs in different brain areas across conditions');
    filename = [mainFolder filesep 'pupilCorrFractionsSignificantNegativeMUAsISViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
  else
    statsFractionsMUAsNegativeISSignificant = [];
  end
  
  
  %% Significant fraction violin plots for negative units comparing wakefulness to anaesthesia: US
  if strcmp(repository, 'uol')
    area1 = 'Th'; area2 = 'VB'; area3 = 'lS1'; area4 = 'lRSC'; area5 = 'CA'; area6 = 'DG';
    areasOIFull = {area1; area2; area3; area4; area5; area6};
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      iCond = 1;
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
      [~, ~, ~, violinData{(iArea*2)-1}] = reduce2unique2(data.rSpearmanUnitsUSThrRecID{iCond}{area},...
        data.rSpearmanUnitsUSThr{iCond}{area}, data.pvalSpearmanUnitsUSThr{iCond}{area});
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(torow(violinData{(iArea*2)-1})');
      iCond = 2;
      violinAreas{iArea*2} = [areasOIFull{iArea} 'Anaest'];
      [~, ~, ~, violinData{iArea*2}] = reduce2unique2(data.rSpearmanUnitsUSThrRecID{iCond}{area},...
        data.rSpearmanUnitsUSThr{iCond}{area}, data.pvalSpearmanUnitsUSThr{iCond}{area});
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(torow(violinData{(iArea*2)})');
      [~, statsFractionsUnitsNegativeUSSignificant{iArea}.p] = ttest2(torow(violinData{(iArea*2)-1})',...
        torow(violinData{(iArea*2)})');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 10.75]);
    end
    set(fH, 'Name','Proportion of US significant negative units in different brain areas across conditions');
    filename = [mainFolder filesep 'pupilCorrFractionsSignificantNegativeUnitsUSViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
  else
    statsFractionsUnitsNegativeUSSignificant = [];
  end
  
  
  %% Significant fraction violin plots for negative MUAs comparing wakefulness to anaesthesia: US
  if strcmp(repository, 'uol')
    area1 = 'Th'; area2 = 'VB'; area3 = 'lS1'; area4 = 'lRSC'; area5 = 'CA'; area6 = 'DG';
    areasOIFull = {area1; area2; area3; area4; area5; area6};
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      iCond = 1;
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
      [~, ~, ~, violinData{(iArea*2)-1}] = reduce2unique2(data.rSpearmanMUAsUSThrRecID{iCond}{area},...
        data.rSpearmanMUAsUSThr{iCond}{area}, data.pvalSpearmanMUAsUSThr{iCond}{area});
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(torow(violinData{(iArea*2)-1})');
      iCond = 2;
      violinAreas{iArea*2} = [areasOIFull{iArea} 'Anaest'];
      [~, ~, ~, violinData{iArea*2}] = reduce2unique2(data.rSpearmanMUAsUSThrRecID{iCond}{area},...
        data.rSpearmanMUAsUSThr{iCond}{area}, data.pvalSpearmanMUAsUSThr{iCond}{area});
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(torow(violinData{(iArea*2)})');
      [~, statsFractionsMUAsNegativeUSSignificant{iArea}.p] = ttest2(torow(violinData{(iArea*2)-1})',...
        torow(violinData{(iArea*2)})');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 10.75]);
    end
    set(fH, 'Name','Proportion of US significant negative MUAs in different brain areas across conditions');
    filename = [mainFolder filesep 'pupilCorrFractionsSignificantNegativeMUAsUSViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
  else
    statsFractionsMUAsNegativeUSSignificant = [];
  end
  
  
  %% Total fraction violin plots for negative units comparing wakefulness to anaesthesia: S
  if strcmp(repository, 'uol')
    area1 = 'Th'; area2 = 'VB'; area3 = 'lS1'; area4 = 'lRSC'; area5 = 'CA'; area6 = 'DG';
    areasOIFull = {area1; area2; area3; area4; area5; area6};
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      iCond = 1;
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
      [violinData{(iArea*2)-1}, ~, ~, ~] = reduce2unique2(data.rSpearmanUnitsSThrRecID{iCond}{area},...
        data.rSpearmanUnitsSThr{iCond}{area}, data.pvalSpearmanUnitsSThr{iCond}{area});
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(torow(violinData{(iArea*2)-1})');
      iCond = 2;
      violinAreas{iArea*2} = [areasOIFull{iArea} 'Anaest'];
      [violinData{iArea*2}, ~, ~, ~] = reduce2unique2(data.rSpearmanUnitsSThrRecID{iCond}{area},...
        data.rSpearmanUnitsSThr{iCond}{area}, data.pvalSpearmanUnitsSThr{iCond}{area});
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(torow(violinData{(iArea*2)})');
      [~, statsFractionsUnitsNegativeS{iArea}.p] = ttest2(torow(violinData{(iArea*2)-1})',...
        torow(violinData{(iArea*2)})');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 10.75]);
    end
    set(fH, 'Name','Proportion of slow positive units in different brain areas across conditions');
    filename = [mainFolder filesep 'pupilCorrFractionsPositiveUnitsSViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
  else
    statsFractionsUnitsNegativeS = [];
  end
  
  
  %% Total fraction violin plots for negative MUAs comparing wakefulness to anaesthesia: S
  if strcmp(repository, 'uol')
    area1 = 'Th'; area2 = 'VB'; area3 = 'lS1'; area4 = 'lRSC'; area5 = 'CA'; area6 = 'DG';
    areasOIFull = {area1; area2; area3; area4; area5; area6};
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      iCond = 1;
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
      [violinData{(iArea*2)-1}, ~, ~, ~] = reduce2unique2(data.rSpearmanMUAsSThrRecID{iCond}{area},...
        data.rSpearmanMUAsSThr{iCond}{area}, data.pvalSpearmanMUAsSThr{iCond}{area});
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(torow(violinData{(iArea*2)-1})');
      iCond = 2;
      violinAreas{iArea*2} = [areasOIFull{iArea} 'Anaest'];
      [violinData{iArea*2}, ~, ~, ~] = reduce2unique2(data.rSpearmanMUAsSThrRecID{iCond}{area},...
        data.rSpearmanMUAsSThr{iCond}{area}, data.pvalSpearmanMUAsSThr{iCond}{area});
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(torow(violinData{(iArea*2)})');
      [~, statsFractionsMUAsNegativeS{iArea}.p] = ttest2(torow(violinData{(iArea*2)-1})',...
        torow(violinData{(iArea*2)})');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 10.75]);
    end
    set(fH, 'Name','Proportion of slow positive MUAs in different brain areas across conditions');
    filename = [mainFolder filesep 'pupilCorrFractionsPositiveMUAsSViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
  else
    statsFractionsMUAsNegativeS = [];
  end
  
  
  %% Total fraction violin plots for negative units comparing wakefulness to anaesthesia: IS
  if strcmp(repository, 'uol')
    area1 = 'Th'; area2 = 'VB'; area3 = 'lS1'; area4 = 'lRSC'; area5 = 'CA'; area6 = 'DG';
    areasOIFull = {area1; area2; area3; area4; area5; area6};
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      iCond = 1;
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
      [violinData{(iArea*2)-1}, ~, ~, ~] = reduce2unique2(data.rSpearmanUnitsISThrRecID{iCond}{area},...
        data.rSpearmanUnitsISThr{iCond}{area}, data.pvalSpearmanUnitsISThr{iCond}{area});
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(torow(violinData{(iArea*2)-1})');
      iCond = 2;
      violinAreas{iArea*2} = [areasOIFull{iArea} 'Anaest'];
      [violinData{iArea*2}, ~, ~, ~] = reduce2unique2(data.rSpearmanUnitsISThrRecID{iCond}{area},...
        data.rSpearmanUnitsISThr{iCond}{area}, data.pvalSpearmanUnitsISThr{iCond}{area});
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(torow(violinData{(iArea*2)})');
      [~, statsFractionsUnitsNegativeIS{iArea}.p] = ttest2(torow(violinData{(iArea*2)-1})',...
        torow(violinData{(iArea*2)})');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 10.75]);
    end
    set(fH, 'Name','Proportion of IS positive units in different brain areas across conditions');
    filename = [mainFolder filesep 'pupilCorrFractionsPositiveUnitsISViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
  else
    statsFractionsUnitsNegativeIS = [];
  end
  
  
  %% Total fraction violin plots for negative MUAs comparing wakefulness to anaesthesia: IS
  if strcmp(repository, 'uol')
    area1 = 'Th'; area2 = 'VB'; area3 = 'lS1'; area4 = 'lRSC'; area5 = 'CA'; area6 = 'DG';
    areasOIFull = {area1; area2; area3; area4; area5; area6};
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      iCond = 1;
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
      [violinData{(iArea*2)-1}, ~, ~, ~] = reduce2unique2(data.rSpearmanMUAsISThrRecID{iCond}{area},...
        data.rSpearmanMUAsISThr{iCond}{area}, data.pvalSpearmanMUAsISThr{iCond}{area});
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(torow(violinData{(iArea*2)-1})');
      iCond = 2;
      violinAreas{iArea*2} = [areasOIFull{iArea} 'Anaest'];
      [violinData{iArea*2}, ~, ~, ~] = reduce2unique2(data.rSpearmanMUAsISThrRecID{iCond}{area},...
        data.rSpearmanMUAsISThr{iCond}{area}, data.pvalSpearmanMUAsISThr{iCond}{area});
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(torow(violinData{(iArea*2)})');
      [~, statsFractionsMUAsNegativeIS{iArea}.p] = ttest2(torow(violinData{(iArea*2)-1})',...
        torow(violinData{(iArea*2)})');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 10.75]);
    end
    set(fH, 'Name','Proportion of IS positive MUAs in different brain areas across conditions');
    filename = [mainFolder filesep 'pupilCorrFractionsPositiveMUAsISViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
  else
    statsFractionsMUAsNegativeIS = [];
  end
  
  
  %% Total fraction violin plots for negative units comparing wakefulness to anaesthesia: US
  if strcmp(repository, 'uol')
    area1 = 'Th'; area2 = 'VB'; area3 = 'lS1'; area4 = 'lRSC'; area5 = 'CA'; area6 = 'DG';
    areasOIFull = {area1; area2; area3; area4; area5; area6};
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      iCond = 1;
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
      [violinData{(iArea*2)-1}, ~, ~, ~] = reduce2unique2(data.rSpearmanUnitsUSThrRecID{iCond}{area},...
        data.rSpearmanUnitsUSThr{iCond}{area}, data.pvalSpearmanUnitsUSThr{iCond}{area});
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(torow(violinData{(iArea*2)-1})');
      iCond = 2;
      violinAreas{iArea*2} = [areasOIFull{iArea} 'Anaest'];
      [violinData{iArea*2}, ~, ~, ~] = reduce2unique2(data.rSpearmanUnitsUSThrRecID{iCond}{area},...
        data.rSpearmanUnitsUSThr{iCond}{area}, data.pvalSpearmanUnitsUSThr{iCond}{area});
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(torow(violinData{(iArea*2)})');
      [~, statsFractionsUnitsNegativeUS{iArea}.p] = ttest2(torow(violinData{(iArea*2)-1})',...
        torow(violinData{(iArea*2)})');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 10.75]);
    end
    set(fH, 'Name','Proportion of US positive units in different brain areas across conditions');
    filename = [mainFolder filesep 'pupilCorrFractionsPositiveUnitsUSViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
  else
    statsFractionsUnitsNegativeUS = [];
  end
  
  
  %% Total fraction violin plots for negative MUAs comparing wakefulness to anaesthesia: US
  if strcmp(repository, 'uol')
    area1 = 'Th'; area2 = 'VB'; area3 = 'lS1'; area4 = 'lRSC'; area5 = 'CA'; area6 = 'DG';
    areasOIFull = {area1; area2; area3; area4; area5; area6};
    dataMean = zeros(1, numel(areasOIFull)*2);
    dataCI95 = zeros(2, numel(areasOIFull)*2);
    for iArea = 1:numel(areasOIFull)
      area = determineArea(areasOIFull{iArea});
      iCond = 1;
      violinAreas{(iArea*2)-1} = [areasOIFull{iArea} 'Awake'];
      [violinData{(iArea*2)-1}, ~, ~, ~] = reduce2unique2(data.rSpearmanMUAsUSThrRecID{iCond}{area},...
        data.rSpearmanMUAsUSThr{iCond}{area}, data.pvalSpearmanMUAsUSThr{iCond}{area});
      [dataMean((iArea*2)-1), dataCI95(:,(iArea*2)-1)] = datamean(torow(violinData{(iArea*2)-1})');
      iCond = 2;
      violinAreas{iArea*2} = [areasOIFull{iArea} 'Anaest'];
      [violinData{iArea*2}, ~, ~, ~] = reduce2unique2(data.rSpearmanMUAsUSThrRecID{iCond}{area},...
        data.rSpearmanMUAsUSThr{iCond}{area}, data.pvalSpearmanMUAsUSThr{iCond}{area});
      [dataMean((iArea*2)), dataCI95(:,(iArea*2))] = datamean(torow(violinData{(iArea*2)})');
      [~, statsFractionsMUAsNegativeUS{iArea}.p] = ttest2(torow(violinData{(iArea*2)-1})',...
        torow(violinData{(iArea*2)})');
    end
    options.yLim = [0 1];
    options.yLabel = 'Fraction';
    options.showNotches = false;
    options.medianPlot = false;
    fH = multiViolinPlots(violinData, violinAreas, dataMean, dataCI95, [], options);
    if strcmp(repository, 'uol')
      xlim([0.25 12.75]);
    elseif strcmp(repository, 'allensdk')
      xlim([0.25 10.75]);
    end
    set(fH, 'Name','Proportion of US positive MUAs in different brain areas across conditions');
    filename = [mainFolder filesep 'pupilCorrFractionsPositiveMUAsUSViolins'];
    savefig(fH, filename, 'compact');
    print(fH, [filename '.png'],'-dpng','-r300');
  else
    statsFractionsMUAsNegativeUS = [];
  end
  
  
  %% Correlation graphs
  
  % Full recordings
  windowSize = [];
  % IS
  [~, rPositiveUnitsIS, pvalPositiveUnitsIS] = corrPlot(areasOI, data, 'positiveUnitsIS', windowSize, mainFolder);
  [~, rPositiveSignificantUnitsIS, pvalPositiveSignificantUnitsIS] = corrPlot(areasOI, data, 'positiveSignificantUnitsIS', windowSize, mainFolder);
  [~, rNegativeUnitsIS, pvalNegativeUnitsIS] = corrPlot(areasOI, data, 'negativeUnitsIS', windowSize, mainFolder);
  [~, rNegativeSignificantUnitsIS, pvalNegativeSignificantUnitsIS] = corrPlot(areasOI, data, 'negativeSignificantUnitsIS', windowSize, mainFolder);
  [~, rPositiveMUAsIS, pvalPositiveMUAsIS] = corrPlot(areasOI, data, 'positiveMUAsIS', windowSize, mainFolder);
  [~, rPositiveSignificantMUAsIS, pvalPositiveSignificantMUAsIS] = corrPlot(areasOI, data, 'positiveSignificantMUAsIS', windowSize, mainFolder);
  [~, rNegativeMUAsIS, pvalNegativeMUAsIS] = corrPlot(areasOI, data, 'negativeMUAsIS', windowSize, mainFolder);
  [~, rNegativeSignificantMUAsIS, pvalNegativeSignificantMUAsIS] = corrPlot(areasOI, data, 'negativeSignificantMUAsIS', windowSize, mainFolder);
  
  % US
  [~, rPositiveUnitsUS, pvalPositiveUnitsUS] = corrPlot(areasOI, data, 'positiveUnitsUS', windowSize, mainFolder);
  [~, rPositiveSignificantUnitsUS, pvalPositiveSignificantUnitsUS] = corrPlot(areasOI, data, 'positiveSignificantUnitsUS', windowSize, mainFolder);
  [~, rNegativeUnitsUS, pvalNegativeUnitsUS] = corrPlot(areasOI, data, 'negativeUnitsUS', windowSize, mainFolder);
  [~, rNegativeSignificantUnitsUS, pvalNegativeSignificantUnitsUS] = corrPlot(areasOI, data, 'negativeSignificantUnitsUS', windowSize, mainFolder);
  [~, rPositiveMUAsUS, pvalPositiveMUAsUS] = corrPlot(areasOI, data, 'positiveMUAsUS', windowSize, mainFolder);
  [~, rPositiveSignificantMUAsUS, pvalPositiveSignificantMUAsUS] = corrPlot(areasOI, data, 'positiveSignificantMUAsUS', windowSize, mainFolder);
  [~, rNegativeMUAsUS, pvalNegativeMUAsUS] = corrPlot(areasOI, data, 'negativeMUAsUS', windowSize, mainFolder);
  [~, rNegativeSignificantMUAsUS, pvalNegativeSignificantMUAsUS] = corrPlot(areasOI, data, 'negativeSignificantMUAsUS', windowSize, mainFolder);
  
  % USIS
  [~, rPositiveUnitsUSIS, pvalPositiveUnitsUSIS] = corrPlot(areasOI, data, 'positiveUnitsUSIS', windowSize, mainFolder);
  [~, rPositiveSignificantUnitsUSIS, pvalPositiveSignificantUnitsUSIS] = corrPlot(areasOI, data, 'positiveSignificantUnitsUSIS', windowSize, mainFolder);
  [~, rNegativeUnitsUSIS, pvalNegativeUnitsUSIS] = corrPlot(areasOI, data, 'negativeUnitsUSIS', windowSize, mainFolder);
  [~, rNegativeSignificantUnitsUSIS, pvalNegativeSignificantUnitsUSIS] = corrPlot(areasOI, data, 'negativeSignificantUnitsUSIS', windowSize, mainFolder);
  [~, rPositiveMUAsUSIS, pvalPositiveMUAsUSIS] = corrPlot(areasOI, data, 'positiveMUAsUSIS', windowSize, mainFolder);
  [~, rPositiveSignificantMUAsUSIS, pvalPositiveSignificantMUAsUSIS] = corrPlot(areasOI, data, 'positiveSignificantMUAsUSIS', windowSize, mainFolder);
  [~, rNegativeMUAsUSIS, pvalNegativeMUAsUSIS] = corrPlot(areasOI, data, 'negativeMUAsUSIS', windowSize, mainFolder);
  [~, rNegativeSignificantMUAsUSIS, pvalNegativeSignificantMUAsUSIS] = corrPlot(areasOI, data, 'negativeSignificantMUAsUSIS', windowSize, mainFolder);
  
  % 10-minute windows
  windowSize = 10;
  % IS
  [~, rPositiveUnitsIS10minWindows, pvalPositiveUnitsIS10minWindows] = corrPlot(areasOI, data, 'positiveUnitsIS', windowSize, mainFolder);
  [~, rPositiveSignificantUnitsIS10minWindows, pvalPositiveSignificantUnitsIS10minWindows] = corrPlot(areasOI, data, 'positiveSignificantUnitsIS', windowSize, mainFolder);
  [~, rNegativeUnitsIS10minWindows, pvalNegativeUnitsIS10minWindows] = corrPlot(areasOI, data, 'negativeUnitsIS', windowSize, mainFolder);
  [~, rNegativeSignificantUnitsIS10minWindows, pvalNegativeSignificantUnitsIS10minWindows] = corrPlot(areasOI, data, 'negativeSignificantUnitsIS', windowSize, mainFolder);
  [~, rPositiveMUAsIS10minWindows, pvalPositiveMUAsIS10minWindows] = corrPlot(areasOI, data, 'positiveMUAsIS', windowSize, mainFolder);
  [~, rPositiveSignificantMUAsIS10minWindows, pvalPositiveSignificantMUAsIS10minWindows] = corrPlot(areasOI, data, 'positiveSignificantMUAsIS', windowSize, mainFolder);
  [~, rNegativeMUAsIS10minWindows, pvalNegativeMUAsIS10minWindows] = corrPlot(areasOI, data, 'negativeMUAsIS', windowSize, mainFolder);
  [~, rNegativeSignificantMUAsIS10minWindows, pvalNegativeSignificantMUAsIS10minWindows] = corrPlot(areasOI, data, 'negativeSignificantMUAsIS', windowSize, mainFolder);
  
  % US
  [~, rPositiveUnitsUS10minWindows, pvalPositiveUnitsUS10minWindows] = corrPlot(areasOI, data, 'positiveUnitsUS', windowSize, mainFolder);
  [~, rPositiveSignificantUnitsUS10minWindows, pvalPositiveSignificantUnitsUS10minWindows] = corrPlot(areasOI, data, 'positiveSignificantUnitsUS', windowSize, mainFolder);
  [~, rNegativeUnitsUS10minWindows, pvalNegativeUnitsUS10minWindows] = corrPlot(areasOI, data, 'negativeUnitsUS', windowSize, mainFolder);
  [~, rNegativeSignificantUnitsUS10minWindows, pvalNegativeSignificantUnitsUS10minWindows] = corrPlot(areasOI, data, 'negativeSignificantUnitsUS', windowSize, mainFolder);
  [~, rPositiveMUAsUS10minWindows, pvalPositiveMUAsUS10minWindows] = corrPlot(areasOI, data, 'positiveMUAsUS', windowSize, mainFolder);
  [~, rPositiveSignificantMUAsUS10minWindows, pvalPositiveSignificantMUAsUS10minWindows] = corrPlot(areasOI, data, 'positiveSignificantMUAsUS', windowSize, mainFolder);
  [~, rNegativeMUAsUS10minWindows, pvalNegativeMUAsUS10minWindows] = corrPlot(areasOI, data, 'negativeMUAsUS', windowSize, mainFolder);
  [~, rNegativeSignificantMUAsUS10minWindows, pvalNegativeSignificantMUAsUS10minWindows] = corrPlot(areasOI, data, 'negativeSignificantMUAsUS', windowSize, mainFolder);
  
  % USIS
  [~, rPositiveUnitsUSIS10minWindows, pvalPositiveUnitsUSIS10minWindows] = corrPlot(areasOI, data, 'positiveUnitsUSIS', windowSize, mainFolder);
  [~, rPositiveSignificantUnitsUSIS10minWindows, pvalPositiveSignificantUnitsUSIS10minWindows] = corrPlot(areasOI, data, 'positiveSignificantUnitsUSIS', windowSize, mainFolder);
  [~, rNegativeUnitsUSIS10minWindows, pvalNegativeUnitsUSIS10minWindows] = corrPlot(areasOI, data, 'negativeUnitsUSIS', windowSize, mainFolder);
  [~, rNegativeSignificantUnitsUSIS10minWindows, pvalNegativeSignificantUnitsUSIS10minWindows] = corrPlot(areasOI, data, 'negativeSignificantUnitsUSIS', windowSize, mainFolder);
  [~, rPositiveMUAsUSIS10minWindows, pvalPositiveMUAsUSIS10minWindows] = corrPlot(areasOI, data, 'positiveMUAsUSIS', windowSize, mainFolder);
  [~, rPositiveSignificantMUAsUSIS10minWindows, pvalPositiveSignificantMUAsUSIS10minWindows] = corrPlot(areasOI, data, 'positiveSignificantMUAsUSIS', windowSize, mainFolder);
  [~, rNegativeMUAsUSIS10minWindows, pvalNegativeMUAsUSIS10minWindows] = corrPlot(areasOI, data, 'negativeMUAsUSIS', windowSize, mainFolder);
  [~, rNegativeSignificantMUAsUSIS10minWindows, pvalNegativeSignificantMUAsUSIS10minWindows] = corrPlot(areasOI, data, 'negativeSignificantMUAsUSIS', windowSize, mainFolder);
  
  % 20-minute windows
  windowSize = 20;
  % IS
  [~, rPositiveUnitsIS20minWindows, pvalPositiveUnitsIS20minWindows] = corrPlot(areasOI, data, 'positiveUnitsIS', windowSize, mainFolder);
  [~, rPositiveSignificantUnitsIS20minWindows, pvalPositiveSignificantUnitsIS20minWindows] = corrPlot(areasOI, data, 'positiveSignificantUnitsIS', windowSize, mainFolder);
  [~, rNegativeUnitsIS20minWindows, pvalNegativeUnitsIS20minWindows] = corrPlot(areasOI, data, 'negativeUnitsIS', windowSize, mainFolder);
  [~, rNegativeSignificantUnitsIS20minWindows, pvalNegativeSignificantUnitsIS20minWindows] = corrPlot(areasOI, data, 'negativeSignificantUnitsIS', windowSize, mainFolder);
  [~, rPositiveMUAsIS20minWindows, pvalPositiveMUAsIS20minWindows] = corrPlot(areasOI, data, 'positiveMUAsIS', windowSize, mainFolder);
  [~, rPositiveSignificantMUAsIS20minWindows, pvalPositiveSignificantMUAsIS20minWindows] = corrPlot(areasOI, data, 'positiveSignificantMUAsIS', windowSize, mainFolder);
  [~, rNegativeMUAsIS20minWindows, pvalNegativeMUAsIS20minWindows] = corrPlot(areasOI, data, 'negativeMUAsIS', windowSize, mainFolder);
  [~, rNegativeSignificantMUAsIS20minWindows, pvalNegativeSignificantMUAsIS20minWindows] = corrPlot(areasOI, data, 'negativeSignificantMUAsIS', windowSize, mainFolder);
  
  % US
  [~, rPositiveUnitsUS20minWindows, pvalPositiveUnitsUS20minWindows] = corrPlot(areasOI, data, 'positiveUnitsUS', windowSize, mainFolder);
  [~, rPositiveSignificantUnitsUS20minWindows, pvalPositiveSignificantUnitsUS20minWindows] = corrPlot(areasOI, data, 'positiveSignificantUnitsUS', windowSize, mainFolder);
  [~, rNegativeUnitsUS20minWindows, pvalNegativeUnitsUS20minWindows] = corrPlot(areasOI, data, 'negativeUnitsUS', windowSize, mainFolder);
  [~, rNegativeSignificantUnitsUS20minWindows, pvalNegativeSignificantUnitsUS20minWindows] = corrPlot(areasOI, data, 'negativeSignificantUnitsUS', windowSize, mainFolder);
  [~, rPositiveMUAsUS20minWindows, pvalPositiveMUAsUS20minWindows] = corrPlot(areasOI, data, 'positiveMUAsUS', windowSize, mainFolder);
  [~, rPositiveSignificantMUAsUS20minWindows, pvalPositiveSignificantMUAsUS20minWindows] = corrPlot(areasOI, data, 'positiveSignificantMUAsUS', windowSize, mainFolder);
  [~, rNegativeMUAsUS20minWindows, pvalNegativeMUAsUS20minWindows] = corrPlot(areasOI, data, 'negativeMUAsUS', windowSize, mainFolder);
  [~, rNegativeSignificantMUAsUS20minWindows, pvalNegativeSignificantMUAsUS20minWindows] = corrPlot(areasOI, data, 'negativeSignificantMUAsUS', windowSize, mainFolder);
  
  % USIS
  [~, rPositiveUnitsUSIS20minWindows, pvalPositiveUnitsUSIS20minWindows] = corrPlot(areasOI, data, 'positiveUnitsUSIS', windowSize, mainFolder);
  [~, rPositiveSignificantUnitsUSIS20minWindows, pvalPositiveSignificantUnitsUSIS20minWindows] = corrPlot(areasOI, data, 'positiveSignificantUnitsUSIS', windowSize, mainFolder);
  [~, rNegativeUnitsUSIS20minWindows, pvalNegativeUnitsUSIS20minWindows] = corrPlot(areasOI, data, 'negativeUnitsUSIS', windowSize, mainFolder);
  [~, rNegativeSignificantUnitsUSIS20minWindows, pvalNegativeSignificantUnitsUSIS20minWindows] = corrPlot(areasOI, data, 'negativeSignificantUnitsUSIS', windowSize, mainFolder);
  [~, rPositiveMUAsUSIS20minWindows, pvalPositiveMUAsUSIS20minWindows] = corrPlot(areasOI, data, 'positiveMUAsUSIS', windowSize, mainFolder);
  [~, rPositiveSignificantMUAsUSIS20minWindows, pvalPositiveSignificantMUAsUSIS20minWindows] = corrPlot(areasOI, data, 'positiveSignificantMUAsUSIS', windowSize, mainFolder);
  [~, rNegativeMUAsUSIS20minWindows, pvalNegativeMUAsUSIS20minWindows] = corrPlot(areasOI, data, 'negativeMUAsUSIS', windowSize, mainFolder);
  [~, rNegativeSignificantMUAsUSIS20minWindows, pvalNegativeSignificantMUAsUSIS20minWindows] = corrPlot(areasOI, data, 'negativeSignificantMUAsUSIS', windowSize, mainFolder);
  
  % 30-minute windows
  windowSize = 30;
  % IS
  [~, rPositiveUnitsIS30minWindows, pvalPositiveUnitsIS30minWindows] = corrPlot(areasOI, data, 'positiveUnitsIS', windowSize, mainFolder);
  [~, rPositiveSignificantUnitsIS30minWindows, pvalPositiveSignificantUnitsIS30minWindows] = corrPlot(areasOI, data, 'positiveSignificantUnitsIS', windowSize, mainFolder);
  [~, rNegativeUnitsIS30minWindows, pvalNegativeUnitsIS30minWindows] = corrPlot(areasOI, data, 'negativeUnitsIS', windowSize, mainFolder);
  [~, rNegativeSignificantUnitsIS30minWindows, pvalNegativeSignificantUnitsIS30minWindows] = corrPlot(areasOI, data, 'negativeSignificantUnitsIS', windowSize, mainFolder);
  [~, rPositiveMUAsIS30minWindows, pvalPositiveMUAsIS30minWindows] = corrPlot(areasOI, data, 'positiveMUAsIS', windowSize, mainFolder);
  [~, rPositiveSignificantMUAsIS30minWindows, pvalPositiveSignificantMUAsIS30minWindows] = corrPlot(areasOI, data, 'positiveSignificantMUAsIS', windowSize, mainFolder);
  [~, rNegativeMUAsIS30minWindows, pvalNegativeMUAsIS30minWindows] = corrPlot(areasOI, data, 'negativeMUAsIS', windowSize, mainFolder);
  [~, rNegativeSignificantMUAsIS30minWindows, pvalNegativeSignificantMUAsIS30minWindows] = corrPlot(areasOI, data, 'negativeSignificantMUAsIS', windowSize, mainFolder);
  
  % US
  [~, rPositiveUnitsUS30minWindows, pvalPositiveUnitsUS30minWindows] = corrPlot(areasOI, data, 'positiveUnitsUS', windowSize, mainFolder);
  [~, rPositiveSignificantUnitsUS30minWindows, pvalPositiveSignificantUnitsUS30minWindows] = corrPlot(areasOI, data, 'positiveSignificantUnitsUS', windowSize, mainFolder);
  [~, rNegativeUnitsUS30minWindows, pvalNegativeUnitsUS30minWindows] = corrPlot(areasOI, data, 'negativeUnitsUS', windowSize, mainFolder);
  [~, rNegativeSignificantUnitsUS30minWindows, pvalNegativeSignificantUnitsUS30minWindows] = corrPlot(areasOI, data, 'negativeSignificantUnitsUS', windowSize, mainFolder);
  [~, rPositiveMUAsUS30minWindows, pvalPositiveMUAsUS30minWindows] = corrPlot(areasOI, data, 'positiveMUAsUS', windowSize, mainFolder);
  [~, rPositiveSignificantMUAsUS30minWindows, pvalPositiveSignificantMUAsUS30minWindows] = corrPlot(areasOI, data, 'positiveSignificantMUAsUS', windowSize, mainFolder);
  [~, rNegativeMUAsUS30minWindows, pvalNegativeMUAsUS30minWindows] = corrPlot(areasOI, data, 'negativeMUAsUS', windowSize, mainFolder);
  [~, rNegativeSignificantMUAsUS30minWindows, pvalNegativeSignificantMUAsUS30minWindows] = corrPlot(areasOI, data, 'negativeSignificantMUAsUS', windowSize, mainFolder);
  
  % USIS
  [~, rPositiveUnitsUSIS30minWindows, pvalPositiveUnitsUSIS30minWindows] = corrPlot(areasOI, data, 'positiveUnitsUSIS', windowSize, mainFolder);
  [~, rPositiveSignificantUnitsUSIS30minWindows, pvalPositiveSignificantUnitsUSIS30minWindows] = corrPlot(areasOI, data, 'positiveSignificantUnitsUSIS', windowSize, mainFolder);
  [~, rNegativeUnitsUSIS30minWindows, pvalNegativeUnitsUSIS30minWindows] = corrPlot(areasOI, data, 'negativeUnitsUSIS', windowSize, mainFolder);
  [~, rNegativeSignificantUnitsUSIS30minWindows, pvalNegativeSignificantUnitsUSIS30minWindows] = corrPlot(areasOI, data, 'negativeSignificantUnitsUSIS', windowSize, mainFolder);
  [~, rPositiveMUAsUSIS30minWindows, pvalPositiveMUAsUSIS30minWindows] = corrPlot(areasOI, data, 'positiveMUAsUSIS', windowSize, mainFolder);
  [~, rPositiveSignificantMUAsUSIS30minWindows, pvalPositiveSignificantMUAsUSIS30minWindows] = corrPlot(areasOI, data, 'positiveSignificantMUAsUSIS', windowSize, mainFolder);
  [~, rNegativeMUAsUSIS30minWindows, pvalNegativeMUAsUSIS30minWindows] = corrPlot(areasOI, data, 'negativeMUAsUSIS', windowSize, mainFolder);
  [~, rNegativeSignificantMUAsUSIS30minWindows, pvalNegativeSignificantMUAsUSIS30minWindows] = corrPlot(areasOI, data, 'negativeSignificantMUAsUSIS', windowSize, mainFolder);
  
  
  %% Save correlation data
  save([mainFolder filesep 'pupilCorrFractionsFilter.mat'],...
    'rPositiveUnitsIS', 'pvalPositiveUnitsIS', 'rPositiveSignificantUnitsIS', 'pvalPositiveSignificantUnitsIS', 'rNegativeUnitsIS', 'pvalNegativeUnitsIS',...
    'rNegativeSignificantUnitsIS', 'pvalNegativeSignificantUnitsIS', 'rPositiveMUAsIS', 'pvalPositiveMUAsIS', 'rPositiveSignificantMUAsIS',...
    'pvalPositiveSignificantMUAsIS', 'rNegativeMUAsIS', 'pvalNegativeMUAsIS', 'rNegativeSignificantMUAsIS', 'pvalNegativeSignificantMUAsIS',...
    'rPositiveUnitsUS', 'pvalPositiveUnitsUS', 'rPositiveSignificantUnitsUS', 'pvalPositiveSignificantUnitsUS', 'rNegativeUnitsUS', 'pvalNegativeUnitsUS',...
    'rNegativeSignificantUnitsUS', 'pvalNegativeSignificantUnitsUS', 'rPositiveMUAsUS', 'pvalPositiveMUAsUS', 'rPositiveSignificantMUAsUS',...
    'pvalPositiveSignificantMUAsUS', 'rNegativeMUAsUS', 'pvalNegativeMUAsUS', 'rNegativeSignificantMUAsUS', 'pvalNegativeSignificantMUAsUS',...
    'rPositiveUnitsUSIS', 'pvalPositiveUnitsUSIS', 'rPositiveSignificantUnitsUSIS', 'pvalPositiveSignificantUnitsUSIS', 'rNegativeUnitsUSIS', 'pvalNegativeUnitsUSIS',...
    'rNegativeSignificantUnitsUSIS', 'pvalNegativeSignificantUnitsUSIS', 'rPositiveMUAsUSIS', 'pvalPositiveMUAsUSIS', 'rPositiveSignificantMUAsUSIS',...
    'pvalPositiveSignificantMUAsUSIS', 'rNegativeMUAsUSIS', 'pvalNegativeMUAsUSIS', 'rNegativeSignificantMUAsUSIS', 'pvalNegativeSignificantMUAsUSIS',...
    'rPositiveUnitsIS10minWindows', 'pvalPositiveUnitsIS10minWindows', 'rPositiveSignificantUnitsIS10minWindows', 'pvalPositiveSignificantUnitsIS10minWindows', 'rNegativeUnitsIS10minWindows', 'pvalNegativeUnitsIS10minWindows',...
    'rNegativeSignificantUnitsIS10minWindows', 'pvalNegativeSignificantUnitsIS10minWindows', 'rPositiveMUAsIS10minWindows', 'pvalPositiveMUAsIS10minWindows', 'rPositiveSignificantMUAsIS10minWindows',...
    'pvalPositiveSignificantMUAsIS10minWindows', 'rNegativeMUAsIS10minWindows', 'pvalNegativeMUAsIS10minWindows', 'rNegativeSignificantMUAsIS10minWindows', 'pvalNegativeSignificantMUAsIS',...
    'rPositiveUnitsUS10minWindows', 'pvalPositiveUnitsUS10minWindows', 'rPositiveSignificantUnitsUS10minWindows', 'pvalPositiveSignificantUnitsUS10minWindows', 'rNegativeUnitsUS10minWindows', 'pvalNegativeUnitsUS10minWindows',...
    'rNegativeSignificantUnitsUS10minWindows', 'pvalNegativeSignificantUnitsUS10minWindows', 'rPositiveMUAsUS10minWindows', 'pvalPositiveMUAsUS10minWindows', 'rPositiveSignificantMUAsUS10minWindows',...
    'pvalPositiveSignificantMUAsUS10minWindows', 'rNegativeMUAsUS10minWindows', 'pvalNegativeMUAsUS10minWindows', 'rNegativeSignificantMUAsUS10minWindows', 'pvalNegativeSignificantMUAsUS',...
    'rPositiveUnitsUSIS10minWindows', 'pvalPositiveUnitsUSIS10minWindows', 'rPositiveSignificantUnitsUSIS10minWindows', 'pvalPositiveSignificantUnitsUSIS10minWindows', 'rNegativeUnitsUSIS10minWindows', 'pvalNegativeUnitsUSIS10minWindows',...
    'rNegativeSignificantUnitsUSIS10minWindows', 'pvalNegativeSignificantUnitsUSIS10minWindows', 'rPositiveMUAsUSIS10minWindows', 'pvalPositiveMUAsUSIS10minWindows', 'rPositiveSignificantMUAsUSIS10minWindows',...
    'pvalPositiveSignificantMUAsUSIS10minWindows', 'rNegativeMUAsUSIS10minWindows', 'pvalNegativeMUAsUSIS10minWindows', 'rNegativeSignificantMUAsUSIS10minWindows', 'pvalNegativeSignificantMUAsUSIS',...
    'rPositiveUnitsIS20minWindows', 'pvalPositiveUnitsIS20minWindows', 'rPositiveSignificantUnitsIS20minWindows', 'pvalPositiveSignificantUnitsIS20minWindows', 'rNegativeUnitsIS20minWindows', 'pvalNegativeUnitsIS20minWindows',...
    'rNegativeSignificantUnitsIS20minWindows', 'pvalNegativeSignificantUnitsIS20minWindows', 'rPositiveMUAsIS20minWindows', 'pvalPositiveMUAsIS20minWindows', 'rPositiveSignificantMUAsIS20minWindows',...
    'pvalPositiveSignificantMUAsIS20minWindows', 'rNegativeMUAsIS20minWindows', 'pvalNegativeMUAsIS20minWindows', 'rNegativeSignificantMUAsIS20minWindows', 'pvalNegativeSignificantMUAsIS',...
    'rPositiveUnitsUS20minWindows', 'pvalPositiveUnitsUS20minWindows', 'rPositiveSignificantUnitsUS20minWindows', 'pvalPositiveSignificantUnitsUS20minWindows', 'rNegativeUnitsUS20minWindows', 'pvalNegativeUnitsUS20minWindows',...
    'rNegativeSignificantUnitsUS20minWindows', 'pvalNegativeSignificantUnitsUS20minWindows', 'rPositiveMUAsUS20minWindows', 'pvalPositiveMUAsUS20minWindows', 'rPositiveSignificantMUAsUS20minWindows',...
    'pvalPositiveSignificantMUAsUS20minWindows', 'rNegativeMUAsUS20minWindows', 'pvalNegativeMUAsUS20minWindows', 'rNegativeSignificantMUAsUS20minWindows', 'pvalNegativeSignificantMUAsUS',...
    'rPositiveUnitsUSIS20minWindows', 'pvalPositiveUnitsUSIS20minWindows', 'rPositiveSignificantUnitsUSIS20minWindows', 'pvalPositiveSignificantUnitsUSIS20minWindows', 'rNegativeUnitsUSIS20minWindows', 'pvalNegativeUnitsUSIS20minWindows',...
    'rNegativeSignificantUnitsUSIS20minWindows', 'pvalNegativeSignificantUnitsUSIS20minWindows', 'rPositiveMUAsUSIS20minWindows', 'pvalPositiveMUAsUSIS20minWindows', 'rPositiveSignificantMUAsUSIS20minWindows',...
    'pvalPositiveSignificantMUAsUSIS20minWindows', 'rNegativeMUAsUSIS20minWindows', 'pvalNegativeMUAsUSIS20minWindows', 'rNegativeSignificantMUAsUSIS20minWindows', 'pvalNegativeSignificantMUAsUSIS',...
    'rPositiveUnitsIS30minWindows', 'pvalPositiveUnitsIS30minWindows', 'rPositiveSignificantUnitsIS30minWindows', 'pvalPositiveSignificantUnitsIS30minWindows', 'rNegativeUnitsIS30minWindows', 'pvalNegativeUnitsIS30minWindows',...
    'rNegativeSignificantUnitsIS30minWindows', 'pvalNegativeSignificantUnitsIS30minWindows', 'rPositiveMUAsIS30minWindows', 'pvalPositiveMUAsIS30minWindows', 'rPositiveSignificantMUAsIS30minWindows',...
    'pvalPositiveSignificantMUAsIS30minWindows', 'rNegativeMUAsIS30minWindows', 'pvalNegativeMUAsIS30minWindows', 'rNegativeSignificantMUAsIS30minWindows', 'pvalNegativeSignificantMUAsIS',...
    'rPositiveUnitsUS30minWindows', 'pvalPositiveUnitsUS30minWindows', 'rPositiveSignificantUnitsUS30minWindows', 'pvalPositiveSignificantUnitsUS30minWindows', 'rNegativeUnitsUS30minWindows', 'pvalNegativeUnitsUS30minWindows',...
    'rNegativeSignificantUnitsUS30minWindows', 'pvalNegativeSignificantUnitsUS30minWindows', 'rPositiveMUAsUS30minWindows', 'pvalPositiveMUAsUS30minWindows', 'rPositiveSignificantMUAsUS30minWindows',...
    'pvalPositiveSignificantMUAsUS30minWindows', 'rNegativeMUAsUS30minWindows', 'pvalNegativeMUAsUS30minWindows', 'rNegativeSignificantMUAsUS30minWindows', 'pvalNegativeSignificantMUAsUS',...
    'rPositiveUnitsUSIS30minWindows', 'pvalPositiveUnitsUSIS30minWindows', 'rPositiveSignificantUnitsUSIS30minWindows', 'pvalPositiveSignificantUnitsUSIS30minWindows', 'rNegativeUnitsUSIS30minWindows', 'pvalNegativeUnitsUSIS30minWindows',...
    'rNegativeSignificantUnitsUSIS30minWindows', 'pvalNegativeSignificantUnitsUSIS30minWindows', 'rPositiveMUAsUSIS30minWindows', 'pvalPositiveMUAsUSIS30minWindows', 'rPositiveSignificantMUAsUSIS30minWindows',...
    'pvalPositiveSignificantMUAsUSIS30minWindows', 'rNegativeMUAsUSIS30minWindows', 'pvalNegativeMUAsUSIS30minWindows', 'rNegativeSignificantMUAsUSIS30minWindows', 'pvalNegativeSignificantMUAsUSIS',...
    'statsFractionsUnitsNegativeSSignificant','statsFractionsMUAsNegativeSSignificant','statsFractionsMUAsNegativeSSignificant','statsFractionsMUAsNegativeSSignificant',...
    'statsFractionsUnitsNegativeISSignificant','statsFractionsMUAsNegativeISSignificant','statsFractionsMUAsNegativeISSignificant','statsFractionsMUAsNegativeISSignificant',...
    'statsFractionsUnitsNegativeUSSignificant','statsFractionsMUAsNegativeUSSignificant','statsFractionsMUAsNegativeUSSignificant','statsFractionsMUAsNegativeUSSignificant', 'areasOI', '-append');
  
else
  load([mainFolder filesep 'pupilCorrFractionsFilter.mat']); %#ok<*UNRCH>
end


if fullRun <= 3
  
  %% Produce tables
  
  % MUAs IS: S1 v RSC
  if strcmp(repository, 'uol')
    inds = findComboEntries('lS1', 'lRSC', areasOI);
    positiveMUAsISTable = [rPositiveMUAsIS10minWindows(inds) rPositiveMUAsIS20minWindows(inds) rPositiveMUAsIS30minWindows(inds) rPositiveMUAsIS(inds)];
    positiveMUAsISTable = [positiveMUAsISTable; positiveMUAsISTable.^2];
    positiveMUAsISTable = [positiveMUAsISTable; [pvalPositiveMUAsIS10minWindows(inds) pvalPositiveMUAsIS20minWindows(inds) pvalPositiveMUAsIS30minWindows(inds) pvalPositiveMUAsIS(inds)]];
    negativeMUAsISTable = [rPositiveMUAsIS10minWindows(inds) rPositiveMUAsIS20minWindows(inds) rPositiveMUAsIS30minWindows(inds) rPositiveMUAsIS(inds)];
    negativeMUAsISTable = [negativeMUAsISTable; negativeMUAsISTable.^2];
    negativeMUAsISTable = [negativeMUAsISTable; [pvalPositiveMUAsIS10minWindows(inds) pvalPositiveMUAsIS20minWindows(inds) pvalPositiveMUAsIS30minWindows(inds) pvalPositiveMUAsIS(inds)]];
    MUAsISTable_S1vRSC = [positiveMUAsISTable; negativeMUAsISTable];
  end
  
  % MUAs IS: S1 v CA
  if strcmp(repository, 'uol')
    inds = findComboEntries('lS1', 'CA', areasOI);
    positiveMUAsISTable = [rPositiveMUAsIS10minWindows(inds) rPositiveMUAsIS20minWindows(inds) rPositiveMUAsIS30minWindows(inds) rPositiveMUAsIS(inds)];
    positiveMUAsISTable = [positiveMUAsISTable; positiveMUAsISTable.^2];
    positiveMUAsISTable = [positiveMUAsISTable; [pvalPositiveMUAsIS10minWindows(inds) pvalPositiveMUAsIS20minWindows(inds) pvalPositiveMUAsIS30minWindows(inds) pvalPositiveMUAsIS(inds)]];
    negativeMUAsISTable = [rPositiveMUAsIS10minWindows(inds) rPositiveMUAsIS20minWindows(inds) rPositiveMUAsIS30minWindows(inds) rPositiveMUAsIS(inds)];
    negativeMUAsISTable = [negativeMUAsISTable; negativeMUAsISTable.^2];
    negativeMUAsISTable = [negativeMUAsISTable; [pvalPositiveMUAsIS10minWindows(inds) pvalPositiveMUAsIS20minWindows(inds) pvalPositiveMUAsIS30minWindows(inds) pvalPositiveMUAsIS(inds)]];
    MUAsISTable_S1vCA = [positiveMUAsISTable; negativeMUAsISTable];
  end
  
  % MUAs IS: RSC v CA
  if strcmp(repository, 'uol')
    inds = findComboEntries('lRSC', 'CA', areasOI);
    positiveMUAsISTable = [rPositiveMUAsIS10minWindows(inds) rPositiveMUAsIS20minWindows(inds) rPositiveMUAsIS30minWindows(inds) rPositiveMUAsIS(inds)];
    positiveMUAsISTable = [positiveMUAsISTable; positiveMUAsISTable.^2];
    positiveMUAsISTable = [positiveMUAsISTable; [pvalPositiveMUAsIS10minWindows(inds) pvalPositiveMUAsIS20minWindows(inds) pvalPositiveMUAsIS30minWindows(inds) pvalPositiveMUAsIS(inds)]];
    negativeMUAsISTable = [rPositiveMUAsIS10minWindows(inds) rPositiveMUAsIS20minWindows(inds) rPositiveMUAsIS30minWindows(inds) rPositiveMUAsIS(inds)];
    negativeMUAsISTable = [negativeMUAsISTable; negativeMUAsISTable.^2];
    negativeMUAsISTable = [negativeMUAsISTable; [pvalPositiveMUAsIS10minWindows(inds) pvalPositiveMUAsIS20minWindows(inds) pvalPositiveMUAsIS30minWindows(inds) pvalPositiveMUAsIS(inds)]];
    MUAsISTable_RSCvCA = [positiveMUAsISTable; negativeMUAsISTable];
  end
  
  % MUAs IS: V1 v CA
  if strcmp(repository, 'allensdk')
    inds = findComboEntries('V1', 'CA', areasOI);
    positiveMUAsISTable = [rPositiveMUAsIS10minWindows(inds) rPositiveMUAsIS20minWindows(inds) rPositiveMUAsIS30minWindows(inds) rPositiveMUAsIS(inds)];
    positiveMUAsISTable = [positiveMUAsISTable; positiveMUAsISTable.^2];
    positiveMUAsISTable = [positiveMUAsISTable; [pvalPositiveMUAsIS10minWindows(inds) pvalPositiveMUAsIS20minWindows(inds) pvalPositiveMUAsIS30minWindows(inds) pvalPositiveMUAsIS(inds)]];
    negativeMUAsISTable = [rPositiveMUAsIS10minWindows(inds) rPositiveMUAsIS20minWindows(inds) rPositiveMUAsIS30minWindows(inds) rPositiveMUAsIS(inds)];
    negativeMUAsISTable = [negativeMUAsISTable; negativeMUAsISTable.^2];
    negativeMUAsISTable = [negativeMUAsISTable; [pvalPositiveMUAsIS10minWindows(inds) pvalPositiveMUAsIS20minWindows(inds) pvalPositiveMUAsIS30minWindows(inds) pvalPositiveMUAsIS(inds)]];
    MUAsISTable_V1vCA = [positiveMUAsISTable; negativeMUAsISTable];
  end
  
  % MUAs IS: S1 v RSC
  if strcmp(repository, 'uol')
    inds = findComboEntries('lS1', 'lRSC', areasOI);
    positiveUnitsISTable = [rPositiveUnitsIS10minWindows(inds) rPositiveUnitsIS20minWindows(inds) rPositiveUnitsIS30minWindows(inds) rPositiveUnitsIS(inds)];
    positiveUnitsISTable = [positiveUnitsISTable; positiveUnitsISTable.^2];
    positiveUnitsISTable = [positiveUnitsISTable; [pvalPositiveUnitsIS10minWindows(inds) pvalPositiveUnitsIS20minWindows(inds) pvalPositiveUnitsIS30minWindows(inds) pvalPositiveUnitsIS(inds)]];
    negativeUnitsISTable = [rPositiveUnitsIS10minWindows(inds) rPositiveUnitsIS20minWindows(inds) rPositiveUnitsIS30minWindows(inds) rPositiveUnitsIS(inds)];
    negativeUnitsISTable = [negativeUnitsISTable; negativeUnitsISTable.^2];
    negativeUnitsISTable = [negativeUnitsISTable; [pvalPositiveUnitsIS10minWindows(inds) pvalPositiveUnitsIS20minWindows(inds) pvalPositiveUnitsIS30minWindows(inds) pvalPositiveUnitsIS(inds)]];
    unitsISTable_S1vRSC = [positiveUnitsISTable; negativeUnitsISTable];
  end
  
  % Units IS: S1 v CA
  if strcmp(repository, 'uol')
    inds = findComboEntries('lS1', 'CA', areasOI);
    positiveUnitsISTable = [rPositiveUnitsIS10minWindows(inds) rPositiveUnitsIS20minWindows(inds) rPositiveUnitsIS30minWindows(inds) rPositiveUnitsIS(inds)];
    positiveUnitsISTable = [positiveUnitsISTable; positiveUnitsISTable.^2];
    positiveUnitsISTable = [positiveUnitsISTable; [pvalPositiveUnitsIS10minWindows(inds) pvalPositiveUnitsIS20minWindows(inds) pvalPositiveUnitsIS30minWindows(inds) pvalPositiveUnitsIS(inds)]];
    negativeUnitsISTable = [rPositiveUnitsIS10minWindows(inds) rPositiveUnitsIS20minWindows(inds) rPositiveUnitsIS30minWindows(inds) rPositiveUnitsIS(inds)];
    negativeUnitsISTable = [negativeUnitsISTable; negativeUnitsISTable.^2];
    negativeUnitsISTable = [negativeUnitsISTable; [pvalPositiveUnitsIS10minWindows(inds) pvalPositiveUnitsIS20minWindows(inds) pvalPositiveUnitsIS30minWindows(inds) pvalPositiveUnitsIS(inds)]];
    unitsISTable_S1vCA = [positiveUnitsISTable; negativeUnitsISTable];
  end
  
  % Units IS: RSC v CA
  if strcmp(repository, 'uol')
    inds = findComboEntries('lRSC', 'CA', areasOI);
    positiveUnitsISTable = [rPositiveUnitsIS10minWindows(inds) rPositiveUnitsIS20minWindows(inds) rPositiveUnitsIS30minWindows(inds) rPositiveUnitsIS(inds)];
    positiveUnitsISTable = [positiveUnitsISTable; positiveUnitsISTable.^2];
    positiveUnitsISTable = [positiveUnitsISTable; [pvalPositiveUnitsIS10minWindows(inds) pvalPositiveUnitsIS20minWindows(inds) pvalPositiveUnitsIS30minWindows(inds) pvalPositiveUnitsIS(inds)]];
    negativeUnitsISTable = [rPositiveUnitsIS10minWindows(inds) rPositiveUnitsIS20minWindows(inds) rPositiveUnitsIS30minWindows(inds) rPositiveUnitsIS(inds)];
    negativeUnitsISTable = [negativeUnitsISTable; negativeUnitsISTable.^2];
    negativeUnitsISTable = [negativeUnitsISTable; [pvalPositiveUnitsIS10minWindows(inds) pvalPositiveUnitsIS20minWindows(inds) pvalPositiveUnitsIS30minWindows(inds) pvalPositiveUnitsIS(inds)]];
    unitsISTable_RSCvCA = [positiveUnitsISTable; negativeUnitsISTable];
  end
  
  % Units IS: V1 v CA
  if strcmp(repository, 'allensdk')
    inds = findComboEntries('V1', 'CA', areasOI);
    positiveUnitsISTable = [rPositiveUnitsIS10minWindows(inds) rPositiveUnitsIS20minWindows(inds) rPositiveUnitsIS30minWindows(inds) rPositiveUnitsIS(inds)];
    positiveUnitsISTable = [positiveUnitsISTable; positiveUnitsISTable.^2];
    positiveUnitsISTable = [positiveUnitsISTable; [pvalPositiveUnitsIS10minWindows(inds) pvalPositiveUnitsIS20minWindows(inds) pvalPositiveUnitsIS30minWindows(inds) pvalPositiveUnitsIS(inds)]];
    negativeUnitsISTable = [rPositiveUnitsIS10minWindows(inds) rPositiveUnitsIS20minWindows(inds) rPositiveUnitsIS30minWindows(inds) rPositiveUnitsIS(inds)];
    negativeUnitsISTable = [negativeUnitsISTable; negativeUnitsISTable.^2];
    negativeUnitsISTable = [negativeUnitsISTable; [pvalPositiveUnitsIS10minWindows(inds) pvalPositiveUnitsIS20minWindows(inds) pvalPositiveUnitsIS30minWindows(inds) pvalPositiveUnitsIS(inds)]];
    unitsISTable_V1vCA = [positiveUnitsISTable; negativeUnitsISTable];
  end
  
  % MUAs US: S1 v RSC
  if strcmp(repository, 'uol')
    inds = findComboEntries('lS1', 'lRSC', areasOI);
    positiveMUAsUSTable = [rPositiveMUAsUS10minWindows(inds) rPositiveMUAsUS20minWindows(inds) rPositiveMUAsUS30minWindows(inds) rPositiveMUAsUS(inds)];
    positiveMUAsUSTable = [positiveMUAsUSTable; positiveMUAsUSTable.^2];
    positiveMUAsUSTable = [positiveMUAsUSTable; [pvalPositiveMUAsUS10minWindows(inds) pvalPositiveMUAsUS20minWindows(inds) pvalPositiveMUAsUS30minWindows(inds) pvalPositiveMUAsUS(inds)]];
    negativeMUAsUSTable = [rPositiveMUAsUS10minWindows(inds) rPositiveMUAsUS20minWindows(inds) rPositiveMUAsUS30minWindows(inds) rPositiveMUAsUS(inds)];
    negativeMUAsUSTable = [negativeMUAsUSTable; negativeMUAsUSTable.^2];
    negativeMUAsUSTable = [negativeMUAsUSTable; [pvalPositiveMUAsUS10minWindows(inds) pvalPositiveMUAsUS20minWindows(inds) pvalPositiveMUAsUS30minWindows(inds) pvalPositiveMUAsUS(inds)]];
    MUAsUSTable_S1vRSC = [positiveMUAsUSTable; negativeMUAsUSTable];
  end
  
  % MUAs US: S1 v CA
  if strcmp(repository, 'uol')
    inds = findComboEntries('lS1', 'CA', areasOI);
    positiveMUAsUSTable = [rPositiveMUAsUS10minWindows(inds) rPositiveMUAsUS20minWindows(inds) rPositiveMUAsUS30minWindows(inds) rPositiveMUAsUS(inds)];
    positiveMUAsUSTable = [positiveMUAsUSTable; positiveMUAsUSTable.^2];
    positiveMUAsUSTable = [positiveMUAsUSTable; [pvalPositiveMUAsUS10minWindows(inds) pvalPositiveMUAsUS20minWindows(inds) pvalPositiveMUAsUS30minWindows(inds) pvalPositiveMUAsUS(inds)]];
    negativeMUAsUSTable = [rPositiveMUAsUS10minWindows(inds) rPositiveMUAsUS20minWindows(inds) rPositiveMUAsUS30minWindows(inds) rPositiveMUAsUS(inds)];
    negativeMUAsUSTable = [negativeMUAsUSTable; negativeMUAsUSTable.^2];
    negativeMUAsUSTable = [negativeMUAsUSTable; [pvalPositiveMUAsUS10minWindows(inds) pvalPositiveMUAsUS20minWindows(inds) pvalPositiveMUAsUS30minWindows(inds) pvalPositiveMUAsUS(inds)]];
    MUAsUSTable_S1vCA = [positiveMUAsUSTable; negativeMUAsUSTable];
  end
  
  % MUAs US: RSC v CA
  if strcmp(repository, 'uol')
    inds = findComboEntries('lRSC', 'CA', areasOI);
    positiveMUAsUSTable = [rPositiveMUAsUS10minWindows(inds) rPositiveMUAsUS20minWindows(inds) rPositiveMUAsUS30minWindows(inds) rPositiveMUAsUS(inds)];
    positiveMUAsUSTable = [positiveMUAsUSTable; positiveMUAsUSTable.^2];
    positiveMUAsUSTable = [positiveMUAsUSTable; [pvalPositiveMUAsUS10minWindows(inds) pvalPositiveMUAsUS20minWindows(inds) pvalPositiveMUAsUS30minWindows(inds) pvalPositiveMUAsUS(inds)]];
    negativeMUAsUSTable = [rPositiveMUAsUS10minWindows(inds) rPositiveMUAsUS20minWindows(inds) rPositiveMUAsUS30minWindows(inds) rPositiveMUAsUS(inds)];
    negativeMUAsUSTable = [negativeMUAsUSTable; negativeMUAsUSTable.^2];
    negativeMUAsUSTable = [negativeMUAsUSTable; [pvalPositiveMUAsUS10minWindows(inds) pvalPositiveMUAsUS20minWindows(inds) pvalPositiveMUAsUS30minWindows(inds) pvalPositiveMUAsUS(inds)]];
    MUAsUSTable_RSCvCA = [positiveMUAsUSTable; negativeMUAsUSTable];
  end
  
  % MUAs US: V1 v CA
  if strcmp(repository, 'allensdk')
    inds = findComboEntries('V1', 'CA', areasOI);
    positiveMUAsUSTable = [rPositiveMUAsUS10minWindows(inds) rPositiveMUAsUS20minWindows(inds) rPositiveMUAsUS30minWindows(inds) rPositiveMUAsUS(inds)];
    positiveMUAsUSTable = [positiveMUAsUSTable; positiveMUAsUSTable.^2];
    positiveMUAsUSTable = [positiveMUAsUSTable; [pvalPositiveMUAsUS10minWindows(inds) pvalPositiveMUAsUS20minWindows(inds) pvalPositiveMUAsUS30minWindows(inds) pvalPositiveMUAsUS(inds)]];
    negativeMUAsUSTable = [rPositiveMUAsUS10minWindows(inds) rPositiveMUAsUS20minWindows(inds) rPositiveMUAsUS30minWindows(inds) rPositiveMUAsUS(inds)];
    negativeMUAsUSTable = [negativeMUAsUSTable; negativeMUAsUSTable.^2];
    negativeMUAsUSTable = [negativeMUAsUSTable; [pvalPositiveMUAsUS10minWindows(inds) pvalPositiveMUAsUS20minWindows(inds) pvalPositiveMUAsUS30minWindows(inds) pvalPositiveMUAsUS(inds)]];
    MUAsUSTable_V1vCA = [positiveMUAsUSTable; negativeMUAsUSTable];
  end
  
  % Units US: S1 v RSC
  if strcmp(repository, 'uol')
    inds = findComboEntries('lS1', 'lRSC', areasOI);
    positiveUnitsUSTable = [rPositiveUnitsUS10minWindows(inds) rPositiveUnitsUS20minWindows(inds) rPositiveUnitsUS30minWindows(inds) rPositiveUnitsUS(inds)];
    positiveUnitsUSTable = [positiveUnitsUSTable; positiveUnitsUSTable.^2];
    positiveUnitsUSTable = [positiveUnitsUSTable; [pvalPositiveUnitsUS10minWindows(inds) pvalPositiveUnitsUS20minWindows(inds) pvalPositiveUnitsUS30minWindows(inds) pvalPositiveUnitsUS(inds)]];
    negativeUnitsUSTable = [rPositiveUnitsUS10minWindows(inds) rPositiveUnitsUS20minWindows(inds) rPositiveUnitsUS30minWindows(inds) rPositiveUnitsUS(inds)];
    negativeUnitsUSTable = [negativeUnitsUSTable; negativeUnitsUSTable.^2];
    negativeUnitsUSTable = [negativeUnitsUSTable; [pvalPositiveUnitsUS10minWindows(inds) pvalPositiveUnitsUS20minWindows(inds) pvalPositiveUnitsUS30minWindows(inds) pvalPositiveUnitsUS(inds)]];
    unitsUSTable_S1vRSC = [positiveUnitsUSTable; negativeUnitsUSTable];
  end
  
  % Units US: S1 v CA
  if strcmp(repository, 'uol')
    inds = findComboEntries('lS1', 'CA', areasOI);
    positiveUnitsUSTable = [rPositiveUnitsUS10minWindows(inds) rPositiveUnitsUS20minWindows(inds) rPositiveUnitsUS30minWindows(inds) rPositiveUnitsUS(inds)];
    positiveUnitsUSTable = [positiveUnitsUSTable; positiveUnitsUSTable.^2];
    positiveUnitsUSTable = [positiveUnitsUSTable; [pvalPositiveUnitsUS10minWindows(inds) pvalPositiveUnitsUS20minWindows(inds) pvalPositiveUnitsUS30minWindows(inds) pvalPositiveUnitsUS(inds)]];
    negativeUnitsUSTable = [rPositiveUnitsUS10minWindows(inds) rPositiveUnitsUS20minWindows(inds) rPositiveUnitsUS30minWindows(inds) rPositiveUnitsUS(inds)];
    negativeUnitsUSTable = [negativeUnitsUSTable; negativeUnitsUSTable.^2];
    negativeUnitsUSTable = [negativeUnitsUSTable; [pvalPositiveUnitsUS10minWindows(inds) pvalPositiveUnitsUS20minWindows(inds) pvalPositiveUnitsUS30minWindows(inds) pvalPositiveUnitsUS(inds)]];
    unitsUSTable_S1vCA = [positiveUnitsUSTable; negativeUnitsUSTable];
  end
  
  % Units US: RSC v CA
  if strcmp(repository, 'uol')
    inds = findComboEntries('lRSC', 'CA', areasOI);
    positiveUnitsUSTable = [rPositiveUnitsUS10minWindows(inds) rPositiveUnitsUS20minWindows(inds) rPositiveUnitsUS30minWindows(inds) rPositiveUnitsUS(inds)];
    positiveUnitsUSTable = [positiveUnitsUSTable; positiveUnitsUSTable.^2];
    positiveUnitsUSTable = [positiveUnitsUSTable; [pvalPositiveUnitsUS10minWindows(inds) pvalPositiveUnitsUS20minWindows(inds) pvalPositiveUnitsUS30minWindows(inds) pvalPositiveUnitsUS(inds)]];
    negativeUnitsUSTable = [rPositiveUnitsUS10minWindows(inds) rPositiveUnitsUS20minWindows(inds) rPositiveUnitsUS30minWindows(inds) rPositiveUnitsUS(inds)];
    negativeUnitsUSTable = [negativeUnitsUSTable; negativeUnitsUSTable.^2];
    negativeUnitsUSTable = [negativeUnitsUSTable; [pvalPositiveUnitsUS10minWindows(inds) pvalPositiveUnitsUS20minWindows(inds) pvalPositiveUnitsUS30minWindows(inds) pvalPositiveUnitsUS(inds)]];
    unitsUSTable_RSCvCA = [positiveUnitsUSTable; negativeUnitsUSTable];
  end
  
  % Units US: V1 v CA
  if strcmp(repository, 'allensdk')
    inds = findComboEntries('V1', 'CA', areasOI);
    positiveUnitsUSTable = [rPositiveUnitsUS10minWindows(inds) rPositiveUnitsUS20minWindows(inds) rPositiveUnitsUS30minWindows(inds) rPositiveUnitsUS(inds)];
    positiveUnitsUSTable = [positiveUnitsUSTable; positiveUnitsUSTable.^2];
    positiveUnitsUSTable = [positiveUnitsUSTable; [pvalPositiveUnitsUS10minWindows(inds) pvalPositiveUnitsUS20minWindows(inds) pvalPositiveUnitsUS30minWindows(inds) pvalPositiveUnitsUS(inds)]];
    negativeUnitsUSTable = [rPositiveUnitsUS10minWindows(inds) rPositiveUnitsUS20minWindows(inds) rPositiveUnitsUS30minWindows(inds) rPositiveUnitsUS(inds)];
    negativeUnitsUSTable = [negativeUnitsUSTable; negativeUnitsUSTable.^2];
    negativeUnitsUSTable = [negativeUnitsUSTable; [pvalPositiveUnitsUS10minWindows(inds) pvalPositiveUnitsUS20minWindows(inds) pvalPositiveUnitsUS30minWindows(inds) pvalPositiveUnitsUS(inds)]];
    unitsUSTable_V1vCA = [positiveUnitsUSTable; negativeUnitsUSTable];
  end
  
  if strcmp(repository, 'uol')
    MUAsTable = [MUAsISTable_S1vRSC MUAsISTable_S1vCA MUAsISTable_RSCvCA MUAsUSTable_S1vRSC MUAsUSTable_S1vCA MUAsUSTable_RSCvCA];
    unitsTable = [unitsISTable_S1vRSC unitsISTable_S1vCA unitsISTable_RSCvCA unitsUSTable_S1vRSC unitsUSTable_S1vCA unitsUSTable_RSCvCA];
  end
  
  
  %% Save the tables
  if strcmp(repository, 'uol')
    save([mainFolder filesep 'pupilCorrFractionsFilter.mat'],...
      'MUAsISTable_S1vRSC', 'MUAsISTable_S1vCA', 'MUAsISTable_RSCvCA',...
      'unitsISTable_S1vRSC', 'unitsISTable_S1vCA', 'unitsISTable_RSCvCA',...
      'MUAsUSTable_S1vRSC', 'MUAsUSTable_S1vCA', 'MUAsUSTable_RSCvCA',...
      'unitsUSTable_S1vRSC', 'unitsUSTable_S1vCA', 'unitsUSTable_RSCvCA',...
      'MUAsTable', 'unitsTable', '-append');
  else
    save([mainFolder filesep 'pupilCorrFractionsFilter.mat'],...
      'MUAsISTable_V1vCA', 'unitsISTable_V1vCA', 'MUAsUSTable_V1vCA', 'unitsUSTable_V1vCA',...
      '-append');
  end
  
else
  load([mainFolder filesep 'pupilCorrFractionsFilter.mat']); %#ok<*UNRCH>
end



%% Local functions
function [dataMean, dataScatter, areaGroup, colourGroup, corrGroup] = barGroup(areaName, colourNumber, barPosition, data, type)

areaCode = determineArea(areaName);
colour = matlabColours(colourNumber);
if strcmp(type, 'IS')
  bar(1+barPosition, data.positiveUnitsFrMeanISThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour, 'FaceAlpha',0.2);
  bar(2+barPosition, data.positiveSignificantUnitsFrMeanISThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
  bar(3+barPosition, data.negativeUnitsFrMeanISThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3), 'FaceAlpha',0.2);
  bar(4+barPosition, data.negativeSignificantUnitsFrMeanISThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [data.positiveUnitsFrMeanISThr{1}{areaCode} data.positiveSignificantUnitsFrMeanISThr{1}{areaCode}...
    data.negativeUnitsFrMeanISThr{1}{areaCode} data.negativeSignificantUnitsFrMeanISThr{1}{areaCode}];
  dataScatter = {data.positiveUnitsFrISThr{1}{areaCode} data.positiveSignificantUnitsFrISThr{1}{areaCode}...
    data.negativeUnitsFrISThr{1}{areaCode} data.negativeSignificantUnitsFrISThr{1}{areaCode}};
elseif strcmp(type, 'US')
  bar(1+barPosition, data.positiveUnitsFrMeanUSThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour, 'FaceAlpha',0.2);
  bar(2+barPosition, data.positiveSignificantUnitsFrMeanUSThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
  bar(3+barPosition, data.negativeUnitsFrMeanUSThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3), 'FaceAlpha',0.2);
  bar(4+barPosition, data.negativeSignificantUnitsFrMeanUSThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [data.positiveUnitsFrMeanUSThr{1}{areaCode} data.positiveSignificantUnitsFrMeanUSThr{1}{areaCode}...
    data.negativeUnitsFrMeanUSThr{1}{areaCode} data.negativeSignificantUnitsFrMeanUSThr{1}{areaCode}];
  dataScatter = {data.positiveUnitsFrUSThr{1}{areaCode} data.positiveSignificantUnitsFrUSThr{1}{areaCode}...
    data.negativeUnitsFrUSThr{1}{areaCode} data.negativeSignificantUnitsFrUSThr{1}{areaCode}};
elseif strcmp(type, 'USIS')
  bar(1+barPosition, data.positiveUnitsFrMeanUSISThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour, 'FaceAlpha',0.2);
  bar(2+barPosition, data.positiveSignificantUnitsFrMeanUSISThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
  bar(3+barPosition, data.negativeUnitsFrMeanUSISThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3), 'FaceAlpha',0.2);
  bar(4+barPosition, data.negativeSignificantUnitsFrMeanUSISThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [data.positiveUnitsFrMeanUSISThr{1}{areaCode} data.positiveSignificantUnitsFrMeanUSISThr{1}{areaCode}...
    data.negativeUnitsFrMeanUSISThr{1}{areaCode} data.negativeSignificantUnitsFrMeanUSISThr{1}{areaCode}];
  dataScatter = {data.positiveUnitsFrUSISThr{1}{areaCode} data.positiveSignificantUnitsFrUSISThr{1}{areaCode}...
    data.negativeUnitsFrUSISThr{1}{areaCode} data.negativeSignificantUnitsFrUSISThr{1}{areaCode}};
end
areaGroup = {areaName, areaName, areaName, areaName};
colourGroup = [colour; colour; colour; colour];
corrGroup = {'positive','positive','negative','negative'};
end


function [dataMean, dataScatter, areaGroup, colourGroup, corrGroup] = barGroupMUAs(areaName, colourNumber, barPosition, data, type)

areaCode = determineArea(areaName);
colour = matlabColours(colourNumber);
if strcmp(type, 'IS')
  bar(1+barPosition, data.positiveMUAsFrMeanISThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour, 'FaceAlpha',0.2);
  bar(2+barPosition, data.positiveSignificantMUAsFrMeanISThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
  bar(3+barPosition, data.negativeMUAsFrMeanISThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3), 'FaceAlpha',0.2);
  bar(4+barPosition, data.negativeSignificantMUAsFrMeanISThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [data.positiveMUAsFrMeanISThr{1}{areaCode} data.positiveSignificantMUAsFrMeanISThr{1}{areaCode}...
    data.negativeMUAsFrMeanISThr{1}{areaCode} data.negativeSignificantMUAsFrMeanISThr{1}{areaCode}];
  dataScatter = {data.positiveMUAsFrISThr{1}{areaCode} data.positiveSignificantMUAsFrISThr{1}{areaCode}...
    data.negativeMUAsFrISThr{1}{areaCode} data.negativeSignificantMUAsFrISThr{1}{areaCode}};
elseif strcmp(type, 'US')
  bar(1+barPosition, data.positiveMUAsFrMeanUSThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour, 'FaceAlpha',0.2);
  bar(2+barPosition, data.positiveSignificantMUAsFrMeanUSThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
  bar(3+barPosition, data.negativeMUAsFrMeanUSThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3), 'FaceAlpha',0.2);
  bar(4+barPosition, data.negativeSignificantMUAsFrMeanUSThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [data.positiveMUAsFrMeanUSThr{1}{areaCode} data.positiveSignificantMUAsFrMeanUSThr{1}{areaCode}...
    data.negativeMUAsFrMeanUSThr{1}{areaCode} data.negativeSignificantMUAsFrMeanUSThr{1}{areaCode}];
  dataScatter = {data.positiveMUAsFrUSThr{1}{areaCode} data.positiveSignificantMUAsFrUSThr{1}{areaCode}...
    data.negativeMUAsFrUSThr{1}{areaCode} data.negativeSignificantMUAsFrUSThr{1}{areaCode}};
elseif strcmp(type, 'USIS')
  bar(1+barPosition, data.positiveMUAsFrMeanUSISThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour, 'FaceAlpha',0.2);
  bar(2+barPosition, data.positiveSignificantMUAsFrMeanUSISThr{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
  bar(3+barPosition, data.negativeMUAsFrMeanUSISThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3), 'FaceAlpha',0.2);
  bar(4+barPosition, data.negativeSignificantMUAsFrMeanUSISThr{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [data.positiveMUAsFrMeanUSISThr{1}{areaCode} data.positiveSignificantMUAsFrMeanUSISThr{1}{areaCode}...
    data.negativeMUAsFrMeanUSISThr{1}{areaCode} data.negativeSignificantMUAsFrMeanUSISThr{1}{areaCode}];
  dataScatter = {data.positiveMUAsFrUSISThr{1}{areaCode} data.positiveSignificantMUAsFrUSISThr{1}{areaCode}...
    data.negativeMUAsFrUSISThr{1}{areaCode} data.negativeSignificantMUAsFrUSISThr{1}{areaCode}};
end
areaGroup = {areaName, areaName, areaName, areaName};
colourGroup = [colour; colour; colour; colour];
corrGroup = {'positive','positive','negative','negative'};
end


function [p, pSig] = ttestGroup(areaName, data, type)

areaCode = determineArea(areaName);
if strcmp(type, 'IS')
  [~,p] = ttest2(data.positiveUnitsFrISThr{1}{areaCode}, data.negativeUnitsFrISThr{1}{areaCode});
  [~,pSig] = ttest2(data.positiveSignificantUnitsFrISThr{1}{areaCode}, data.negativeSignificantUnitsFrISThr{1}{areaCode});
elseif strcmp(type, 'US')
  [~,p] = ttest2(data.positiveUnitsFrUSThr{1}{areaCode}, data.negativeUnitsFrUSThr{1}{areaCode});
  [~,pSig] = ttest2(data.positiveSignificantUnitsFrUSThr{1}{areaCode}, data.negativeSignificantUnitsFrUSThr{1}{areaCode});
elseif strcmp(type, 'USIS')
  [~,p] = ttest2(data.positiveUnitsFrUSISThr{1}{areaCode}, data.negativeUnitsFrUSISThr{1}{areaCode});
  [~,pSig] = ttest2(data.positiveSignificantUnitsFrUSISThr{1}{areaCode}, data.negativeSignificantUnitsFrUSISThr{1}{areaCode});
end
end


function [p, pSig] = ttestGroupMUAs(areaName, data, type)

areaCode = determineArea(areaName);
if strcmp(type, 'IS')
  [~,p] = ttest2(data.positiveMUAsFrISThr{1}{areaCode}, data.negativeMUAsFrISThr{1}{areaCode});
  [~,pSig] = ttest2(data.positiveSignificantMUAsFrISThr{1}{areaCode}, data.negativeSignificantMUAsFrISThr{1}{areaCode});
elseif strcmp(type, 'US')
  [~,p] = ttest2(data.positiveMUAsFrUSThr{1}{areaCode}, data.negativeMUAsFrUSThr{1}{areaCode});
  [~,pSig] = ttest2(data.positiveSignificantMUAsFrUSThr{1}{areaCode}, data.negativeSignificantMUAsFrUSThr{1}{areaCode});
elseif strcmp(type, 'USIS')
  [~,p] = ttest2(data.positiveMUAsFrUSISThr{1}{areaCode}, data.negativeMUAsFrUSISThr{1}{areaCode});
  [~,pSig] = ttest2(data.positiveSignificantMUAsFrUSISThr{1}{areaCode}, data.negativeSignificantMUAsFrUSISThr{1}{areaCode});
end
end


function [fH, scatterGroups, areaGroups, colourGroups, colorCode, corrGroups] = barPlotUnits(data, area1, area2, area3, area4, yLim, type)

[pVB, pVBSig] = ttestGroup(area1, data, type);
[plS1, plS1Sig] = ttestGroup(area2, data, type);
[plRSC, plRSCSig] = ttestGroup(area3, data, type);
[pCA, pCASig] = ttestGroup(area4, data, type);


%% Draw the bar graphs
fH = figProperties('Bar plot for units', 'normalized', [0, .005, .97, .90], 'w', 'on');
hold on


% Data bars
gap = 1;
nBarsPerGroup = 4;
[barGroup1, scatterGroup1, areaGroup1, colourGroup1, corrGroup1] = barGroup(area1, 1, 0*(nBarsPerGroup+gap), data, type);
[barGroup2, scatterGroup2, areaGroup2, colourGroup2, corrGroup2] = barGroup(area2, 2, 1*(nBarsPerGroup+gap), data, type);
[barGroup3, scatterGroup3, areaGroup3, colourGroup3, corrGroup3] = barGroup(area3, 3, 2*(nBarsPerGroup+gap), data, type);
[barGroup4, scatterGroup4, areaGroup4, colourGroup4, corrGroup4] = barGroup(area4, 4, 3*(nBarsPerGroup+gap), data, type);


% Error bars
gaps = 5:5:15;
bars = sort([1 1+gaps 2 2+gaps 3 3+gaps 4 4+gaps]);
dataMean = [barGroup1 barGroup2 barGroup3 barGroup4];
if strcmp(type, 'IS')
  err = [data.positiveUnitsFrISCI95Thr{1}{determineArea(area1)} data.positiveSignificantUnitsFrISCI95Thr{1}{determineArea(area1)}...
    data.negativeUnitsFrISCI95Thr{1}{determineArea(area1)} data.negativeSignificantUnitsFrISCI95Thr{1}{determineArea(area1)}...
    data.positiveUnitsFrISCI95Thr{1}{determineArea(area2)} data.positiveSignificantUnitsFrISCI95Thr{1}{determineArea(area2)}...
    data.negativeUnitsFrISCI95Thr{1}{determineArea(area2)} data.negativeSignificantUnitsFrISCI95Thr{1}{determineArea(area2)}...
    data.positiveUnitsFrISCI95Thr{1}{determineArea(area3)} data.positiveSignificantUnitsFrISCI95Thr{1}{determineArea(area3)}...
    data.negativeUnitsFrISCI95Thr{1}{determineArea(area3)} data.negativeSignificantUnitsFrISCI95Thr{1}{determineArea(area3)}...
    data.positiveUnitsFrISCI95Thr{1}{determineArea(area4)} data.positiveSignificantUnitsFrISCI95Thr{1}{determineArea(area4)}...
    data.negativeUnitsFrISCI95Thr{1}{determineArea(area4)} data.negativeSignificantUnitsFrISCI95Thr{1}{determineArea(area4)}];
elseif strcmp(type, 'US')
  err = [data.positiveUnitsFrUSCI95Thr{1}{determineArea(area1)} data.positiveSignificantUnitsFrUSCI95Thr{1}{determineArea(area1)}...
    data.negativeUnitsFrUSCI95Thr{1}{determineArea(area1)} data.negativeSignificantUnitsFrUSCI95Thr{1}{determineArea(area1)}...
    data.positiveUnitsFrUSCI95Thr{1}{determineArea(area2)} data.positiveSignificantUnitsFrUSCI95Thr{1}{determineArea(area2)}...
    data.negativeUnitsFrUSCI95Thr{1}{determineArea(area2)} data.negativeSignificantUnitsFrUSCI95Thr{1}{determineArea(area2)}...
    data.positiveUnitsFrUSCI95Thr{1}{determineArea(area3)} data.positiveSignificantUnitsFrUSCI95Thr{1}{determineArea(area3)}...
    data.negativeUnitsFrUSCI95Thr{1}{determineArea(area3)} data.negativeSignificantUnitsFrUSCI95Thr{1}{determineArea(area3)}...
    data.positiveUnitsFrUSCI95Thr{1}{determineArea(area4)} data.positiveSignificantUnitsFrUSCI95Thr{1}{determineArea(area4)}...
    data.negativeUnitsFrUSCI95Thr{1}{determineArea(area4)} data.negativeSignificantUnitsFrUSCI95Thr{1}{determineArea(area4)}];
elseif strcmp(type, 'USIS')
  err = [data.positiveUnitsFrUSISCI95Thr{1}{determineArea(area1)} data.positiveSignificantUnitsFrUSISCI95Thr{1}{determineArea(area1)}...
    data.negativeUnitsFrUSISCI95Thr{1}{determineArea(area1)} data.negativeSignificantUnitsFrUSISCI95Thr{1}{determineArea(area1)}...
    data.positiveUnitsFrUSISCI95Thr{1}{determineArea(area2)} data.positiveSignificantUnitsFrUSISCI95Thr{1}{determineArea(area2)}...
    data.negativeUnitsFrUSISCI95Thr{1}{determineArea(area2)} data.negativeSignificantUnitsFrUSISCI95Thr{1}{determineArea(area2)}...
    data.positiveUnitsFrUSISCI95Thr{1}{determineArea(area3)} data.positiveSignificantUnitsFrUSISCI95Thr{1}{determineArea(area3)}...
    data.negativeUnitsFrUSISCI95Thr{1}{determineArea(area3)} data.negativeSignificantUnitsFrUSISCI95Thr{1}{determineArea(area3)}...
    data.positiveUnitsFrUSISCI95Thr{1}{determineArea(area4)} data.positiveSignificantUnitsFrUSISCI95Thr{1}{determineArea(area4)}...
    data.negativeUnitsFrUSISCI95Thr{1}{determineArea(area4)} data.negativeSignificantUnitsFrUSISCI95Thr{1}{determineArea(area4)}];
end

er = errorbar(bars,dataMean,err(2,:),err(1,:));
er.Color = [0 0 0];
er.LineStyle = 'none';


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
ax1.XTickLabel = {area1,area2,area3,area4};

xLim = xlim;
xAxisLength = xLim(2)-xLim(1);
yAxisLength = yLim(2)-yLim(1);
textStr = ['t-test +v- p-val for all units: ' num2str(pVB) '(VB) ' num2str(plS1) '(lS1) ' num2str(plRSC) '(lRSC) ' num2str(pCA) '(CA)'];
text(xLim(2)-xAxisLength*0.85, yLim(2)-yAxisLength*0.05, textStr, 'FontSize',20);
textStr = ['t-test +v- p-val for significant units only: ' num2str(pVBSig) '(VB) ' num2str(plS1Sig) '(lS1) ' num2str(plRSCSig) '(lRSC) ' num2str(pCASig) '(CA)'];
text(xLim(2)-xAxisLength*0.85, yLim(2)-yAxisLength*0.01, textStr, 'FontSize',20);
text(1-0.2,0.03*yLim(2), 'All', 'FontSize',20)
text(2-0.25,0.03*yLim(2), 'Sig', 'FontSize',20)
for iBar = bars([1 2 5 6 9 10 13 14])
  text(iBar-0.075,-0.01*yLim(2), '+', 'FontSize',12)
  text(iBar+2-0.075,-0.01*yLim(2), '-', 'FontSize',12)
end
end


function [fH, scatterGroups, areaGroups, colourGroups, colorCode, corrGroups] = barPlotMUAs(data, area1, area2, area3, area4, yLim, type)

[pVB, pVBSig] = ttestGroupMUAs(area1, data, type);
[plS1, plS1Sig] = ttestGroupMUAs(area2, data, type);
[plRSC, plRSCSig] = ttestGroupMUAs(area3, data, type);
[pCA, pCASig] = ttestGroupMUAs(area4, data, type);


%% Draw the bar graphs
fH = figProperties('Bar plot for MUAs', 'normalized', [0, .005, .97, .90], 'w', 'on');
hold on


% Data bars
gap = 1;
nBarsPerGroup = 4;
[barGroup1, scatterGroup1, areaGroup1, colourGroup1, corrGroup1] = barGroupMUAs(area1, 1, 0*(nBarsPerGroup+gap), data, type);
[barGroup2, scatterGroup2, areaGroup2, colourGroup2, corrGroup2] = barGroupMUAs(area2, 2, 1*(nBarsPerGroup+gap), data, type);
[barGroup3, scatterGroup3, areaGroup3, colourGroup3, corrGroup3] = barGroupMUAs(area3, 3, 2*(nBarsPerGroup+gap), data, type);
[barGroup4, scatterGroup4, areaGroup4, colourGroup4, corrGroup4] = barGroupMUAs(area4, 4, 3*(nBarsPerGroup+gap), data, type);


% Error bars
gaps = 5:5:15;
bars = sort([1 1+gaps 2 2+gaps 3 3+gaps 4 4+gaps]);
dataMean = [barGroup1 barGroup2 barGroup3 barGroup4];
if strcmp(type, 'IS')
  err = [data.positiveMUAsFrISCI95Thr{1}{determineArea(area1)} data.positiveSignificantMUAsFrISCI95Thr{1}{determineArea(area1)}...
    data.negativeMUAsFrISCI95Thr{1}{determineArea(area1)} data.negativeSignificantMUAsFrISCI95Thr{1}{determineArea(area1)}...
    data.positiveMUAsFrISCI95Thr{1}{determineArea(area2)} data.positiveSignificantMUAsFrISCI95Thr{1}{determineArea(area2)}...
    data.negativeMUAsFrISCI95Thr{1}{determineArea(area2)} data.negativeSignificantMUAsFrISCI95Thr{1}{determineArea(area2)}...
    data.positiveMUAsFrISCI95Thr{1}{determineArea(area3)} data.positiveSignificantMUAsFrISCI95Thr{1}{determineArea(area3)}...
    data.negativeMUAsFrISCI95Thr{1}{determineArea(area3)} data.negativeSignificantMUAsFrISCI95Thr{1}{determineArea(area3)}...
    data.positiveMUAsFrISCI95Thr{1}{determineArea(area4)} data.positiveSignificantMUAsFrISCI95Thr{1}{determineArea(area4)}...
    data.negativeMUAsFrISCI95Thr{1}{determineArea(area4)} data.negativeSignificantMUAsFrISCI95Thr{1}{determineArea(area4)}];
elseif strcmp(type, 'US')
  err = [data.positiveMUAsFrUSCI95Thr{1}{determineArea(area1)} data.positiveSignificantMUAsFrUSCI95Thr{1}{determineArea(area1)}...
    data.negativeMUAsFrUSCI95Thr{1}{determineArea(area1)} data.negativeSignificantMUAsFrUSCI95Thr{1}{determineArea(area1)}...
    data.positiveMUAsFrUSCI95Thr{1}{determineArea(area2)} data.positiveSignificantMUAsFrUSCI95Thr{1}{determineArea(area2)}...
    data.negativeMUAsFrUSCI95Thr{1}{determineArea(area2)} data.negativeSignificantMUAsFrUSCI95Thr{1}{determineArea(area2)}...
    data.positiveMUAsFrUSCI95Thr{1}{determineArea(area3)} data.positiveSignificantMUAsFrUSCI95Thr{1}{determineArea(area3)}...
    data.negativeMUAsFrUSCI95Thr{1}{determineArea(area3)} data.negativeSignificantMUAsFrUSCI95Thr{1}{determineArea(area3)}...
    data.positiveMUAsFrUSCI95Thr{1}{determineArea(area4)} data.positiveSignificantMUAsFrUSCI95Thr{1}{determineArea(area4)}...
    data.negativeMUAsFrUSCI95Thr{1}{determineArea(area4)} data.negativeSignificantMUAsFrUSCI95Thr{1}{determineArea(area4)}];
elseif strcmp(type, 'USIS')
  err = [data.positiveMUAsFrUSISCI95Thr{1}{determineArea(area1)} data.positiveSignificantMUAsFrUSISCI95Thr{1}{determineArea(area1)}...
    data.negativeMUAsFrUSISCI95Thr{1}{determineArea(area1)} data.negativeSignificantMUAsFrUSISCI95Thr{1}{determineArea(area1)}...
    data.positiveMUAsFrUSISCI95Thr{1}{determineArea(area2)} data.positiveSignificantMUAsFrUSISCI95Thr{1}{determineArea(area2)}...
    data.negativeMUAsFrUSISCI95Thr{1}{determineArea(area2)} data.negativeSignificantMUAsFrUSISCI95Thr{1}{determineArea(area2)}...
    data.positiveMUAsFrUSISCI95Thr{1}{determineArea(area3)} data.positiveSignificantMUAsFrUSISCI95Thr{1}{determineArea(area3)}...
    data.negativeMUAsFrUSISCI95Thr{1}{determineArea(area3)} data.negativeSignificantMUAsFrUSISCI95Thr{1}{determineArea(area3)}...
    data.positiveMUAsFrUSISCI95Thr{1}{determineArea(area4)} data.positiveSignificantMUAsFrUSISCI95Thr{1}{determineArea(area4)}...
    data.negativeMUAsFrUSISCI95Thr{1}{determineArea(area4)} data.negativeSignificantMUAsFrUSISCI95Thr{1}{determineArea(area4)}];
end

er = errorbar(bars,dataMean,err(2,:),err(1,:));
er.Color = [0 0 0];
er.LineStyle = 'none';


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
ax1.XTickLabel = {area1,area2,area3,area4};

xLim = xlim;
xAxisLength = xLim(2)-xLim(1);
yAxisLength = yLim(2)-yLim(1);
textStr = ['t-test +v- p-val for all units+MUAs: ' num2str(pVB) '(VB) ' num2str(plS1) '(lS1) ' num2str(plRSC) '(lRSC) ' num2str(pCA) '(CA)'];
text(xLim(2)-xAxisLength*0.85, yLim(2)-yAxisLength*0.05, textStr, 'FontSize',20);
textStr = ['t-test +v- p-val for significant units+MUAs only: ' num2str(pVBSig) '(VB) ' num2str(plS1Sig) '(lS1) ' num2str(plRSCSig) '(lRSC) ' num2str(pCASig) '(CA)'];
text(xLim(2)-xAxisLength*0.85, yLim(2)-yAxisLength*0.1, textStr, 'FontSize',20);
text(1-0.2,0.03*yLim(2), 'All', 'FontSize',20)
text(2-0.25,0.03*yLim(2), 'Sig', 'FontSize',20)
for iBar = bars([1 2 5 6 9 10 13 14])
  text(iBar-0.075,-0.01*yLim(2), '+', 'FontSize',12)
  text(iBar+2-0.075,-0.01*yLim(2), '-', 'FontSize',12)
end
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
  if strcmpi(type, 'positiveUnitsIS')
    areaInd1 = ismember(data.(['positiveUnitsFrISThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveUnitsFrISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveUnitsFrISThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveUnitsFrISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveUnitsFrISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveUnitsFrISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveUnitsFrISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveUnitsFrISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['IS positive unit fraction in ' areaName1], ['IS positive unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['IS positive unit fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['IS positive unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveSignificantUnitsIS')
    areaInd1 = ismember(data.(['positiveSignificantUnitsFrISThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveSignificantUnitsFrISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveSignificantUnitsFrISThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveSignificantUnitsFrISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveSignificantUnitsFrISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveSignificantUnitsFrISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveSignificantUnitsFrISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveSignificantUnitsFrISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['IS positive significant unit fraction in ' areaName1], ['IS positive significant unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['IS positive significant unit fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['IS positive significant unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeUnitsIS')
    areaInd1 = ismember(data.(['negativeUnitsFrISThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeUnitsFrISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeUnitsFrISThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeUnitsFrISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeUnitsFrISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeUnitsFrISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeUnitsFrISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeUnitsFrISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['IS negative unit fraction in ' areaName1], ['IS negative unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['IS negative unit fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['IS negative unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeSignificantUnitsIS')
    areaInd1 = ismember(data.(['negativeSignificantUnitsFrISThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeSignificantUnitsFrISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeSignificantUnitsFrISThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeSignificantUnitsFrISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeSignificantUnitsFrISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeSignificantUnitsFrISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeSignificantUnitsFrISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeSignificantUnitsFrISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['IS negative significant unit fraction in ' areaName1], ['IS negative significant unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['IS negative significant unit fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['IS negative significant unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveMUAsIS')
    areaInd1 = ismember(data.(['positiveMUAsFrISThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveMUAsFrISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveMUAsFrISThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveMUAsFrISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveMUAsFrISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveMUAsFrISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveMUAsFrISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveMUAsFrISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['IS positive unit+MUA fraction in ' areaName1], ['IS positive unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['IS positive unit and MUA fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['IS positive unit and MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveSignificantMUAsIS')
    areaInd1 = ismember(data.(['positiveSignificantMUAsFrISThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveSignificantMUAsFrISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveSignificantMUAsFrISThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveSignificantMUAsFrISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveSignificantMUAsFrISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveSignificantMUAsFrISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveSignificantMUAsFrISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveSignificantMUAsFrISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1, areaVec2, 'Spearman');
    fH(iCombo) = xyPlot(areaVec1', areaVec2', 1, ['IS positive significant unit+MUA fraction in ' areaName1], ['IS positive significant unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['IS positive singificant unit and MUA fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['IS positive singificant unit and MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeMUAsIS')
    areaInd1 = ismember(data.(['negativeMUAsFrISThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeMUAsFrISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeMUAsFrISThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeMUAsFrISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeMUAsFrISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeMUAsFrISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeMUAsFrISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeMUAsFrISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['IS negative unit+MUA fraction in ' areaName1], ['IS negative unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['IS negative unit and MUA fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['IS negative unit and MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeSignificantMUAsIS')
    areaInd1 = ismember(data.(['negativeSignificantMUAsFrISThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeSignificantMUAsFrISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeSignificantMUAsFrISThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeSignificantMUAsFrISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeSignificantMUAsFrISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeSignificantMUAsFrISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeSignificantMUAsFrISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeSignificantMUAsFrISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['IS negative significant unit+MUA fraction in ' areaName1], ['IS negative significant unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['IS negative significant unit and MUA fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['IS negative significant unit and MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveUnitsUS')
    areaInd1 = ismember(data.(['positiveUnitsFrUSThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveUnitsFrUSThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveUnitsFrUSThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveUnitsFrUSThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveUnitsFrUSThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveUnitsFrUSThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveUnitsFrUSThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveUnitsFrUSThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['US positive unit fraction in ' areaName1], ['US positive unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['US positive unit fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['US positive unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveSignificantUnitsUS')
    areaInd1 = ismember(data.(['positiveSignificantUnitsFrUSThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveSignificantUnitsFrUSThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveSignificantUnitsFrUSThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveSignificantUnitsFrUSThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveSignificantUnitsFrUSThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveSignificantUnitsFrUSThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveSignificantUnitsFrUSThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveSignificantUnitsFrUSThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['US positive significant unit fraction in ' areaName1], ['US positive significant unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['US positive significant unit fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['US positive significant unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeUnitsUS')
    areaInd1 = ismember(data.(['negativeUnitsFrUSThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeUnitsFrUSThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeUnitsFrUSThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeUnitsFrUSThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeUnitsFrUSThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeUnitsFrUSThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeUnitsFrUSThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeUnitsFrUSThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['US negative unit fraction in ' areaName1], ['US negative unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['US negative unit fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['US negative unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeSignificantUnitsUS')
    areaInd1 = ismember(data.(['negativeSignificantUnitsFrUSThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeSignificantUnitsFrUSThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeSignificantUnitsFrUSThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeSignificantUnitsFrUSThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeSignificantUnitsFrUSThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeSignificantUnitsFrUSThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeSignificantUnitsFrUSThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeSignificantUnitsFrUSThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['US negative significant unit fraction in ' areaName1], ['US negative significant unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['US negative significant unit fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['US negative significant unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveMUAsUS')
    areaInd1 = ismember(data.(['positiveMUAsFrUSThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveMUAsFrUSThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveMUAsFrUSThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveMUAsFrUSThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveMUAsFrUSThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveMUAsFrUSThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveMUAsFrUSThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveMUAsFrUSThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['US positive unit+MUA fraction in ' areaName1], ['US positive unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['US positive unit and MUA fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['US positive unit and MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveSignificantMUAsUS')
    areaInd1 = ismember(data.(['positiveSignificantMUAsFrUSThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveSignificantMUAsFrUSThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveSignificantMUAsFrUSThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveSignificantMUAsFrUSThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveSignificantMUAsFrUSThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveSignificantMUAsFrUSThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveSignificantMUAsFrUSThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveSignificantMUAsFrUSThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1, areaVec2, 'Spearman');
    fH(iCombo) = xyPlot(areaVec1', areaVec2', 1, ['US positive significant unit+MUA fraction in ' areaName1], ['US positive significant unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['US positive singificant unit and MUA fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['US positive singificant unit and MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeMUAsUS')
    areaInd1 = ismember(data.(['negativeMUAsFrUSThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeMUAsFrUSThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeMUAsFrUSThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeMUAsFrUSThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeMUAsFrUSThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeMUAsFrUSThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeMUAsFrUSThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeMUAsFrUSThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['US negative unit+MUA fraction in ' areaName1], ['US negative unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['US negative unit and MUA fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['US negative unit and MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeSignificantMUAsUS')
    areaInd1 = ismember(data.(['negativeSignificantMUAsFrUSThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeSignificantMUAsFrUSThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeSignificantMUAsFrUSThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeSignificantMUAsFrUSThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeSignificantMUAsFrUSThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeSignificantMUAsFrUSThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeSignificantMUAsFrUSThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeSignificantMUAsFrUSThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['US negative significant unit+MUA fraction in ' areaName1], ['US negative significant unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['US negative significant unit and MUA fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['US negative significant unit and MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveUnitsUSIS')
    areaInd1 = ismember(data.(['positiveUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveUnitsFrUSISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveUnitsFrUSISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Positive US_I_S/positive US unit fraction in ' areaName1], ['Positive US_I_S/positive US unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Positive US_I_S/positive US unit fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['Positive US_I_S/positive US unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveSignificantUnitsUSIS')
    areaInd1 = ismember(data.(['positiveSignificantUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveSignificantUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveSignificantUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveSignificantUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveSignificantUnitsFrUSISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveSignificantUnitsFrUSISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveSignificantUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveSignificantUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Positive US_I_S/positive US significant unit fraction in ' areaName1], ['Positive US_I_S/positive US significant unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Positive US_I_S/positive US significant unit fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['Positive US_I_S/positive US significant unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeUnitsUSIS')
    areaInd1 = ismember(data.(['negativeUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeUnitsFrUSISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeUnitsFrUSISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Negative US_I_S/negative US unit fraction in ' areaName1], ['Negative US_I_S/negative US unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Negative US_I_S/negative US unit fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['Negative US_I_S/negative US unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeSignificantUnitsUSIS')
    areaInd1 = ismember(data.(['negativeSignificantUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeSignificantUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeSignificantUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeSignificantUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeSignificantUnitsFrUSISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeSignificantUnitsFrUSISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeSignificantUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeSignificantUnitsFrUSISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Negative US_I_S/negative US significant unit fraction in ' areaName1], ['Negative US_I_S/negative US significant unit fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Negative US_I_S/negative US significant unit fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['Negative US_I_S/negative US significant unit fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveMUAsUSIS')
    areaInd1 = ismember(data.(['positiveMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveMUAsFrUSISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveMUAsFrUSISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Positive US_I_S/positive US unit+MUA fraction in ' areaName1], ['Positive US_I_S/positive US unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Positive US_I_S/positive US unit and MUA fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['Positive US_I_S/positive US unit and MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'positiveSignificantMUAsUSIS')
    areaInd1 = ismember(data.(['positiveSignificantMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode1}, data.(['positiveSignificantMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['positiveSignificantMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode2}, data.(['positiveSignificantMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['positiveSignificantMUAsFrUSISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['positiveSignificantMUAsFrUSISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['positiveSignificantMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['positiveSignificantMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1, areaVec2, 'Spearman');
    fH(iCombo) = xyPlot(areaVec1', areaVec2', 1, ['Positive US_I_S/positive US significant unit+MUA fraction in ' areaName1], ['Positive US_I_S/positive US significant unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Positive US_I_S/positive US singificant unit and MUA fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['Positive US_I_S/positive US singificant unit and MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeMUAsUSIS')
    areaInd1 = ismember(data.(['negativeMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeMUAsFrUSISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeMUAsFrUSISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Negative US_I_S/negative US unit+MUA fraction in ' areaName1], ['Negative US_I_S/negative US unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Negative US_I_S/negative US unit and MUA fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['Negative US_I_S/negative US unit and MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
    end
  elseif strcmpi(type, 'negativeSignificantMUAsUSIS')
    areaInd1 = ismember(data.(['negativeSignificantMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode1}, data.(['negativeSignificantMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode2});
    areaInd2 = ismember(data.(['negativeSignificantMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode2}, data.(['negativeSignificantMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode1});
    areaVec1 = data.(['negativeSignificantMUAsFrUSISThr' windowSize]){1}{areaCode1}(areaInd1);
    areaVec2 = data.(['negativeSignificantMUAsFrUSISThr' windowSize]){1}{areaCode2}(areaInd2);
    areaVec1 = reduce2unique(data.(['negativeSignificantMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode1}(areaInd1), areaVec1);
    areaVec2 = reduce2unique(data.(['negativeSignificantMUAsFrUSISThr' windowSize 'RecID']){1}{areaCode2}(areaInd2), areaVec2);
    [r(iCombo), pval(iCombo)] = corrSimple(areaVec1', areaVec2', 'Spearman');
    fH(iCombo) = xyPlot(areaVec1, areaVec2, 1, ['Negative US_I_S/negative US significant unit+MUA fraction in ' areaName1], ['Negative US_I_S/negative US significant unit+MUA fraction in ' areaName2]);
    if isempty(windowSize)
      figName = ['Negative US_I_S/negative US significant unit and MUA fraction in ' areaName1 ' v ' areaName2];
    else
      figName = ['Negative US_I_S/negative US significant unit and MUA fraction in ' areaName1 ' v ' areaName2 ': ' windowSize];
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
  
  title(figName);
  figName = strrep(figName, 'I_S', 'IS');
  set(fH(iCombo), 'Name',figName);
  if isempty(areaInd1) && isempty(areaInd2)
    continue
  end
  figName = strrep(figName, ':', '_');
  filename = [outputDir filesep figName];
  filename = strrep(filename, ' ', '_');
  filename = strrep(filename, '/', '_among_');
  savefig(fH(iCombo), filename, 'compact');
  print(fH(iCombo), [filename '.png'],'-dpng','-r300');
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


function [uniqueVecPos, uniqueVecPosSignificant, uniqueVecNeg, uniqueVecNegSignificant] = reduce2unique2(nameVec, corrVec, significanceVec)

[~, uniqueElements] = unique(nameVec, 'stable');
uniqueVecPos = zeros(numel(uniqueElements),1);
for i = 1:numel(uniqueElements)
  if i == numel(uniqueElements)
    corrVeci = corrVec(uniqueElements(i):end);
    significanceVeci = significanceVec(uniqueElements(i):end);
    uniqueVecPos(i) = sum(corrVeci >= 0)/numel(corrVeci);
    uniqueVecPosSignificant(i) = sum(corrVeci >= 0 & significanceVeci < 0.05)/numel(corrVeci);
    uniqueVecNeg(i) = sum(corrVeci < 0)/numel(corrVeci);
    uniqueVecNegSignificant(i) = sum(corrVeci < 0 & significanceVeci < 0.05)/numel(corrVeci);
  else
    corrVeci = corrVec(uniqueElements(i):uniqueElements(i+1)-1);
    significanceVeci = significanceVec(uniqueElements(i):uniqueElements(i+1)-1);
    uniqueVecPos(i) = sum(corrVeci >= 0)/numel(corrVeci);
    uniqueVecPosSignificant(i) = sum(corrVeci >= 0 & significanceVeci < 0.05)/numel(corrVeci);
    uniqueVecNeg(i) = sum(corrVeci < 0)/numel(corrVeci);
    uniqueVecNegSignificant(i) = sum(corrVeci < 0 & significanceVeci < 0.05)/numel(corrVeci);
  end
end
end


function data = assignValues2Windows(dbStruct, windowDuration, data, iCondPlusAll, iAreaPlusAll, seriesName, unitThr, alpha) %#ok<*DEFNU>

filterType = {'IS', 'US', 'USIS'};
for iType = 1:numel(filterType)
  
  % units
  type = [filterType{iType} num2str(windowDuration) 'minWindows'];
  shankIDs = fieldnames(dbStruct.shankData);
  nUnitsWindows = zeros(numel(dbStruct.popData.(['meanPupilArea' num2str(windowDuration) 'minWindows'])),1);
  for iWindow = 1:numel(dbStruct.popData.(['meanPupilArea' num2str(windowDuration) 'minWindows']))
    rSpearmanWindowIS = [];
    pvalSpearmanWindowIS = [];
    rSpearmanWindowUS = [];
    pvalSpearmanWindowUS = [];
    for sh = 1:numel(shankIDs)
      rSpearmanWindowIS = [rSpearmanWindowIS; dbStruct.shankData.(shankIDs{sh}).(['rSpearmanISIS' num2str(windowDuration) 'minWindows']){iWindow}'];
      pvalSpearmanWindowIS = [pvalSpearmanWindowIS; dbStruct.shankData.(shankIDs{sh}).(['pvalSpearmanISIS' num2str(windowDuration) 'minWindows']){iWindow}'];
      rSpearmanWindowUS = [rSpearmanWindowUS; dbStruct.shankData.(shankIDs{sh}).(['rSpearmanUSUS' num2str(windowDuration) 'minWindows']){iWindow}'];
      pvalSpearmanWindowUS = [pvalSpearmanWindowUS; dbStruct.shankData.(shankIDs{sh}).(['pvalSpearmanUSUS' num2str(windowDuration) 'minWindows']){iWindow}'];
    end
    if strcmp(filterType{iType}, 'IS')
      rSpearmanWindow = rSpearmanWindowIS;
      pvalSpearmanWindow = pvalSpearmanWindowIS;
      positiveDenominator = numel(rSpearmanWindow);
      negativeDenominator = numel(rSpearmanWindow);
    elseif strcmp(filterType{iType}, 'US')
      rSpearmanWindow = rSpearmanWindowUS;
      pvalSpearmanWindow = pvalSpearmanWindowUS;
      positiveDenominator = numel(rSpearmanWindow);
      negativeDenominator = numel(rSpearmanWindow);
    elseif strcmp(filterType{iType}, 'USIS')
      rSpearmanWindow = NaN(size(rSpearmanWindowIS));
      rSpearmanWindow(rSpearmanWindowIS >= 0 & rSpearmanWindowUS >= 0) = 1;
      rSpearmanWindow(rSpearmanWindowIS < 0 & rSpearmanWindowUS < 0) = -1;
      pvalSpearmanWindow = NaN(size(pvalSpearmanWindowIS));
      pvalSpearmanWindow(pvalSpearmanWindowIS < alpha & pvalSpearmanWindowUS < alpha) = 0;
      pvalSpearmanWindow(~(pvalSpearmanWindowIS < alpha & pvalSpearmanWindowUS < alpha)) = 1;
      positiveDenominator = sum(rSpearmanWindowUS >= 0);
      negativeDenominator = sum(rSpearmanWindowUS < 0);
    end
    nUnitsWindows(iWindow) = numel(rSpearmanWindow);
    data.(['positiveUnitsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['positiveUnitsFr' type]){iCondPlusAll}{iAreaPlusAll};...
      sum(rSpearmanWindow >= 0)/positiveDenominator];
    data.(['positiveUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['positiveUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    data.(['positiveSignificantUnitsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['positiveSignificantUnitsFr' type]){iCondPlusAll}{iAreaPlusAll};...
      sum(rSpearmanWindow >= 0 & pvalSpearmanWindow < alpha)/positiveDenominator];
    data.(['positiveSignificantUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['positiveSignificantUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    data.(['negativeUnitsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['negativeUnitsFr' type]){iCondPlusAll}{iAreaPlusAll};...
      sum(rSpearmanWindow < 0)/negativeDenominator];
    data.(['negativeUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['negativeUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    data.(['negativeSignificantUnitsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['negativeSignificantUnitsFr' type]){iCondPlusAll}{iAreaPlusAll};...
      sum(rSpearmanWindow < 0 & pvalSpearmanWindow < alpha)/negativeDenominator];
    data.(['negativeSignificantUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['negativeSignificantUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    if nUnitsWindows(iWindow) >= unitThr
      type = [filterType{iType} 'Thr' num2str(windowDuration) 'minWindows'];
      data.(['positiveUnitsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['positiveUnitsFr' type]){iCondPlusAll}{iAreaPlusAll};...
        sum(rSpearmanWindow >= 0)/positiveDenominator];
      data.(['positiveUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['positiveUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
        [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
      data.(['positiveSignificantUnitsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['positiveSignificantUnitsFr' type]){iCondPlusAll}{iAreaPlusAll};...
        sum(rSpearmanWindow >= 0 & pvalSpearmanWindow < alpha)/positiveDenominator];
      data.(['positiveSignificantUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['positiveSignificantUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
        [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
      data.(['negativeUnitsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['negativeUnitsFr' type]){iCondPlusAll}{iAreaPlusAll};...
        sum(rSpearmanWindow < 0)/negativeDenominator];
      data.(['negativeUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['negativeUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
        [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
      data.(['negativeSignificantUnitsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['negativeSignificantUnitsFr' type]){iCondPlusAll}{iAreaPlusAll};...
        sum(rSpearmanWindow < 0 & pvalSpearmanWindow < alpha)/negativeDenominator];
      data.(['negativeSignificantUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['negativeSignificantUnitsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
        [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    end
  end
  
  % MUAs
  type = [filterType{iType} num2str(windowDuration) 'minWindows'];
  for iWindow = 1:numel(dbStruct.popData.(['meanPupilArea' num2str(windowDuration) 'minWindows']))
    rSpearmanWindowIS = dbStruct.popData.(['rSpearmanISIS' num2str(windowDuration) 'minWindows']){iWindow}';
    pvalSpearmanWindowIS = dbStruct.popData.(['pvalSpearmanISIS' num2str(windowDuration) 'minWindows']){iWindow}';
    rSpearmanWindowUS = dbStruct.popData.(['rSpearmanUSUS' num2str(windowDuration) 'minWindows']){iWindow}';
    pvalSpearmanWindowUS = dbStruct.popData.(['pvalSpearmanUSUS' num2str(windowDuration) 'minWindows']){iWindow}';
    if strcmp(filterType{iType}, 'IS')
      rSpearmanWindow = rSpearmanWindowIS;
      pvalSpearmanWindow = pvalSpearmanWindowIS;
      positiveDenominator = numel(rSpearmanWindow);
      negativeDenominator = numel(rSpearmanWindow);
    elseif strcmp(filterType{iType}, 'US')
      rSpearmanWindow = rSpearmanWindowUS;
      pvalSpearmanWindow = pvalSpearmanWindowUS;
      positiveDenominator = numel(rSpearmanWindow);
      negativeDenominator = numel(rSpearmanWindow);
    elseif strcmp(filterType{iType}, 'USIS')
      rSpearmanWindow = NaN(size(rSpearmanWindowIS));
      rSpearmanWindow(rSpearmanWindowIS >= 0 & rSpearmanWindowUS >= 0) = 1;
      rSpearmanWindow(rSpearmanWindowIS < 0 & rSpearmanWindowUS < 0) = -1;
      pvalSpearmanWindow = NaN(size(pvalSpearmanWindowIS));
      pvalSpearmanWindow(pvalSpearmanWindowIS < alpha & pvalSpearmanWindowUS < alpha) = 0;
      pvalSpearmanWindow(~(pvalSpearmanWindowIS < alpha & pvalSpearmanWindowUS < alpha)) = 1;
      positiveDenominator = sum(rSpearmanWindowUS >= 0);
      negativeDenominator = sum(rSpearmanWindowUS < 0);
    end
    data.(['positiveMUAsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['positiveMUAsFr' type]){iCondPlusAll}{iAreaPlusAll};...
      sum(rSpearmanWindow >= 0)/positiveDenominator];
    data.(['positiveMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['positiveMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    data.(['positiveSignificantMUAsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['positiveSignificantMUAsFr' type]){iCondPlusAll}{iAreaPlusAll};...
      sum(rSpearmanWindow >= 0 & pvalSpearmanWindow < alpha)/positiveDenominator];
    data.(['positiveSignificantMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['positiveSignificantMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    data.(['negativeMUAsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['negativeMUAsFr' type]){iCondPlusAll}{iAreaPlusAll};
      sum(rSpearmanWindow < 0)/negativeDenominator];
    data.(['negativeMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['negativeMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    data.(['negativeSignificantMUAsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['negativeSignificantMUAsFr' type]){iCondPlusAll}{iAreaPlusAll};...
      sum(rSpearmanWindow < 0 & pvalSpearmanWindow < alpha)/negativeDenominator];
    data.(['negativeSignificantMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
      [data.(['negativeSignificantMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
      [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    if nUnitsWindows(iWindow) >= unitThr
      type = [filterType{iType} 'Thr' num2str(windowDuration) 'minWindows'];
      data.(['positiveMUAsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['positiveMUAsFr' type]){iCondPlusAll}{iAreaPlusAll};...
        sum(rSpearmanWindow >= 0)/positiveDenominator];
      data.(['positiveMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['positiveMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
        [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
      data.(['positiveSignificantMUAsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['positiveSignificantMUAsFr' type]){iCondPlusAll}{iAreaPlusAll};...
        sum(rSpearmanWindow >= 0 & pvalSpearmanWindow < alpha)/positiveDenominator];
      data.(['positiveSignificantMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['positiveSignificantMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
        [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
      data.(['negativeMUAsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['negativeMUAsFr' type]){iCondPlusAll}{iAreaPlusAll};...
        sum(rSpearmanWindow < 0)/negativeDenominator];
      data.(['negativeMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['negativeMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
        [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
      data.(['negativeSignificantMUAsFr' type]){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['negativeSignificantMUAsFr' type]){iCondPlusAll}{iAreaPlusAll};...
        sum(rSpearmanWindow < 0 & pvalSpearmanWindow < alpha)/negativeDenominator];
      data.(['negativeSignificantMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll} =...
        [data.(['negativeSignificantMUAsFr' type 'RecID']){iCondPlusAll}{iAreaPlusAll};...
        [seriesName(1:min([numel(seriesName) 14])) '_window' num2str(iWindow)]];
    end
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