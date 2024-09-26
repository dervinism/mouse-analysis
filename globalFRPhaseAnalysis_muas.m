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
end

significantUnits = false;
crossValidation = false;
crossValidationType = '50percent';
predictionType = 1; % 1 - binned pre-fit, 2 - smooth post-fit
distroType = 'skewnormal'; % 'gamma' or 'skewnormal'
exclude = true;

drawDistros = true;
drawMeans = [true true];
drawVars = true;
drawMeansVars = true;
drawModulationDistros = true;
drawCorrelations = true;
drawCorrelationsLog = true;
drawCorrelationsLog2 = true;
drawPredictions = true;
drawCombinedDistros = true;
drawCombinedDistrosSplit = true;
drawSmoothCombinedDistros = true;
drawSmoothCombinedDistrosSplit = true;

if ~exclude
  rootFolder = [rootFolder filesep 'allRecs']; %#ok<*UNRCH>
  if strcmp(subpop, 'all')
    rootFolderPositive = [rootFolderPositive filesep 'allRecs'];
    rootFolderNegative = [rootFolderNegative filesep 'allRecs'];
  end
end
rootFolder = [rootFolder filesep MUAsFolder];
if strcmp(subpop, 'all')
  rootFolderPositive = [rootFolderPositive filesep MUAsFolder];
  rootFolderNegative = [rootFolderNegative filesep MUAsFolder];
end
mainFolder = [rootFolder filesep firingRatesSubfolder]; %#ok<*NASGU>
if strcmp(subpop, 'all')
  mainFolderPositive = [rootFolderPositive filesep firingRatesSubfolder];
  mainFolderNegative = [rootFolderNegative filesep firingRatesSubfolder];
end
if significantUnits
  if strcmpi(subpop, 'all')
    error('Significant units can only be analysed for positive or negative subpopulations.');
  end
  rootFolder = [mainFolder filesep significantMUAsFolder];
  mainFolder = rootFolder;
  if strcmp(subpop, 'all')
    rootFolderPositive = [mainFolderPositive filesep significantMUAsFolder];
    mainFolderPositive = rootFolderPositive;
    rootFolderNegative = [mainFolderNegative filesep significantMUAsFolder];
    mainFolderNegative = rootFolderNegative;
  end
end
if crossValidation
  rootFolder = [mainFolder filesep 'crossValidation' crossValidationType];
  mainFolder = rootFolder;
  if strcmp(subpop, 'all')
    rootFolderPositive = [mainFolderPositive filesep 'crossValidation' crossValidationType];
    mainFolderPositive = rootFolderPositive;
    rootFolderNegative = [mainFolderNegative filesep 'crossValidation' crossValidationType];
    mainFolderNegative = rootFolderNegative;
  end
end


%% COMPUTE VARIABLES AND INITIALISE FIGURES NEEDED FOR DISPLAYING FIRING RATE PLOTS
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
    srData = dataStruct.seriesData.(fnsData{1}).conf.samplingParams.srData;
    if animal == 1 || ~exist('areaFRPercentileIndividual_50', 'var')
      arearSpearman5secIndividual = {};
      areapvalSpearman5secIndividual = {};
      areaRecIDSpearman5secIndividual = {};
      areaFRIndividual = {};
      areaFRPercentileIndividual_50 = {};
      areaFRPercentileIndividual_25 = {};
      areaFRPercentileIndividual_12p5 = {};
      areaFRPercentileIndividual_third = {};
      areaFRRiseDecayIndividual = {};
      areaFRFiltLP0p001HzIndividual_50 = {};
      areaFRFiltLP0p001HzIndividual_25 = {};
      areaFRFiltLP0p001HzIndividual_12p5 = {};
      areaFRFiltLP0p001HzIndividual_third = {};
      areaFRFiltLP0p01HzIndividual_50 = {};
      areaFRFiltLP0p1HzIndividual_50 = {};
      areaFRFiltBP0p01to0p05HzIndividual_pi = {};
      areaFRFiltBP0p01to0p05HzIndividual_2piOver3 = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver2 = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver3 = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver4 = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver6 = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver8 = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver12 = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver16 = {};
      areaFRFiltBP0p01to0p05HzIndividual_piShifted = {};
      areaFRFiltBP0p01to0p05HzIndividual_2piOver3Shifted = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver2Shifted = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver3Shifted = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver4Shifted = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver6Shifted = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver8Shifted = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver12Shifted = {};
      areaFRFiltBP0p01to0p05HzIndividual_piOver16Shifted = {};
      areaFRFiltBP0p1to0p5HzIndividual_pi = {};
      areaFRFiltBP0p1to0p5HzIndividual_2piOver3 = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver2 = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver3 = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver4 = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver6 = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver8 = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver12 = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver16 = {};
      areaFRFiltBP0p1to0p5HzIndividual_piShifted = {};
      areaFRFiltBP0p1to0p5HzIndividual_2piOver3Shifted = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver2Shifted = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver3Shifted = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver4Shifted = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver6Shifted = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver8Shifted = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver12Shifted = {};
      areaFRFiltBP0p1to0p5HzIndividual_piOver16Shifted = {};
      for iCond = 1:numel(conditions)
        arearSpearman5secIndividualCond = {};
        areapvalSpearman5secIndividualCond = {};
        areaRecIDSpearman5secIndividualCond = {};
        areaFRIndividualCond = {};
        areaFRPercentileIndividual_50Cond = {};
        areaFRPercentileIndividual_25Cond = {};
        areaFRPercentileIndividual_12p5Cond = {};
        areaFRPercentileIndividual_thirdCond = {};
        areaFRRiseDecayIndividualCond = {};
        areaFRFiltLP0p001HzIndividual_50Cond = {};
        areaFRFiltLP0p001HzIndividual_25Cond = {};
        areaFRFiltLP0p001HzIndividual_12p5Cond = {};
        areaFRFiltLP0p001HzIndividual_thirdCond = {};
        areaFRFiltLP0p01HzIndividual_50Cond = {};
        areaFRFiltLP0p1HzIndividual_50Cond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piCond = {};
        areaFRFiltBP0p01to0p05HzIndividual_2piOver3Cond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver2Cond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver3Cond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver4Cond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver6Cond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver8Cond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver12Cond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver16Cond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piShiftedCond = {};
        areaFRFiltBP0p01to0p05HzIndividual_2piOver3ShiftedCond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver2ShiftedCond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver3ShiftedCond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver4ShiftedCond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver6ShiftedCond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver8ShiftedCond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver12ShiftedCond = {};
        areaFRFiltBP0p01to0p05HzIndividual_piOver16ShiftedCond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piCond = {};
        areaFRFiltBP0p1to0p5HzIndividual_2piOver3Cond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver2Cond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver3Cond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver4Cond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver6Cond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver8Cond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver12Cond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver16Cond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piShiftedCond = {};
        areaFRFiltBP0p1to0p5HzIndividual_2piOver3ShiftedCond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver2ShiftedCond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver3ShiftedCond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver4ShiftedCond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver6ShiftedCond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver8ShiftedCond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver12ShiftedCond = {};
        areaFRFiltBP0p1to0p5HzIndividual_piOver16ShiftedCond = {};
        for iArea = 1:numel(areas)
          arearSpearman5secIndividualCond{iArea} = [];
          areapvalSpearman5secIndividualCond{iArea} = [];
          areaRecIDSpearman5secIndividualCond{iArea} = {};
          areaFRIndividualCond{iArea} = []; %#ok<*SAGROW>
          areaFRPercentileIndividual_50Cond{iArea} = [];
          areaFRPercentileIndividual_25Cond{iArea} = [];
          areaFRPercentileIndividual_12p5Cond{iArea} = [];
          areaFRPercentileIndividual_thirdCond{iArea} = [];
          areaFRRiseDecayIndividualCond{iArea} = [];
          areaFRFiltLP0p001HzIndividual_50Cond{iArea} = [];
          areaFRFiltLP0p001HzIndividual_25Cond{iArea} = [];
          areaFRFiltLP0p001HzIndividual_12p5Cond{iArea} = [];
          areaFRFiltLP0p001HzIndividual_thirdCond{iArea} = [];
          areaFRFiltLP0p01HzIndividual_50Cond{iArea} = [];
          areaFRFiltLP0p1HzIndividual_50Cond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piCond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_2piOver3Cond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver2Cond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver3Cond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver4Cond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver6Cond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver8Cond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver12Cond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver16Cond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piShiftedCond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_2piOver3ShiftedCond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver2ShiftedCond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver3ShiftedCond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver4ShiftedCond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver6ShiftedCond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver8ShiftedCond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver12ShiftedCond{iArea} = [];
          areaFRFiltBP0p01to0p05HzIndividual_piOver16ShiftedCond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piCond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_2piOver3Cond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver2Cond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver3Cond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver4Cond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver6Cond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver8Cond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver12Cond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver16Cond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piShiftedCond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_2piOver3ShiftedCond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver2ShiftedCond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver3ShiftedCond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver4ShiftedCond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver6ShiftedCond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver8ShiftedCond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver12ShiftedCond{iArea} = [];
          areaFRFiltBP0p1to0p5HzIndividual_piOver16ShiftedCond{iArea} = [];
        end
        arearSpearman5secIndividual{iCond} = arearSpearman5secIndividualCond;
        areapvalSpearman5secIndividual{iCond} = areapvalSpearman5secIndividualCond;
        areaRecIDSpearman5secIndividual{iCond} = areaRecIDSpearman5secIndividualCond;
        areaFRIndividual{iCond} = areaFRIndividualCond;
        areaFRPercentileIndividual_50{iCond} = areaFRPercentileIndividual_50Cond;
        areaFRPercentileIndividual_25{iCond} = areaFRPercentileIndividual_25Cond;
        areaFRPercentileIndividual_12p5{iCond} = areaFRPercentileIndividual_12p5Cond;
        areaFRPercentileIndividual_third{iCond} = areaFRPercentileIndividual_thirdCond;
        areaFRRiseDecayIndividual{iCond} = areaFRRiseDecayIndividualCond;
        areaFRFiltLP0p001HzIndividual_50{iCond} = areaFRFiltLP0p001HzIndividual_50Cond;
        areaFRFiltLP0p001HzIndividual_25{iCond} = areaFRFiltLP0p001HzIndividual_25Cond;
        areaFRFiltLP0p001HzIndividual_12p5{iCond} = areaFRFiltLP0p001HzIndividual_12p5Cond;
        areaFRFiltLP0p001HzIndividual_third{iCond} = areaFRFiltLP0p001HzIndividual_thirdCond;
        areaFRFiltLP0p01HzIndividual_50{iCond} = areaFRFiltLP0p01HzIndividual_50Cond;
        areaFRFiltLP0p1HzIndividual_50{iCond} = areaFRFiltLP0p1HzIndividual_50Cond;
        areaFRFiltBP0p01to0p05HzIndividual_pi{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piCond;
        areaFRFiltBP0p01to0p05HzIndividual_2piOver3{iCond} = areaFRFiltBP0p01to0p05HzIndividual_2piOver3Cond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver2{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver2Cond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver3{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver3Cond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver4{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver4Cond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver6{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver6Cond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver8{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver8Cond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver12{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver12Cond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver16{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver16Cond;
        areaFRFiltBP0p01to0p05HzIndividual_piShifted{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piShiftedCond;
        areaFRFiltBP0p01to0p05HzIndividual_2piOver3Shifted{iCond} = areaFRFiltBP0p01to0p05HzIndividual_2piOver3ShiftedCond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver2Shifted{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver2ShiftedCond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver3Shifted{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver3ShiftedCond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver4Shifted{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver4ShiftedCond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver6Shifted{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver6ShiftedCond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver8Shifted{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver8ShiftedCond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver12Shifted{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver12ShiftedCond;
        areaFRFiltBP0p01to0p05HzIndividual_piOver16Shifted{iCond} = areaFRFiltBP0p01to0p05HzIndividual_piOver16ShiftedCond;
        areaFRFiltBP0p1to0p5HzIndividual_pi{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piCond;
        areaFRFiltBP0p1to0p5HzIndividual_2piOver3{iCond} = areaFRFiltBP0p1to0p5HzIndividual_2piOver3Cond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver2{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver2Cond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver3{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver3Cond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver4{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver4Cond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver6{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver6Cond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver8{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver8Cond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver12{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver12Cond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver16{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver16Cond;
        areaFRFiltBP0p1to0p5HzIndividual_piShifted{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piShiftedCond;
        areaFRFiltBP0p1to0p5HzIndividual_2piOver3Shifted{iCond} = areaFRFiltBP0p1to0p5HzIndividual_2piOver3ShiftedCond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver2Shifted{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver2ShiftedCond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver3Shifted{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver3ShiftedCond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver4Shifted{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver4ShiftedCond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver6Shifted{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver6ShiftedCond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver8Shifted{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver8ShiftedCond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver12Shifted{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver12ShiftedCond;
        areaFRFiltBP0p1to0p5HzIndividual_piOver16Shifted{iCond} = areaFRFiltBP0p1to0p5HzIndividual_piOver16ShiftedCond;
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
      if exclude && exceptionTest(except, seriesName)
        continue
      elseif ~exclude && exceptionTest(exceptFR, seriesName)
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
      if ~isfield(dbStruct.popData, 'pupil')
        continue
      end
      units = dbStruct.popData.spkDB_units;
      if crossValidation % not implemented yet to work with the latest +/- splitting based on phase
        if strcmpi(crossValidationType, '10percent')
          rSpearman = dbStruct.popData.rSpearman10percent(1,:)';
          pvalSpearman = dbStruct.popData.pvalSpearman10percent(1,:)';
        elseif strcmpi(crossValidationType, '25percent')
          rSpearman = dbStruct.popData.rSpearman25percent(1,:)';
          pvalSpearman = dbStruct.popData.pvalSpearman25percent(1,:)';
        elseif strcmpi(crossValidationType, '33percent')
          rSpearman = dbStruct.popData.rSpearman33percent(1,:)';
          pvalSpearman = dbStruct.popData.pvalSpearman33percent(1,:)';
        elseif strcmpi(crossValidationType, '50percent')
          rSpearman = dbStruct.popData.rSpearman50percent(1,:)';
          pvalSpearman = dbStruct.popData.pvalSpearman50percent(1,:)';
        elseif strcmpi(crossValidationType, '10percentBack')
          rSpearman = dbStruct.popData.rSpearman10percent(end,:)';
          pvalSpearman = dbStruct.popData.pvalSpearman10percent(end,:)';
        elseif strcmpi(crossValidationType, '25percentBack')
          rSpearman = dbStruct.popData.rSpearman25percent(end,:)';
          pvalSpearman = dbStruct.popData.pvalSpearman25percent(end,:)';
        elseif strcmpi(crossValidationType, '33percentBack')
          rSpearman = dbStruct.popData.rSpearman33percent(end,:)';
          pvalSpearman = dbStruct.popData.pvalSpearman33percent(end,:)';
        elseif strcmpi(crossValidationType, '50percentBack')
          rSpearman = dbStruct.popData.rSpearman50percent(end,:)';
          pvalSpearman = dbStruct.popData.pvalSpearman50percent(end,:)';
        end
      else
        if iCond == 1
          if ~isfield(dbStruct.popData, 'pupil')
            continue
          else
            phase = spkPhase(dbStruct.popData.pupil.phaseCoh.unitData, fRef)';
            phaseSignificant = ~isnan(spkPhaseSignificant(dbStruct.popData.pupil.phaseCoh.unitData, fRef)');
          end
        elseif iCond == 2
          if ~isfield(dbStruct.popData, 'rSpearman')
            continue
          else
            rSpearman = dbStruct.popData.rSpearman';
            pvalSpearman = dbStruct.popData.pvalSpearman';
          end
        end
      end
      if significantUnits
        if iCond == 1
          significantInd = phaseSignificant;
        elseif iCond == 2
          significantInd = pvalSpearman < alpha;
        end
      else
        if iCond == 1
          significantInd = true(size(phaseSignificant));
        elseif iCond == 2
          significantInd = pvalSpearman <= 1;
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
      
      if numel(conditions) > 1
        condLoop = [iCond numel(conditions)];
      else
        condLoop = iCond;
      end
      for iCondPlusAll = condLoop % Loop through the main and pooled conditions
        for iAreaPlusAll = area % Loop through the main and pooled areas
          close all
          
          % Get unit spiking pupil correlation data
          spk = full(dbStruct.popData.spkDB);
          if iCond == 1
            if ~isfield(dbStruct.popData,'pupil')
              continue
            end
            phase = spkPhase(dbStruct.popData.pupil.phaseCoh.unitData, fRef)';
            negativeInd = false(size(phase));
            modeBoundaries = recentrePhase(modeBoundaries, modeBoundaries(3));
            phase = recentrePhase(phase, modeBoundaries(3));
            negativeInd(phase > modeBoundaries(2) & phase <= modeBoundaries(end)) = true;
            rSpearman = spkCoh(dbStruct.popData.pupil.phaseCoh.unitData, fRef);
            rSpearman(negativeInd) = -rSpearman(negativeInd); % signed coherence rather than Spearman rho, so remember not to confuse
            pvalSpearman = ~isnan(spkPhaseSignificant(dbStruct.popData.pupil.phaseCoh.unitData, fRef)');
          elseif iCond == 2
            if ~isfield(dbStruct.popData,'rSpearman')...
                || isempty(dbStruct.popData.rSpearman)
              continue
            end
            rSpearman = dbStruct.popData.rSpearman';
            pvalSpearman = dbStruct.popData.pvalSpearman';
          end
          arearSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll} = [arearSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll}; torow(rSpearman(correlatedInd))'];
          areapvalSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll} = [areapvalSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll}; torow(pvalSpearman(correlatedInd))'];
          recIDs = [];
          [recIDs{1:numel(rSpearman(correlatedInd))}] = deal(fnsData{dbCount});
          areaRecIDSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll} = [areaRecIDSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll}; recIDs'];
          period = dbStruct.db(dbCount).period;
          assert(numel(recIDs) == numel(rSpearman(correlatedInd)));
          assert(numel(areaRecIDSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll}) == numel(arearSpearman5secIndividual{iCondPlusAll}{iAreaPlusAll}));
          
          % Get eye data
          eyeData = dataStruct.eyeData.([animals{animal} '_s' seriesName(1:min([14 numel(seriesName)]))]);
          eyeData.pupilArea = double(eyeData.pupilAreaFilt.pupilAreaFiltHP0p01Hz);
          eyeData.frameTimes = eyeData.pupilAreaFilt.timesFiltStart:eyeData.pupilAreaFilt.timesFiltStep:eyeData.pupilAreaFilt.timesFiltStop;
          
          % Combine periods
          if exclude || (~exclude && ismember(animals{animal},animalsNeuronexusUOL))
            commonPeriod = combinePeriods(period, eyeData.period, srData);
          else
            commonPeriod = combinePeriods(eyeData.period, eyeData.period, srData);
          end
          
          % Interpolate and filter pupil area data
          if isempty(commonPeriod)
            continue
          end
          [pupilArea, interpTimes, interpInds] = pupilFilt(eyeData, srData, sum(spk,1), 0, commonPeriod, srData);
          if crossValidation
            if strcmpi(crossValidationType, '10percent') || strcmpi(crossValidationType, '10percentBack')
              lengthFractionalDuration = round(0.1*numel(pupilArea));
            elseif strcmpi(crossValidationType, '25percent') || strcmpi(crossValidationType, '25percentBack')
              lengthFractionalDuration = round(0.25*numel(pupilArea));
            elseif strcmpi(crossValidationType, '33percent') || strcmpi(crossValidationType, '33percentBack')
              lengthFractionalDuration = round((1/3)*numel(pupilArea));
            elseif strcmpi(crossValidationType, '50percent') || strcmpi(crossValidationType, '50percentBack')
              lengthFractionalDuration = round(0.5*numel(pupilArea));
            else
              error('Only 10, 25, 33, and 50 percent forward and backward cross-validation types are supported.');
            end
            if strcmpi(crossValidationType, '10percent') || strcmpi(crossValidationType, '25percent') || strcmpi(crossValidationType, '33percent') || strcmpi(crossValidationType, '50percent')
              pupilArea = pupilArea(lengthFractionalDuration+1:end);
              interpTimes = interpTimes(lengthFractionalDuration+1:end);
              interpInds = interpInds(lengthFractionalDuration+1:end);
            elseif strcmpi(crossValidationType, '10percentBack') || strcmpi(crossValidationType, '25percentBack') || strcmpi(crossValidationType, '33percentBack') || strcmpi(crossValidationType, '50percentBack')
              pupilArea = pupilArea(1:end-lengthFractionalDuration);
              interpTimes = interpTimes(1:end-lengthFractionalDuration);
              interpInds = interpInds(1:end-lengthFractionalDuration);
            end
          end
          pupilAreaDiff = [diff(pupilArea) 0];
          pupilAreaFiltTimes = eyeData.pupilAreaFilt.timesFiltStart:eyeData.pupilAreaFilt.timesFiltStep:eyeData.pupilAreaFilt.timesFiltStop;
          pupilAreaFiltLP0p001Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.pupilAreaFiltLP0p001Hz, interpTimes);
          pupilAreaFiltLP0p01Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.pupilAreaFiltLP0p01Hz, interpTimes);
          pupilAreaFiltLP0p1Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.pupilAreaFiltLP0p1Hz, interpTimes);
          pupilAreaFiltBP0p01to0p05Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.pupilAreaFiltBP0p01to0p05Hz, interpTimes);
          pupilAreaFiltBP0p1to0p5Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.pupilAreaFiltBP0p1to0p5Hz, interpTimes);
          phaseFiltBP0p01to0p05Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.phaseFiltBP0p01to0p05Hz, interpTimes);
          phaseFiltBP0p1to0p5Hz = interp1(pupilAreaFiltTimes, eyeData.pupilAreaFilt.phaseFiltBP0p1to0p5Hz, interpTimes);
          
          % Truncate spiking data
          spk = spk(correlatedInd, interpInds);
          
          % Calculate firing rates
          % Full rates
          areaFRIndividual{iCondPlusAll}{iAreaPlusAll} = [areaFRIndividual{iCondPlusAll}{iAreaPlusAll}; (sum(spk(:,:),2).*srData)./size(spk(:,:),2)];
          
          % Percentiles
          percentiles = [12.5 25 100/3 37.5 50 62.5 200/3 75 87.5];
          percentileValues = prctile(pupilArea,percentiles);
          
          areaFRPercentileIndividual_50DB = zeros(size(spk,1),2);
          areaFRPercentileIndividual_50DB(:,1) = (sum(spk(:,pupilArea<=percentileValues(5)),2).*srData)./size(spk(:,pupilArea<=percentileValues(5)),2);
          areaFRPercentileIndividual_50DB(:,2) = (sum(spk(:,pupilArea>percentileValues(5)),2).*srData)./size(spk(:,pupilArea>percentileValues(5)),2);
          areaFRPercentileIndividual_50{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_50{iCondPlusAll}{iAreaPlusAll}; areaFRPercentileIndividual_50DB];
          
          areaFRPercentileIndividual_thirdDB = zeros(size(spk,1),3);
          areaFRPercentileIndividual_thirdDB(:,1) = (sum(spk(:,pupilArea<=percentileValues(3)),2).*srData)./size(spk(:,pupilArea<=percentileValues(3)),2);
          areaFRPercentileIndividual_thirdDB(:,2) = (sum(spk(:,pupilArea>percentileValues(3) & pupilArea<=percentileValues(7)),2).*srData)./size(spk(:,pupilArea>percentileValues(3) & pupilArea<=percentileValues(7)),2);
          areaFRPercentileIndividual_thirdDB(:,3) = (sum(spk(:,pupilArea>percentileValues(7)),2).*srData)./size(spk(:,pupilArea>percentileValues(7)),2);
          areaFRPercentileIndividual_third{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_third{iCondPlusAll}{iAreaPlusAll}; areaFRPercentileIndividual_thirdDB];
          
          areaFRPercentileIndividual_25DB = zeros(size(spk,1),4);
          areaFRPercentileIndividual_25DB(:,1) = (sum(spk(:,pupilArea<=percentileValues(2)),2).*srData)./size(spk(:,pupilArea<=percentileValues(2)),2);
          areaFRPercentileIndividual_25DB(:,2) = (sum(spk(:,pupilArea>percentileValues(2) & pupilArea<=percentileValues(5)),2).*srData)./size(spk(:,pupilArea>percentileValues(2) & pupilArea<=percentileValues(5)),2);
          areaFRPercentileIndividual_25DB(:,3) = (sum(spk(:,pupilArea>percentileValues(5) & pupilArea<=percentileValues(8)),2).*srData)./size(spk(:,pupilArea>percentileValues(5) & pupilArea<=percentileValues(8)),2);
          areaFRPercentileIndividual_25DB(:,4) = (sum(spk(:,pupilArea>percentileValues(8)),2).*srData)./size(spk(:,pupilArea>percentileValues(8)),2);
          areaFRPercentileIndividual_25{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_25{iCondPlusAll}{iAreaPlusAll}; areaFRPercentileIndividual_25DB];
          
          areaFRPercentileIndividual_12p5DB = zeros(size(spk,1),8);
          areaFRPercentileIndividual_12p5DB(:,1) = (sum(spk(:,pupilArea<=percentileValues(1)),2).*srData)./size(spk(:,pupilArea<=percentileValues(1)),2);
          areaFRPercentileIndividual_12p5DB(:,2) = (sum(spk(:,pupilArea>percentileValues(1) & pupilArea<=percentileValues(2)),2).*srData)./size(spk(:,pupilArea>percentileValues(1) & pupilArea<=percentileValues(2)),2);
          areaFRPercentileIndividual_12p5DB(:,3) = (sum(spk(:,pupilArea>percentileValues(2) & pupilArea<=percentileValues(4)),2).*srData)./size(spk(:,pupilArea>percentileValues(2) & pupilArea<=percentileValues(4)),2);
          areaFRPercentileIndividual_12p5DB(:,4) = (sum(spk(:,pupilArea>percentileValues(4) & pupilArea<=percentileValues(5)),2).*srData)./size(spk(:,pupilArea>percentileValues(4) & pupilArea<=percentileValues(5)),2);
          areaFRPercentileIndividual_12p5DB(:,5) = (sum(spk(:,pupilArea>percentileValues(5) & pupilArea<=percentileValues(6)),2).*srData)./size(spk(:,pupilArea>percentileValues(5) & pupilArea<=percentileValues(6)),2);
          areaFRPercentileIndividual_12p5DB(:,6) = (sum(spk(:,pupilArea>percentileValues(6) & pupilArea<=percentileValues(8)),2).*srData)./size(spk(:,pupilArea>percentileValues(6) & pupilArea<=percentileValues(8)),2);
          areaFRPercentileIndividual_12p5DB(:,7) = (sum(spk(:,pupilArea>percentileValues(8) & pupilArea<=percentileValues(9)),2).*srData)./size(spk(:,pupilArea>percentileValues(8) & pupilArea<=percentileValues(9)),2);
          areaFRPercentileIndividual_12p5DB(:,8) = (sum(spk(:,pupilArea>percentileValues(9)),2).*srData)./size(spk(:,pupilArea>percentileValues(9)),2);
          areaFRPercentileIndividual_12p5{iCondPlusAll}{iAreaPlusAll} = [areaFRPercentileIndividual_12p5{iCondPlusAll}{iAreaPlusAll}; areaFRPercentileIndividual_12p5DB];
          
          % Rise and decay
          areaFRRiseDecayIndividualDB = zeros(size(spk,1),2);
          areaFRRiseDecayIndividualDB(:,1) = (sum(spk(:,pupilAreaDiff>0),2).*srData)./size(spk(:,pupilAreaDiff>0),2);
          areaFRRiseDecayIndividualDB(:,2) = (sum(spk(:,pupilAreaDiff<0),2).*srData)./size(spk(:,pupilAreaDiff<0),2);
          areaFRRiseDecayIndividual{iCondPlusAll}{iAreaPlusAll} = [areaFRRiseDecayIndividual{iCondPlusAll}{iAreaPlusAll}; areaFRRiseDecayIndividualDB];
          
          % Filtered percentiles
          percentileValues = percentileValues - percentileValues(5);
          percentileValuesFilt = zeros(numel(percentileValues),numel(pupilAreaFiltLP0p001Hz));
          for iCent = 1:numel(percentileValues)
            percentileValuesFilt(iCent,:) = pupilAreaFiltLP0p001Hz + percentileValues(iCent);
          end
          
          areaFRFiltLP0p001HzIndividual_50DB = zeros(size(spk,1),2);
          areaFRFiltLP0p001HzIndividual_50DB(:,1) = (sum(spk(:,pupilArea<=percentileValuesFilt(5,:)),2).*srData)./size(spk(:,pupilArea<=percentileValuesFilt(5,:)),2);
          areaFRFiltLP0p001HzIndividual_50DB(:,2) = (sum(spk(:,pupilArea>percentileValuesFilt(5,:)),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt(5,:)),2);
          areaFRFiltLP0p001HzIndividual_50{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltLP0p001HzIndividual_50{iCondPlusAll}{iAreaPlusAll}; areaFRFiltLP0p001HzIndividual_50DB];
          
          areaFRFiltLP0p001HzIndividual_thirdDB = zeros(size(spk,1),3);
          areaFRFiltLP0p001HzIndividual_thirdDB(:,1) = (sum(spk(:,pupilArea<=percentileValuesFilt(3,:)),2).*srData)./size(spk(:,pupilArea<=percentileValuesFilt(3,:)),2);
          areaFRFiltLP0p001HzIndividual_thirdDB(:,2) = (sum(spk(:,pupilArea>percentileValuesFilt(3,:) & pupilArea<=percentileValuesFilt(7,:)),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt(3,:) & pupilArea<=percentileValuesFilt(7,:)),2);
          areaFRFiltLP0p001HzIndividual_thirdDB(:,3) = (sum(spk(:,pupilArea>percentileValuesFilt(7,:)),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt(7,:)),2);
          areaFRFiltLP0p001HzIndividual_third{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltLP0p001HzIndividual_third{iCondPlusAll}{iAreaPlusAll}; areaFRFiltLP0p001HzIndividual_thirdDB];
          
          areaFRFiltLP0p001HzIndividual_25DB = zeros(size(spk,1),4);
          areaFRFiltLP0p001HzIndividual_25DB(:,1) = (sum(spk(:,pupilArea<=percentileValuesFilt(2,:)),2).*srData)./size(spk(:,pupilArea<=percentileValuesFilt(2,:)),2);
          areaFRFiltLP0p001HzIndividual_25DB(:,2) = (sum(spk(:,pupilArea>percentileValuesFilt(2,:) & pupilArea<=percentileValuesFilt(5,:)),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt(2,:) & pupilArea<=percentileValuesFilt(5,:)),2);
          areaFRFiltLP0p001HzIndividual_25DB(:,3) = (sum(spk(:,pupilArea>percentileValuesFilt(5,:) & pupilArea<=percentileValuesFilt(8,:)),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt(5,:) & pupilArea<=percentileValuesFilt(8,:)),2);
          areaFRFiltLP0p001HzIndividual_25DB(:,4) = (sum(spk(:,pupilArea>percentileValuesFilt(8,:)),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt(8,:)),2);
          areaFRFiltLP0p001HzIndividual_25{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltLP0p001HzIndividual_25{iCondPlusAll}{iAreaPlusAll}; areaFRFiltLP0p001HzIndividual_25DB];
          
          areaFRFiltLP0p001HzIndividual_12p5DB = zeros(size(spk,1),8);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,1) = (sum(spk(:,pupilArea<=percentileValuesFilt(1,:)),2).*srData)./size(spk(:,pupilArea<=percentileValuesFilt(1,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,2) = (sum(spk(:,pupilArea>percentileValuesFilt(1,:) & pupilArea<=percentileValuesFilt(2,:)),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt(1,:) & pupilArea<=percentileValuesFilt(2,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,3) = (sum(spk(:,pupilArea>percentileValuesFilt(2,:) & pupilArea<=percentileValuesFilt(4,:)),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt(2,:) & pupilArea<=percentileValuesFilt(4,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,4) = (sum(spk(:,pupilArea>percentileValuesFilt(4,:) & pupilArea<=percentileValuesFilt(5,:)),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt(4,:) & pupilArea<=percentileValuesFilt(5,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,5) = (sum(spk(:,pupilArea>percentileValuesFilt(5,:) & pupilArea<=percentileValuesFilt(6,:)),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt(5,:) & pupilArea<=percentileValuesFilt(6,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,6) = (sum(spk(:,pupilArea>percentileValuesFilt(6,:) & pupilArea<=percentileValuesFilt(8,:)),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt(6,:) & pupilArea<=percentileValuesFilt(8,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,7) = (sum(spk(:,pupilArea>percentileValuesFilt(8,:) & pupilArea<=percentileValuesFilt(9,:)),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt(8,:) & pupilArea<=percentileValuesFilt(9,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5DB(:,8) = (sum(spk(:,pupilArea>percentileValuesFilt(9,:)),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt(9,:)),2);
          areaFRFiltLP0p001HzIndividual_12p5{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltLP0p001HzIndividual_12p5{iCondPlusAll}{iAreaPlusAll}; areaFRFiltLP0p001HzIndividual_12p5DB];
          
          percentileValuesFilt = pupilAreaFiltLP0p01Hz;
          areaFRFiltLP0p01HzIndividual_50DB = zeros(size(spk,1),2);
          areaFRFiltLP0p01HzIndividual_50DB(:,1) = (sum(spk(:,pupilArea<=percentileValuesFilt),2).*srData)./size(spk(:,pupilArea<=percentileValuesFilt),2);
          areaFRFiltLP0p01HzIndividual_50DB(:,2) = (sum(spk(:,pupilArea>percentileValuesFilt),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt),2);
          areaFRFiltLP0p01HzIndividual_50{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltLP0p01HzIndividual_50{iCondPlusAll}{iAreaPlusAll}; areaFRFiltLP0p01HzIndividual_50DB];
          
          percentileValuesFilt = pupilAreaFiltLP0p1Hz;
          areaFRFiltLP0p1HzIndividual_50DB = zeros(size(spk,1),2);
          areaFRFiltLP0p1HzIndividual_50DB(:,1) = (sum(spk(:,pupilArea<=percentileValuesFilt),2).*srData)./size(spk(:,pupilArea<=percentileValuesFilt),2);
          areaFRFiltLP0p1HzIndividual_50DB(:,2) = (sum(spk(:,pupilArea>percentileValuesFilt),2).*srData)./size(spk(:,pupilArea>percentileValuesFilt),2);
          areaFRFiltLP0p1HzIndividual_50{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltLP0p1HzIndividual_50{iCondPlusAll}{iAreaPlusAll}; areaFRFiltLP0p1HzIndividual_50DB];
          
          % Phases
          areaFRFiltBP0p01to0p05HzIndividual_pi{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_pi{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi, 0)];
          areaFRFiltBP0p01to0p05HzIndividual_2piOver3{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_2piOver3{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, 2*pi/3, 0)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver2{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver2{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/2, 0)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver3{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver3{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/3, 0)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver4{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver4{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/4, 0)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver6{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver6{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/6, 0)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver8{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver8{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/8, 0)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver12{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver12{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/12, 0)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver16{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver16{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/16, 0)];
          
          areaFRFiltBP0p1to0p5HzIndividual_pi{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_pi{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi, 0)];
          areaFRFiltBP0p1to0p5HzIndividual_2piOver3{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_2piOver3{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, 2*pi/3, 0)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver2{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver2{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/2, 0)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver3{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver3{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/3, 0)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver4{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver4{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/4, 0)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver6{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver6{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/6, 0)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver8{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver8{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/8, 0)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver12{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver12{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/12, 0)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver16{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver16{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/16, 0)];
          
          % Shifted phases
          areaFRFiltBP0p01to0p05HzIndividual_piShifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piShifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi, pi/2)];
          areaFRFiltBP0p01to0p05HzIndividual_2piOver3Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_2piOver3Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, 2*pi/3, pi/3)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver2Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver2Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/2, pi/4)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver3Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver3Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/3, pi/6)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver4Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver4Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/4, pi/8)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver6Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver6Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/6, pi/12)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver8Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver8Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/8, pi/16)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver12Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver12Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/12, pi/24)];
          areaFRFiltBP0p01to0p05HzIndividual_piOver16Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p01to0p05HzIndividual_piOver16Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p01to0p05Hz, pi/16, pi/32)];
          
          areaFRFiltBP0p1to0p5HzIndividual_piShifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piShifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi, pi/2)];
          areaFRFiltBP0p1to0p5HzIndividual_2piOver3Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_2piOver3Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, 2*pi/3, pi/3)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver2Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver2Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/2, pi/4)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver3Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver3Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/3, pi/6)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver4Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver4Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/4, pi/8)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver6Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver6Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/6, pi/12)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver8Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver8Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/8, pi/16)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver12Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver12Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/12, pi/24)];
          areaFRFiltBP0p1to0p5HzIndividual_piOver16Shifted{iCondPlusAll}{iAreaPlusAll} = [areaFRFiltBP0p1to0p5HzIndividual_piOver16Shifted{iCondPlusAll}{iAreaPlusAll}; firingRateByPhase(spk, srData, phaseFiltBP0p1to0p5Hz, pi/16, pi/32)];
        end
      end
    end
  end
end

% Determine the file name and either save or load the data
filename = [rootFolder filesep 'globalMUAs.mat'];
if strcmp(subpop, 'all')
  filenamePositive = [rootFolderPositive filesep 'globalMUAs.mat'];
  filenameNegative = [rootFolderNegative filesep 'globalMUAs.mat'];
end
if fullRun
  if ~exist(rootFolder, 'file')
    mkdir(rootFolder);
  end
  if ~exist(mainFolder, 'file')
    mkdir(mainFolder);
  end
  if ~exist(filename, 'file')
    save(filename, 'conditions','areas','arearSpearman5secIndividual','areapvalSpearman5secIndividual','areaRecIDSpearman5secIndividual','areaFRIndividual',...
      'areaFRPercentileIndividual_50','areaFRPercentileIndividual_25','areaFRPercentileIndividual_12p5','areaFRPercentileIndividual_third','areaFRRiseDecayIndividual',...
      'areaFRFiltLP0p001HzIndividual_50','areaFRFiltLP0p001HzIndividual_25','areaFRFiltLP0p001HzIndividual_12p5','areaFRFiltLP0p001HzIndividual_third',...
      'areaFRFiltLP0p01HzIndividual_50','areaFRFiltLP0p1HzIndividual_50','areaFRFiltBP0p01to0p05HzIndividual_pi','areaFRFiltBP0p01to0p05HzIndividual_piOver2','areaFRFiltBP0p01to0p05HzIndividual_2piOver3',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver3','areaFRFiltBP0p01to0p05HzIndividual_piOver4','areaFRFiltBP0p01to0p05HzIndividual_piOver6',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver8','areaFRFiltBP0p01to0p05HzIndividual_piOver12','areaFRFiltBP0p01to0p05HzIndividual_piOver16',...
      'areaFRFiltBP0p01to0p05HzIndividual_piShifted','areaFRFiltBP0p01to0p05HzIndividual_2piOver3Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver2Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver3Shifted',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver4Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver6Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver8Shifted',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver12Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver16Shifted','areaFRFiltBP0p1to0p5HzIndividual_pi','areaFRFiltBP0p1to0p5HzIndividual_2piOver3',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver2','areaFRFiltBP0p1to0p5HzIndividual_piOver3','areaFRFiltBP0p1to0p5HzIndividual_piOver4',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver6','areaFRFiltBP0p1to0p5HzIndividual_piOver8','areaFRFiltBP0p1to0p5HzIndividual_piOver12',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver16','areaFRFiltBP0p1to0p5HzIndividual_piShifted','areaFRFiltBP0p1to0p5HzIndividual_piOver2Shifted','areaFRFiltBP0p1to0p5HzIndividual_2piOver3Shifted',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver3Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver4Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver6Shifted',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver8Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver12Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver16Shifted');
  else
    save(filename, 'conditions','areas','arearSpearman5secIndividual','areapvalSpearman5secIndividual','areaRecIDSpearman5secIndividual','areaFRIndividual',...
      'areaFRPercentileIndividual_50','areaFRPercentileIndividual_25','areaFRPercentileIndividual_12p5','areaFRPercentileIndividual_third','areaFRRiseDecayIndividual',...
      'areaFRFiltLP0p001HzIndividual_50','areaFRFiltLP0p001HzIndividual_25','areaFRFiltLP0p001HzIndividual_12p5','areaFRFiltLP0p001HzIndividual_third',...
      'areaFRFiltLP0p01HzIndividual_50','areaFRFiltLP0p1HzIndividual_50','areaFRFiltBP0p01to0p05HzIndividual_pi','areaFRFiltBP0p01to0p05HzIndividual_piOver2','areaFRFiltBP0p01to0p05HzIndividual_2piOver3',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver3','areaFRFiltBP0p01to0p05HzIndividual_piOver4','areaFRFiltBP0p01to0p05HzIndividual_piOver6',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver8','areaFRFiltBP0p01to0p05HzIndividual_piOver12','areaFRFiltBP0p01to0p05HzIndividual_piOver16',...
      'areaFRFiltBP0p01to0p05HzIndividual_piShifted','areaFRFiltBP0p01to0p05HzIndividual_2piOver3Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver2Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver3Shifted',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver4Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver6Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver8Shifted',...
      'areaFRFiltBP0p01to0p05HzIndividual_piOver12Shifted','areaFRFiltBP0p01to0p05HzIndividual_piOver16Shifted','areaFRFiltBP0p1to0p5HzIndividual_pi','areaFRFiltBP0p1to0p5HzIndividual_2piOver3',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver2','areaFRFiltBP0p1to0p5HzIndividual_piOver3','areaFRFiltBP0p1to0p5HzIndividual_piOver4',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver6','areaFRFiltBP0p1to0p5HzIndividual_piOver8','areaFRFiltBP0p1to0p5HzIndividual_piOver12',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver16','areaFRFiltBP0p1to0p5HzIndividual_piShifted','areaFRFiltBP0p1to0p5HzIndividual_piOver2Shifted','areaFRFiltBP0p1to0p5HzIndividual_2piOver3Shifted',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver3Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver4Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver6Shifted',...
      'areaFRFiltBP0p1to0p5HzIndividual_piOver8Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver12Shifted','areaFRFiltBP0p1to0p5HzIndividual_piOver16Shifted', '-append');
  end
else
  load(filename);
end


%% PLOT THE LOG FIRING RATE DISTRIBUTIONS
if drawDistros
  edges = -3:0.2:2;
  options = struct();
  options.xLabel = 'Log_{10}(firing rate)';
  if strcmpi(subpop, 'all')
    options.yLabel = 'Probability';
  elseif strcmpi(subpop, 'positive')
    %options.yLabel = 'Positive unit count';
    options.yLabel = 'Probability';
  elseif strcmpi(subpop, 'negative')
    %options.yLabel = 'Negative unit count';
    options.yLabel = 'Probability';
  end
  options.figFolder = mainFolder;
  options.figSize = 18;
  options.saveFig = true;
  options.pdf = true;
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      options.figName = ['FiringRateCentile50_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      options.legendLabels = {'0:50 cent','50:100 cent'};
      if strcmpi(subpop, 'all')
        options.figName2 = ['FiringRateCentile50Split_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        [statsMeanFRPercentile50{iCond}{iAreasOI(iArea)}, statsVarFRPercentile50{iCond}{iAreasOI(iArea)}, statsMeanFRPercentile50Split{iCond}{iAreasOI(iArea)},...
          statsVarFRPercentile50Split{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)},...
          arearSpearman5secIndividual{iCond}{iAreasOI(iArea)}, options);
      else
        [statsMeanFRPercentile50{iCond}{iAreasOI(iArea)}, statsVarFRPercentile50{iCond}{iAreasOI(iArea)}, statsMeanFRPercentile50Split{iCond}{iAreasOI(iArea)},...
          statsVarFRPercentile50Split{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}, [], options);
      end
      
      options.figName = ['FiringRateCentile25_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      options.legendLabels = {'0:25 cent','25:50 cent','50:75 cent','75:100 cent'};
      if strcmpi(subpop, 'all')
        options.figName2 = ['FiringRateCentile25Split_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        [statsMeanFRPercentile25{iCond}{iAreasOI(iArea)}, statsVarFRPercentile25{iCond}{iAreasOI(iArea)}, statsMeanFRPercentile25Split{iCond}{iAreasOI(iArea)},...
          statsVarFRPercentile25Split{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)},...
          arearSpearman5secIndividual{iCond}{iAreasOI(iArea)}, options);
      else
        [statsMeanFRPercentile25{iCond}{iAreasOI(iArea)}, statsVarFRPercentile25{iCond}{iAreasOI(iArea)}, statsMeanFRPercentile25Split{iCond}{iAreasOI(iArea)},...
          statsVarFRPercentile25Split{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}, [], options);
      end
      
      options.figName = ['FiringRateCentile12p5_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      options.legendLabels = {'0:12.5 cent','12.5:25 cent','25:37.5 cent','37.5:50 cent','50:62.5 cent','62.5:75 cent','75:87.5 cent','87.5:100 cent'};
      if strcmpi(subpop, 'all')
        options.figName2 = ['FiringRateCentile12p5Split_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        [statsMeanFRPercentile12p5{iCond}{iAreasOI(iArea)}, statsVarFRPercentile12p5{iCond}{iAreasOI(iArea)}, statsMeanFRPercentile12p5Split{iCond}{iAreasOI(iArea)},...
          statsVarFRPercentile12p5Split{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)},...
          arearSpearman5secIndividual{iCond}{iAreasOI(iArea)}, options);
      else
        [statsMeanFRPercentile12p5{iCond}{iAreasOI(iArea)}, statsVarFRPercentile12p5{iCond}{iAreasOI(iArea)}, statsMeanFRPercentile12p5Split{iCond}{iAreasOI(iArea)},...
          statsVarFRPercentile12p5Split{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}, [], options);
      end
      
      options.figName = ['FiringRateCentileThird_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      options.legendLabels = {'0:33.3 cent','33.3:66.7 cent','67.7:100 cent'};
      if strcmpi(subpop, 'all')
        options.figName2 = ['FiringRateCentileThirdSplit_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        [statsMeanFRPercentileThird{iCond}{iAreasOI(iArea)}, statsVarFRPercentileThird{iCond}{iAreasOI(iArea)}, statsMeanFRPercentileThirdSplit{iCond}{iAreasOI(iArea)},...
          statsVarFRPercentileThirdSplit{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)},...
          arearSpearman5secIndividual{iCond}{iAreasOI(iArea)}, options);
      else
        [statsMeanFRPercentileThird{iCond}{iAreasOI(iArea)}, statsVarFRPercentileThird{iCond}{iAreasOI(iArea)}, statsMeanFRPercentileThirdSplit{iCond}{iAreasOI(iArea)},...
          statsVarFRPercentileThirdSplit{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}, [], options);
      end
      
      options.plotType = 'extremes';
      if ~isempty(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)})
        options.figName = ['FiringRateCentile25Grey_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        options.legendLabels = {'0:25 cent','75:100 cent'};
        logRates = getLog(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)});
        logRates(isinf(logRates)) = -9;
        options.stats = varTest(logRates(:,[1 end]));
        histPlotFR(edges, logRates, options);
      end
      
      if ~isempty(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)})
        options.figName = ['FiringRateCentile12p5Grey_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        options.legendLabels = {'0:12.5 cent','87.5:100 cent'};
        logRates = getLog(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)});
        logRates(isinf(logRates)) = -9;
        options.stats = varTest(logRates(:,[1 end]));
        histPlotFR(edges, logRates, options);
      end
      
      if ~isempty(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)})
        options.figName = ['FiringRateCentileThirdGrey_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        options.legendLabels = {'0:33.3 cent','67.7:100 cent'};
        logRates = getLog(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)});
        logRates(isinf(logRates)) = -9;
        options.stats = varTest(logRates(:,[1 end]));
        histPlotFR(edges, logRates, options);
      end
      
      options.plotType = 'regular';
      options.figName = ['FiringRateRiseDecay_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      options.legendLabels = {'Dilating','Constricting'};
      [statsMeanFRRiseDecay{iCond}{iAreasOI(iArea)}, statsVarFRRiseDecay{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRRiseDecayIndividual{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltLP0p001Hz50_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      options.legendLabels = {'0:50 cent','50:100 cent'};
      [statsMeanFRFiltLP0p001Hz50{iCond}{iAreasOI(iArea)}, statsVarFRFiltLP0p001Hz50{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltLP0p001HzIndividual_50{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltLP0p001Hz25_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      options.legendLabels = {'0:25 cent','25:50 cent','50:75 cent','75:100 cent'};
      [statsMeanFRFiltLP0p001Hz25{iCond}{iAreasOI(iArea)}, statsVarFRFiltLP0p001Hz25{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltLP0p001HzIndividual_25{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltLP0p001Hz12p5_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      options.legendLabels = {'0:12.5 cent','12.5:25 cent','25:37.5 cent','37.5:50 cent','50:62.5 cent','62.5:75 cent','75:87.5 cent','87.5:100 cent'};
      [statsMeanFRFiltLP0p001Hz12p5{iCond}{iAreasOI(iArea)}, statsVarFRFiltLP0p001Hz12p5{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltLP0p001HzIndividual_12p5{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltLP0p001HzThird_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      options.legendLabels = {'0:33.3 cent','33.3:66.7 cent','67.7:100 cent'};
      [statsMeanFRFiltLP0p001HzThird{iCond}{iAreasOI(iArea)}, statsVarFRFiltLP0p001HzThird{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltLP0p001HzIndividual_third{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltLP0p01Hz50_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      options.legendLabels = {'0:50 cent','50:100 cent'};
      [statsMeanFRFiltLP0p01Hz50{iCond}{iAreasOI(iArea)}, statsVarFRFiltLP0p01Hz50{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltLP0p01HzIndividual_50{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltLP0p1Hz50_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      options.legendLabels = {'0:50 cent','50:100 cent'};
      [statsMeanFRFiltLP0p1Hz50{iCond}{iAreasOI(iArea)}, statsVarFRFiltLP0p1Hz50{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltLP0p1HzIndividual_50{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05HzPi_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:\pi','\pi:2\pi'};
      options.legendLabels = {'\pi:0','0:-\pi'};
      [statsMeanFRFiltBP0p01to0p05HzPi{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05HzPi{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_pi{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05Hz2PiOver3_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:2\pi/3','2\pi/3:4\pi/3','4\pi/3:2\pi'};
      options.legendLabels = {'\pi:\pi/3','\pi/3:-\pi/3','-\pi/3:-\pi'};
      [statsMeanFRFiltBP0p01to0p05Hz2piOver3{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05Hz2piOver3{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_2piOver3{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05HzPiOver2_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:\pi/2','\pi/2:\pi','\pi:3\pi/2','3\pi/2:2\pi'};
      options.legendLabels = {'\pi:\pi/2','\pi/2:0','0:-\pi/2','-\pi/2:-pi'};
      [statsMeanFRFiltBP0p01to0p05HzPiOver2{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05HzPiOver2{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_piOver2{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05HzPiOver3_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:\pi/3','\pi/3:2\pi/3','2\pi/3:\pi','\pi:4\pi/3','4\pi/3:5\pi/3','5\pi/3:2\pi'};
      options.legendLabels = {'\pi:2\pi/3','2\pi/3:\pi/3','\pi/3:0','0:-\pi/3','-\pi/3:-2\pi/3','-2\pi/3:-pi'};
      [statsMeanFRFiltBP0p01to0p05HzPiOver3{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05HzPiOver3{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_piOver3{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05HzPiOver4_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:\pi/4','\pi/4:2\pi/4','2\pi/4:3\pi/4','3\pi/4:\pi','\pi:5\pi/4','5\pi/4:6\pi/4','6\pi/4:7\pi/4','7\pi/4:2\pi'};
      options.legendLabels = {'\pi:3\pi/4','3\pi/4:\pi/2','\pi/2:\pi/4','\pi/4:0','0:-\pi/4','-\pi/4:-\pi/2','-\pi/2:-3\pi/4','-3\pi/4:-\pi'};
      [statsMeanFRFiltBP0p01to0p05HzPiOver4{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05HzPiOver4{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_piOver4{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05HzPiOver6_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:\pi/6','\pi/6:2\pi/6','2\pi/6:3\pi/6','3\pi/6:4\pi/6','4\pi/6:5\pi/6','5\pi/6:\pi','\pi:7\pi/6','7\pi/6:8\pi/6','8\pi/6:9\pi/6','9\pi/6:10\pi/6','10\pi/6:11\pi/6','11\pi/6:2\pi'};
      options.legendLabels = {'\pi:5\pi/6','5\pi/6:4\pi/6','4\pi/6:3\pi/6','3\pi/6:2\pi/6','2\pi/6:\pi/6','\pi/6:0','0:-\pi/6','-\pi/6:-2\pi/6','-2\pi/6:-3\pi/6','-3\pi/6:-4\pi/6','-4\pi/6:-5\pi/6','-5\pi/6:-\pi'};
      [statsMeanFRFiltBP0p01to0p05HzPiOver6{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05HzPiOver6{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_piOver6{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05HzPiOver8_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:\pi/8','\pi/8:2\pi/8','2\pi/8:3\pi/8','3\pi/8:4\pi/8','4\pi/8:5\pi/8','5\pi/8:6\pi/8','6\pi/8:7\pi/8','7\pi/8:\pi','\pi:9\pi/8','9\pi/8:10\pi/8','10\pi/8:11\pi/8','11\pi/8:12\pi/8','12\pi/8:13\pi/8','13\pi/8:14\pi/8','14\pi/8:15\pi/8','15\pi/8:2\pi'};
      options.legendLabels = {'\pi:7\pi/8','7\pi/8:6\pi/8','6\pi/8:5\pi/8','5\pi/8:4\pi/8','4\pi/8:3\pi/8','3\pi/8:2\pi/8','2\pi/8:\pi/8','\pi/8:0','0:-\pi/8','-\pi/8:-2\pi/8','-2\pi/8:-3\pi/8','-3\pi/8:-4\pi/8','-4\pi/8:-5\pi/8','-5\pi/8:-6\pi/8','-6\pi/8:-7\pi/8','-7\pi/8:-\pi'};
      [statsMeanFRFiltBP0p01to0p05HzPiOver8{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05HzPiOver8{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_piOver8{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05HzPiShifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:\pi]-\pi/2','[\pi:2\pi]-\pi/2'};
      options.legendLabels = {'[\pi:0]-\pi/2','[0:-\pi]-\pi/2'};
      [statsMeanFRFiltBP0p01to0p05HzPiShifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05HzPiShifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_piShifted{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05Hz2PiOver3Shifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:2\pi/3]-\pi/3','[2\pi/3:4\pi/3]-\pi/3','[4\pi/3:2\pi]-\pi/3'};
      options.legendLabels = {'[\pi:\pi/3]-\pi/3','[\pi/3:-\pi/3]-\pi/3','[-\pi/3:-\pi]-\pi/3'};
      [statsMeanFRFiltBP0p01to0p05Hz2piOver3Shifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05Hz2piOver3Shifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_2piOver3Shifted{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05HzPiOver2Shifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:\pi/2]-\pi/4','[\pi/2:\pi]-\pi/4','[\pi:3\pi/2]-\pi/4','[3\pi/2:2\pi]-\pi/4'};
      options.legendLabels = {'[\pi:\pi/2]-\pi/4','[\pi/2:0]-\pi/4','[0:-\pi/2]-\pi/4','[-\pi/2:-pi]-\pi/4'};
      [statsMeanFRFiltBP0p01to0p05HzPiOver2Shifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05HzPiOver2Shifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_piOver2Shifted{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05HzPiOver3Shifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:\pi/3]-\pi/6','[\pi/3:2\pi/3]-\pi/6','[2\pi/3:\pi]-\pi/6','[\pi:4\pi/3]-\pi/6','[4\pi/3:5\pi/3]-\pi/6','[5\pi/3:2\pi]-\pi/6'};
      options.legendLabels = {'[\pi:2\pi/3]-\pi/6','[2\pi/3:\pi/3]-\pi/6','[\pi/3:0]-\pi/6','[0:-\pi/3]-\pi/6','[-\pi/3:-2\pi/3]-\pi/6','[-2\pi/3:-pi]-\pi/6'};
      [statsMeanFRFiltBP0p01to0p05HzPiOver3Shifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05HzPiOver3Shifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_piOver3Shifted{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05HzPiOver4Shifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:\pi/4]-\pi/8','[\pi/4:2\pi/4]-\pi/8','[2\pi/4:3\pi/4]-\pi/8','[3\pi/4:\pi]-\pi/8','[\pi:5\pi/4]-\pi/8','[5\pi/4:6\pi/4]-\pi/8','[6\pi/4:7\pi/4]-\pi/8','[7\pi/4:2\pi]-\pi/8'};
      options.legendLabels = {'[\pi:3\pi/4]-\pi/8','[3\pi/4:\pi/2]-\pi/8','[\pi/2:\pi/4]-\pi/8','[\pi/4:0]-\pi/8','[0:-\pi/4]-\pi/8','[-\pi/4:-\pi/2]-\pi/8','[-\pi/2:-3\pi/4]-\pi/8','[-3\pi/4:-\pi]-\pi/8'};
      [statsMeanFRFiltBP0p01to0p05HzPiOver4Shifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05HzPiOver4Shifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_piOver4Shifted{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05HzPiOver6Shifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:\pi/6]-\pi/12','[\pi/6:2\pi/6]-\pi/12','[2\pi/6:3\pi/6]-\pi/12','[3\pi/6:4\pi/6]-\pi/12','[4\pi/6:5\pi/6]-\pi/12','[5\pi/6:\pi]-\pi/12','[\pi:7\pi/6]-\pi/12','[7\pi/6:8\pi/6]-\pi/12','[8\pi/6:9\pi/6]-\pi/12','[9\pi/6:10\pi/6]-\pi/12','[10\pi/6:11\pi/6]-\pi/12','[11\pi/6:2\pi]-\pi/12'};
      options.legendLabels = {'[\pi:5\pi/6]-\pi/12','[5\pi/6:4\pi/6]-\pi/12','[4\pi/6:3\pi/6]-\pi/12','[3\pi/6:2\pi/6]-\pi/12','[2\pi/6:\pi/6]-\pi/12','[\pi/6:0]-\pi/12','[0:-\pi/6]-\pi/12','[-\pi/6:-2\pi/6]-\pi/12','[-2\pi/6:-3\pi/6]-\pi/12','[-3\pi/6:-4\pi/6]-\pi/12','[-4\pi/6:-5\pi/6]-\pi/12','[-5\pi/6:-\pi]-\pi/12'};
      [statsMeanFRFiltBP0p01to0p05HzPiOver6Shifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05HzPiOver6Shifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_piOver6Shifted{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p01to0p05HzPiOver8Shifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:\pi/8]-\pi/16','[\pi/8:2\pi/8]-\pi/16','[2\pi/8:3\pi/8]-\pi/16','[3\pi/8:4\pi/8]-\pi/16','[4\pi/8:5\pi/8]-\pi/16','[5\pi/8:6\pi/8]-\pi/16','[6\pi/8:7\pi/8]-\pi/16','[7\pi/8:\pi]-\pi/16','[\pi:9\pi/8]-\pi/16','[9\pi/8:10\pi/8]-\pi/16','[10\pi/8:11\pi/8]-\pi/16','[11\pi/8:12\pi/8]-\pi/16','[12\pi/8:13\pi/8]-\pi/16','[13\pi/8:14\pi/8]-\pi/16','[14\pi/8:15\pi/8]-\pi/16','[15\pi/8:2\pi]-\pi/16'};
      options.legendLabels = {'[\pi:7\pi/8]-\pi/16','[7\pi/8:6\pi/8]-\pi/16','[6\pi/8:5\pi/8]-\pi/16','[5\pi/8:4\pi/8]-\pi/16','[4\pi/8:3\pi/8]-\pi/16','[3\pi/8:2\pi/8]-\pi/16','[2\pi/8:\pi/8]-\pi/16','[\pi/8:0]-\pi/16','[0:-\pi/8]-\pi/16','[-\pi/8:-2\pi/8]-\pi/16','[-2\pi/8:-3\pi/8]-\pi/16','[-3\pi/8:-4\pi/8]-\pi/16','[-4\pi/8:-5\pi/8]-\pi/16','[-5\pi/8:-6\pi/8]-\pi/16','[-6\pi/8:-7\pi/8]-\pi/16','[-7\pi/8:-\pi]-\pi/16'};
      [statsMeanFRFiltBP0p01to0p05HzPiOver8Shifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p01to0p05HzPiOver8Shifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p01to0p05HzIndividual_piOver8Shifted{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5HzPi_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:\pi','\pi:2\pi'};
      options.legendLabels = {'\pi:0','0:-\pi'};
      [statsMeanFRFiltBP0p1to0p5HzPi{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5HzPi{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_pi{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5Hz2PiOver3_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:2\pi/3','2\pi/3:4\pi/3','4\pi/3:2\pi'};
      options.legendLabels = {'\pi:\pi/3','\pi/3:-\pi/3','-\pi/3:-\pi'};
      [statsMeanFRFiltBP0p1to0p5Hz2piOver3{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5Hz2piOver3{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_2piOver3{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5HzPiOver2_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:\pi/2','\pi/2:\pi','\pi:3\pi/2','3\pi/2:2\pi'};
      options.legendLabels = {'\pi:\pi/2','\pi/2:0','0:-\pi/2','-\pi/2:-pi'};
      [statsMeanFRFiltBP0p1to0p5HzPiOver2{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5HzPiOver2{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_piOver2{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5HzPiOver3_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:\pi/3','\pi/3:2\pi/3','2\pi/3:\pi','\pi:4\pi/3','4\pi/3:5\pi/3','5\pi/3:2\pi'};
      options.legendLabels = {'\pi:2\pi/3','2\pi/3:\pi/3','\pi/3:0','0:-\pi/3','-\pi/3:-2\pi/3','-2\pi/3:-pi'};
      [statsMeanFRFiltBP0p1to0p5HzPiOver3{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5HzPiOver3{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_piOver3{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5HzPiOver4_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:\pi/4','\pi/4:2\pi/4','2\pi/4:3\pi/4','3\pi/4:\pi','\pi:5\pi/4','5\pi/4:6\pi/4','6\pi/4:7\pi/4','7\pi/4:2\pi'};
      options.legendLabels = {'\pi:3\pi/4','3\pi/4:\pi/2','\pi/2:\pi/4','\pi/4:0','0:-\pi/4','-\pi/4:-\pi/2','-\pi/2:-3\pi/4','-3\pi/4:-\pi'};
      [statsMeanFRFiltBP0p1to0p5HzPiOver4{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5HzPiOver4{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_piOver4{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5HzPiOver6_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:\pi/6','\pi/6:2\pi/6','2\pi/6:3\pi/6','3\pi/6:4\pi/6','4\pi/6:5\pi/6','5\pi/6:\pi','\pi:7\pi/6','7\pi/6:8\pi/6','8\pi/6:9\pi/6','9\pi/6:10\pi/6','10\pi/6:11\pi/6','11\pi/6:2\pi'};
      options.legendLabels = {'\pi:5\pi/6','5\pi/6:4\pi/6','4\pi/6:3\pi/6','3\pi/6:2\pi/6','2\pi/6:\pi/6','\pi/6:0','0:-\pi/6','-\pi/6:-2\pi/6','-2\pi/6:-3\pi/6','-3\pi/6:-4\pi/6','-4\pi/6:-5\pi/6','-5\pi/6:-\pi'};
      [statsMeanFRFiltBP0p1to0p5HzPiOver6{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5HzPiOver6{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_piOver6{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5HzPiOver8_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'0:\pi/8','\pi/8:2\pi/8','2\pi/8:3\pi/8','3\pi/8:4\pi/8','4\pi/8:5\pi/8','5\pi/8:6\pi/8','6\pi/8:7\pi/8','7\pi/8:\pi','\pi:9\pi/8','9\pi/8:10\pi/8','10\pi/8:11\pi/8','11\pi/8:12\pi/8','12\pi/8:13\pi/8','13\pi/8:14\pi/8','14\pi/8:15\pi/8','15\pi/8:2\pi'};
      options.legendLabels = {'\pi:7\pi/8','7\pi/8:6\pi/8','6\pi/8:5\pi/8','5\pi/8:4\pi/8','4\pi/8:3\pi/8','3\pi/8:2\pi/8','2\pi/8:\pi/8','\pi/8:0','0:-\pi/8','-\pi/8:-2\pi/8','-2\pi/8:-3\pi/8','-3\pi/8:-4\pi/8','-4\pi/8:-5\pi/8','-5\pi/8:-6\pi/8','-6\pi/8:-7\pi/8','-7\pi/8:-\pi'};
      [statsMeanFRFiltBP0p1to0p5HzPiOver8{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5HzPiOver8{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_piOver8{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5HzPiShifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:\pi]-\pi/2','[\pi:2\pi]-\pi/2'};
      options.legendLabels = {'[\pi:0]-\pi/2','[0:-\pi]-\pi/2'};
      [statsMeanFRFiltBP0p1to0p5HzPiShifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5HzPiShifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_piShifted{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5Hz2PiOver3Shifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:2\pi/3]-\pi/3','[2\pi/3:4\pi/3]-\pi/3','[4\pi/3:2\pi]-\pi/3'};
      options.legendLabels = {'[\pi:\pi/3]-\pi/3','[\pi/3:-\pi/3]-\pi/3','[-\pi/3:-\pi]-\pi/3'};
      [statsMeanFRFiltBP0p1to0p5Hz2piOver3Shifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5Hz2piOver3Shifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_2piOver3Shifted{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5HzPiOver2Shifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:\pi/2]-\pi/4','[\pi/2:\pi]-\pi/4','[\pi:3\pi/2]-\pi/4','[3\pi/2:2\pi]-\pi/4'};
      options.legendLabels = {'[\pi:\pi/2]-\pi/4','[\pi/2:0]-\pi/4','[0:-\pi/2]-\pi/4','[-\pi/2:-pi]-\pi/4'};
      [statsMeanFRFiltBP0p1to0p5HzPiOver2Shifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5HzPiOver2Shifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_piOver2Shifted{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5HzPiOver3Shifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:\pi/3]-\pi/6','[\pi/3:2\pi/3]-\pi/6','[2\pi/3:\pi]-\pi/6','[\pi:4\pi/3]-\pi/6','[4\pi/3:5\pi/3]-\pi/6','[5\pi/3:2\pi]-\pi/6'};
      options.legendLabels = {'[\pi:2\pi/3]-\pi/6','[2\pi/3:\pi/3]-\pi/6','[\pi/3:0]-\pi/6','[0:-\pi/3]-\pi/6','[-\pi/3:-2\pi/3]-\pi/6','[-2\pi/3:-pi]-\pi/6'};
      [statsMeanFRFiltBP0p1to0p5HzPiOver3Shifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5HzPiOver3Shifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_piOver3Shifted{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5HzPiOver4Shifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:\pi/4]-\pi/8','[\pi/4:2\pi/4]-\pi/8','[2\pi/4:3\pi/4]-\pi/8','[3\pi/4:\pi]-\pi/8','[\pi:5\pi/4]-\pi/8','[5\pi/4:6\pi/4]-\pi/8','[6\pi/4:7\pi/4]-\pi/8','[7\pi/4:2\pi]-\pi/8'};
      options.legendLabels = {'[\pi:3\pi/4]-\pi/8','[3\pi/4:\pi/2]-\pi/8','[\pi/2:\pi/4]-\pi/8','[\pi/4:0]-\pi/8','[0:-\pi/4]-\pi/8','[-\pi/4:-\pi/2]-\pi/8','[-\pi/2:-3\pi/4]-\pi/8','[-3\pi/4:-\pi]-\pi/8'};
      [statsMeanFRFiltBP0p1to0p5HzPiOver4Shifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5HzPiOver4Shifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_piOver4Shifted{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5HzPiOver6Shifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:\pi/6]-\pi/12','[\pi/6:2\pi/6]-\pi/12','[2\pi/6:3\pi/6]-\pi/12','[3\pi/6:4\pi/6]-\pi/12','[4\pi/6:5\pi/6]-\pi/12','[5\pi/6:\pi]-\pi/12','[\pi:7\pi/6]-\pi/12','[7\pi/6:8\pi/6]-\pi/12','[8\pi/6:9\pi/6]-\pi/12','[9\pi/6:10\pi/6]-\pi/12','[10\pi/6:11\pi/6]-\pi/12','[11\pi/6:2\pi]-\pi/12'};
      options.legendLabels = {'[\pi:5\pi/6]-\pi/12','[5\pi/6:4\pi/6]-\pi/12','[4\pi/6:3\pi/6]-\pi/12','[3\pi/6:2\pi/6]-\pi/12','[2\pi/6:\pi/6]-\pi/12','[\pi/6:0]-\pi/12','[0:-\pi/6]-\pi/12','[-\pi/6:-2\pi/6]-\pi/12','[-2\pi/6:-3\pi/6]-\pi/12','[-3\pi/6:-4\pi/6]-\pi/12','[-4\pi/6:-5\pi/6]-\pi/12','[-5\pi/6:-\pi]-\pi/12'};
      [statsMeanFRFiltBP0p1to0p5HzPiOver6Shifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5HzPiOver6Shifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_piOver6Shifted{iCond}{iAreasOI(iArea)}, [], options);
      
      options.figName = ['FiringRateFiltBP0p1to0p5HzPiOver8Shifted_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
      %options.legendLabels = {'[0:\pi/8]-\pi/16','[\pi/8:2\pi/8]-\pi/16','[2\pi/8:3\pi/8]-\pi/16','[3\pi/8:4\pi/8]-\pi/16','[4\pi/8:5\pi/8]-\pi/16','[5\pi/8:6\pi/8]-\pi/16','[6\pi/8:7\pi/8]-\pi/16','[7\pi/8:\pi]-\pi/16','[\pi:9\pi/8]-\pi/16','[9\pi/8:10\pi/8]-\pi/16','[10\pi/8:11\pi/8]-\pi/16','[11\pi/8:12\pi/8]-\pi/16','[12\pi/8:13\pi/8]-\pi/16','[13\pi/8:14\pi/8]-\pi/16','[14\pi/8:15\pi/8]-\pi/16','[15\pi/8:2\pi]-\pi/16'};
      options.legendLabels = {'[\pi:7\pi/8]-\pi/16','[7\pi/8:6\pi/8]-\pi/16','[6\pi/8:5\pi/8]-\pi/16','[5\pi/8:4\pi/8]-\pi/16','[4\pi/8:3\pi/8]-\pi/16','[3\pi/8:2\pi/8]-\pi/16','[2\pi/8:\pi/8]-\pi/16','[\pi/8:0]-\pi/16','[0:-\pi/8]-\pi/16','[-\pi/8:-2\pi/8]-\pi/16','[-2\pi/8:-3\pi/8]-\pi/16','[-3\pi/8:-4\pi/8]-\pi/16','[-4\pi/8:-5\pi/8]-\pi/16','[-5\pi/8:-6\pi/8]-\pi/16','[-6\pi/8:-7\pi/8]-\pi/16','[-7\pi/8:-\pi]-\pi/16'};
      [statsMeanFRFiltBP0p1to0p5HzPiOver8Shifted{iCond}{iAreasOI(iArea)}, statsVarFRFiltBP0p1to0p5HzPiOver8Shifted{iCond}{iAreasOI(iArea)}] = firingRateTestAndPlot(edges, areaFRFiltBP0p1to0p5HzIndividual_piOver8Shifted{iCond}{iAreasOI(iArea)}, [], options);
    end
  end
  
  % Save variance data
  save(filename, 'statsMeanFRPercentile50','statsMeanFRPercentile25','statsMeanFRPercentile12p5','statsMeanFRPercentileThird','statsMeanFRRiseDecay',...
    'statsMeanFRFiltLP0p001Hz50','statsMeanFRFiltLP0p001Hz25','statsMeanFRFiltLP0p001Hz12p5','statsMeanFRFiltLP0p001HzThird',...
    'statsMeanFRFiltLP0p01Hz50','statsMeanFRFiltLP0p1Hz50',...
    'statsMeanFRFiltBP0p01to0p05HzPi','statsMeanFRFiltBP0p01to0p05Hz2piOver3','statsMeanFRFiltBP0p01to0p05HzPiOver2','statsMeanFRFiltBP0p01to0p05HzPiOver3',...
    'statsMeanFRFiltBP0p01to0p05HzPiOver4','statsMeanFRFiltBP0p01to0p05HzPiOver6','statsMeanFRFiltBP0p01to0p05HzPiOver8',...
    'statsMeanFRFiltBP0p01to0p05HzPiShifted','statsMeanFRFiltBP0p01to0p05Hz2piOver3Shifted','statsMeanFRFiltBP0p01to0p05HzPiOver2Shifted',...
    'statsMeanFRFiltBP0p01to0p05HzPiOver3Shifted','statsMeanFRFiltBP0p01to0p05HzPiOver4Shifted','statsMeanFRFiltBP0p01to0p05HzPiOver6Shifted',...
    'statsMeanFRFiltBP0p01to0p05HzPiOver8Shifted','statsMeanFRFiltBP0p1to0p5HzPi','statsMeanFRFiltBP0p1to0p5Hz2piOver3','statsMeanFRFiltBP0p1to0p5HzPiOver2',...
    'statsMeanFRFiltBP0p1to0p5HzPiOver3','statsMeanFRFiltBP0p1to0p5HzPiOver4','statsMeanFRFiltBP0p1to0p5HzPiOver6','statsMeanFRFiltBP0p1to0p5HzPiOver8',...
    'statsMeanFRFiltBP0p1to0p5HzPiShifted','statsMeanFRFiltBP0p1to0p5Hz2piOver3Shifted','statsMeanFRFiltBP0p1to0p5HzPiOver2Shifted',...
    'statsMeanFRFiltBP0p1to0p5HzPiOver3Shifted','statsMeanFRFiltBP0p1to0p5HzPiOver4Shifted','statsMeanFRFiltBP0p1to0p5HzPiOver6Shifted',...
    'statsMeanFRFiltBP0p1to0p5HzPiOver8Shifted',...
    'statsVarFRPercentile50','statsVarFRPercentile25','statsVarFRPercentile12p5','statsVarFRPercentileThird','statsVarFRRiseDecay',...
    'statsVarFRFiltLP0p001Hz50','statsVarFRFiltLP0p001Hz25','statsVarFRFiltLP0p001Hz12p5','statsVarFRFiltLP0p001HzThird',...
    'statsVarFRFiltLP0p01Hz50','statsVarFRFiltLP0p1Hz50',...
    'statsVarFRFiltBP0p01to0p05HzPi','statsVarFRFiltBP0p01to0p05Hz2piOver3','statsVarFRFiltBP0p01to0p05HzPiOver2','statsVarFRFiltBP0p01to0p05HzPiOver3',...
    'statsVarFRFiltBP0p01to0p05HzPiOver4','statsVarFRFiltBP0p01to0p05HzPiOver6','statsVarFRFiltBP0p01to0p05HzPiOver8',...
    'statsVarFRFiltBP0p01to0p05HzPiShifted','statsVarFRFiltBP0p01to0p05Hz2piOver3Shifted','statsVarFRFiltBP0p01to0p05HzPiOver2Shifted',...
    'statsVarFRFiltBP0p01to0p05HzPiOver3Shifted','statsVarFRFiltBP0p01to0p05HzPiOver4Shifted','statsVarFRFiltBP0p01to0p05HzPiOver6Shifted',...
    'statsVarFRFiltBP0p01to0p05HzPiOver8Shifted','statsVarFRFiltBP0p1to0p5HzPi','statsVarFRFiltBP0p1to0p5Hz2piOver3','statsVarFRFiltBP0p1to0p5HzPiOver2',...
    'statsVarFRFiltBP0p1to0p5HzPiOver3','statsVarFRFiltBP0p1to0p5HzPiOver4','statsVarFRFiltBP0p1to0p5HzPiOver6','statsVarFRFiltBP0p1to0p5HzPiOver8',...
    'statsVarFRFiltBP0p1to0p5HzPiShifted','statsVarFRFiltBP0p1to0p5Hz2piOver3Shifted','statsVarFRFiltBP0p1to0p5HzPiOver2Shifted',...
    'statsVarFRFiltBP0p1to0p5HzPiOver3Shifted','statsVarFRFiltBP0p1to0p5HzPiOver4Shifted','statsVarFRFiltBP0p1to0p5HzPiOver6Shifted',...
    'statsVarFRFiltBP0p1to0p5HzPiOver8Shifted',...
    'statsMeanFRPercentile50Split','statsMeanFRPercentile25Split','statsMeanFRPercentile12p5Split','statsMeanFRPercentileThirdSplit',...
    'statsVarFRPercentile50Split','statsVarFRPercentile25Split','statsVarFRPercentile12p5Split','statsVarFRPercentileThirdSplit', '-append');
end


%% PLOT THE LOG FIRING RATE DISTRIBUTION MEANS
if drawMeans(1) || drawMeans(2)
  options = struct();
  options.xLabel = 'Phase range (rad)';
  if strcmpi(subpop, 'all')
    options.yLabel = 'Log_{10}(firing rate)';
  elseif strcmpi(subpop, 'positive')
    %options.yLabel = 'Positive log_{10}(firing rate)';
    options.yLabel = 'Log_{10}(firing rate)';
  elseif strcmpi(subpop, 'negative')
    %options.yLabel = 'Negative log_{10}(firing rate)';
    options.yLabel = 'Log_{10}(firing rate)';
  end
  options.figFolder = mainFolder;
  options.saveFig = true;
  for iCond = 1:min([2 numel(conditions)])
    if strcmpi(repository, 'uol')
      if iCond == 1
        options.yLim = [-0.5 2]; %[0 1.3];
      elseif iCond == 2
        options.yLim = [-1 1.5]; %[-0.5 0.5];
      else
        options.yLim = [-1 2]; %[-0.5 1.3];
      end
    elseif strcmpi(repository, 'allensdk')
      options.yLim = [-0.5 2]; %[0.3 1.3];
    end
    for iArea = 1:numel(iAreasOI)
      if size(areaFRFiltBP0p01to0p05HzIndividual_pi{iCond}{iAreasOI(iArea)},1) > 1
        if drawMeans(1)
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi
          %xTickLabels = {'0:\pi','\pi:2\pi'};
          xTickLabels = {'\pi:0','0:-\pi'};
          firingRateViolinPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_pi{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05HzPi{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiMeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond} '.fig'], options);
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: 2*pi/3
          %xTickLabels = {'0:2\pi/3','2\pi/3:4\pi/3','4\pi/3:2\pi'};
          xTickLabels = {'\pi:\pi/3','\pi/3:-\pi/3','-\pi/3:-\pi'};
          firingRateViolinPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_2piOver3{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05Hz2piOver3{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05Hz2PiOver3MeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/2
          %xTickLabels = {'0:\pi/2','\pi/2:\pi','\pi:3\pi/2','3\pi/2:2\pi'};
          xTickLabels = {'\pi:\pi/2','\pi/2:0','0:-\pi/2','-\pi/2:-pi'};
          firingRateViolinPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver2{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05HzPiOver2{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver2MeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/3
          %xTickLabels = {'0:\pi/3','\pi/3:2\pi/3','2\pi/3:\pi','\pi:4\pi/3','4\pi/3:5\pi/3','5\pi/3:2\pi'};
          xTickLabels = {'\pi:2\pi/3','2\pi/3:\pi/3','\pi/3:0','0:-\pi/3','-\pi/3:-2\pi/3','-2\pi/3:-pi'};
          firingRateViolinPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver3{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05HzPiOver3{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver3MeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/4
          %xTickLabels = {'0:\pi/4','\pi/4:2\pi/4','2\pi/4:3\pi/4','3\pi/4:\pi','\pi:5\pi/4','5\pi/4:6\pi/4','6\pi/4:7\pi/4','7\pi/4:2\pi'};
          xTickLabels = {'\pi:3\pi/4','3\pi/4:\pi/2','\pi/2:\pi/4','\pi/4:0','0:-\pi/4','-\pi/4:-\pi/2','-\pi/2:-3\pi/4','-3\pi/4:-\pi'};
          firingRateViolinPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver4{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05HzPiOver4{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver4MeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/6
          %xTickLabels = {'0:\pi/6','\pi/6:2\pi/6','2\pi/6:3\pi/6','3\pi/6:4\pi/6','4\pi/6:5\pi/6','5\pi/6:\pi','\pi:7\pi/6','7\pi/6:8\pi/6','8\pi/6:9\pi/6','9\pi/6:10\pi/6','10\pi/6:11\pi/6','11\pi/6:2\pi'};
          xTickLabels = {'\pi:5\pi/6','5\pi/6:4\pi/6','4\pi/6:3\pi/6','3\pi/6:2\pi/6','2\pi/6:\pi/6','\pi/6:0','0:-\pi/6','-\pi/6:-2\pi/6','-2\pi/6:-3\pi/6','-3\pi/6:-4\pi/6','-4\pi/6:-5\pi/6','-5\pi/6:-\pi'};
          firingRateViolinPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver6{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05HzPiOver6{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver6MeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/8
          %xTickLabels = {'0:\pi/8','\pi/8:2\pi/8','2\pi/8:3\pi/8','3\pi/8:4\pi/8','4\pi/8:5\pi/8','5\pi/8:6\pi/8','6\pi/8:7\pi/8','7\pi/8:\pi','\pi:9\pi/8','9\pi/8:10\pi/8','10\pi/8:11\pi/8','11\pi/8:12\pi/8','12\pi/8:13\pi/8','13\pi/8:14\pi/8','14\pi/8:15\pi/8','15\pi/8:2\pi'};
          xTickLabels = {'\pi:7\pi/8','7\pi/8:6\pi/8','6\pi/8:5\pi/8','5\pi/8:4\pi/8','4\pi/8:3\pi/8','3\pi/8:2\pi/8','2\pi/8:\pi/8','\pi/8:0','0:-\pi/8','-\pi/8:-2\pi/8','-2\pi/8:-3\pi/8','-3\pi/8:-4\pi/8','-4\pi/8:-5\pi/8','-5\pi/8:-6\pi/8','-6\pi/8:-7\pi/8','-7\pi/8:-\pi'};
          firingRateViolinPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver8{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05HzPiOver8{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver8MeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi
          %xTickLabels = {'0:\pi','\pi:2\pi'};
          xTickLabels = {'\pi:0','0:-\pi'};
          firingRateViolinPlot({'0:\pi','\pi:2\pi'}, getLog(areaFRFiltBP0p1to0p5HzIndividual_pi{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5HzPi{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiMeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: 2*pi/3
          %xTickLabels = {'0:2\pi/3','2\pi/3:4\pi/3','4\pi/3:2\pi'};
          xTickLabels = {'\pi:\pi/3','\pi/3:-\pi/3','-\pi/3:-\pi'};
          firingRateViolinPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_2piOver3{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5Hz2piOver3{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5Hz2PiOver3MeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/2
          %xTickLabels = {'0:\pi/2','\pi/2:\pi','\pi:3\pi/2','3\pi/2:2\pi'};
          xTickLabels = {'\pi:\pi/2','\pi/2:0','0:-\pi/2','-\pi/2:-pi'};
          firingRateViolinPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver2{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5HzPiOver2{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver2MeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/3
          %xTickLabels = {'0:\pi/3','\pi/3:2\pi/3','2\pi/3:\pi','\pi:4\pi/3','4\pi/3:5\pi/3','5\pi/3:2\pi'};
          xTickLabels = {'\pi:2\pi/3','2\pi/3:\pi/3','\pi/3:0','0:-\pi/3','-\pi/3:-2\pi/3','-2\pi/3:-pi'};
          firingRateViolinPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver3{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5HzPiOver3{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver3MeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/4
          %xTickLabels = {'0:\pi/4','\pi/4:2\pi/4','2\pi/4:3\pi/4','3\pi/4:\pi','\pi:5\pi/4','5\pi/4:6\pi/4','6\pi/4:7\pi/4','7\pi/4:2\pi'};
          xTickLabels = {'\pi:3\pi/4','3\pi/4:\pi/2','\pi/2:\pi/4','\pi/4:0','0:-\pi/4','-\pi/4:-\pi/2','-\pi/2:-3\pi/4','-3\pi/4:-\pi'};
          firingRateViolinPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver4{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5HzPiOver4{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver4MeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/6
          %xTickLabels = {'0:\pi/6','\pi/6:2\pi/6','2\pi/6:3\pi/6','3\pi/6:4\pi/6','4\pi/6:5\pi/6','5\pi/6:\pi','\pi:7\pi/6','7\pi/6:8\pi/6','8\pi/6:9\pi/6','9\pi/6:10\pi/6','10\pi/6:11\pi/6','11\pi/6:2\pi'};
          xTickLabels = {'\pi:5\pi/6','5\pi/6:4\pi/6','4\pi/6:3\pi/6','3\pi/6:2\pi/6','2\pi/6:\pi/6','\pi/6:0','0:-\pi/6','-\pi/6:-2\pi/6','-2\pi/6:-3\pi/6','-3\pi/6:-4\pi/6','-4\pi/6:-5\pi/6','-5\pi/6:-\pi'};
          firingRateViolinPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver6{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5HzPiOver6{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver6MeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/8
          %xTickLabels = {'0:\pi/8','\pi/8:2\pi/8','2\pi/8:3\pi/8','3\pi/8:4\pi/8','4\pi/8:5\pi/8','5\pi/8:6\pi/8','6\pi/8:7\pi/8','7\pi/8:\pi','\pi:9\pi/8','9\pi/8:10\pi/8','10\pi/8:11\pi/8','11\pi/8:12\pi/8','12\pi/8:13\pi/8','13\pi/8:14\pi/8','14\pi/8:15\pi/8','15\pi/8:2\pi'};
          xTickLabels = {'\pi:7\pi/8','7\pi/8:6\pi/8','6\pi/8:5\pi/8','5\pi/8:4\pi/8','4\pi/8:3\pi/8','3\pi/8:2\pi/8','2\pi/8:\pi/8','\pi/8:0','0:-\pi/8','-\pi/8:-2\pi/8','-2\pi/8:-3\pi/8','-3\pi/8:-4\pi/8','-4\pi/8:-5\pi/8','-5\pi/8:-6\pi/8','-6\pi/8:-7\pi/8','-7\pi/8:-\pi'};
          firingRateViolinPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver8{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5HzPiOver8{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver8MeanViolins_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        end
        if drawMeans(2)
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi
          %xTickLabels = {'0:\pi','\pi:2\pi'};
          xTickLabels = {'\pi:0','0:-\pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_pi{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05HzPi{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiMeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond} '.fig'], options);
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: 2*pi/3
          %xTickLabels = {'0:2\pi/3','2\pi/3:4\pi/3','4\pi/3:2\pi'};
          xTickLabels = {'\pi:\pi/3','\pi/3:-\pi/3','-\pi/3:-\pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_2piOver3{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05Hz2piOver3{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05Hz2PiOver3MeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/2
          %xTickLabels = {'0:\pi/2','\pi/2:\pi','\pi:3\pi/2','3\pi/2:2\pi'};
          xTickLabels = {'\pi:\pi/2','\pi/2:0','0:-\pi/2','-\pi/2:-pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver2{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05HzPiOver2{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver2MeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/3
          %xTickLabels = {'0:\pi/3','\pi/3:2\pi/3','2\pi/3:\pi','\pi:4\pi/3','4\pi/3:5\pi/3','5\pi/3:2\pi'};
          xTickLabels = {'\pi:2\pi/3','2\pi/3:\pi/3','\pi/3:0','0:-\pi/3','-\pi/3:-2\pi/3','-2\pi/3:-pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver3{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05HzPiOver3{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver3MeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/4
          %xTickLabels = {'0:\pi/4','\pi/4:2\pi/4','2\pi/4:3\pi/4','3\pi/4:\pi','\pi:5\pi/4','5\pi/4:6\pi/4','6\pi/4:7\pi/4','7\pi/4:2\pi'};
          xTickLabels = {'\pi:3\pi/4','3\pi/4:\pi/2','\pi/2:\pi/4','\pi/4:0','0:-\pi/4','-\pi/4:-\pi/2','-\pi/2:-3\pi/4','-3\pi/4:-\pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver4{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05HzPiOver4{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver4MeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/6
          %xTickLabels = {'0:\pi/6','\pi/6:2\pi/6','2\pi/6:3\pi/6','3\pi/6:4\pi/6','4\pi/6:5\pi/6','5\pi/6:\pi','\pi:7\pi/6','7\pi/6:8\pi/6','8\pi/6:9\pi/6','9\pi/6:10\pi/6','10\pi/6:11\pi/6','11\pi/6:2\pi'};
          xTickLabels = {'\pi:5\pi/6','5\pi/6:4\pi/6','4\pi/6:3\pi/6','3\pi/6:2\pi/6','2\pi/6:\pi/6','\pi/6:0','0:-\pi/6','-\pi/6:-2\pi/6','-2\pi/6:-3\pi/6','-3\pi/6:-4\pi/6','-4\pi/6:-5\pi/6','-5\pi/6:-\pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver6{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05HzPiOver6{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver6MeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/8
          %xTickLabels = {'0:\pi/8','\pi/8:2\pi/8','2\pi/8:3\pi/8','3\pi/8:4\pi/8','4\pi/8:5\pi/8','5\pi/8:6\pi/8','6\pi/8:7\pi/8','7\pi/8:\pi','\pi:9\pi/8','9\pi/8:10\pi/8','10\pi/8:11\pi/8','11\pi/8:12\pi/8','12\pi/8:13\pi/8','13\pi/8:14\pi/8','14\pi/8:15\pi/8','15\pi/8:2\pi'};
          xTickLabels = {'\pi:7\pi/8','7\pi/8:6\pi/8','6\pi/8:5\pi/8','5\pi/8:4\pi/8','4\pi/8:3\pi/8','3\pi/8:2\pi/8','2\pi/8:\pi/8','\pi/8:0','0:-\pi/8','-\pi/8:-2\pi/8','-2\pi/8:-3\pi/8','-3\pi/8:-4\pi/8','-4\pi/8:-5\pi/8','-5\pi/8:-6\pi/8','-6\pi/8:-7\pi/8','-7\pi/8:-\pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver8{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p01to0p05HzPiOver8{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver8MeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi
          %xTickLabels = {'0:\pi','\pi:2\pi'};
          xTickLabels = {'\pi:0','0:-\pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_pi{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5HzPi{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiMeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: 2*pi/3
          %xTickLabels = {'0:2\pi/3','2\pi/3:4\pi/3','4\pi/3:2\pi'};
          xTickLabels = {'\pi:\pi/3','\pi/3:-\pi/3','-\pi/3:-\pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_2piOver3{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5Hz2piOver3{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5Hz2PiOver3MeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/2
          %xTickLabels = {'0:\pi/2','\pi/2:\pi','\pi:3\pi/2','3\pi/2:2\pi'};
          xTickLabels = {'\pi:\pi/2','\pi/2:0','0:-\pi/2','-\pi/2:-pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver2{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5HzPiOver2{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver2MeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/3
          %xTickLabels = {'0:\pi/3','\pi/3:2\pi/3','2\pi/3:\pi','\pi:4\pi/3','4\pi/3:5\pi/3','5\pi/3:2\pi'};
          xTickLabels = {'\pi:2\pi/3','2\pi/3:\pi/3','\pi/3:0','0:-\pi/3','-\pi/3:-2\pi/3','-2\pi/3:-pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver3{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5HzPiOver3{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver3MeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/4
          %xTickLabels = {'0:\pi/4','\pi/4:2\pi/4','2\pi/4:3\pi/4','3\pi/4:\pi','\pi:5\pi/4','5\pi/4:6\pi/4','6\pi/4:7\pi/4','7\pi/4:2\pi'};
          xTickLabels = {'\pi:3\pi/4','3\pi/4:\pi/2','\pi/2:\pi/4','\pi/4:0','0:-\pi/4','-\pi/4:-\pi/2','-\pi/2:-3\pi/4','-3\pi/4:-\pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver4{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5HzPiOver4{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver4MeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/6
          %xTickLabels = {'0:\pi/6','\pi/6:2\pi/6','2\pi/6:3\pi/6','3\pi/6:4\pi/6','4\pi/6:5\pi/6','5\pi/6:\pi','\pi:7\pi/6','7\pi/6:8\pi/6','8\pi/6:9\pi/6','9\pi/6:10\pi/6','10\pi/6:11\pi/6','11\pi/6:2\pi'};
          xTickLabels = {'\pi:5\pi/6','5\pi/6:4\pi/6','4\pi/6:3\pi/6','3\pi/6:2\pi/6','2\pi/6:\pi/6','\pi/6:0','0:-\pi/6','-\pi/6:-2\pi/6','-2\pi/6:-3\pi/6','-3\pi/6:-4\pi/6','-4\pi/6:-5\pi/6','-5\pi/6:-\pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver6{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5HzPiOver6{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver6MeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
          
          % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/8
          %xTickLabels = {'0:\pi/8','\pi/8:2\pi/8','2\pi/8:3\pi/8','3\pi/8:4\pi/8','4\pi/8:5\pi/8','5\pi/8:6\pi/8','6\pi/8:7\pi/8','7\pi/8:\pi','\pi:9\pi/8','9\pi/8:10\pi/8','10\pi/8:11\pi/8','11\pi/8:12\pi/8','12\pi/8:13\pi/8','13\pi/8:14\pi/8','14\pi/8:15\pi/8','15\pi/8:2\pi'};
          xTickLabels = {'\pi:7\pi/8','7\pi/8:6\pi/8','6\pi/8:5\pi/8','5\pi/8:4\pi/8','4\pi/8:3\pi/8','3\pi/8:2\pi/8','2\pi/8:\pi/8','\pi/8:0','0:-\pi/8','-\pi/8:-2\pi/8','-2\pi/8:-3\pi/8','-3\pi/8:-4\pi/8','-4\pi/8:-5\pi/8','-5\pi/8:-6\pi/8','-6\pi/8:-7\pi/8','-7\pi/8:-\pi'};
          firingRateDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver8{iCond}{iAreasOI(iArea)}),...
            statsMeanFRFiltBP0p1to0p5HzPiOver8{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver8MeanDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        end
      end
    end
  end
end


%% PLOT THE LOG FIRING RATE DISTRIBUTION VARIANCES
if drawVars
  options = struct();
  options.yLim = [];
  options.xLabel = 'Phase range (rad)';
  if strcmpi(subpop, 'all')
    options.yLabel = 'Variance of log_{10}(firing rate)';
  elseif strcmpi(subpop, 'positive')
    %options.yLabel = 'Positive variance of log_{10}(firing rate)';
    options.yLabel = 'Variance of log_{10}(firing rate)';
  elseif strcmpi(subpop, 'negative')
    %options.yLabel = 'Negative variance of log_{10}(firing rate)';
    options.yLabel = 'Variance of log_{10}(firing rate)';
  end
  options.figFolder = mainFolder;
  options.transparent = true;
  options.coloured = true;
  options.figSize = 18;
  options.saveFig = true;
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if size(areaFRFiltBP0p01to0p05HzIndividual_pi{iCond}{iAreasOI(iArea)},1) > 1
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi
        %xTickLabels = {'0:\pi','\pi:2\pi'};
        xTickLabels = {'\pi:0','0:-\pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_pi{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p01to0p05HzPi{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiVarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: 2*pi/3
        %xTickLabels = {'0:2\pi/3','2\pi/3:4\pi/3','4\pi/3:2\pi'};
        xTickLabels = {'\pi:\pi/3','\pi/3:-\pi/3','-\pi/3:-\pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_2piOver3{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p01to0p05Hz2piOver3{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05Hz2PiOver3VarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/2
        %xTickLabels = {'0:\pi/2','\pi/2:\pi','\pi:3\pi/2','3\pi/2:2\pi'};
        xTickLabels = {'\pi:\pi/2','\pi/2:0','0:-\pi/2','-\pi/2:-pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver2{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p01to0p05HzPiOver2{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver2VarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/3
        %xTickLabels = {'0:\pi/3','\pi/3:2\pi/3','2\pi/3:\pi','\pi:4\pi/3','4\pi/3:5\pi/3','5\pi/3:2\pi'};
        xTickLabels = {'\pi:2\pi/3','2\pi/3:\pi/3','\pi/3:0','0:-\pi/3','-\pi/3:-2\pi/3','-2\pi/3:-pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver3{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p01to0p05HzPiOver3{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver3VarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/4
        %xTickLabels = {'0:\pi/4','\pi/4:2\pi/4','2\pi/4:3\pi/4','3\pi/4:\pi','\pi:5\pi/4','5\pi/4:6\pi/4','6\pi/4:7\pi/4','7\pi/4:2\pi'};
        xTickLabels = {'\pi:3\pi/4','3\pi/4:\pi/2','\pi/2:\pi/4','\pi/4:0','0:-\pi/4','-\pi/4:-\pi/2','-\pi/2:-3\pi/4','-3\pi/4:-\pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver4{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p01to0p05HzPiOver4{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver4VarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/6
        %xTickLabels = {'0:\pi/6','\pi/6:2\pi/6','2\pi/6:3\pi/6','3\pi/6:4\pi/6','4\pi/6:5\pi/6','5\pi/6:\pi','\pi:7\pi/6','7\pi/6:8\pi/6','8\pi/6:9\pi/6','9\pi/6:10\pi/6','10\pi/6:11\pi/6','11\pi/6:2\pi'};
        xTickLabels = {'\pi:5\pi/6','5\pi/6:4\pi/6','4\pi/6:3\pi/6','3\pi/6:2\pi/6','2\pi/6:\pi/6','\pi/6:0','0:-\pi/6','-\pi/6:-2\pi/6','-2\pi/6:-3\pi/6','-3\pi/6:-4\pi/6','-4\pi/6:-5\pi/6','-5\pi/6:-\pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver6{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p01to0p05HzPiOver6{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver6VarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/8
        %xTickLabels = {'0:\pi/8','\pi/8:2\pi/8','2\pi/8:3\pi/8','3\pi/8:4\pi/8','4\pi/8:5\pi/8','5\pi/8:6\pi/8','6\pi/8:7\pi/8','7\pi/8:\pi','\pi:9\pi/8','9\pi/8:10\pi/8','10\pi/8:11\pi/8','11\pi/8:12\pi/8','12\pi/8:13\pi/8','13\pi/8:14\pi/8','14\pi/8:15\pi/8','15\pi/8:2\pi'};
        xTickLabels = {'\pi:7\pi/8','7\pi/8:6\pi/8','6\pi/8:5\pi/8','5\pi/8:4\pi/8','4\pi/8:3\pi/8','3\pi/8:2\pi/8','2\pi/8:\pi/8','\pi/8:0','0:-\pi/8','-\pi/8:-2\pi/8','-2\pi/8:-3\pi/8','-3\pi/8:-4\pi/8','-4\pi/8:-5\pi/8','-5\pi/8:-6\pi/8','-6\pi/8:-7\pi/8','-7\pi/8:-\pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver8{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p01to0p05HzPiOver8{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p01to0p05HzPiOver8VarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi
        %xTickLabels = {'0:\pi','\pi:2\pi'};
        xTickLabels = {'\pi:0','0:-\pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_pi{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p1to0p5HzPi{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiVarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: 2*pi/3
        %xTickLabels = {'0:2\pi/3','2\pi/3:4\pi/3','4\pi/3:2\pi'};
        xTickLabels = {'\pi:\pi/3','\pi/3:-\pi/3','-\pi/3:-\pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_2piOver3{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p1to0p5Hz2piOver3{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5Hz2PiOver3VarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/2
        %xTickLabels = {'0:\pi/2','\pi/2:\pi','\pi:3\pi/2','3\pi/2:2\pi'};
        xTickLabels = {'\pi:\pi/2','\pi/2:0','0:-\pi/2','-\pi/2:-pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver2{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p1to0p5HzPiOver2{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver2VarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/3
        %xTickLabels = {'0:\pi/3','\pi/3:2\pi/3','2\pi/3:\pi','\pi:4\pi/3','4\pi/3:5\pi/3','5\pi/3:2\pi'};
        xTickLabels = {'\pi:2\pi/3','2\pi/3:\pi/3','\pi/3:0','0:-\pi/3','-\pi/3:-2\pi/3','-2\pi/3:-pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver3{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p1to0p5HzPiOver3{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver3VarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/4
        %xTickLabels = {'0:\pi/4','\pi/4:2\pi/4','2\pi/4:3\pi/4','3\pi/4:\pi','\pi:5\pi/4','5\pi/4:6\pi/4','6\pi/4:7\pi/4','7\pi/4:2\pi'};
        xTickLabels = {'\pi:3\pi/4','3\pi/4:\pi/2','\pi/2:\pi/4','\pi/4:0','0:-\pi/4','-\pi/4:-\pi/2','-\pi/2:-3\pi/4','-3\pi/4:-\pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver4{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p1to0p5HzPiOver4{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver4VarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/6
        %xTickLabels = {'0:\pi/6','\pi/6:2\pi/6','2\pi/6:3\pi/6','3\pi/6:4\pi/6','4\pi/6:5\pi/6','5\pi/6:\pi','\pi:7\pi/6','7\pi/6:8\pi/6','8\pi/6:9\pi/6','9\pi/6:10\pi/6','10\pi/6:11\pi/6','11\pi/6:2\pi'};
        xTickLabels = {'\pi:5\pi/6','5\pi/6:4\pi/6','4\pi/6:3\pi/6','3\pi/6:2\pi/6','2\pi/6:\pi/6','\pi/6:0','0:-\pi/6','-\pi/6:-2\pi/6','-2\pi/6:-3\pi/6','-3\pi/6:-4\pi/6','-4\pi/6:-5\pi/6','-5\pi/6:-\pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver6{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p1to0p5HzPiOver6{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver6VarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/8
        %xTickLabels = {'0:\pi/8','\pi/8:2\pi/8','2\pi/8:3\pi/8','3\pi/8:4\pi/8','4\pi/8:5\pi/8','5\pi/8:6\pi/8','6\pi/8:7\pi/8','7\pi/8:\pi','\pi:9\pi/8','9\pi/8:10\pi/8','10\pi/8:11\pi/8','11\pi/8:12\pi/8','12\pi/8:13\pi/8','13\pi/8:14\pi/8','14\pi/8:15\pi/8','15\pi/8:2\pi'};
        xTickLabels = {'\pi:7\pi/8','7\pi/8:6\pi/8','6\pi/8:5\pi/8','5\pi/8:4\pi/8','4\pi/8:3\pi/8','3\pi/8:2\pi/8','2\pi/8:\pi/8','\pi/8:0','0:-\pi/8','-\pi/8:-2\pi/8','-2\pi/8:-3\pi/8','-3\pi/8:-4\pi/8','-4\pi/8:-5\pi/8','-5\pi/8:-6\pi/8','-6\pi/8:-7\pi/8','-7\pi/8:-\pi'};
        firingRateBarPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver8{iCond}{iAreasOI(iArea)}),...
          statsVarFRFiltBP0p1to0p5HzPiOver8{iCond}{iAreasOI(iArea)}, ['FiringRateFiltBP0p1to0p5HzPiOver8VarBars_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options);
      end
    end
  end
end


%% PLOT THE LOG FIRING RATE DISTRIBUTION MEANS AND VARIANCES COMBINED
if drawMeansVars
  options = struct();
  options.xLabel = 'Phase (rad)';
  if strcmpi(subpop, 'all')
    options.yLabel = 'Log_{10}(firing rate)';
    options.colour = {'k','k'};
  elseif strcmpi(subpop, 'positive')
    %options.yLabel = 'Positive log_{10}(firing rate)';
    options.yLabel = 'Log_{10}(firing rate)';
    options.colour = {matlabColours(10), matlabColours(3)};
  elseif strcmpi(subpop, 'negative')
    %options.yLabel = 'Negative log_{10}(firing rate)';
    options.yLabel = 'Log_{10}(firing rate)';
    options.colour = {matlabColours(8), matlabColours(6)};
  end
  options.figFolder = mainFolder;
  options.figSize = [24 18];
  options.saveFig = true;
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if size(areaFRFiltBP0p01to0p05HzIndividual_pi{iCond}{iAreasOI(iArea)},1) > 1
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi
        %xTickLabels = {'0:\pi','\pi:2\pi'};
        xTickLabels = {'\pi:0','0:-\pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_pi{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p01to0p05HzPiMeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: 2*pi/3
        %xTickLabels = {'0:2\pi/3','2\pi/3:4\pi/3','4\pi/3:2\pi'};
        xTickLabels = {'\pi:\pi/3','\pi/3:-\pi/3','-\pi/3:-\pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_2piOver3{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p01to0p05Hz2piOver3MeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/2
        %xTickLabels = {'0:\pi/2','\pi/2:\pi','\pi:3\pi/2','3\pi/2:2\pi'};
        xTickLabels = {'\pi:\pi/2','\pi/2:0','0:-\pi/2','-\pi/2:-pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver2{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p01to0p05HzPiOver2MeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/3
        %xTickLabels = {'0:\pi/3','\pi/3:2\pi/3','2\pi/3:\pi','\pi:4\pi/3','4\pi/3:5\pi/3','5\pi/3:2\pi'};
        xTickLabels = {'\pi:2\pi/3','2\pi/3:\pi/3','\pi/3:0','0:-\pi/3','-\pi/3:-2\pi/3','-2\pi/3:-pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver3{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p01to0p05HzPiOver3MeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/4
        %xTickLabels = {'0:\pi/4','\pi/4:2\pi/4','2\pi/4:3\pi/4','3\pi/4:\pi','\pi:5\pi/4','5\pi/4:6\pi/4','6\pi/4:7\pi/4','7\pi/4:2\pi'};
        xTickLabels = {'\pi:3\pi/4','3\pi/4:\pi/2','\pi/2:\pi/4','\pi/4:0','0:-\pi/4','-\pi/4:-\pi/2','-\pi/2:-3\pi/4','-3\pi/4:-\pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver4{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p01to0p05HzPiOver4MeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/6
        %xTickLabels = {'0:\pi/6','\pi/6:2\pi/6','2\pi/6:3\pi/6','3\pi/6:4\pi/6','4\pi/6:5\pi/6','5\pi/6:\pi','\pi:7\pi/6','7\pi/6:8\pi/6','8\pi/6:9\pi/6','9\pi/6:10\pi/6','10\pi/6:11\pi/6','11\pi/6:2\pi'};
        xTickLabels = {'\pi:5\pi/6','5\pi/6:4\pi/6','4\pi/6:3\pi/6','3\pi/6:2\pi/6','2\pi/6:\pi/6','\pi/6:0','0:-\pi/6','-\pi/6:-2\pi/6','-2\pi/6:-3\pi/6','-3\pi/6:-4\pi/6','-4\pi/6:-5\pi/6','-5\pi/6:-\pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver6{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p01to0p05HzPiOver6MeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
        
        % Band-pass filtered with a bandwidth of 0.01-0.05 Hz: pi/8
        %xTickLabels = {'0:\pi/8','\pi/8:2\pi/8','2\pi/8:3\pi/8','3\pi/8:4\pi/8','4\pi/8:5\pi/8','5\pi/8:6\pi/8','6\pi/8:7\pi/8','7\pi/8:\pi','\pi:9\pi/8','9\pi/8:10\pi/8','10\pi/8:11\pi/8','11\pi/8:12\pi/8','12\pi/8:13\pi/8','13\pi/8:14\pi/8','14\pi/8:15\pi/8','15\pi/8:2\pi'};
        xTickLabels = {'\pi:7\pi/8','7\pi/8:6\pi/8','6\pi/8:5\pi/8','5\pi/8:4\pi/8','4\pi/8:3\pi/8','3\pi/8:2\pi/8','2\pi/8:\pi/8','\pi/8:0','0:-\pi/8','-\pi/8:-2\pi/8','-2\pi/8:-3\pi/8','-3\pi/8:-4\pi/8','-4\pi/8:-5\pi/8','-5\pi/8:-6\pi/8','-6\pi/8:-7\pi/8','-7\pi/8:-\pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p01to0p05HzIndividual_piOver8{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p01to0p05HzPiOver8MeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi
        %xTickLabels = {'0:\pi','\pi:2\pi'};
        xTickLabels = {'\pi:0','0:-\pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_pi{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p1to0p5HzPiMeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: 2*pi/3
        %xTickLabels = {'0:2\pi/3','2\pi/3:4\pi/3','4\pi/3:2\pi'};
        xTickLabels = {'\pi:\pi/3','\pi/3:-\pi/3','-\pi/3:-\pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_2piOver3{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p1to0p5Hz2piOver3MeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/2
        %xTickLabels = {'0:\pi/2','\pi/2:\pi','\pi:3\pi/2','3\pi/2:2\pi'};
        xTickLabels = {'\pi:\pi/2','\pi/2:0','0:-\pi/2','-\pi/2:-pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver2{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p1to0p5HzPiOver2MeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/3
        %xTickLabels = {'0:\pi/3','\pi/3:2\pi/3','2\pi/3:\pi','\pi:4\pi/3','4\pi/3:5\pi/3','5\pi/3:2\pi'};
        xTickLabels = {'\pi:2\pi/3','2\pi/3:\pi/3','\pi/3:0','0:-\pi/3','-\pi/3:-2\pi/3','-2\pi/3:-pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver3{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p1to0p5HzPiOver3MeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/4
        %xTickLabels = {'0:\pi/4','\pi/4:2\pi/4','2\pi/4:3\pi/4','3\pi/4:\pi','\pi:5\pi/4','5\pi/4:6\pi/4','6\pi/4:7\pi/4','7\pi/4:2\pi'};
        xTickLabels = {'\pi:3\pi/4','3\pi/4:\pi/2','\pi/2:\pi/4','\pi/4:0','0:-\pi/4','-\pi/4:-\pi/2','-\pi/2:-3\pi/4','-3\pi/4:-\pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver4{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p1to0p5HzPiOver4MeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/6
        %xTickLabels = {'0:\pi/6','\pi/6:2\pi/6','2\pi/6:3\pi/6','3\pi/6:4\pi/6','4\pi/6:5\pi/6','5\pi/6:\pi','\pi:7\pi/6','7\pi/6:8\pi/6','8\pi/6:9\pi/6','9\pi/6:10\pi/6','10\pi/6:11\pi/6','11\pi/6:2\pi'};
        xTickLabels = {'\pi:5\pi/6','5\pi/6:4\pi/6','4\pi/6:3\pi/6','3\pi/6:2\pi/6','2\pi/6:\pi/6','\pi/6:0','0:-\pi/6','-\pi/6:-2\pi/6','-2\pi/6:-3\pi/6','-3\pi/6:-4\pi/6','-4\pi/6:-5\pi/6','-5\pi/6:-\pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver6{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p1to0p5HzPiOver6MeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
        
        % Band-pass filtered with a bandwidth of 0.1-0.5 Hz: pi/8
        %xTickLabels = {'0:\pi/8','\pi/8:2\pi/8','2\pi/8:3\pi/8','3\pi/8:4\pi/8','4\pi/8:5\pi/8','5\pi/8:6\pi/8','6\pi/8:7\pi/8','7\pi/8:\pi','\pi:9\pi/8','9\pi/8:10\pi/8','10\pi/8:11\pi/8','11\pi/8:12\pi/8','12\pi/8:13\pi/8','13\pi/8:14\pi/8','14\pi/8:15\pi/8','15\pi/8:2\pi'};
        xTickLabels = {'\pi:7\pi/8','7\pi/8:6\pi/8','6\pi/8:5\pi/8','5\pi/8:4\pi/8','4\pi/8:3\pi/8','3\pi/8:2\pi/8','2\pi/8:\pi/8','\pi/8:0','0:-\pi/8','-\pi/8:-2\pi/8','-2\pi/8:-3\pi/8','-3\pi/8:-4\pi/8','-4\pi/8:-5\pi/8','-5\pi/8:-6\pi/8','-6\pi/8:-7\pi/8','-7\pi/8:-\pi'};
        firingRateCombinedDotPlot(xTickLabels, getLog(areaFRFiltBP0p1to0p5HzIndividual_piOver8{iCond}{iAreasOI(iArea)}),...
          ['FiringRateFiltBP0p1to0p5HzPiOver8MeanVarDots_' areas{iAreasOI(iArea)} '_' conditions{iCond}], options)
      end
    end
  end
end


%% FIRING RATE MODULATION DISTRIBUTIONS
if drawModulationDistros
  options = struct();
  options.pdf = true;
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      
      % 50 cent
      if ~isempty(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)})
        
        % Signed modulations
        frModulations = areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,1);
        nBins = 96;
        edges = -30:60/nBins:30;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['signedModulationDistro50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = 'Unit firing rate change (APs/s)';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
        
        % Unsigned modulations
        frModulations = abs(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,1));
        nBins = nBins/2;
        edges = 0:30/nBins:30;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['unsignedModulationDistro50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = '|Unit firing rate change| (APs/s)';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
        
        % Signed normalised modulations
        frModulations = (areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,1))./areaFRIndividual{iCond}{iAreasOI(iArea)};
        nBins = 96;
        edges = -3:6/nBins:3;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['signedNormalisedModulationDistro50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = 'Normalised unit firing rate change';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
        
        % Unsigned normalised modulations
        frModulations = abs((areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,1))./areaFRIndividual{iCond}{iAreasOI(iArea)});
        nBins = nBins/2;
        edges = 0:3/nBins:3;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['unsignedNormalisedModulationDistro50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = '|Normalised unit firing rate change|';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
      end
      
      % 25 cent
      if ~isempty(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)})
        
        % Signed modulations
        frModulations = areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,1);
        nBins = 96;
        edges = -30:60/nBins:30;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['signedModulationDistro25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = 'Unit firing rate change (APs/s)';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
        
        % Unsigned modulations
        frModulations = abs(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,1));
        nBins = nBins/2;
        edges = 0:30/nBins:30;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['unsignedModulationDistro25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = '|Unit firing rate change| (APs/s)';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
        
        % Signed normalised modulations
        frModulations = (areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,1))./areaFRIndividual{iCond}{iAreasOI(iArea)};
        nBins = 96;
        edges = -3:6/nBins:3;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['signedNormalisedModulationDistro25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = 'Normalised unit firing rate change';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
        
        % Unsigned modulations
        frModulations = abs((areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,1))./areaFRIndividual{iCond}{iAreasOI(iArea)});
        nBins = nBins/2;
        edges = 0:3/nBins:3;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['unsignedNormalisedModulationDistro25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = '|Normalised unit firing rate change|';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
      end
      
      % 12.5 cent
      if ~isempty(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)})
        
        % Signed modulations
        frModulations = areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,1);
        nBins = 96;
        edges = -30:60/nBins:30;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['signedModulationDistro12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = 'Unit firing rate change (APs/s)';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
        
        % Unsigned modulations
        frModulations = abs(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,1));
        nBins = nBins/2;
        edges = 0:30/nBins:30;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['unsignedModulationDistro12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = '|Unit firing rate change| (APs/s)';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
        
        % Signed normalised modulations
        frModulations = (areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,1))./areaFRIndividual{iCond}{iAreasOI(iArea)};
        nBins = 96;
        edges = -3:6/nBins:3;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['signedNormalisedModulationDistro12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = 'Normalised unit firing rate change';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
        
        % Unsigned modulations
        frModulations = abs((areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,1))./areaFRIndividual{iCond}{iAreasOI(iArea)});
        nBins = nBins/2;
        edges = 0:3/nBins:3;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['unsignedNormalisedModulationDistro12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = '|Normalised unit firing rate change|';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
      end
      
      % 33.33 cent
      if ~isempty(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)})
        
        % Signed modulations
        frModulations = areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,1);
        nBins = 96;
        edges = -30:60/nBins:30;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['signedModulationDistroThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = 'Unit firing rate change (APs/s)';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
        
        % Unsigned modulations
        frModulations = abs(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,1));
        nBins = nBins/2;
        edges = 0:30/nBins:30;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['unsignedModulationDistroThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = '|Unit firing rate change| (APs/s)';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
        
        % Signed normalised modulations
        frModulations = (areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,1))./areaFRIndividual{iCond}{iAreasOI(iArea)};
        nBins = 96;
        edges = -3:6/nBins:3;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['signedNormalisedModulationDistroThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = 'Normalised unit firing rate change';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
        
        % Unsigned modulations
        frModulations = abs((areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,1))./areaFRIndividual{iCond}{iAreasOI(iArea)});
        nBins = nBins/2;
        edges = 0:3/nBins:3;
        frModulations = histcounts(frModulations, edges);
        options.figName = ['unsignedNormalisedModulationDistroThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options.figFolder = mainFolder;
        options.xLabel = '|Normalised unit firing rate change|';
        options.yLabel = 'Probability';
        options.figSize = figSize;
        options.saveFig = true;
        options.splitModes = false;
        options.dataScatter = [];
        histPlot(edges, frModulations, options);
      end
    end
  end
end


%% FIRING RATE MODULATION CORRELATIONS
options = struct();
options.corrType = 'Spearman';
options.figFolder = mainFolder;
options.diagonal = true;
options.fitLine = true;
options.fitLineDisplay = 'on';
options.xLim = [];
options.xTicks = [];
options.yLim = [];
options.yTicks = [];
options.figSize = figSize;
options.saveFig = true; %#ok<*STRNU>
if drawCorrelations
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      fr = areaFRIndividual{iCond}{iAreasOI(iArea)};
      rSpearman = arearSpearman5secIndividual{iCond}{iAreasOI(iArea)};
      pvalSpearman = areapvalSpearman5secIndividual{iCond}{iAreasOI(iArea)};
      options2 = options;
      options2.colourVector = zeros(size(rSpearman));
      options2.colourVector(rSpearman >= 0 & pvalSpearman < 0.05) = 1;
      options2.colourVector(rSpearman >= 0 & pvalSpearman >= 0.05) = 2;
      options2.colourVector(rSpearman < 0 & pvalSpearman < 0.05) = 3;
      options2.colourVector(rSpearman < 0 & pvalSpearman >= 0.05) = 4;
      inds = logical(options2.colourVector);
      options2.colourVector = options2.colourVector(inds);
      options2.colourCodes = zeros(numel(unique(options2.colourVector)), 3);
      options2.colourCodes(1,:) = matlabColours(10); % Red
      options2.colourCodes(2,:) = [255 128 128]./255; % Light red
      options2.colourCodes(3,:) = matlabColours(8); % Blue
      options2.colourCodes(4,:) = [128 128 255]./255; % Light blue
      
      % 50 cent
      if ~isempty(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)})
        
        % Signed modulations: Firing rate
        frModulations = areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,1);
        options2.figTitle = ['Unit firing rate vs firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrSignedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = 'Unit firing rate change (APs/s)';
        options2.axesType = 'semilogx';
        [~, rCoefFRVsFRSignedModulation50{iCond}{iAreasOI(iArea)}, pvalFRVsFRSignedModulation50{iCond}{iAreasOI(iArea)},...
          modelFRVsFRSignedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrSignedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Unit firing rate change (APs/s)';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRSignedModulation50{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRSignedModulation50{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRSignedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit firing rate vs absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrUnsignedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = '|Unit firing rate change| (APs/s)';
        options2.axesType = 'loglog';
        [~, rCoefFRVsFRUnsignedModulation50{iCond}{iAreasOI(iArea)}, pvalFRVsFRUnsignedModulation50{iCond}{iAreasOI(iArea)},...
          modelFRVsFRUnsignedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Unsigned modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrUnsignedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Unit firing rate change| (APs/s)';
        options2.axesType = 'semilogy';
        [~, rCoefrSpearman5secVsFRUnsignedModulation50{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRUnsignedModulation50{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRUnsignedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Signed normalised modulations: Firing rate
        frModulations = (areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,1))./areaFRIndividual{iCond}{iAreasOI(iArea)};
        options2.figTitle = ['Unit firing rate vs normalised firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrSignedNormalisedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = 'Normalised unit firing rate change';
        options2.axesType = 'semilogx';
        [~, rCoefFRVsFRSignedNormalisedModulation50{iCond}{iAreasOI(iArea)}, pvalFRVsFRSignedNormalisedModulation50{iCond}{iAreasOI(iArea)},...
          modelFRVsFRSignedNormalisedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrSignedNormalisedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Normalised unit firing rate change';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRSignedNormalisedModulation50{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRSignedNormalisedModulation50{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRSignedNormalisedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned normalised modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit firing rate vs normalised absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrUnsignedNormalisedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = '|Normalised unit firing rate change|';
        options2.axesType = 'loglog';
        [~, rCoefFRVsFRUnsignedNormalisedModulation50{iCond}{iAreasOI(iArea)}, pvalFRVsFRUnsignedNormalisedModulation50{iCond}{iAreasOI(iArea)},...
          modelFRVsFRUnsignedNormalisedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Unsigned normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrUnsignedNormalisedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Normalised unit firing rate change|';
        options2.axesType = 'semilogy';
        [~, rCoefrSpearman5secVsFRUnsignedNormalisedModulation50{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRUnsignedNormalisedModulation50{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRUnsignedNormalisedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Firing rate: constricted vs dilated pupil
        options2.figTitle = ['Unit firing rate associated to constricted pupil phase vs unit firing rate associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frConstrictedVsfrDilated50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Constricted firing rate';
        options2.yLabel = 'Dilated firing rate';
        options2.axesType = 'loglog';
        [~, rCoefFRConstrictedVsFRDilated50{iCond}{iAreasOI(iArea)}, pvalFRConstrictedVsFRDilated50{iCond}{iAreasOI(iArea)},...
          modelFRConstrictedVsFRDilated50{iCond}{iAreasOI(iArea)}] = corrPlot(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(inds,1), areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(inds,end), options2);
      end
      
      % 25 cent
      if ~isempty(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)})
        
        % Signed modulations: Firing rate
        frModulations = areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,1);
        options2.figTitle = ['Unit firing rate vs firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrSignedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = 'Unit firing rate change (APs/s)';
        options2.axesType = 'semilogx';
        [~, rCoefFRVsFRSignedModulation25{iCond}{iAreasOI(iArea)}, pvalFRVsFRSignedModulation25{iCond}{iAreasOI(iArea)},...
          modelFRVsFRSignedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrSignedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Unit firing rate change (APs/s)';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRSignedModulation25{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRSignedModulation25{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRSignedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit firing rate vs absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrUnsignedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = '|Unit firing rate change| (APs/s)';
        options2.axesType = 'loglog';
        [~, rCoefFRVsFRUnsignedModulation25{iCond}{iAreasOI(iArea)}, pvalFRVsFRUnsignedModulation25{iCond}{iAreasOI(iArea)},...
          modelFRVsFRUnsignedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Unsigned modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrUnsignedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Unit firing rate change| (APs/s)';
        options2.axesType = 'semilogy';
        [~, rCoefrSpearman5secVsFRUnsignedModulation25{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRUnsignedModulation25{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRUnsignedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Signed normalised modulations: Firing rate
        frModulations = (areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,1))./areaFRIndividual{iCond}{iAreasOI(iArea)};
        options2.figTitle = ['Unit firing rate vs normalised firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrSignedNormalisedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = 'Normalised unit firing rate change';
        options2.axesType = 'semilogx';
        [~, rCoefFRVsFRSignedNormalisedModulation25{iCond}{iAreasOI(iArea)}, pvalFRVsFRSignedNormalisedModulation25{iCond}{iAreasOI(iArea)},...
          modelFRVsFRSignedNormalisedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrSignedNormalisedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Normalised unit firing rate change';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRSignedNormalisedModulation25{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRSignedNormalisedModulation25{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRSignedNormalisedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned normalised modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit firing rate vs normalised absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrUnsignedNormalisedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = '|Normalised unit firing rate change|';
        options2.axesType = 'loglog';
        [~, rCoefFRVsFRUnsignedNormalisedModulation25{iCond}{iAreasOI(iArea)}, pvalFRVsFRUnsignedNormalisedModulation25{iCond}{iAreasOI(iArea)},...
          modelFRVsFRUnsignedNormalisedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Unsigned normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrUnsignedNormalisedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Normalised unit firing rate change|';
        options2.axesType = 'semilogy';
        [~, rCoefrSpearman5secVsFRUnsignedNormalisedModulation25{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRUnsignedNormalisedModulation25{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRUnsignedNormalisedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Firing rate: constricted vs dilated pupil
        options2.figTitle = ['Unit firing rate associated to constricted pupil phase vs unit firing rate associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frConstrictedVsfrDilated25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Constricted firing rate';
        options2.yLabel = 'Dilated firing rate';
        options2.axesType = 'loglog';
        [~, rCoefFRConstrictedVsFRDilated25{iCond}{iAreasOI(iArea)}, pvalFRConstrictedVsFRDilated25{iCond}{iAreasOI(iArea)},...
          modelFRConstrictedVsFRDilated25{iCond}{iAreasOI(iArea)}] = corrPlot(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(inds,1), areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(inds,end), options2);
      end
      
      % 12.5 cent
      if ~isempty(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)})
        
        % Signed modulations: Firing rate
        frModulations = areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,1);
        options2.figTitle = ['Unit firing rate vs firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrSignedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = 'Unit firing rate change (APs/s)';
        options2.axesType = 'semilogx';
        [~, rCoefFRVsFRSignedModulation12p5{iCond}{iAreasOI(iArea)}, pvalFRVsFRSignedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelFRVsFRSignedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrSignedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Unit firing rate change (APs/s)';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRSignedModulation12p5{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRSignedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRSignedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit firing rate vs absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrUnsignedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = '|Unit firing rate change| (APs/s)';
        options2.axesType = 'loglog';
        [~, rCoefFRVsFRUnsignedModulation12p5{iCond}{iAreasOI(iArea)}, pvalFRVsFRUnsignedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelFRVsFRUnsignedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Unsigned modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrUnsignedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Unit firing rate change| (APs/s)';
        options2.axesType = 'semilogy';
        [~, rCoefrSpearman5secVsFRUnsignedModulation12p5{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRUnsignedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRUnsignedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Signed normalised modulations: Firing rate
        frModulations = (areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,1))./areaFRIndividual{iCond}{iAreasOI(iArea)};
        options2.figTitle = ['Unit firing rate vs normalised firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrSignedNormalisedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = 'Normalised unit firing rate change';
        options2.axesType = 'semilogx';
        [~, rCoefFRVsFRSignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}, pvalFRVsFRSignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelFRVsFRSignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrSignedNormalisedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Normalised unit firing rate change';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRSignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRSignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRSignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned normalised modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit firing rate vs normalised absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrUnsignedNormalisedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = '|Normalised unit firing rate change|';
        options2.axesType = 'loglog';
        [~, rCoefFRVsFRUnsignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}, pvalFRVsFRUnsignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelFRVsFRUnsignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Unsigned normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrUnsignedNormalisedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Normalised unit firing rate change|';
        options2.axesType = 'semilogy';
        [~, rCoefrSpearman5secVsFRUnsignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRUnsignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRUnsignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Firing rate: constricted vs dilated pupil
        options2.figTitle = ['Unit firing rate associated to constricted pupil phase vs unit firing rate associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frConstrictedVsfrDilated12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Constricted firing rate';
        options2.yLabel = 'Dilated firing rate';
        options2.axesType = 'loglog';
        [~, rCoefFRConstrictedVsFRDilated12p5{iCond}{iAreasOI(iArea)}, pvalFRConstrictedVsFRDilated12p5{iCond}{iAreasOI(iArea)},...
          modelFRConstrictedVsFRDilated12p5{iCond}{iAreasOI(iArea)}] = corrPlot(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(inds,1), areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(inds,end), options2);
      end
      
      % 33.33 cent
      if ~isempty(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)})
        
        % Signed modulations: Firing rate
        frModulations = areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,1);
        options2.figTitle = ['Unit firing rate vs firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrSignedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = 'Unit firing rate change (APs/s)';
        options2.axesType = 'semilogx';
        [~, rCoefFRVsFRSignedModulationThird{iCond}{iAreasOI(iArea)}, pvalFRVsFRSignedModulationThird{iCond}{iAreasOI(iArea)},...
          modelFRVsFRSignedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrSignedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Unit firing rate change (APs/s)';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRSignedModulationThird{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRSignedModulationThird{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRSignedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit firing rate vs absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrUnsignedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = '|Unit firing rate change| (APs/s)';
        options2.axesType = 'loglog';
        [~, rCoefFRVsFRUnsignedModulationThird{iCond}{iAreasOI(iArea)}, pvalFRVsFRUnsignedModulationThird{iCond}{iAreasOI(iArea)},...
          modelFRVsFRUnsignedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Unsigned modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrUnsignedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Unit firing rate change| (APs/s)';
        options2.axesType = 'semilogy';
        [~, rCoefrSpearman5secVsFRUnsignedModulationThird{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRUnsignedModulationThird{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRUnsignedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Signed normalised modulations: Firing rate
        frModulations = (areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,1))./areaFRIndividual{iCond}{iAreasOI(iArea)};
        options2.figTitle = ['Unit firing rate vs normalised firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrSignedNormalisedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = 'Normalised unit firing rate change';
        options2.axesType = 'semilogx';
        [~, rCoefFRVsFRSignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}, pvalFRVsFRSignedNormalisedModulationThird{iCond}{iAreasOI(iArea)},...
          modelFRVsFRSignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrSignedNormalisedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Normalised unit firing rate change';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRSignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRSignedNormalisedModulationThird{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRSignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned normalised modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit firing rate vs normalised absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frVsfrUnsignedNormalisedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit firing rate (APs/s)';
        options2.yLabel = '|Normalised unit firing rate change|';
        options2.axesType = 'loglog';
        [~, rCoefFRVsFRUnsignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}, pvalFRVsFRUnsignedNormalisedModulationThird{iCond}{iAreasOI(iArea)},...
          modelFRVsFRUnsignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Unsigned normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised absolute firing rate modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrUnsignedNormalisedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Normalised unit firing rate change|';
        options2.axesType = 'semilogy';
        [~, rCoefrSpearman5secVsFRUnsignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRUnsignedNormalisedModulationThird{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRUnsignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Firing rate: constricted vs dilated pupil
        options2.figTitle = ['Unit firing rate associated to constricted pupil phase vs unit firing rate associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frConstrictedVsfrDilatedThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Constricted firing rate';
        options2.yLabel = 'Dilated firing rate';
        options2.axesType = 'loglog';
        [~, rCoefFRConstrictedVsFRDilatedThird{iCond}{iAreasOI(iArea)}, pvalFRConstrictedVsFRDilatedThird{iCond}{iAreasOI(iArea)},...
          modelFRConstrictedVsFRDilatedThird{iCond}{iAreasOI(iArea)}] = corrPlot(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(inds,1), areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(inds,end), options2);
      end
      
      % Firing rate: r
      options2.figTitle = ['Unit rCoef over 5sec bins vs firing rate in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
      options2.figName = ['rSpearman5secVsfr in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
      options2.xLabel = 'r coef';
      options2.yLabel = 'Unit firing rate (APs/s)';
      options2.axesType = 'semilogy';
      [~, rCoefrSpearman5secVsFR{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFR{iCond}{iAreasOI(iArea)},...
        modelrSpearman5secVsFR{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), fr(inds), options2);
      
      % Firing rate: unsigned r
      options2.figTitle = ['Unit rCoef over 5sec bins vs firing rate in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
      options2.figName = ['rSpearmanUnsigned5secVsfr in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
      options2.xLabel = '|r coef|';
      options2.yLabel = 'Unit firing rate (APs/s)';
      options2.axesType = 'semilogy';
      [~, rCoefrSpearmanUnsigned5secVsFR{iCond}{iAreasOI(iArea)}, pvalrSpearmanUnsigned5secVsFR{iCond}{iAreasOI(iArea)},...
        modelrSpearmanUnsigned5secVsFR{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), fr(inds), options2);
    end
  end
  
  % Save correlation data
  save(filename, 'rCoefFRVsFRSignedModulation50','pvalFRVsFRSignedModulation50','rCoefrSpearman5secVsFRSignedModulation50','pvalrSpearman5secVsFRSignedModulation50',...
    'rCoefFRVsFRUnsignedModulation50','pvalFRVsFRUnsignedModulation50','rCoefrSpearman5secVsFRUnsignedModulation50','pvalrSpearman5secVsFRUnsignedModulation50',...
    'rCoefFRVsFRSignedNormalisedModulation50','pvalFRVsFRSignedNormalisedModulation50','rCoefrSpearman5secVsFRSignedNormalisedModulation50','pvalrSpearman5secVsFRSignedNormalisedModulation50',...
    'rCoefFRVsFRUnsignedNormalisedModulation50','pvalFRVsFRUnsignedNormalisedModulation50','rCoefrSpearman5secVsFRUnsignedNormalisedModulation50','pvalrSpearman5secVsFRUnsignedNormalisedModulation50',...
    'rCoefFRVsFRSignedModulation25','pvalFRVsFRSignedModulation25','rCoefrSpearman5secVsFRSignedModulation25','pvalrSpearman5secVsFRSignedModulation25',...
    'rCoefFRVsFRUnsignedModulation25','pvalFRVsFRUnsignedModulation25','rCoefrSpearman5secVsFRUnsignedModulation25','pvalrSpearman5secVsFRUnsignedModulation25',...
    'rCoefFRVsFRSignedNormalisedModulation25','pvalFRVsFRSignedNormalisedModulation25','rCoefrSpearman5secVsFRSignedNormalisedModulation25','pvalrSpearman5secVsFRSignedNormalisedModulation25',...
    'rCoefFRVsFRUnsignedNormalisedModulation25','pvalFRVsFRUnsignedNormalisedModulation25','rCoefrSpearman5secVsFRUnsignedNormalisedModulation25','pvalrSpearman5secVsFRUnsignedNormalisedModulation25',...
    'rCoefFRVsFRSignedModulation12p5','pvalFRVsFRSignedModulation12p5','rCoefrSpearman5secVsFRSignedModulation12p5','pvalrSpearman5secVsFRSignedModulation12p5',...
    'rCoefFRVsFRUnsignedModulation12p5','pvalFRVsFRUnsignedModulation12p5','rCoefrSpearman5secVsFRUnsignedModulation12p5','pvalrSpearman5secVsFRUnsignedModulation12p5',...
    'rCoefFRVsFRSignedNormalisedModulation12p5','pvalFRVsFRSignedNormalisedModulation12p5','rCoefrSpearman5secVsFRSignedNormalisedModulation12p5','pvalrSpearman5secVsFRSignedNormalisedModulation12p5',...
    'rCoefFRVsFRUnsignedNormalisedModulation12p5','pvalFRVsFRUnsignedNormalisedModulation12p5','rCoefrSpearman5secVsFRUnsignedNormalisedModulation12p5','pvalrSpearman5secVsFRUnsignedNormalisedModulation12p5',...
    'rCoefFRVsFRSignedModulationThird','pvalFRVsFRSignedModulationThird','rCoefrSpearman5secVsFRSignedModulationThird','pvalrSpearman5secVsFRSignedModulationThird',...
    'rCoefFRVsFRUnsignedModulationThird','pvalFRVsFRUnsignedModulationThird','rCoefrSpearman5secVsFRUnsignedModulationThird','pvalrSpearman5secVsFRUnsignedModulationThird',...
    'rCoefFRVsFRSignedNormalisedModulationThird','pvalFRVsFRSignedNormalisedModulationThird','rCoefrSpearman5secVsFRSignedNormalisedModulationThird','pvalrSpearman5secVsFRSignedNormalisedModulationThird',...
    'rCoefFRVsFRUnsignedNormalisedModulationThird','pvalFRVsFRUnsignedNormalisedModulationThird','rCoefrSpearman5secVsFRUnsignedNormalisedModulationThird','pvalrSpearman5secVsFRUnsignedNormalisedModulationThird',...
    'rCoefrSpearman5secVsFR','pvalrSpearman5secVsFR','rCoefrSpearmanUnsigned5secVsFR','pvalrSpearmanUnsigned5secVsFR',...
    'modelFRVsFRSignedModulation50','modelrSpearman5secVsFRSignedModulation50','modelFRVsFRUnsignedModulation50','modelrSpearman5secVsFRUnsignedModulation50','modelFRVsFRSignedNormalisedModulation50',...
    'modelrSpearman5secVsFRSignedNormalisedModulation50','modelFRVsFRUnsignedNormalisedModulation50','modelrSpearman5secVsFRUnsignedNormalisedModulation50','modelFRVsFRSignedModulation25',...
    'modelrSpearman5secVsFRSignedModulation25','modelFRVsFRUnsignedModulation25','modelrSpearman5secVsFRUnsignedModulation25','modelFRVsFRSignedNormalisedModulation25',...
    'modelrSpearman5secVsFRSignedNormalisedModulation25','modelrSpearman5secVsFRSignedNormalisedModulation25','modelFRVsFRUnsignedNormalisedModulation25','modelrSpearman5secVsFRUnsignedNormalisedModulation25',...
    'modelFRVsFRSignedModulation12p5','modelrSpearman5secVsFRSignedModulation12p5','modelFRVsFRUnsignedModulation12p5','modelrSpearman5secVsFRUnsignedModulation12p5','modelFRVsFRSignedNormalisedModulation12p5',...
    'modelrSpearman5secVsFRSignedNormalisedModulation12p5','modelFRVsFRUnsignedNormalisedModulation12p5','modelrSpearman5secVsFRUnsignedNormalisedModulation12p5','modelFRVsFRSignedModulationThird',...
    'modelrSpearman5secVsFRSignedModulationThird','modelFRVsFRUnsignedModulationThird','modelrSpearman5secVsFRUnsignedModulationThird','modelFRVsFRSignedNormalisedModulationThird',...
    'modelrSpearman5secVsFRSignedNormalisedModulationThird','modelFRVsFRUnsignedNormalisedModulationThird','modelrSpearman5secVsFRUnsignedNormalisedModulationThird',...
    'modelrSpearman5secVsFR','modelrSpearmanUnsigned5secVsFR',...
    'rCoefFRConstrictedVsFRDilated50','pvalFRConstrictedVsFRDilated50','modelFRConstrictedVsFRDilated50','rCoefFRConstrictedVsFRDilated25','pvalFRConstrictedVsFRDilated25',...
    'modelFRConstrictedVsFRDilated25','rCoefFRConstrictedVsFRDilated12p5','pvalFRConstrictedVsFRDilated12p5','modelFRConstrictedVsFRDilated12p5','rCoefFRConstrictedVsFRDilatedThird',...
    'pvalFRConstrictedVsFRDilatedThird','modelFRConstrictedVsFRDilatedThird', '-append');
end


%% FIRING RATE MODULATION CORRELATIONS: LOG
options = struct();
options.corrType = 'Spearman';
options.figFolder = mainFolder;
options.diagonal = true;
options.fitLine = true;
options.fitLineDisplay = 'on';
options.xLim = [];
options.xTicks = [];
options.yLim = [];
options.yTicks = [];
options.figSize = figSize;
options.saveFig = true;
if drawCorrelationsLog
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      fr = getLog(areaFRIndividual{iCond}{iAreasOI(iArea)});
      rSpearman = arearSpearman5secIndividual{iCond}{iAreasOI(iArea)};
      pvalSpearman = areapvalSpearman5secIndividual{iCond}{iAreasOI(iArea)};
      options2 = options;
      options2.colourVector = zeros(size(rSpearman));
      options2.colourVector(rSpearman >= 0 & pvalSpearman < 0.05) = 1;
      options2.colourVector(rSpearman >= 0 & pvalSpearman >= 0.05) = 2;
      options2.colourVector(rSpearman < 0 & pvalSpearman < 0.05) = 3;
      options2.colourVector(rSpearman < 0 & pvalSpearman >= 0.05) = 4;
      inds = logical(options2.colourVector);
      options2.colourVector = options2.colourVector(inds);
      options2.colourCodes = zeros(numel(unique(options2.colourVector)), 3);
      options2.colourCodes(1,:) = matlabColours(10); % Red
      options2.colourCodes(2,:) = [255 128 128]./255; % Light red
      options2.colourCodes(3,:) = matlabColours(8); % Blue
      options2.colourCodes(4,:) = [128 128 255]./255; % Light blue
      
      % 50 cent
      if ~isempty(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)})
        
        % Signed modulations: Firing rate
        frModulations = getLog(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,end)) - getLog(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,1));
        options2.figTitle = ['Unit log10(firing rate) vs log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogSignedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogSignedModulation50{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogSignedModulation50{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogSignedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogSignedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogSignedModulation50{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogSignedModulation50{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogSignedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit log10(firing rate) vs absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = '|Unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedModulation50{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedModulation50{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        frModulations2 = getLog(abs(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,1)));
        options2.figTitle = ['Unit log10(firing rate) vs absolute log10(firing rate modulation) in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedModulation50_2 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Unit log_{10}(|firing rate change|)';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedModulation50_2{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedModulation50_2{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedModulation50_2{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations2(inds), options2);
        
        % Unsigned modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogUnsignedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogUnsignedModulation50{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogUnsignedModulation50{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogUnsignedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Signed normalised modulations: Firing rate
        frModulations = (getLog(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,end)) - getLog(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(:,1)))./getLog(areaFRIndividual{iCond}{iAreasOI(iArea)});
        options2.figTitle = ['Unit log10(firing rate) vs normalised log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogSignedNormalisedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Normalised unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogSignedNormalisedModulation50{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogSignedNormalisedModulation50{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogSignedNormalisedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogSignedNormalisedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Normalised unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogSignedNormalisedModulation50{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogSignedNormalisedModulation50{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogSignedNormalisedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned normalised modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit log10(firing rate) vs normalised absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedNormalisedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = '|Normalised unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedNormalisedModulation50{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedNormalisedModulation50{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedNormalisedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        frModulations2 = frModulations2./getLog(areaFRIndividual{iCond}{iAreasOI(iArea)});
        options2.figTitle = ['Unit log10(firing rate) vs normalised absolute log10(firing rate modulation) in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedNormalisedModulation50_2 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Normalised unit log_{10}(|firing rate change|)';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedNormalisedModulation50_2{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedNormalisedModulation50_2{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedNormalisedModulation50_2{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations2(inds), options2);
        
        % Unsigned normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogUnsignedNormalisedModulation50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Normalised unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogUnsignedNormalisedModulation50{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogUnsignedNormalisedModulation50{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogUnsignedNormalisedModulation50{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Firing rate: constricted vs dilated pupil
        options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogConstrictedVsfrLogDilated50 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Constricted log_{10}(firing rate)';
        options2.yLabel = 'Dilated log_{10}(firing rate)';
        options2.axesType = 'regular';
        [~, rCoefFRLogConstrictedVsFRLogDilated50{iCond}{iAreasOI(iArea)}, pvalFRLogConstrictedVsFRLogDilated50{iCond}{iAreasOI(iArea)},...
          modelFRLogConstrictedVsFRLogDilated50{iCond}{iAreasOI(iArea)}] = corrPlot(getLog(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(inds,1)), getLog(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(inds,end)), options2);
      end
      
      % 25 cent
      if ~isempty(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)})
        
        % Signed modulations: Firing rate
        frModulations = getLog(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,end)) - getLog(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,1));
        options2.figTitle = ['Unit log10(firing rate) vs log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogSignedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogSignedModulation25{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogSignedModulation25{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogSignedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogSignedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogSignedModulation25{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogSignedModulation25{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogSignedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit log10(firing rate) vs absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = '|Unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedModulation25{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedModulation25{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        frModulations2 = getLog(abs(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,1)));
        options2.figTitle = ['Unit log10(firing rate) vs absolute log10(firing rate modulation) in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedModulation25_2 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Unit log_{10}(|firing rate change|)';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedModulation25_2{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedModulation25_2{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedModulation25_2{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations2(inds), options2);
        
        % Unsigned modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogUnsignedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogUnsignedModulation25{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogUnsignedModulation25{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogUnsignedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Signed normalised modulations: Firing rate
        frModulations = (getLog(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,end)) - getLog(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(:,1)))./getLog(areaFRIndividual{iCond}{iAreasOI(iArea)});
        options2.figTitle = ['Unit log10(firing rate) vs normalised log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogSignedNormalisedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Normalised unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogSignedNormalisedModulation25{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogSignedNormalisedModulation25{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogSignedNormalisedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogSignedNormalisedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Normalised unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogSignedNormalisedModulation25{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogSignedNormalisedModulation25{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogSignedNormalisedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned normalised modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit log10(firing rate) vs normalised absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedNormalisedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = '|Normalised unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedNormalisedModulation25{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedNormalisedModulation25{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedNormalisedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        frModulations2 = frModulations2./getLog(areaFRIndividual{iCond}{iAreasOI(iArea)});
        options2.figTitle = ['Unit log10(firing rate) vs normalised absolute log10(firing rate modulation) in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedNormalisedModulation25_2 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Normalised unit log_{10}(|firing rate change|)';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedNormalisedModulation25_2{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedNormalisedModulation25_2{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedNormalisedModulation25_2{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations2(inds), options2);
        
        % Unsigned normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogUnsignedNormalisedModulation25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Normalised unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogUnsignedNormalisedModulation25{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogUnsignedNormalisedModulation25{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogUnsignedNormalisedModulation25{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Firing rate: constricted vs dilated pupil
        options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogConstrictedVsfrLogDilated25 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Constricted log_{10}(firing rate)';
        options2.yLabel = 'Dilated log_{10}(firing rate)';
        options2.axesType = 'regular';
        [~, rCoefFRLogConstrictedVsFRLogDilated25{iCond}{iAreasOI(iArea)}, pvalFRLogConstrictedVsFRLogDilated25{iCond}{iAreasOI(iArea)},...
          modelFRLogConstrictedVsFRLogDilated25{iCond}{iAreasOI(iArea)}] = corrPlot(getLog(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(inds,1)), getLog(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(inds,end)), options2);
      end
      
      % 12.5 cent
      if ~isempty(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)})
        
        % Signed modulations: Firing rate
        frModulations = getLog(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,end)) - getLog(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,1));
        options2.figTitle = ['Unit log10(firing rate) vs log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogSignedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogSignedModulation12p5{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogSignedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogSignedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogSignedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogSignedModulation12p5{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogSignedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogSignedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit log10(firing rate) vs absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = '|Unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedModulation12p5{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        frModulations2 = getLog(abs(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,1)));
        options2.figTitle = ['Unit log10(firing rate) vs absolute log10(firing rate modulation) in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedModulation12p5_2 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Unit log_{10}(|firing rate change|)';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedModulation12p5_2{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedModulation12p5_2{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedModulation12p5_2{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations2(inds), options2);
        
        % Unsigned modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogUnsignedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogUnsignedModulation12p5{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogUnsignedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogUnsignedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Signed normalised modulations: Firing rate
        frModulations = (getLog(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,end)) - getLog(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(:,1)))./getLog(areaFRIndividual{iCond}{iAreasOI(iArea)});
        options2.figTitle = ['Unit log10(firing rate) vs normalised log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogSignedNormalisedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Normalised unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogSignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogSignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogSignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogSignedNormalisedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Normalised unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogSignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogSignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogSignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned normalised modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit log10(firing rate) vs normalised absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedNormalisedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = '|Normalised unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        frModulations2 = frModulations2./getLog(areaFRIndividual{iCond}{iAreasOI(iArea)});
        options2.figTitle = ['Unit log10(firing rate) vs normalised absolute log10(firing rate modulation) in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedNormalisedModulation12p5_2 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Normalised unit log_{10}(|firing rate change|)';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedNormalisedModulation12p5_2{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedNormalisedModulation12p5_2{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedNormalisedModulation12p5_2{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations2(inds), options2);
        
        % Unsigned normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogUnsignedNormalisedModulation12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Normalised unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogUnsignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogUnsignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogUnsignedNormalisedModulation12p5{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Firing rate: constricted vs dilated pupil
        options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogConstrictedVsfrLogDilated12p5 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Constricted log_{10}(firing rate)';
        options2.yLabel = 'Dilated log_{10}(firing rate)';
        options2.axesType = 'regular';
        [~, rCoefFRLogConstrictedVsFRLogDilated12p5{iCond}{iAreasOI(iArea)}, pvalFRLogConstrictedVsFRLogDilated12p5{iCond}{iAreasOI(iArea)},...
          modelFRLogConstrictedVsFRLogDilated12p5{iCond}{iAreasOI(iArea)}] = corrPlot(getLog(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(inds,1)), getLog(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(inds,end)), options2);
      end
      
      % 33.33 cent
      if ~isempty(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)})
        
        % Signed modulations: Firing rate
        frModulations = getLog(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,end)) - getLog(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,1));
        options2.figTitle = ['Unit log10(firing rate) vs log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogSignedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogSignedModulationThird{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogSignedModulationThird{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogSignedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogSignedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogSignedModulationThird{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogSignedModulationThird{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogSignedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit log10(firing rate) vs absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = '|Unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedModulationThird{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedModulationThird{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        frModulations2 = getLog(abs(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,end) - areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,1)));
        options2.figTitle = ['Unit log10(firing rate) vs absolute log10(firing rate modulation) in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedModulationThird_2 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Unit log_{10}(|firing rate change|)';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedModulationThird_2{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedModulationThird_2{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedModulationThird_2{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations2(inds), options2);
        
        % Unsigned modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogUnsignedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogUnsignedModulationThird{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogUnsignedModulationThird{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogUnsignedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Signed normalised modulations: Firing rate
        frModulations = (getLog(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,end)) - getLog(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(:,1)))./getLog(areaFRIndividual{iCond}{iAreasOI(iArea)});
        options2.figTitle = ['Unit log10(firing rate) vs normalised log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogSignedNormalisedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Normalised unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogSignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogSignedNormalisedModulationThird{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogSignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        % Signed normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogSignedNormalisedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'r coef';
        options2.yLabel = 'Normalised unit log_{10}(firing rate) change';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogSignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogSignedNormalisedModulationThird{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogSignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), frModulations(inds), options2);
        
        % Unsigned normalised modulations: Firing rate
        frModulations = abs(frModulations);
        options2.figTitle = ['Unit log10(firing rate) vs normalised absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedNormalisedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = '|Normalised unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedNormalisedModulationThird{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations(inds), options2);
        
        frModulations2 = frModulations2./getLog(areaFRIndividual{iCond}{iAreasOI(iArea)});
        options2.figTitle = ['Unit log10(firing rate) vs normalised absolute log10(firing rate modulation) in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogVsfrLogUnsignedNormalisedModulationThird_2 in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Unit log_{10}(firing rate)';
        options2.yLabel = 'Normalised unit log_{10}(|firing rate change|)';
        options2.axesType = 'regular';
        [~, rCoefFRLogVsFRLogUnsignedNormalisedModulationThird_2{iCond}{iAreasOI(iArea)}, pvalFRLogVsFRLogUnsignedNormalisedModulationThird_2{iCond}{iAreasOI(iArea)},...
          modelFRLogVsFRLogUnsignedNormalisedModulationThird_2{iCond}{iAreasOI(iArea)}] = corrPlot(fr(inds), frModulations2(inds), options2);
        
        % Unsigned normalised modulations: r
        options2.figTitle = ['Unit rCoef over 5sec bins vs normalised absolute log10(firing rate) modulation in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['rSpearman5secVsfrLogUnsignedNormalisedModulationThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = '|r coef|';
        options2.yLabel = '|Normalised unit log_{10}(firing rate) change|';
        options2.axesType = 'regular';
        [~, rCoefrSpearman5secVsFRLogUnsignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLogUnsignedNormalisedModulationThird{iCond}{iAreasOI(iArea)},...
          modelrSpearman5secVsFRLogUnsignedNormalisedModulationThird{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), frModulations(inds), options2);
        
        % Firing rate: constricted vs dilated pupil
        options2.figTitle = ['Unit log10(firing rate) associated to constricted pupil phase vs unit log10(firing rate) associated to dilated pupil phase in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.figName = ['frLogConstrictedVsfrLogDilatedThird in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
        options2.xLabel = 'Constricted log_{10}(firing rate)';
        options2.yLabel = 'Dilated log_{10}(firing rate)';
        options2.axesType = 'regular';
        [~, rCoefFRLogConstrictedVsFRLogDilatedThird{iCond}{iAreasOI(iArea)}, pvalFRLogConstrictedVsFRLogDilatedThird{iCond}{iAreasOI(iArea)},...
          modelFRLogConstrictedVsFRLogDilatedThird{iCond}{iAreasOI(iArea)}] = corrPlot(getLog(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(inds,1)), getLog(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(inds,end)), options2);
      end
      
      % Firing rate: r
      options2.figTitle = ['Unit rCoef over 5sec bins vs log10(firing rate) in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
      options2.figName = ['rSpearman5secVsfrLog in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
      options2.xLabel = 'r coef';
      options2.yLabel = 'Unit log_{10}(firing rate)';
      options2.axesType = 'regular';
      [~, rCoefrSpearman5secVsFRLog{iCond}{iAreasOI(iArea)}, pvalrSpearman5secVsFRLog{iCond}{iAreasOI(iArea)},...
        modelrSpearman5secVsFRLog{iCond}{iAreasOI(iArea)}] = corrPlot(rSpearman(inds), fr(inds), options2);
      
      % Firing rate: unsigned r
      options2.figTitle = ['Unit rCoef over 5sec bins vs log10(firing rate) in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
      options2.figName = ['rSpearmanUnsigned5secVsfrLog in ' areas{iAreasOI(iArea)} ' during ' conditions{iCond}];
      options2.xLabel = '|r coef|';
      options2.yLabel = 'Unit log_{10}(firing rate)';
      options2.axesType = 'regular';
      [~, rCoefrSpearmanUnsigned5secVsFRLog{iCond}{iAreasOI(iArea)}, pvalrSpearmanUnsigned5secVsFRLog{iCond}{iAreasOI(iArea)},...
        modelrSpearmanUnsigned5secVsFRLog{iCond}{iAreasOI(iArea)}] = corrPlot(abs(rSpearman(inds)), fr(inds), options2);
    end
  end
  
  % Save correlation data
  save(filename, 'rCoefFRLogVsFRLogSignedModulation50','pvalFRLogVsFRLogSignedModulation50','rCoefrSpearman5secVsFRLogSignedModulation50','pvalrSpearman5secVsFRLogSignedModulation50',...
    'rCoefFRLogVsFRLogUnsignedModulation50','pvalFRLogVsFRLogUnsignedModulation50','rCoefrSpearman5secVsFRLogUnsignedModulation50','pvalrSpearman5secVsFRLogUnsignedModulation50',...
    'rCoefFRLogVsFRLogSignedNormalisedModulation50','pvalFRLogVsFRLogSignedNormalisedModulation50','rCoefrSpearman5secVsFRLogSignedNormalisedModulation50','pvalrSpearman5secVsFRLogSignedNormalisedModulation50',...
    'rCoefFRLogVsFRLogUnsignedNormalisedModulation50','pvalFRLogVsFRLogUnsignedNormalisedModulation50','rCoefrSpearman5secVsFRLogUnsignedNormalisedModulation50','pvalrSpearman5secVsFRLogUnsignedNormalisedModulation50',...
    'rCoefFRLogVsFRLogSignedModulation25','pvalFRLogVsFRLogSignedModulation25','rCoefrSpearman5secVsFRLogSignedModulation25','pvalrSpearman5secVsFRLogSignedModulation25',...
    'rCoefFRLogVsFRLogUnsignedModulation25','pvalFRLogVsFRLogUnsignedModulation25','rCoefrSpearman5secVsFRLogUnsignedModulation25','pvalrSpearman5secVsFRLogUnsignedModulation25',...
    'rCoefFRLogVsFRLogSignedNormalisedModulation25','pvalFRLogVsFRLogSignedNormalisedModulation25','rCoefrSpearman5secVsFRLogSignedNormalisedModulation25','pvalrSpearman5secVsFRLogSignedNormalisedModulation25',...
    'rCoefFRLogVsFRLogUnsignedNormalisedModulation25','pvalFRLogVsFRLogUnsignedNormalisedModulation25','rCoefrSpearman5secVsFRLogUnsignedNormalisedModulation25','pvalrSpearman5secVsFRLogUnsignedNormalisedModulation25',...
    'rCoefFRLogVsFRLogSignedModulation12p5','pvalFRLogVsFRLogSignedModulation12p5','rCoefrSpearman5secVsFRLogSignedModulation12p5','pvalrSpearman5secVsFRLogSignedModulation12p5',...
    'rCoefFRLogVsFRLogUnsignedModulation12p5','pvalFRLogVsFRLogUnsignedModulation12p5','rCoefrSpearman5secVsFRLogUnsignedModulation12p5','pvalrSpearman5secVsFRLogUnsignedModulation12p5',...
    'rCoefFRLogVsFRLogSignedNormalisedModulation12p5','pvalFRLogVsFRLogSignedNormalisedModulation12p5','rCoefrSpearman5secVsFRLogSignedNormalisedModulation12p5','pvalrSpearman5secVsFRLogSignedNormalisedModulation12p5',...
    'rCoefFRLogVsFRLogUnsignedNormalisedModulation12p5','pvalFRLogVsFRLogUnsignedNormalisedModulation12p5','rCoefrSpearman5secVsFRLogUnsignedNormalisedModulation12p5','pvalrSpearman5secVsFRLogUnsignedNormalisedModulation12p5',...
    'rCoefFRLogVsFRLogSignedModulationThird','pvalFRLogVsFRLogSignedModulationThird','rCoefrSpearman5secVsFRLogSignedModulationThird','pvalrSpearman5secVsFRLogSignedModulationThird',...
    'rCoefFRLogVsFRLogUnsignedModulationThird','pvalFRLogVsFRLogUnsignedModulationThird','rCoefrSpearman5secVsFRLogUnsignedModulationThird','pvalrSpearman5secVsFRLogUnsignedModulationThird',...
    'rCoefFRLogVsFRLogSignedNormalisedModulationThird','pvalFRLogVsFRLogSignedNormalisedModulationThird','rCoefrSpearman5secVsFRLogSignedNormalisedModulationThird','pvalrSpearman5secVsFRLogSignedNormalisedModulationThird',...
    'rCoefFRLogVsFRLogUnsignedNormalisedModulationThird','pvalFRLogVsFRLogUnsignedNormalisedModulationThird','rCoefrSpearman5secVsFRLogUnsignedNormalisedModulationThird','pvalrSpearman5secVsFRLogUnsignedNormalisedModulationThird',...
    'rCoefrSpearman5secVsFRLog','pvalrSpearman5secVsFRLog','rCoefrSpearmanUnsigned5secVsFRLog','pvalrSpearmanUnsigned5secVsFRLog',...
    'modelFRLogVsFRLogSignedModulation50','modelrSpearman5secVsFRLogSignedModulation50','modelFRLogVsFRLogUnsignedModulation50','modelrSpearman5secVsFRLogUnsignedModulation50','modelFRLogVsFRLogSignedNormalisedModulation50',...
    'modelrSpearman5secVsFRLogSignedNormalisedModulation50','modelFRLogVsFRLogUnsignedNormalisedModulation50','modelrSpearman5secVsFRLogUnsignedNormalisedModulation50','modelFRLogVsFRLogSignedModulation25',...
    'modelrSpearman5secVsFRLogSignedModulation25','modelFRLogVsFRLogUnsignedModulation25','modelrSpearman5secVsFRLogUnsignedModulation25','modelFRLogVsFRLogSignedNormalisedModulation25',...
    'modelrSpearman5secVsFRLogSignedNormalisedModulation25','modelrSpearman5secVsFRLogSignedNormalisedModulation25','modelFRLogVsFRLogUnsignedNormalisedModulation25','modelrSpearman5secVsFRLogUnsignedNormalisedModulation25',...
    'modelFRLogVsFRLogSignedModulation12p5','modelrSpearman5secVsFRLogSignedModulation12p5','modelFRLogVsFRLogUnsignedModulation12p5','modelrSpearman5secVsFRLogUnsignedModulation12p5','modelFRLogVsFRLogSignedNormalisedModulation12p5',...
    'modelrSpearman5secVsFRLogSignedNormalisedModulation12p5','modelFRLogVsFRLogUnsignedNormalisedModulation12p5','modelrSpearman5secVsFRLogUnsignedNormalisedModulation12p5','modelFRLogVsFRLogSignedModulationThird',...
    'modelrSpearman5secVsFRLogSignedModulationThird','modelFRLogVsFRLogUnsignedModulationThird','modelrSpearman5secVsFRLogUnsignedModulationThird','modelFRLogVsFRLogSignedNormalisedModulationThird',...
    'modelrSpearman5secVsFRLogSignedNormalisedModulationThird','modelFRLogVsFRLogUnsignedNormalisedModulationThird','modelrSpearman5secVsFRLogUnsignedNormalisedModulationThird',...
    'modelrSpearman5secVsFRLog','modelrSpearmanUnsigned5secVsFRLog',...
    'rCoefFRLogConstrictedVsFRLogDilated50','pvalFRLogConstrictedVsFRLogDilated50','modelFRLogConstrictedVsFRLogDilated50','rCoefFRLogConstrictedVsFRLogDilated25','pvalFRLogConstrictedVsFRLogDilated25',...
    'modelFRLogConstrictedVsFRLogDilated25','rCoefFRLogConstrictedVsFRLogDilated12p5','pvalFRLogConstrictedVsFRLogDilated12p5','modelFRLogConstrictedVsFRLogDilated12p5','rCoefFRLogConstrictedVsFRLogDilatedThird',...
    'pvalFRLogConstrictedVsFRLogDilatedThird','modelFRLogConstrictedVsFRLogDilatedThird',...
    'rCoefFRLogVsFRLogUnsignedModulation50_2','pvalFRLogVsFRLogUnsignedModulation50_2','modelFRLogVsFRLogUnsignedModulation50_2',...
    'rCoefFRLogVsFRLogUnsignedModulation25_2','pvalFRLogVsFRLogUnsignedModulation25_2','modelFRLogVsFRLogUnsignedModulation25_2',...
    'rCoefFRLogVsFRLogUnsignedModulation12p5_2','pvalFRLogVsFRLogUnsignedModulation12p5_2','modelFRLogVsFRLogUnsignedModulation12p5_2',...
    'rCoefFRLogVsFRLogUnsignedModulationThird_2','pvalFRLogVsFRLogUnsignedModulationThird_2','modelFRLogVsFRLogUnsignedModulationThird_2', '-append');
end


%% FIRING RATE MEAN AND VARIANCE MODULATION CORRELATIONS: LOG
options = struct();
options.corrType = 'Spearman';
options.figFolder = mainFolder;
options.diagonal = true;
options.fitLine = true;
options.fitLineDisplay = 'on';
options.xLim = [];
options.xTicks = [];
options.yLim = [];
options.yTicks = [];
options.figSize = figSize;
options.saveFig = true;
if drawCorrelationsLog2 && strcmp(subpop, 'all')
  for iCond = 1:min([2 numel(conditions)])
    fractions = [];
    colourCount = 0;
    colourVector = [];
    colourCodes = [];
    meanModulations50 = [];
    varModulations50 = [];
    meanModulations25 = [];
    varModulations25 = [];
    meanModulations12p5 = [];
    varModulations12p5 = [];
    meanModulationsThird = [];
    varModulationsThird = [];
    for iArea = 1:numel(iAreasOI)
      recIDSpearman = areaRecIDSpearman5secIndividual{iCond}{iAreasOI(iArea)};
      if ~isempty(recIDSpearman)
        %[~, ~, areaName] = determineArea(seriesFromEntry(recIDSpearman{1}));
        areaName = areas{iAreasOI(iArea)};
        if strcmpi(areaName, 'Th') || strcmpi(areaName, 'VB') || strcmpi(areaName, 'Po') ||...
            strcmpi(areaName, 'LGN') || strcmpi(areaName, 'LP') || strcmpi(areaName, 'S1') ||...
            strcmpi(areaName, 'RSC') || strcmpi(areaName, 'V1') || strcmpi(areaName, 'V2') ||...
            strcmpi(areaName, 'VIS') || strcmpi(areaName, 'CA') || strcmpi(areaName, 'DG') || strcmpi(areaName, 'Cx')
          colourCount = colourCount + 1;
          rSpearman = arearSpearman5secIndividual{iCond}{iAreasOI(iArea)};
          pvalSpearman = areapvalSpearman5secIndividual{iCond}{iAreasOI(iArea)};
          uniqueRecs = unique(recIDSpearman);
          for iRec = 1:numel(uniqueRecs)
            unitInds = ismember(recIDSpearman, uniqueRecs{iRec});
            nUnits = sum(unitInds);
            if nUnits >= 10
              fractions = [fractions; sum(unitInds & rSpearman >= 0)/sum(unitInds)];
              colourVector = [colourVector; colourCount];
              colourCodes = [colourCodes; areaColours2(areaName)];
              
              % 50 cent
              dilatedFR = getLog(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(unitInds,end));
              constrictedFR = getLog(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}(unitInds,1));
              meanModulations50 = [meanModulations50; mean(dilatedFR - constrictedFR)];
              varModulations50 = [varModulations50; var(dilatedFR) - var(constrictedFR)];
              
              % 25 cent
              dilatedFR = getLog(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(unitInds,end));
              constrictedFR = getLog(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}(unitInds,1));
              meanModulations25 = [meanModulations25; mean(dilatedFR - constrictedFR)];
              varModulations25 = [varModulations25; var(dilatedFR) - var(constrictedFR)];
              
              % 12.5 cent
              dilatedFR = getLog(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(unitInds,end));
              constrictedFR = getLog(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}(unitInds,1));
              meanModulations12p5 = [meanModulations12p5; mean(dilatedFR - constrictedFR)];
              varModulations12p5 = [varModulations12p5; var(dilatedFR) - var(constrictedFR)];
              
              % 33.33 cent
              dilatedFR = getLog(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(unitInds,end));
              constrictedFR = getLog(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}(unitInds,1));
              meanModulationsThird = [meanModulationsThird; mean(dilatedFR - constrictedFR)];
              varModulationsThird = [varModulationsThird; var(dilatedFR) - var(constrictedFR)];
            end
          end
        end
      end
    end
    options.colourVector = colourVector;
    options.colourCodes = unique(colourCodes,'rows');
    options2 = options;
    
    % 50 cent
    options2.figTitle = ['Fraction of positive units vs log10(firing rate) mean modulation during ' conditions{iCond}];
    options2.figName = ['positiveFractionVsfrLogSignedMeanModulation50 during ' conditions{iCond}];
    options2.xLabel = 'Fraction of positive units';
    options2.yLabel = 'Unit log_{10}(firing rate) change';
    options2.axesType = 'regular';
    [~, rCoefPositiveFractionVsFRLogSignedMeanModulation50{iCond}, pvalPositiveFractionVsFRLogSignedMeanModulation50{iCond},...
      modelPositiveFractionVsFRLogSignedMeanModulation50{iCond}] = corrPlot(fractions, meanModulations50, options2);
    
    options2.figTitle = ['Fraction of positive units vs log10(firing rate) variance modulation during ' conditions{iCond}];
    options2.figName = ['positiveFractionVsfrLogSignedVarModulation50 during ' conditions{iCond}];
    options2.xLabel = 'Fraction of positive units';
    options2.yLabel = 'Unit log_{10}(firing rate) variance change';
    options2.axesType = 'regular';
    [~, rCoefPositiveFractionVsFRLogSignedVarModulation50{iCond}, pvalPositiveFractionVsFRLogSignedVarModulation50{iCond},...
      modelPositiveFractionVsFRLogSignedVarModulation50{iCond}] = corrPlot(fractions, varModulations50, options2);
    
    % 25 cent
    options2.figTitle = ['Fraction of positive units vs log10(firing rate) mean modulation during ' conditions{iCond}];
    options2.figName = ['positiveFractionVsfrLogSignedMeanModulation25 during ' conditions{iCond}];
    options2.xLabel = 'Fraction of positive units';
    options2.yLabel = 'Unit log_{10}(firing rate) change';
    options2.axesType = 'regular';
    [~, rCoefPositiveFractionVsFRLogSignedMeanModulation25{iCond}, pvalPositiveFractionVsFRLogSignedMeanModulation25{iCond},...
      modelPositiveFractionVsFRLogSignedMeanModulation25{iCond}] = corrPlot(fractions, meanModulations25, options2);
    
    options2.figTitle = ['Fraction of positive units vs log10(firing rate) variance modulation during ' conditions{iCond}];
    options2.figName = ['positiveFractionVsfrLogSignedVarModulation25 during ' conditions{iCond}];
    options2.xLabel = 'Fraction of positive units';
    options2.yLabel = 'Unit log_{10}(firing rate) variance change';
    options2.axesType = 'regular';
    [~, rCoefPositiveFractionVsFRLogSignedVarModulation25{iCond}, pvalPositiveFractionVsFRLogSignedVarModulation25{iCond},...
      modelPositiveFractionVsFRLogSignedVarModulation25{iCond}] = corrPlot(fractions, varModulations25, options2);
    
    % 12.5 cent
    options2.figTitle = ['Fraction of positive units vs log10(firing rate) mean modulation during ' conditions{iCond}];
    options2.figName = ['positiveFractionVsfrLogSignedMeanModulation12p5 during ' conditions{iCond}];
    options2.xLabel = 'Fraction of positive units';
    options2.yLabel = 'Unit log_{10}(firing rate) change';
    options2.axesType = 'regular';
    [~, rCoefPositiveFractionVsFRLogSignedMeanModulation12p5{iCond}, pvalPositiveFractionVsFRLogSignedMeanModulation12p5{iCond},...
      modelPositiveFractionVsFRLogSignedMeanModulation12p5{iCond}] = corrPlot(fractions, meanModulations12p5, options2);
    
    options2.figTitle = ['Fraction of positive units vs log10(firing rate) variance modulation during ' conditions{iCond}];
    options2.figName = ['positiveFractionVsfrLogSignedVarModulation12p5 during ' conditions{iCond}];
    options2.xLabel = 'Fraction of positive units';
    options2.yLabel = 'Unit log_{10}(firing rate) variance change';
    options2.axesType = 'regular';
    [~, rCoefPositiveFractionVsFRLogSignedVarModulation12p5{iCond}, pvalPositiveFractionVsFRLogSignedVarModulation12p5{iCond},...
      modelPositiveFractionVsFRLogSignedVarModulation12p5{iCond}] = corrPlot(fractions, varModulations12p5, options2);
    
    % 33.33 cent
    options2.figTitle = ['Fraction of positive units vs log10(firing rate) mean modulation during ' conditions{iCond}];
    options2.figName = ['positiveFractionVsfrLogSignedMeanModulationThird during ' conditions{iCond}];
    options2.xLabel = 'Fraction of positive units';
    options2.yLabel = 'Unit log_{10}(firing rate) change';
    options2.axesType = 'regular';
    [~, rCoefPositiveFractionVsFRLogSignedMeanModulationThird{iCond}, pvalPositiveFractionVsFRLogSignedMeanModulationThird{iCond},...
      modelPositiveFractionVsFRLogSignedMeanModulationThird{iCond}] = corrPlot(fractions, meanModulationsThird, options2);
    
    options2.figTitle = ['Fraction of positive units vs log10(firing rate) variance modulation during ' conditions{iCond}];
    options2.figName = ['positiveFractionVsfrLogSignedVarModulationThird during ' conditions{iCond}];
    options2.xLabel = 'Fraction of positive units';
    options2.yLabel = 'Unit log_{10}(firing rate) variance change';
    options2.axesType = 'regular';
    [~, rCoefPositiveFractionVsFRLogSignedVarModulationThird{iCond}, pvalPositiveFractionVsFRLogSignedVarModulationThird{iCond},...
      modelPositiveFractionVsFRLogSignedVarModulationThird{iCond}] = corrPlot(fractions, varModulationsThird, options2);
  end
  
  % Save correlation data
  save(filename, 'rCoefPositiveFractionVsFRLogSignedMeanModulation50','pvalPositiveFractionVsFRLogSignedMeanModulation50','modelPositiveFractionVsFRLogSignedMeanModulation50',...
    'rCoefPositiveFractionVsFRLogSignedMeanModulation25','pvalPositiveFractionVsFRLogSignedMeanModulation25','modelPositiveFractionVsFRLogSignedMeanModulation25',...
    'rCoefPositiveFractionVsFRLogSignedMeanModulation12p5','pvalPositiveFractionVsFRLogSignedMeanModulation12p5','modelPositiveFractionVsFRLogSignedMeanModulation12p5',...
    'rCoefPositiveFractionVsFRLogSignedMeanModulationThird','pvalPositiveFractionVsFRLogSignedMeanModulationThird','modelPositiveFractionVsFRLogSignedMeanModulationThird', '-append');
end


%% PREDICTING FIRING RATE MODULATION EFFECTS ASSOCIATED WITH THE PUPIL AREA SIZE
edges = -3:0.2:3;
centres = (edges(1:end-1)+0.1);
x = centres(1):0.001:centres(end);
options = struct();
options.xLabel = 'Log_{10}(firing rate)';
options.xLim = [edges(1) 2];
if strcmpi(subpop, 'all')
  options.yLabel = 'Probability';
elseif strcmpi(subpop, 'positive')
  %options.yLabel = 'Positive unit count';
  options.yLabel = 'Probability';
elseif strcmpi(subpop, 'negative')
  %options.yLabel = 'Negative unit count';
  options.yLabel = 'Probability';
end
options.distroType = distroType;
options.figFolder = mainFolder;
options.figSize = 18;
options.saveFig = true;
options.pdf = true;
if drawPredictions
  if predictionType == 1
    for iCond = 1:min([2 numel(conditions)])
      for iArea = 1:numel(iAreasOI)
        if ~isempty(areaFRIndividual{iCond}{iAreasOI(iArea)})
          frData = areaFRIndividual{iCond}{iAreasOI(iArea)};
          
          % 50 cent
          options.figName = ['FiringRateCentile50Predicted_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          options.legendLabels = {'0:50 cent predicted','50:100 cent predicted'};
          model = modelFRLogVsFRLogSignedModulation50{iCond}{iAreasOI(iArea)};
          [~, predictions50{iCond}{iAreasOI(iArea)}] = firingRatePredictPlot(edges, centres, x, frData, model, options);
          
          % 25 cent
          options.figName = ['FiringRateCentile25Predicted_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          options.legendLabels = {'0:25 cent predicted','75:100 cent predicted'};
          model = modelFRLogVsFRLogSignedModulation25{iCond}{iAreasOI(iArea)};
          [~, predictions25{iCond}{iAreasOI(iArea)}] = firingRatePredictPlot(edges, centres, x, frData, model, options);
          
          % 12.5 cent
          options.figName = ['FiringRateCentile12p5Predicted_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          options.legendLabels = {'0:12.5 cent predicted','87.5:100 cent predicted'};
          model = modelFRLogVsFRLogSignedModulation12p5{iCond}{iAreasOI(iArea)};
          [~, predictions12p5{iCond}{iAreasOI(iArea)}] = firingRatePredictPlot(edges, centres, x, frData, model, options);
          
          % 33.33 cent
          options.figName = ['FiringRateCentileThirdPredicted_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          options.legendLabels = {'0:33.3 cent predicted','66.7:100 cent predicted'};
          model = modelFRLogVsFRLogSignedModulationThird{iCond}{iAreasOI(iArea)};
          [~, predictionsThird{iCond}{iAreasOI(iArea)}] = firingRatePredictPlot(edges, centres, x, frData, model, options);
        end
      end
    end
    save(filename, 'predictions50','predictions25','predictions12p5','predictionsThird', '-append');
  elseif predictionType == 2
    for iCond = 1:min([2 numel(conditions)])
      for iArea = 1:numel(iAreasOI)
        if ~isempty(areaFRIndividual{iCond}{iAreasOI(iArea)})
          frData = areaFRIndividual{iCond}{iAreasOI(iArea)};
          
          % 50 cent
          options.figName = ['FiringRateCentile50Predicted2_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          options.legendLabels = {'0:50 cent predicted','50:100 cent predicted'};
          model = modelFRLogVsFRLogSignedModulation50{iCond}{iAreasOI(iArea)};
          [~, predictions50_2{iCond}{iAreasOI(iArea)}] = firingRatePredictPlot2(edges, centres, x, frData, model, options);
          
          % 25 cent
          options.figName = ['FiringRateCentile25Predicted2_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          options.legendLabels = {'0:25 cent predicted','75:100 cent predicted'};
          model = modelFRLogVsFRLogSignedModulation25{iCond}{iAreasOI(iArea)};
          [~, predictions25_2{iCond}{iAreasOI(iArea)}] = firingRatePredictPlot2(edges, centres, x, frData, model, options);
          
          % 12.5 cent
          options.figName = ['FiringRateCentile12p5Predicted2_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          options.legendLabels = {'0:12.5 cent predicted','87.5:100 cent predicted'};
          model = modelFRLogVsFRLogSignedModulation12p5{iCond}{iAreasOI(iArea)};
          [~, predictions12p5_2{iCond}{iAreasOI(iArea)}] = firingRatePredictPlot2(edges, centres, x, frData, model, options);
          
          % 33.33 cent
          options.figName = ['FiringRateCentileThirdPredicted2_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          options.legendLabels = {'0:33.3 cent predicted','66.7:100 cent predicted'};
          model = modelFRLogVsFRLogSignedModulationThird{iCond}{iAreasOI(iArea)};
          [~, predictionsThird_2{iCond}{iAreasOI(iArea)}] = firingRatePredictPlot2(edges, centres, x, frData, model, options);
        end
      end
    end
    save(filename, 'predictions50_2','predictions25_2','predictions12p5_2','predictionsThird_2', '-append');
  end
end


%% PLOT THE REAL AND PREDICTED LOG FIRING RATE DISTRIBUTIONS
if drawCombinedDistros && strcmpi(subpop, 'all')
  edges = -3:0.2:3;
  centres = (edges(1:end-1)+0.1);
  x = centres(1):0.001:centres(end);
  options = struct();
  options.xLabel = 'Log_{10}(firing rate)';
  options.xLim = [edges(1)+0.5 2];
  if strcmpi(subpop, 'all')
    options.yLabel = 'Probability';
  elseif strcmpi(subpop, 'positive')
    %options.yLabel = 'Positive unit count';
    options.yLabel = 'Probability';
  elseif strcmpi(subpop, 'negative')
    %options.yLabel = 'Negative unit count';
    options.yLabel = 'Probability';
  end
  options.distroType = distroType;
  options.figFolder = mainFolder;
  options.figSize = 18;
  options.saveFig = false;
  options.pdf = true;
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if ~isempty(areaFRIndividual{iCond}{iAreasOI(iArea)})
        if predictionType == 1
          frData = areaFRIndividual{iCond}{iAreasOI(iArea)};
          
          % 50 cent
          options.figName = ['FiringRateCentile50Combined_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          legendLabels1 = {'0:50% data','50:100% data'};
          legendLabels2 = {'0:50% prediction','50:100% prediction'};
          firingRateCombinedPlot(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}, edges, centres, x, frData,...
            modelFRLogVsFRLogSignedModulation50{iCond}{iAreasOI(iArea)}, legendLabels1, legendLabels2, options)
          
          % 25 cent
          options.figName = ['FiringRateCentile25Combined_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          legendLabels1 = {'0:25% data','75:100% data'};
          legendLabels2 = {'0:25% prediction','75:100% prediction'};
          firingRateCombinedPlot(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}, edges, centres, x, frData,...
            modelFRLogVsFRLogSignedModulation25{iCond}{iAreasOI(iArea)}, legendLabels1, legendLabels2, options)
          
          % 12.5 cent
          options.figName = ['FiringRateCentile12p5Combined_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          legendLabels1 = {'0:12.5% data','87.5:100% data'};
          legendLabels2 = {'0:12.5% prediction','87.5:100% prediction'};
          firingRateCombinedPlot(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}, edges, centres, x, frData,...
            modelFRLogVsFRLogSignedModulation12p5{iCond}{iAreasOI(iArea)}, legendLabels1, legendLabels2, options)
          
          % 33.33 cent
          options.figName = ['FiringRateCentileThirdCombined_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          legendLabels1 = {'0:33.3% data','66.7:100% data'};
          legendLabels2 = {'0:33.3% prediction','66.7:100% prediction'};
          firingRateCombinedPlot(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}, edges, centres, x, frData,...
            modelFRLogVsFRLogSignedModulationThird{iCond}{iAreasOI(iArea)}, legendLabels1, legendLabels2, options)
        elseif predictionType == 2
          frData = areaFRIndividual{iCond}{iAreasOI(iArea)};
          
          % 50 cent
          options.figName = ['FiringRateCentile50Combined2_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          legendLabels1 = {'0:50% data','50:100% data'};
          legendLabels2 = {'0:50% prediction','50:100% prediction'};
          firingRateCombinedPlot2(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}, edges, centres, x, frData,...
            modelFRLogVsFRLogSignedModulation50{iCond}{iAreasOI(iArea)}, legendLabels1, legendLabels2, options)
          
          % 25 cent
          options.figName = ['FiringRateCentile25Combined2_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          legendLabels1 = {'0:25% data','75:100% data'};
          legendLabels2 = {'0:25% prediction','75:100% prediction'};
          firingRateCombinedPlot2(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}, edges, centres, x, frData,...
            modelFRLogVsFRLogSignedModulation25{iCond}{iAreasOI(iArea)}, legendLabels1, legendLabels2, options)
          
          % 12.5 cent
          options.figName = ['FiringRateCentile12p5Combined2_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          legendLabels1 = {'0:12.5% data','87.5:100% data'};
          legendLabels2 = {'0:12.5% prediction','87.5:100% prediction'};
          firingRateCombinedPlot2(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}, edges, centres, x, frData,...
            modelFRLogVsFRLogSignedModulation12p5{iCond}{iAreasOI(iArea)}, legendLabels1, legendLabels2, options)
          
          % 33.33 cent
          options.figName = ['FiringRateCentileThirdCombined2_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
          legendLabels1 = {'0:33.3% data','66.7:100% data'};
          legendLabels2 = {'0:33.3% prediction','66.7:100% prediction'};
          firingRateCombinedPlot2(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}, edges, centres, x, frData,...
            modelFRLogVsFRLogSignedModulationThird{iCond}{iAreasOI(iArea)}, legendLabels1, legendLabels2, options)
        end
        
      end
    end
  end
end


%% PLOT THE REAL AND PREDICTED LOG FIRING RATE DISTRIBUTIONS WITH POPULATIONS SPLIT
if drawCombinedDistrosSplit && strcmpi(subpop, 'all')
  positiveData = load(filenamePositive);
  negativeData = load(filenameNegative);
  edges = -3:0.2:3;
  centres = (edges(1:end-1)+0.1);
  x = centres(1):0.001:centres(end);
  options = struct();
  options.xLabel = 'Log_{10}(firing rate)';
  options.xLim = [edges(1)+0.5 2];
  options.yLabel = 'Probability';
  options.figFolder = mainFolder;
  options.figSize = 18;
  options.saveFig = false;
  options.pdf = true;
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if ~isempty(areaFRIndividual{iCond}{iAreasOI(iArea)})
        
        % 50 cent
        options.figName = ['FiringRateCentile50CombinedSplit_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        legendLabels1 = {'0:50% data','50:100% data'};
        legendLabels2 = {'0:50% prediction','50:100% prediction'};
        predictedDistros = positiveData.predictions50{iCond}{iAreasOI(iArea)} + negativeData.predictions50{iCond}{iAreasOI(iArea)};
        firingRateCombinedSplitPlot(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}, edges, predictedDistros,...
          legendLabels1, legendLabels2, options)
        
        % 25 cent
        options.figName = ['FiringRateCentile25CombinedSplit_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        legendLabels1 = {'0:25% data','75:100% data'};
        legendLabels2 = {'0:25% prediction','75:100% prediction'};
        predictedDistros = positiveData.predictions25{iCond}{iAreasOI(iArea)} + negativeData.predictions25{iCond}{iAreasOI(iArea)};
        firingRateCombinedSplitPlot(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}, edges, predictedDistros,...
          legendLabels1, legendLabels2, options)
        
        % 12.5 cent
        options.figName = ['FiringRateCentile12p5CombinedSplit_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        legendLabels1 = {'0:12.5% data','87.5:100% data'};
        legendLabels2 = {'0:12.5% prediction','87.5:100% prediction'};
        predictedDistros = positiveData.predictions12p5{iCond}{iAreasOI(iArea)} + negativeData.predictions12p5{iCond}{iAreasOI(iArea)};
        firingRateCombinedSplitPlot(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}, edges, predictedDistros,...
          legendLabels1, legendLabels2, options)
        
        % 33.33 cent
        options.figName = ['FiringRateCentileThirdCombinedSplit_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        legendLabels1 = {'0:33.3% data','66.7:100% data'};
        legendLabels2 = {'0:33.3% prediction','66.7:100% prediction'};
        predictedDistros = positiveData.predictionsThird{iCond}{iAreasOI(iArea)} + negativeData.predictionsThird{iCond}{iAreasOI(iArea)};
        firingRateCombinedSplitPlot(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}, edges, predictedDistros,...
          legendLabels1, legendLabels2, options)
        
      end
    end
  end
end


%% PLOT THE GAMMA-FITTED FIRING RATE DISTRIBUTIONS
if drawSmoothCombinedDistros && strcmpi(subpop, 'all')
  binSize = 0.02;
  edges = -3:binSize:3;
  centres = (edges(1:end-1)+(binSize/2));
  x = centres(1):0.001:centres(end);
  options = struct();
  options.xLabel = 'Log_{10}(firing rate)';
  options.xLim = [edges(1)+0.5 2];
  if strcmpi(subpop, 'all')
    options.yLabel = 'Probability';
  elseif strcmpi(subpop, 'positive')
    %options.yLabel = 'Positive unit count';
    options.yLabel = 'Probability';
  elseif strcmpi(subpop, 'negative')
    %options.yLabel = 'Negative unit count';
    options.yLabel = 'Probability';
  end
  options.distroType = distroType;
  options.figFolder = mainFolder;
  options.figSize = 18;
  options.saveFig = false;
  options.pdf = true;
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if ~isempty(areaFRIndividual{iCond}{iAreasOI(iArea)})
        % 50 cent
        options.figName = ['FiringRateCentile50Smoothed_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        legendLabels = {'0:50% pupil area','50:100% pupil area'};
        firingRateSmoothedPlot(areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}, edges, legendLabels, options)
        
        % 25 cent
        options.figName = ['FiringRateCentile25Smoothed_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        legendLabels = {'0:25% pupil area','75:100% pupil area'};
        firingRateSmoothedPlot(areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}, edges, legendLabels, options)
        
        % 12.5 cent
        options.figName = ['FiringRateCentile12p5Smoothed_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        legendLabels = {'0:12.5% pupil area','87.5:100% pupil area'};
        firingRateSmoothedPlot(areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}, edges, legendLabels, options)
        
        % 33.33 cent
        options.figName = ['FiringRateCentileThirdSmoothed_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        legendLabels = {'0:33.3% pupil area','66.7:100% pupil area'};
        firingRateSmoothedPlot(areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}, edges, legendLabels, options)
        
      end
    end
  end
end


%% PLOT THE SPLIT GAMMA-FITTED FIRING RATE DISTRIBUTIONS
if drawSmoothCombinedDistrosSplit && strcmpi(subpop, 'all')
  positiveData = load(filenamePositive);
  negativeData = load(filenameNegative);
  binSize = 0.02;
  edges = -3:binSize:3;
  centres = (edges(1:end-1)+(binSize/2));
  x = centres(1):0.001:centres(end);
  options = struct();
  options.xLabel = 'Log_{10}(firing rate)';
  options.xLim = [edges(1)+0.5 2];
  if strcmpi(subpop, 'all')
    options.yLabel = 'Probability';
  elseif strcmpi(subpop, 'positive')
    %options.yLabel = 'Positive unit count';
    options.yLabel = 'Probability';
  elseif strcmpi(subpop, 'negative')
    %options.yLabel = 'Negative unit count';
    options.yLabel = 'Probability';
  end
  options.distroType = distroType;
  options.figFolder = mainFolder;
  options.figSize = 18;
  options.saveFig = false;
  options.pdf = true;
  for iCond = 1:min([2 numel(conditions)])
    for iArea = 1:numel(iAreasOI)
      if ~isempty(areaFRIndividual{iCond}{iAreasOI(iArea)})
        % 50 cent
        options.figName = ['FiringRateCentile50SplitSmoothed_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        firingRateSplitSmoothedPlot(positiveData.areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)},...
          negativeData.areaFRPercentileIndividual_50{iCond}{iAreasOI(iArea)}, edges, options)
        
        % 25 cent
        options.figName = ['FiringRateCentile25SplitSmoothed_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        firingRateSplitSmoothedPlot(positiveData.areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)},...
          negativeData.areaFRPercentileIndividual_25{iCond}{iAreasOI(iArea)}, edges, options)
        
        % 12.5 cent
        options.figName = ['FiringRateCentile12p5SplitSmoothed_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        firingRateSplitSmoothedPlot(positiveData.areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)},...
          negativeData.areaFRPercentileIndividual_12p5{iCond}{iAreasOI(iArea)}, edges, options)
        
        % 33.33 cent
        options.figName = ['FiringRateCentileThirdSplitSmoothed_' distroType '_' areas{iAreasOI(iArea)} '_' conditions{iCond}];
        firingRateSplitSmoothedPlot(positiveData.areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)},...
          negativeData.areaFRPercentileIndividual_third{iCond}{iAreasOI(iArea)}, edges, options)
        
      end
    end
  end
end



%% Local functions
function unitFiringRates = firingRateByPhase(spk, srSpk, phase, phaseStep, phaseShift)

phaseValues = (pi:-phaseStep:-pi) + phaseShift;
phase = recentrePhase(phase, phaseShift);
unitFiringRates = zeros(size(spk,1),numel(phaseValues)-1);
for iPhase = 1:numel(phaseValues)-1
  unitFiringRates(:,iPhase) = (sum(spk(:,phase<phaseValues(iPhase) & phase>=phaseValues(iPhase+1)),2).*srSpk)./size(spk(:,phase<phaseValues(iPhase) & phase>=phaseValues(iPhase+1)),2);
end
end


function [statsMean, statsVar, statsMeanSplit, statsVarSplit] = firingRateTestAndPlot(edges, data, rSpearman, options) %#ok<*DEFNU>

% Combined histograms
logRates = getLog(data);
statsVar = varTest(logRates);
options.stats = statsVar;
histPlotFR(edges, logRates, options);
statsMean = meanTest(logRates, 'ANOVARM');

% Split histograms
if ~isempty(rSpearman)
  if min(size(rSpearman)) == 1
    options.colours = [matlabColours(1); matlabColours(2); matlabColours(1); matlabColours(2)];
    options.lineStyles = {'--','--','-','-'};
    options.markerStyles = {'v','v','^','^'};
    options.legendLabels = {'Constricted negative','Constricted positive','Dilated negative','Dilated positive'};
    options.stats = [];
    options.saveFig = false;
    logRates2 = {logRates(rSpearman < 0, 1), logRates(rSpearman >= 0, 1), logRates(rSpearman < 0, end), logRates(rSpearman >= 0, end)};
    h = histPlotFR(edges, logRates2, options);
    l = legend;
    l.Position = l.Position - [0.08 0.08 0 0];
    
    % Display var test stats
    logRates3 = NaN(size(data,1),4);
    for col = 1:numel(logRates2)
      logRates3(1:numel(logRates2{col}),col) = logRates2{col};
    end
    statsVarSplit = varTest(logRates3);
    xLim = xlim;
    xAxisLength = xLim(2)-xLim(1);
    yLim = ylim;
    yAxisLength = yLim(2)-yLim(1);
    stats = struct2table(statsVarSplit, 'AsArray',true);
    iComp1 = stats.iCol1 == 1 & stats.iCol2 == 3;
    iComp2 = stats.iCol1 == 2 & stats.iCol2 == 4;
    iComp3 = stats.iCol1 == 1 & stats.iCol2 == 2;
    iComp4 = stats.iCol1 == 3 & stats.iCol2 == 4;
    textStr = ['CN(\sigma^2=' num2str(stats.var1(iComp1)) ')vDN(\sigma^2=' num2str(stats.var2(iComp1)) ') p=' num2str(stats.p(iComp1)) ', '...
      'CP(\sigma^2=' num2str(stats.var1(iComp2)) ')vDP(\sigma^2=' num2str(stats.var2(iComp2)) ') p=' num2str(stats.p(iComp2))];
    text(xLim(2)-xAxisLength*0.99, yLim(2)-yAxisLength*0.045, textStr, 'FontSize',10);
    textStr = ['CN(\sigma^2=' num2str(stats.var1(iComp3)) ')vCP(\sigma^2=' num2str(stats.var2(iComp3)) ') p=' num2str(stats.p(iComp3)) ', '...
      'DN(\sigma^2=' num2str(stats.var1(iComp4)) ')vDP(\sigma^2=' num2str(stats.var2(iComp4)) ') p=' num2str(stats.p(iComp4))];
    text(xLim(2)-xAxisLength*0.99, yLim(2)-yAxisLength*0.08, textStr, 'FontSize',10);
    
    % Display means test stats
%     [statsMeanSplit.pNegative, ~, statsMeanSplit.statsNegative] = signrank(logRates2{1}, logRates2{3});
%     [statsMeanSplit.pPositive, ~, statsMeanSplit.statsPositive] = signrank(logRates2{2}, logRates2{4});
%     [statsMeanSplit.pConstricted, ~, statsMeanSplit.statsConstricted] = ranksum(logRates2{1}, logRates2{2});
%     [statsMeanSplit.pDilated, ~, statsMeanSplit.statsDilated] = ranksum(logRates2{3}, logRates2{4});
    [~, statsMeanSplit.pNegative, ~, statsMeanSplit.statsNegative] = ttest(logRates2{1}, logRates2{3});
    [~, statsMeanSplit.pPositive, ~, statsMeanSplit.statsPositive] = ttest(logRates2{2}, logRates2{4});
    [~, statsMeanSplit.pConstricted, ~, statsMeanSplit.statsConstricted] = ttest2(logRates2{1}, logRates2{2});
    [~, statsMeanSplit.pDilated, ~, statsMeanSplit.statsDilated] = ttest2(logRates2{3}, logRates2{4});
    statsMeanSplit.meansNegative = [mean(logRates2{1}), mean(logRates2{3})];
    statsMeanSplit.meansPositive = [mean(logRates2{2}), mean(logRates2{4})];
    statsMeanSplit.meansConstricted = [mean(logRates2{1}), mean(logRates2{2})];
    statsMeanSplit.meansDilated = [mean(logRates2{3}), mean(logRates2{4})];
    [~, ci1] = datamean(logRates2{1}); [~, ci2] = datamean(logRates2{3}); statsMeanSplit.ciNegative = [ci1 ci2];
    [~, ci1] = datamean(logRates2{2}); [~, ci2] = datamean(logRates2{4}); statsMeanSplit.ciPositive = [ci1 ci2];
    [~, ci1] = datamean(logRates2{1}); [~, ci2] = datamean(logRates2{2}); statsMeanSplit.ciConstricted = [ci1 ci2];
    [~, ci1] = datamean(logRates2{3}); [~, ci2] = datamean(logRates2{4}); statsMeanSplit.ciDilated = [ci1 ci2];
    stats = statsMeanSplit;
    textStr = ['CN(\mu=' num2str(mean(logRates2{1},'omitnan')) ')vDN(\mu=' num2str(mean(logRates2{3},'omitnan')) ') p=' num2str(stats.pNegative) ', '...
      'CP(\mu=' num2str(mean(logRates2{2},'omitnan')) ')vDP(\mu=' num2str(mean(logRates2{4},'omitnan')) ') p=' num2str(stats.pPositive)];
    text(xLim(2)-xAxisLength*0.99, yLim(2)+yAxisLength*0.025, textStr, 'FontSize',10);
    textStr = ['CN(\mu=' num2str(mean(logRates2{1},'omitnan')) ')vCP(\mu=' num2str(mean(logRates2{2},'omitnan')) ') p=' num2str(stats.pConstricted) ', '...
      'DN(\mu=' num2str(mean(logRates2{3},'omitnan')) ')vDP(\mu=' num2str(mean(logRates2{4},'omitnan')) ') p=' num2str(stats.pDilated)];
    text(xLim(2)-xAxisLength*0.99, yLim(2)-yAxisLength*0.01, textStr, 'FontSize',10);
    
    % Save the figure
    label = [3.5 3];
    margin = [0.55 0.75];
    width = options.figSize-label(1)-margin(1);
    height = options.figSize-label(2)-margin(2);
    paperSize = resizeFig(h, gca, width, height, label, margin, 0);
    hgsave(h, [options.figFolder filesep options.figName2 '.fig']);
    exportFig(h, [options.figFolder filesep options.figName2 '.png'],'-dpng','-r300', paperSize);
    close(h);
  elseif min(size(rSpearman)) == 2
    if size(rSpearman,1) < size(rSpearman,2)
      rSpearmanPositive = rSpearman(1,:);
      rSpearmanNegative = rSpearman(2,:);
    else
      rSpearmanPositive = rSpearman(:,1);
      rSpearmanNegative = rSpearman(:,2);
    end
    
    % Positive cells
    options.yLabel = 'Positive unit count';
    options.figName = options.figName2{1};
    statsVarSplit.positive = varTest(logRates(rSpearmanPositive >= 0));
    options.stats = statsVarSplit.positive;
    histPlotFR(edges, logRates(rSpearmanPositive >= 0), options);
    statsMeanSplit.positive = meanTest(logRates(rSpearmanPositive >= 0), 'ANOVARM');
    
    % Negative cells
    options.yLabel = 'Negative unit count';
    options.figName = options.figName2{2};
    statsVarSplit.negative = varTest(logRates(rSpearmanNegative < 0));
    options.stats = statsVarSplit.negative;
    histPlotFR(edges, logRates(rSpearmanNegative < 0), options);
    statsMeanSplit.negative = meanTest(logRates(rSpearmanNegative < 0), 'ANOVARM');
  else
    statsMeanSplit = []; statsVarSplit = [];
  end
else
  statsMeanSplit = []; statsVarSplit = [];
end
end


function firingRateViolinPlot(tickLabels, data, statsInit, filenameFig, options)

if isempty(data)
  return
end

data(isinf(data)) = -3;
[dataMean, dataCI95] = datamean(data);
nCombos = size(statsInit,2);
if numel(dataMean) > 2
  stats.pval = zeros(nCombos-1,1);
  stats.area1 = cell(nCombos-1,1);
  stats.area2 = cell(nCombos-1,1);
else
  stats.pval = zeros(nCombos,1);
  stats.area1 = cell(nCombos,1);
  stats.area2 = cell(nCombos,1);
end
for combo = 1:nCombos
  if numel(dataMean) > 2
    if combo == nCombos
      break
    end
    stats.pval(combo) = statsInit(1,combo+1).p;
    stats.area1{combo} = ['band' num2str(statsInit(1,combo+1).iCol1)];
    stats.area2{combo} = ['band' num2str(statsInit(1,combo+1).iCol2)];
  else
    stats.pval(combo) = statsInit(1,combo).p;
    stats.area1{combo} = ['band' num2str(statsInit(1,combo).iCol1)];
    stats.area2{combo} = ['band' num2str(statsInit(1,combo).iCol2)];
  end
end
phaseBands = cell(1,numel(dataMean));
textStr = 'Means: ';
for iCol = 1:numel(dataMean)
  phaseBands{iCol} = ['band' num2str(iCol)];
  textStr = [textStr num2str(dataMean(iCol)) ' '];
end
fH = multiViolinPlots(data, phaseBands, dataMean, dataCI95, stats, options);
ax1 = gca;
ax1.XTickLabel = tickLabels;
xlabel(options.xLabel);
hold on
x = -pi:0.001:pi;
if ~isempty(options.yLim)
  y = cos(x).*((options.yLim(2)-options.yLim(1))/2)+mean(options.yLim);
else
  y = cos(x);
end
p = plot(1-0.5:(numel(dataMean))/(numel(x)-1):numel(dataMean)+0.5,y, 'k:');
uistack(p, 'bottom');
hold off
xLim = xlim;
xAxisLength = xLim(2)-xLim(1);
if ~isempty(options.yLim)
  ylim(options.yLim);
end
yLim = ylim;
yAxisLength = yLim(2)-yLim(1);
text(xLim(2)-xAxisLength*0.95, yLim(2)-yAxisLength*0.1, textStr, 'FontSize',20);
set(fH, 'Name','Log Firing Rate Distributions across Phase Bands');
filenameFig = [options.figFolder filesep filenameFig];
if options.saveFig
  savefig(fH, filenameFig, 'compact');
  print(fH, filenameFig,'-dpng','-r300');
  close(fH);
end
end


function firingRateDotPlot(tickLabels, data, statsInit, filenameFig, options)

if isempty(data)
  return
end

data(isinf(data)) = -3;
[dataMean, dataCI95] = datamean(data);
nCombos = size(statsInit,2);
if numel(dataMean) > 2
  stats.pval = zeros(nCombos-1,1);
  stats.area1 = cell(nCombos-1,1);
  stats.area2 = cell(nCombos-1,1);
else
  stats.pval = zeros(nCombos,1);
  stats.area1 = cell(nCombos,1);
  stats.area2 = cell(nCombos,1);
end
for combo = 1:nCombos
  if numel(dataMean) > 2
    if combo == nCombos
      break
    end
    stats.pval(combo) = statsInit(1,combo+1).p;
    stats.area1{combo} = ['band' num2str(statsInit(1,combo+1).iCol1)];
    stats.area2{combo} = ['band' num2str(statsInit(1,combo+1).iCol2)];
  else
    stats.pval(combo) = statsInit(1,combo).p;
    stats.area1{combo} = ['band' num2str(statsInit(1,combo).iCol1)];
    stats.area2{combo} = ['band' num2str(statsInit(1,combo).iCol2)];
  end
end
phaseBands = cell(1,numel(dataMean));
textStr = 'Means: ';
for iCol = 1:numel(dataMean)
  phaseBands{iCol} = ['band' num2str(iCol)];
  textStr = [textStr num2str(dataMean(iCol)) ' '];
end
if isfield(options, 'fH')
  figure(options.fH); hold on
else
  fH = figProperties('Dot plot for units', 'normalized', [0, .005, .97, .90], 'w', 'on');
end
for iDot = 1:numel(dataMean)
  errorbar(iDot, dataMean(iDot), abs(dataCI95(1,iDot)), abs(dataCI95(2,iDot)), 'o', 'MarkerSize',20,...
    'Color',matlabColours(iDot), 'MarkerEdgeColor',matlabColours(iDot), 'MarkerFaceColor',matlabColours(iDot));
  if iDot == 1
    hold on
  end
end
xLim = xlim + [-1 1];
xlim(xLim);
yLim = [min(dataMean+dataCI95(1,:))*0.999 max(dataMean+dataCI95(2,:))*1.001];
ylim(yLim);
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Arial', 25, 4/3, 2, [0.003 0.003], 'out',...
  'on', 'k', {options.xLabel}, xLim, 1:numel(dataMean),...
  'on', 'k', {options.yLabel}, yLim, yticks);
ax1.XTickLabel = tickLabels;
x = -pi:0.001:pi;
y = cos(x).*((yLim(2)-yLim(1))/2)+mean(yLim);
p = plot(1-0.5:(numel(dataMean))/(numel(x)-1):numel(dataMean)+0.5,y, 'k:');
uistack(p, 'bottom');
hold off
xLim = xlim;
xAxisLength = xLim(2)-xLim(1);
ylim(yLim);
yAxisLength = yLim(2)-yLim(1);
text(xLim(2)-xAxisLength*0.95, yLim(2)-yAxisLength*0.1, textStr, 'FontSize',20);
if ~isempty(stats) && (~isfield(options, 'textStr') || isempty(options.textStr))
  textStr = [];
  for iBand = 1:numel(phaseBands)
    textStr = [textStr 'v' num2str(iBand)]; %#ok<*AGROW>
  end
  textStr = ['comp ' textStr(2:end) 'v1 p-val: '];
  for iBand = 1:numel(phaseBands)-1
    textStr = [textStr num2str(stats.pval(ismember(stats.area1,phaseBands{iBand}) & ismember(stats.area2,phaseBands{iBand+1}))) ' '];
  end
  textStr = [textStr num2str(stats.pval(ismember(stats.area1,phaseBands{1}) & ismember(stats.area2,phaseBands{iBand+1})))];
  text(xLim(2)-xAxisLength*0.95, yLim(2)-yAxisLength*0.05, textStr, 'FontSize',20);
elseif isfield(options, 'textStr') && ~isempty(options.textStr)
  text(xLim(2)-xAxisLength*0.95, yLim(2)-yAxisLength*0.05, options.textStr, 'FontSize',20);
end
set(fH, 'Name','Log Firing Rate Means across Phase Bands');
filenameFig = [options.figFolder filesep filenameFig];
if options.saveFig
  savefig(fH, filenameFig, 'compact');
  print(fH, filenameFig,'-dpng','-r300');
  close(fH);
end
end


function firingRateBarPlot(tickLabels, data, statsInit, filenameFig, options)

if isempty(data)
  return
end

data(isinf(data)) = -3;
dataVar = var(data, 'omitnan');
nCombos = size(statsInit,2);
if numel(dataVar) > 2
  stats.pval = zeros(nCombos-1,1);
  stats.area1 = cell(nCombos-1,1);
  stats.area2 = cell(nCombos-1,1);
else
  stats.pval = zeros(nCombos,1);
  stats.area1 = cell(nCombos,1);
  stats.area2 = cell(nCombos,1);
end
for combo = 1:nCombos
  if numel(dataVar) > 2
    if combo == nCombos
      break
    end
    stats.pval(combo) = statsInit(1,combo+1).p;
    stats.area1{combo} = ['band' num2str(statsInit(1,combo+1).iCol1)];
    stats.area2{combo} = ['band' num2str(statsInit(1,combo+1).iCol2)];
  else
    stats.pval(combo) = statsInit(1,combo).p;
    stats.area1{combo} = ['band' num2str(statsInit(1,combo).iCol1)];
    stats.area2{combo} = ['band' num2str(statsInit(1,combo).iCol2)];
  end
end
phaseBands = cell(1,numel(dataVar));
textStr = 's^2: ';
for iCol = 1:numel(dataVar)
  phaseBands{iCol} = ['band' num2str(iCol)];
  textStr = [textStr num2str(dataVar(iCol)) ' '];
end
saveFig = options.saveFig;
options.saveFig = false;
fH = histPlot(0:numel(dataVar), dataVar, options);
xticks(1:numel(dataVar));
ax1 = gca;
ax1.XTickLabel = tickLabels;
xlabel(options.xLabel);
yTick = ceil(max(dataVar)*100)/200;
yticks([0 yTick 2*yTick]);
hold on
x = -pi:0.001:pi;
if ~isempty(options.yLim)
  yLim = options.yLim;
else
  yLim = [min(dataVar)*0.9 max(dataVar)*1.1];
end
y = cos(x).*((yLim(2)-yLim(1))/2)+mean(yLim);
p = plot(1-0.5:(numel(dataVar))/(numel(x)-1):numel(dataVar)+0.5,y, 'k:');
uistack(p, 'bottom');
hold off
xLim = xlim;
xAxisLength = xLim(2)-xLim(1);
ylim(yLim);
yAxisLength = yLim(2)-yLim(1);
text(xLim(2)-xAxisLength*0.95, yLim(2)-yAxisLength*0.1, textStr, 'FontSize',20);
textStr = [];
for iArea = 1:numel(dataVar)
  textStr = [textStr 'v' num2str(iArea)]; %#ok<*AGROW>
end
textStr = ['comp ' textStr(2:end) 'v1 p-val: '];
for iArea = 1:numel(dataVar)-1
  textStr = [textStr num2str(stats.pval(ismember(stats.area1,phaseBands{iArea}) & ismember(stats.area2,phaseBands{iArea+1}))) ' '];
end
textStr = [textStr num2str(stats.pval(ismember(stats.area1,phaseBands{1}) & ismember(stats.area2,phaseBands{iArea+1})))];
text(xLim(2)-xAxisLength*0.95, yLim(2)-yAxisLength*0.05, textStr, 'FontSize',20);
yticks(yLim(1):yAxisLength/5:yLim(2));
set(fH, 'Name','Log Firing Rate Distributions across Phase Bands');
filenameFig = [options.figFolder filesep filenameFig];
if saveFig
  savefig(fH, filenameFig, 'compact');
  print(fH, filenameFig,'-dpng','-r300');
  close(fH);
end
end


function firingRateCombinedDotPlot(tickLabels, data, filenameFig, options)

if isempty(data)
  return
end

% Descriptive stats
[dataMean, dataCI95] = datamean(data);
dataVar = var(data, 'omitnan');
alpha = 0.05;
df = size(data,1)-1;
dataVarCI95 = [(df.*var(data,1,'omitnan'))./chi2inv(1-alpha/2,df); (df.*var(data,1,'omitnan'))./chi2inv(alpha/2,df)] - dataVar;

% Figure:
fH = figProperties('Dot plot for units', 'normalized', [0, .005, .97, .90], 'w', 'on');
% Mean
yyaxis left
xInterp = 1:0.001:numel(dataMean);
plot(xInterp, interp1(1:numel(dataMean),dataMean,xInterp,'cubic'), 'Color',options.colour{1});
hold on
for iDot = 1:numel(dataMean)
  errorbar(iDot, dataMean(iDot), abs(dataCI95(1,iDot)), abs(dataCI95(2,iDot)), 'o', 'MarkerSize',10,...
    'Color',options.colour{1}, 'MarkerEdgeColor',options.colour{1}, 'MarkerFaceColor',options.colour{1});
end
yLim = [min(dataMean+dataCI95(1,:))*0.999 max(dataMean+dataCI95(2,:))*1.001];
ylim(yLim);
xLim = xlim + [-1 1];
xlim(xLim);
axesProperties({}, 1, 'normal', 'off', 'w', 'Arial', 25, 4/3, 2, [0.003 0.003], 'out',...
  'on', 'k', {options.xLabel}, xLim, 1:numel(dataMean)/4:numel(dataMean),...
  'on', options.colour{1}, {options.yLabel}, yLim, yticks);
hold off
% Variance
yyaxis right
plot(xInterp, interp1(1:numel(dataVar),dataVar,xInterp,'cubic'), 'Color',options.colour{2});
hold on
for iDot = 1:numel(dataVar)
  errorbar(iDot, dataVar(iDot), abs(dataVarCI95(1,iDot)), abs(dataVarCI95(2,iDot)), 's', 'MarkerSize',10,...
    'Color',options.colour{2}, 'MarkerEdgeColor',options.colour{2}, 'MarkerFaceColor',options.colour{2});
end
yLim = [min(dataVar+dataVarCI95(1,:))*0.999 max(dataVar+dataVarCI95(2,:))*1.001];
ylim(yLim);
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Arial', 25, 4/3, 2, [0.003 0.003], 'out',...
  'on', 'k', {options.xLabel}, xLim, 0.5:numel(dataMean)/4:numel(dataMean)+0.5,...
  'on', options.colour{2}, {options.yLabel}, yLim, yticks);
% Background phase sinusoid
x = -pi:0.001:pi;
y = cos(x).*((yLim(2)-yLim(1))/2)+mean(yLim);
p = plot(1-0.5:(numel(dataMean))/(numel(x)-1):numel(dataMean)+0.5,y, 'k:'); %#ok<*NASGU>
%uistack(p, 'bottom');
hold off
%ax1.XTickLabel = tickLabels;
%ax1.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
ax1.XTickLabel = {'\pi','\pi/2','0','-\pi/2','-\pi'};

% Stats display
% Var test
stats = struct2table(varTest(data), 'AsArray',true);
varMat = [stats.var1 stats.var2];
iComp = min(varMat,[],2) == min(min(varMat)) & max(varMat,[],2) == max(max(varMat));
if sum(iComp) > 1
  iComp0 = find(iComp,1);
  iComp = false(size(iComp));
  iComp(iComp0) = true;
end
dist1 = tickLabels{stats.iCol1(iComp)};
dist2 = tickLabels{stats.iCol2(iComp)};
var1 = round(stats.var1(iComp),2);
var2 = round(stats.var2(iComp),2);
pval = stats.p(iComp);

xAxisLength = xLim(2)-xLim(1);
yAxisLength = yLim(2)-yLim(1);
textStr = [dist1 ' (\sigma^2=' num2str(var1) ') vs ' dist2 ' (\sigma^2=' num2str(var2) ') p=' num2str(pval)];
text(xLim(2)-xAxisLength*0.85, yLim(2)-yAxisLength*0.01, textStr, 'FontSize',12);

% Means t-test
[~, minCol] = min(dataMean);
[~, maxCol] = max(dataMean);
[~, pval, ~, stats] = ttest(data(:,minCol), data(:,maxCol));
stats.p = pval;
dist1 = tickLabels{minCol};
dist2 = tickLabels{maxCol};
mean1 = dataMean(minCol);
mean2 = dataMean(maxCol);

textStr = [dist1 ' (\mu=' num2str(mean1) ') vs ' dist2 ' (\mu=' num2str(mean2) ') p=' num2str(pval)];
text(xLim(2)-xAxisLength*0.85, yLim(2)+yAxisLength*0.025, textStr, 'FontSize',12);

% Save the figure
label = [4.25 3];
margin = [4.25 0.55];
width = options.figSize(1)-label(1)-margin(1);
height = options.figSize(2)-label(2)-margin(2);
resizeFig(fH, gca, width, height, label, margin, 0);
set(fH, 'Name','Log Firing Rate Means and Variances across Phase Bands');
if options.saveFig
  savefig(fH, [options.figFolder filesep filenameFig '.fig'], 'compact');
  print(fH, [options.figFolder filesep filenameFig '.png'],'-dpng','-r300');
  close(fH);
end
end


function [h, predictions] = firingRatePredictPlot(edges, centres, x, frData, model, options)

% Preprocess data
exclude = frData(:,1) < 1e-9;
frData = frData(~exclude,:);
frDistro = histcounts(getLog(frData), edges);

% Fit the mean firing rate distribution
if strcmpi(options.distroType, 'skewnormal')
  frDistroInterp = interp1(centres, frDistro, x);
  [mode, iMode] = max(frDistro);
  ft = fittype('skewedgaussian(x, alpha, xi, omega, omega2)');
  f = fit(x', frDistroInterp', ft, 'StartPoint',[-centres(iMode), mode, 1, centres(iMode)]);
  fittedDistro = skewedgaussian(x, f.alpha, f.xi, f.omega, f.omega2);
  fittedDistroDS = interp1(x, fittedDistro, centres);
elseif strcmpi(options.distroType, 'gamma')
  fitParams = gamfit(frData);
  xNotLog = 10^edges(1):0.01:10^edges(end);
  fittedDistro = gampdf(xNotLog,fitParams(1),fitParams(2));
end

% Simulate firing rates based on the fitted distribution
if strcmpi(options.distroType, 'skewnormal')
  nSamples = 1000000;
  fittedCDF = cumsum(fittedDistro)./sum(fittedDistro);
  R = rand(nSamples,1);
  p = @(r) find(r<fittedCDF,1,'first');
  rR = arrayfun(p,R);
  %hist(rR,1:numel(fittedDistro));
  simulatedLogFR = x(rR);
elseif strcmpi(options.distroType, 'gamma')
  nSamples = 1000000;
  fittedCDF = cumsum(fittedDistro)./sum(fittedDistro);
  R = rand(nSamples,1);
  p = @(r) find(r<fittedCDF,1,'first');
  rR = arrayfun(p,R);
  simulatedLogFR = getLog(xNotLog(rR));
  %hist(simulatedLogFR,edges);
end

% Predict constricted and dilated firing rates
if strcmpi(options.distroType, 'skewnormal')
  predictedModulationSize = simulatedLogFR.*model(1) + model(2);
  predictedModulationSizeRelativeToMean = predictedModulationSize/2;
  predictedFRConsitricted = simulatedLogFR - predictedModulationSizeRelativeToMean;
  %predictedFRConsitrictedDistro = histcounts(predictedFRConsitricted, [x(1)-(x(2)-x(1)) x+(x(2)-x(1))]).*(sum(fittedDistro)/numel(predictedFRConsitricted));
  predictedFRConsitrictedDistroDS = histcounts(predictedFRConsitricted, edges).*(sum(fittedDistroDS)/nSamples);
  predictedFRDilated = simulatedLogFR + predictedModulationSizeRelativeToMean;
  %predictedFRDilatedDistro = histcounts(predictedFRDilated, [x(1)-(x(2)-x(1)) x+(x(2)-x(1))]).*(sum(fittedDistro)/numel(predictedFRDilated));
  predictedFRDilatedDistroDS = histcounts(predictedFRDilated, edges).*(sum(fittedDistroDS)/nSamples);
  predictions = [predictedFRConsitrictedDistroDS; predictedFRDilatedDistroDS];
elseif strcmpi(options.distroType, 'gamma')
  predictedModulationSize = simulatedLogFR.*model(1) + model(2);
  predictedModulationSizeRelativeToMean = predictedModulationSize/2;
  predictedFRConsitricted = simulatedLogFR - predictedModulationSizeRelativeToMean;
  %predictedFRConsitrictedDistro = histcounts(predictedFRConsitricted, [x(1)-(x(2)-x(1)) x+(x(2)-x(1))]).*(sum(fittedDistro)/numel(predictedFRConsitricted));
  predictedFRConsitrictedDistroDS = histcounts(predictedFRConsitricted, edges).*(size(frData,1)/nSamples);
  predictedFRDilated = simulatedLogFR + predictedModulationSizeRelativeToMean;
  %predictedFRDilatedDistro = histcounts(predictedFRDilated, [x(1)-(x(2)-x(1)) x+(x(2)-x(1))]).*(sum(fittedDistro)/numel(predictedFRDilated));
  predictedFRDilatedDistroDS = histcounts(predictedFRDilated, edges).*(size(frData,1)/nSamples);
  predictions = [predictedFRConsitrictedDistroDS; predictedFRDilatedDistroDS];
end

% Plot the figure
optionsIntermediate = options;
optionsIntermediate.saveFig = false;
optionsIntermediate.binnedData = predictions;
h = histPlotFR(edges, [predictedFRConsitricted' predictedFRDilated'], optionsIntermediate); hold on
if ~isfield(options,'fittedMean') || (isfield(options,'fittedMean') && options.fittedMean)
  plot(x, fittedDistro, 'k:'); hold off
  l = legend;
  l.String = [options.legendLabels 'Fitted mean'];
end

% Save the figure
label = [3.5 3];
margin = [0.55 0.55];
width = options.figSize-label(1)-margin(1);
height = options.figSize-label(2)-margin(2);
paperSize = resizeFig(h, gca, width, height, label, margin, 0);
if options.saveFig
  hgsave(h, [options.figFolder filesep options.figName '.fig']);
  exportFig(h, [options.figFolder filesep options.figName '.png'],'-dpng','-r300', paperSize);
  close(h);
end
end


function [h, predictions] = firingRatePredictPlot2(edges, centres, x, frData, model, options)

% Preprocess data
exclude1 = frData(:,1) < 1e-9;
exclude2 = frData(:,2) < 1e-9;
frData = frData(~exclude1 & ~exclude2,:);
frDistro = histcounts(getLog(frData), edges);

% Fit the mean firing rate distribution
if ~isfield(options,'fittedMean') || (isfield(options,'fittedMean') && options.fittedMean)
  ft = fittype('skewedgaussian(x, alpha, xi, omega, omega2)');
  frDistroInterp = interp1(centres, frDistro, x);
  [mode, iMode] = max(frDistro);
  f = fit(x', frDistroInterp', ft, 'StartPoint',[-centres(iMode), mode, 1, centres(iMode)]);
  fittedDistro = skewedgaussian(x, f.alpha, f.xi, f.omega, f.omega2);
  fittedDistroDS = interp1(x, fittedDistro, centres);
end

% Predict constricted and dilated firing rates
predictedModulationSize = frData.*model(1) + model(2);
predictedModulationSizeRelativeToMean = predictedModulationSize/2;
predictedFRConsitricted = frData - predictedModulationSizeRelativeToMean;
predictedFRConsitrictedDistroDS = histcounts(predictedFRConsitricted, edges);
predictedFRDilated = frData + predictedModulationSizeRelativeToMean;
predictedFRDilatedDistroDS = histcounts(predictedFRDilated, edges);
predictions = [predictedFRConsitrictedDistroDS; predictedFRDilatedDistroDS];

% Fit the predicted distributions
ft = fittype('skewedgaussian(x, alpha, xi, omega, omega2)');
predictedFRConsitrictedDistroDS = interp1(centres, predictedFRConsitrictedDistroDS, x);
[mode, iMode] = max(predictedFRConsitrictedDistroDS);
f = fit(x', predictedFRConsitrictedDistroDS', ft, 'StartPoint',[-x(iMode), mode, 1, x(iMode)]);
predictedFRConsitrictedDistroDS_fit = skewedgaussian(x, f.alpha, f.xi, f.omega, f.omega2);
%predictedFRConsitrictedDistroDS_fit = interp1(x, fittedDistro, centres);

ft = fittype('skewedgaussian(x, alpha, xi, omega, omega2)');
predictedFRDilatedDistroDS = interp1(centres, predictedFRDilatedDistroDS, x);
[mode, iMode] = max(predictedFRDilatedDistroDS);
f = fit(x', predictedFRDilatedDistroDS', ft, 'StartPoint',[-x(iMode), mode, 1, x(iMode)]);
predictedFRDilatedDistroDS_fit = skewedgaussian(x, f.alpha, f.xi, f.omega, f.omega2);
%predictedFRDilatedDistroDS_fit = interp1(x, fittedDistro, centres);
predictionFits = [predictedFRConsitrictedDistroDS_fit; predictedFRDilatedDistroDS_fit];

% Plot the figure
optionsIntermediate = options;
optionsIntermediate.saveFig = false;
optionsIntermediate.binnedData = predictionFits;
h = histPlotFR([edges(1) x], [predictedFRConsitrictedDistroDS_fit' predictedFRDilatedDistroDS_fit'], optionsIntermediate); hold on
%optionsIntermediate.binnedData = predictions;
%h = histPlotFR(edges, [predictedFRConsitrictedDistroDS' predictedFRDilatedDistroDS'], optionsIntermediate); hold on
if ~isfield(options,'fittedMean') || (isfield(options,'fittedMean') && options.fittedMean)
  plot(x, fittedDistro, 'k:'); hold off
  l = legend;
  l.String = [options.legendLabels 'Fitted mean'];
end

% Save the figure
label = [3.5 3];
margin = [0.55 0.55];
width = options.figSize-label(1)-margin(1);
height = options.figSize-label(2)-margin(2);
paperSize = resizeFig(h, gca, width, height, label, margin, 0);
if options.saveFig
  hgsave(h, [options.figFolder filesep options.figName '.fig']);
  exportFig(h, [options.figFolder filesep options.figName '.png'],'-dpng','-r300', paperSize);
  close(h);
end
end


function firingRateCombinedPlot(data, edges, centres, x, frData, model, legendLabels1, legendLabels2, options)

% Plot data distributions
logRates = getLog(data(:,[1 end]));
options.legendLabels = legendLabels1;
options.lineStyles = {':',':'};
options.lineWidths = [1.5 1.5];
options.markerStyles = {'.','o'};
options.displayMeans = true;
options.displayVars = true;
options.markerVPos = 1.15*max([max(histcounts(logRates(:,1), edges)) max(histcounts(logRates(:,2), edges))]);
options.fH = [];
options.fH = histPlotFR(edges, logRates, options);

% Plot predicted distributions
options.legendLabels = legendLabels2;
options.lineStyles = {'-','-'};
options.lineWidths = [2 2];
options.markerStyles = {'.','o'};
options.displayMeans = false;
options.displayVars = false;
options.fittedMean = false;
options.markerVPos = 1.15*max([max(histcounts(logRates(:,1), edges)) max(histcounts(logRates(:,2), edges))]);
options.fH = firingRatePredictPlot(edges, centres, x, frData, model, options);

% Display stats
digit = 3;
statsMean = meanTest(logRates, 'ANOVARM');
statsVar = varTest(logRates);
text(max(datamean(logRates))+0.2,options.markerVPos, ['p=' num2str(round(statsMean.p,digit,'significant'))], 'FontSize',20);
text(max(datamean(logRates))-0.6,0.125*max([max(histcounts(logRates(:,1), edges)) max(histcounts(logRates(:,2), edges))]),...
  ['p=' num2str(round(statsVar.p,digit,'significant'))], 'FontSize',20);

% Save the figure
label = [3.5 3];
margin = [0.6 0.55];
width = options.figSize-label(1)-margin(1);
height = options.figSize-label(2)-margin(2);
paperSize = resizeFig(options.fH, gca, width, height, label, margin, 0);
hgsave(options.fH, [options.figFolder filesep options.figName '.fig']);
exportFig(options.fH, [options.figFolder filesep options.figName '.png'],'-dpng','-r300', paperSize);
close(options.fH);
end


function firingRateCombinedPlot2(data, edges, centres, x, frData, model, legendLabels1, legendLabels2, options)

% Plot data distributions
logRates = getLog(data(:,[1 end]));
options.legendLabels = legendLabels1;
options.lineStyles = {':',':'};
options.lineWidths = [1.5 1.5];
options.markerStyles = {'.','o'};
options.displayMeans = true;
options.displayVars = true;
options.markerVPos = 1.15*max([max(histcounts(logRates(:,1), edges)) max(histcounts(logRates(:,2), edges))]);
options.fH = [];
options.fH = histPlotFR(edges, logRates, options);

% Plot predicted distributions
options.legendLabels = legendLabels2;
options.lineStyles = {'-','-'};
options.lineWidths = [2 2];
options.markerStyles = {'.','o'};
options.displayMeans = false;
options.displayVars = false;
options.fittedMean = false;
options.markerVPos = 1.15*max([max(histcounts(logRates(:,1), edges)) max(histcounts(logRates(:,2), edges))]);
options.fH = firingRatePredictPlot2(edges, centres, x, frData, model, options);

% Display stats
digit = 3;
statsMean = meanTest(logRates, 'ANOVARM');
statsVar = varTest(logRates);
text(max(datamean(logRates))+0.2,options.markerVPos, ['p=' num2str(round(statsMean.p,digit,'significant'))], 'FontSize',20);
text(max(datamean(logRates))-0.6,0.125*max([max(histcounts(logRates(:,1), edges)) max(histcounts(logRates(:,2), edges))]),...
  ['p=' num2str(round(statsVar.p,digit,'significant'))], 'FontSize',20);

% Save the figure
label = [3.5 3];
margin = [0.6 0.55];
width = options.figSize-label(1)-margin(1);
height = options.figSize-label(2)-margin(2);
paperSize = resizeFig(options.fH, gca, width, height, label, margin, 0);
hgsave(options.fH, [options.figFolder filesep options.figName '.fig']);
exportFig(options.fH, [options.figFolder filesep options.figName '.png'],'-dpng','-r300', paperSize);
close(options.fH);
end


function firingRateCombinedSplitPlot(data, edges, dataDistro, legendLabels1, legendLabels2, options)

% Plot data distributions
logRates = getLog(data(:,[1 end]));
options.legendLabels = legendLabels1;
options.lineStyles = {':',':'};
options.lineWidths = [1.5 1.5];
options.markerStyles = {'.','o'};
options.displayMeans = true;
options.displayVars = true;
options.markerVPos = 1.15*max([max(histcounts(logRates(:,1), edges)) max(histcounts(logRates(:,2), edges))]);
options.fH = [];
options.fH = histPlotFR(edges, logRates, options);

% Plot predicted distributions
options.legendLabels = legendLabels2;
options.lineStyles = {'-','-'};
options.lineWidths = [2 2];
options.markerStyles = {'.','o'};
options.displayMeans = false;
options.displayVars = false;
options.fittedMean = false;
options.markerVPos = 1.15*max([max(histcounts(logRates(:,1), edges)) max(histcounts(logRates(:,2), edges))]);
options.binnedData = dataDistro;
options.fH = histPlotFR(edges, logRates, options);

% Display stats
digit = 3;
statsMean = meanTest(logRates, 'ANOVARM');
statsVar = varTest(logRates);
text(max(datamean(logRates))+0.2,options.markerVPos, ['p=' num2str(round(statsMean.p,digit,'significant'))], 'FontSize',20);
text(max(datamean(logRates))-0.6,0.125*max([max(histcounts(logRates(:,1), edges)) max(histcounts(logRates(:,2), edges))]),...
  ['p=' num2str(round(statsVar.p,digit,'significant'))], 'FontSize',20);

% Save the figure
label = [3.5 3];
margin = [0.6 0.55];
width = options.figSize-label(1)-margin(1);
height = options.figSize-label(2)-margin(2);
paperSize = resizeFig(options.fH, gca, width, height, label, margin, 0);
hgsave(options.fH, [options.figFolder filesep options.figName '.fig']);
exportFig(options.fH, [options.figFolder filesep options.figName '.png'],'-dpng','-r300', paperSize);
close(options.fH);
end


function firingRateSmoothedPlot(data, edges, legendLabels, options)

options.lineStyles = {'-','-'};
options.displayMeans = true;
options.displayVars = true;
options.fH = [];
%[options, logRates] = firingRateSimulatePlot(data, edges, legendLabels, options);

options.binnedData = [smoothdata(histcounts(getLog(data(:,1)), edges),'gaussian',31);...
  smoothdata(histcounts(getLog(data(:,end)), edges),'gaussian',31)];
logRates = [getLog(data(:,1)) getLog(data(:,end))];
options.legendLabels = legendLabels;
options.lineWidths = [1.5 1.5];
options.markerStyles = {'.','o'};
% options.markerVPos = 1.15*max([max(options.binnedData(1,:))/sum(options.binnedData(1,:))...
%   max(options.binnedData(2,:))/sum(options.binnedData(2,:))]);
[fH, ~, options] = histPlotFR(edges, logRates, options);
options.fH = fH;

% Display stats
digit = 3;
statsMean = meanTest(logRates, 'ANOVARM');
statsVar = varTest(logRates);
text(max(datamean(logRates))+0.2,options.markerVPos, ['p=' num2str(round(statsMean.p,digit,'significant'))], 'FontSize',20);
text(max(datamean(logRates))-0.6,0.125*options.markerVPos, ['p=' num2str(round(statsVar.p,digit,'significant'))], 'FontSize',20);

% Save the figure
label = [3.5 3];
margin = [1.5 0.55];
width = options.figSize-label(1)-margin(1);
height = options.figSize-label(2)-margin(2);
paperSize = resizeFig(options.fH, gca, width, height, label, margin, 0);
hgsave(options.fH, [options.figFolder filesep options.figName '.fig']);
exportFig(options.fH, [options.figFolder filesep options.figName '.png'],'-dpng','-r300', paperSize);
close(options.fH);
end


function firingRateSplitSmoothedPlot(positiveData, negativeData, edges, options)

% options.fH = [];
% options.displayMeans = false;
% options.displayVars = false;
% [options, positiveLogRates] = firingRateSimulatePlot(positiveData, edges, legendLabels(1:2), options);
% [options, negativeLogRates] = firingRateSimulatePlot(negativeData, edges, legendLabels(3:4), options);
positiveLogRates = getLog(positiveData);
negativeLogRates = getLog(negativeData);

options.colours = [matlabColours(1); matlabColours(2); matlabColours(1); matlabColours(2)];
options.lineStyles = {'--','--','-','-'};
options.markerStyles = {'v','v','^','^'};
options.legendLabels = {'Constricted negative','Constricted positive','Dilated negative','Dilated positive'};
options.stats = [];
options.saveFig = false;
options.binnedData = [smoothdata(histcounts(negativeLogRates(:,1), edges),'gaussian',31);...
  smoothdata(histcounts(positiveLogRates(:,1), edges),'gaussian',31);...
  smoothdata(histcounts(negativeLogRates(:,end), edges),'gaussian',31);...
  smoothdata(histcounts(positiveLogRates(:,end), edges),'gaussian',31)];
logRates = {negativeLogRates(:,1), positiveLogRates(:,1), negativeLogRates(:,end), positiveLogRates(:,end)};
h = histPlotFR(edges, logRates, options);
l = legend;
l.Position = l.Position - [0.08 0.08 0 0];

% Display stats
logRates = NaN(max([size(positiveLogRates,1) size(negativeLogRates,1)]),4);
logRates(1:size(negativeLogRates,1),1) = negativeLogRates(:,1);
logRates(1:size(positiveLogRates,1),2) = positiveLogRates(:,1);
logRates(1:size(negativeLogRates,1),3) = negativeLogRates(:,end);
logRates(1:size(positiveLogRates,1),4) = positiveLogRates(:,end);
statsVarSplit = varTest(logRates);
[statsMeanSplit.pNegative, ~, statsMeanSplit.statsNegative] = signrank(negativeLogRates(:,1), negativeLogRates(:,end));
[statsMeanSplit.pPositive, ~, statsMeanSplit.statsPositive] = signrank(positiveLogRates(:,1), positiveLogRates(:,end));
[statsMeanSplit.pConstricted, ~, statsMeanSplit.statsConstricted] = ranksum(negativeLogRates(:,1), positiveLogRates(:,1));
[statsMeanSplit.pDilated, ~, statsMeanSplit.statsDilated] = ranksum(negativeLogRates(:,end), positiveLogRates(:,end));

% Display var test stats
statsVarSplit = varTest(logRates);
xLim = xlim;
xAxisLength = xLim(2)-xLim(1);
yLim = ylim;
yAxisLength = yLim(2)-yLim(1);
stats = struct2table(statsVarSplit, 'AsArray',true);
iComp1 = stats.iCol1 == 1 & stats.iCol2 == 3;
iComp2 = stats.iCol1 == 2 & stats.iCol2 == 4;
iComp3 = stats.iCol1 == 1 & stats.iCol2 == 2;
iComp4 = stats.iCol1 == 3 & stats.iCol2 == 4;
textStr = ['CN(\sigma^2=' num2str(stats.var1(iComp1)) ')vDN(\sigma^2=' num2str(stats.var2(iComp1)) ') p=' num2str(stats.p(iComp1)) ', '...
  'CP(\sigma^2=' num2str(stats.var1(iComp2)) ')vDP(\sigma^2=' num2str(stats.var2(iComp2)) ') p=' num2str(stats.p(iComp2))];
text(xLim(2)-xAxisLength*0.99, yLim(2)-yAxisLength*0.045, textStr, 'FontSize',10);
textStr = ['CN(\sigma^2=' num2str(stats.var1(iComp3)) ')vCP(\sigma^2=' num2str(stats.var2(iComp3)) ') p=' num2str(stats.p(iComp3)) ', '...
  'DN(\sigma^2=' num2str(stats.var1(iComp4)) ')vDP(\sigma^2=' num2str(stats.var2(iComp4)) ') p=' num2str(stats.p(iComp4))];
text(xLim(2)-xAxisLength*0.99, yLim(2)-yAxisLength*0.08, textStr, 'FontSize',10);

% Display means test stats
%     [statsMeanSplit.pNegative, ~, statsMeanSplit.statsNegative] = signrank(logRates(:,1), logRates(:,3));
%     [statsMeanSplit.pPositive, ~, statsMeanSplit.statsPositive] = signrank(logRates(:,2), logRates(:,4));
%     [statsMeanSplit.pConstricted, ~, statsMeanSplit.statsConstricted] = ranksum(logRates(:,1), logRates(:,2));
%     [statsMeanSplit.pDilated, ~, statsMeanSplit.statsDilated] = ranksum(logRates(:,3), logRates(:,4));
[~, statsMeanSplit.pNegative, ~, statsMeanSplit.statsNegative] = ttest(logRates(:,1), logRates(:,3));
[~, statsMeanSplit.pPositive, ~, statsMeanSplit.statsPositive] = ttest(logRates(:,2), logRates(:,4));
[~, statsMeanSplit.pConstricted, ~, statsMeanSplit.statsConstricted] = ttest2(logRates(:,1), logRates(:,2));
[~, statsMeanSplit.pDilated, ~, statsMeanSplit.statsDilated] = ttest2(logRates(:,3), logRates(:,4));
statsMeanSplit.meansNegative = [mean(logRates(:,1)), mean(logRates(:,3))];
statsMeanSplit.meansPositive = [mean(logRates(:,2)), mean(logRates(:,4))];
statsMeanSplit.meansConstricted = [mean(logRates(:,1)), mean(logRates(:,2))];
statsMeanSplit.meansDilated = [mean(logRates(:,3)), mean(logRates(:,4))];
[~, ci1] = datamean(logRates(:,1)); [~, ci2] = datamean(logRates(:,3)); statsMeanSplit.ciNegative = [ci1 ci2];
[~, ci1] = datamean(logRates(:,2)); [~, ci2] = datamean(logRates(:,4)); statsMeanSplit.ciPositive = [ci1 ci2];
[~, ci1] = datamean(logRates(:,1)); [~, ci2] = datamean(logRates(:,2)); statsMeanSplit.ciConstricted = [ci1 ci2];
[~, ci1] = datamean(logRates(:,3)); [~, ci2] = datamean(logRates(:,4)); statsMeanSplit.ciDilated = [ci1 ci2];
stats = statsMeanSplit;
textStr = ['CN(\mu=' num2str(mean(logRates(:,1),'omitnan')) ')vDN(\mu=' num2str(mean(logRates(:,3),'omitnan')) ') p=' num2str(stats.pNegative) ', '...
  'CP(\mu=' num2str(mean(logRates(:,2),'omitnan')) ')vDP(\mu=' num2str(mean(logRates(:,4),'omitnan')) ') p=' num2str(stats.pPositive)];
text(xLim(2)-xAxisLength*0.99, yLim(2)+yAxisLength*0.025, textStr, 'FontSize',10);
textStr = ['CN(\mu=' num2str(mean(logRates(:,1),'omitnan')) ')vCP(\mu=' num2str(mean(logRates(:,2),'omitnan')) ') p=' num2str(stats.pConstricted) ', '...
  'DN(\mu=' num2str(mean(logRates(:,3),'omitnan')) ')vDP(\mu=' num2str(mean(logRates(:,4),'omitnan')) ') p=' num2str(stats.pDilated)];
text(xLim(2)-xAxisLength*0.99, yLim(2)-yAxisLength*0.01, textStr, 'FontSize',10);

% Save the figure
label = [3.5 3];
margin = [0.6 0.55];
width = options.figSize-label(1)-margin(1);
height = options.figSize-label(2)-margin(2);
paperSize = resizeFig(h, gca, width, height, label, margin, 0);
hgsave(h, [options.figFolder filesep options.figName '.fig']);
exportFig(h, [options.figFolder filesep options.figName '.png'],'-dpng','-r300', paperSize);
close(h);
end


function [options, logRates, scaleFactorSmall, scaleFactorLarge] = firingRateSimulatePlot(data, edges, legendLabels, options)

% Preprocess data
dataSmallPupil = data(:,1);
dataLargePupil = data(:,2);
exclude = dataSmallPupil < 1e-9 | dataLargePupil < 1e-9;
dataSmallPupil = dataSmallPupil(~exclude,:);
dataLargePupil = dataLargePupil(~exclude,:);
logRates = [getLog(dataSmallPupil) getLog(dataLargePupil)]; 
frDistroSmallPupil = histcounts(getLog(dataSmallPupil), edges);
frDistroLargePupil = histcounts(getLog(dataLargePupil), edges);

% Fit the distributions
fitParams = gamfit(dataSmallPupil);
xNotLog = 10^edges(1):0.01:10^edges(end);
frDistroSmallPupilFitted = gampdf(xNotLog,fitParams(1),fitParams(2));

fitParams = gamfit(dataLargePupil);
xNotLog = 10^edges(1):0.01:10^edges(end);
frDistroLargePupilFitted = gampdf(xNotLog,fitParams(1),fitParams(2));

% Simulate firing rates based on the fitted distribution
nSamples = 1000000;
fittedCDF = cumsum(frDistroSmallPupilFitted)./sum(frDistroSmallPupilFitted);
R = rand(nSamples,1);
p = @(r) find(r<fittedCDF,1,'first');
rR = arrayfun(p,R);
simulatedLogFRSmallPupil = getLog(xNotLog(rR));

fittedCDF = cumsum(frDistroLargePupilFitted)./sum(frDistroLargePupilFitted);
R = rand(nSamples,1);
p = @(r) find(r<fittedCDF,1,'first');
rR = arrayfun(p,R);
simulatedLogFRLargePupil = getLog(xNotLog(rR));

% Plot data distributions
dataFitted = [simulatedLogFRSmallPupil' simulatedLogFRLargePupil'];
scaleFactorSmall = length(dataSmallPupil)/length(simulatedLogFRSmallPupil);
scaleFactorLarge = length(dataLargePupil)/length(simulatedLogFRLargePupil);
options.binnedData = [histcounts(dataFitted(:,1), edges).*scaleFactorSmall;...
  histcounts(dataFitted(:,2), edges).*scaleFactorLarge];
options.legendLabels = legendLabels;
options.lineWidths = [1.5 1.5];
options.markerStyles = {'.','o'};
options.markerVPos = 1.15*max([max(histcounts(dataFitted(:,1), edges))*scaleFactorSmall...
  max(histcounts(dataFitted(:,2), edges))*scaleFactorLarge]);
options.fH = histPlotFR(edges, dataFitted, options);
end