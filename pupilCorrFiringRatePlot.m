clearvars -except allData
params
lists

repository = 'allensdk';
alpha = 0.05; % significance level
fullRun = false;

if strcmp(repository,'all')
  mainFolder = [dataDir filesep paDir];
  animals = animalsOI;
elseif strcmp(repository,'uol')
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
  if fullRun
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
if fullRun
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
      firingRatePositiveUnits = {};
      firingRatePositiveSignificantUnits = {};
      firingRatePositiveNonsignificantUnits = {};
      firingRatePositiveMUAs = {};
      firingRatePositiveSignificantMUAs = {};
      firingRatePositiveNonsignificantMUAs = {};
      firingRateNegativeUnits = {};
      firingRateNegativeSignificantUnits = {};
      firingRateNegativeNonsignificantUnits = {};
      firingRateNegativeMUAs = {};
      firingRateNegativeSignificantMUAs = {};
      firingRateNegativeNonsignificantMUAs = {};
      
      for iCond = 1:numel(conditions)
        firingRatePositiveUnitsCond = {};
        firingRatePositiveSignificantUnitsCond = {};
        firingRatePositiveNonsignificantUnitsCond = {};
        firingRatePositiveMUAsCond = {};
        firingRatePositiveSignificantMUAsCond = {};
        firingRatePositiveNonsignificantMUAsCond = {};
        firingRateNegativeUnitsCond = {};
        firingRateNegativeSignificantUnitsCond = {};
        firingRateNegativeNonsignificantUnitsCond = {};
        firingRateNegativeMUAsCond = {};
        firingRateNegativeSignificantMUAsCond = {};
        firingRateNegativeNonsignificantMUAsCond = {};
        
        for iArea = 1:numel(areas)
          firingRatePositiveUnitsCond{iArea} = []; %#ok<*SAGROW>
          firingRatePositiveSignificantUnitsCond{iArea} = [];
          firingRatePositiveNonsignificantUnitsCond{iArea} = [];
          firingRatePositiveMUAsCond{iArea} = [];
          firingRatePositiveSignificantMUAsCond{iArea} = [];
          firingRatePositiveNonsignificantMUAsCond{iArea} = [];
          firingRateNegativeUnitsCond{iArea} = [];
          firingRateNegativeSignificantUnitsCond{iArea} = [];
          firingRateNegativeNonsignificantUnitsCond{iArea} = [];
          firingRateNegativeMUAsCond{iArea} = [];
          firingRateNegativeSignificantMUAsCond{iArea} = [];
          firingRateNegativeNonsignificantMUAsCond{iArea} = [];
        end
        
        firingRatePositiveUnits{iCond} = firingRatePositiveUnitsCond;
        firingRatePositiveSignificantUnits{iCond} = firingRatePositiveSignificantUnitsCond;
        firingRatePositiveNonsignificantUnits{iCond} = firingRatePositiveNonsignificantUnitsCond;
        firingRatePositiveMUAs{iCond} = firingRatePositiveMUAsCond;
        firingRatePositiveSignificantMUAs{iCond} = firingRatePositiveSignificantMUAsCond;
        firingRatePositiveNonsignificantMUAs{iCond} = firingRatePositiveNonsignificantMUAsCond;
        firingRateNegativeUnits{iCond} = firingRateNegativeUnitsCond;
        firingRateNegativeSignificantUnits{iCond} = firingRateNegativeSignificantUnitsCond;
        firingRateNegativeNonsignificantUnits{iCond} = firingRateNegativeNonsignificantUnitsCond;
        firingRateNegativeMUAs{iCond} = firingRateNegativeMUAsCond;
        firingRateNegativeSignificantMUAs{iCond} = firingRateNegativeSignificantMUAsCond;
        firingRateNegativeNonsignificantMUAs{iCond} = firingRateNegativeNonsignificantMUAsCond;
      end
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
      
      for iAreaPlusAll = area % Loop through the main and pooled areas
        if numel(conditions) > 1
          condLoop = [iCond numel(conditions)];
        else
          condLoop = iCond;
        end
        for iCondPlusAll = condLoop % Loop through the main and pooled conditions
          
          if ~isempty(dbStructPositive)
            srData = dbStructPositive.conf.samplingParams.srData;
            spk = dbStructPositive.shankData.shank1.spk;
            unitMetadataPositive = dbStructPositive.shankData.shank1.unitMetadata;
            if ~isempty(unitMetadataPositive)
              significantUnitInd = unitMetadataPositive(:,11) < alpha;
              firingRatePositiveUnits{iCondPlusAll}{iAreaPlusAll} =...
                [firingRatePositiveUnits{iCondPlusAll}{iAreaPlusAll};...
                (sum(spk,2).*srData)./size(spk,2)];
              firingRatePositiveSignificantUnits{iCondPlusAll}{iAreaPlusAll} =...
                [firingRatePositiveSignificantUnits{iCondPlusAll}{iAreaPlusAll};...
                (sum(spk(significantUnitInd,:),2).*srData)./size(spk(significantUnitInd,:),2)];
              firingRatePositiveNonsignificantUnits{iCondPlusAll}{iAreaPlusAll} =...
                [firingRatePositiveNonsignificantUnits{iCondPlusAll}{iAreaPlusAll};...
                (sum(spk(~significantUnitInd,:),2).*srData)./size(spk(~significantUnitInd,:),2)];
            end
            spkDB = dbStructPositive.popData.spkDB;
            if ~isempty(spkDB)
              significantMUAInd = dbStructPositive.popData.pvalSpearman < alpha;
              firingRatePositiveMUAs{iCondPlusAll}{iAreaPlusAll} =...
                [firingRatePositiveMUAs{iCondPlusAll}{iAreaPlusAll};...
                (sum(spkDB,2).*srData)./size(spkDB,2)];
              firingRatePositiveSignificantMUAs{iCondPlusAll}{iAreaPlusAll} =...
                [firingRatePositiveSignificantMUAs{iCondPlusAll}{iAreaPlusAll};...
                (sum(spkDB(significantMUAInd,:),2).*srData)./size(spkDB(significantMUAInd,:),2)];
              firingRatePositiveNonsignificantMUAs{iCondPlusAll}{iAreaPlusAll} =...
                [firingRatePositiveNonsignificantMUAs{iCondPlusAll}{iAreaPlusAll};...
                (sum(spkDB(~significantMUAInd,:),2).*srData)./size(spkDB(~significantMUAInd,:),2)];
            end
          end
          
          if ~isempty(dbStructNegative)
            srData = dbStructNegative.conf.samplingParams.srData;
            spk = dbStructNegative.shankData.shank1.spk;
            unitMetadataNegative = dbStructNegative.shankData.shank1.unitMetadata;
            if ~isempty(unitMetadataNegative)
              significantUnitInd = unitMetadataNegative(:,11) < alpha;
              firingRateNegativeUnits{iCondPlusAll}{iAreaPlusAll} =...
                [firingRateNegativeUnits{iCondPlusAll}{iAreaPlusAll};...
                (sum(spk,2).*srData)./size(spk,2)];
              firingRateNegativeSignificantUnits{iCondPlusAll}{iAreaPlusAll} =...
                [firingRateNegativeSignificantUnits{iCondPlusAll}{iAreaPlusAll};...
                (sum(spk(significantUnitInd,:),2).*srData)./size(spk(significantUnitInd,:),2)];
              firingRateNegativeNonsignificantUnits{iCondPlusAll}{iAreaPlusAll} =...
                [firingRateNegativeNonsignificantUnits{iCondPlusAll}{iAreaPlusAll};...
                (sum(spk(~significantUnitInd,:),2).*srData)./size(spk(~significantUnitInd,:),2)];
            end
            spkDB = dbStructNegative.popData.spkDB;
            if ~isempty(spkDB)
              significantMUAInd = dbStructNegative.popData.pvalSpearman < alpha;
              firingRateNegativeMUAs{iCondPlusAll}{iAreaPlusAll} =...
                [firingRateNegativeMUAs{iCondPlusAll}{iAreaPlusAll};...
                (sum(spkDB,2).*srData)./size(spkDB,2)];
              firingRateNegativeSignificantMUAs{iCondPlusAll}{iAreaPlusAll} =...
                [firingRateNegativeSignificantMUAs{iCondPlusAll}{iAreaPlusAll};...
                (sum(spkDB(significantMUAInd,:),2).*srData)./size(spkDB(significantMUAInd,:),2)];
              firingRateNegativeNonsignificantMUAs{iCondPlusAll}{iAreaPlusAll} =...
                [firingRateNegativeNonsignificantMUAs{iCondPlusAll}{iAreaPlusAll};...
                (sum(spkDB(~significantMUAInd,:),2).*srData)./size(spkDB(~significantMUAInd,:),2)];
            end
          end
          
        end
      end
    end
  end
  
  % Calculate mean firing rates and 95% confidence limits on mean rates
  for iCond = 1:numel(conditions)
    for iArea = 1:numel(areas)
      [firingRatePositiveUnitsCI95{iCond}{iArea}, ~, ~, firingRatePositiveUnitsMean{iCond}{iArea}] = calc95CI(firingRatePositiveUnits{iCond}{iArea});
      [firingRatePositiveSignificantUnitsCI95{iCond}{iArea}, ~, ~, firingRatePositiveSignificantUnitsMean{iCond}{iArea}] = calc95CI(firingRatePositiveSignificantUnits{iCond}{iArea});
      [firingRatePositiveNonsignificantUnitsCI95{iCond}{iArea}, ~, ~, firingRatePositiveNonsignificantUnitsMean{iCond}{iArea}] = calc95CI(firingRatePositiveNonsignificantUnits{iCond}{iArea});
      [firingRateNegativeUnitsCI95{iCond}{iArea}, ~, ~, firingRateNegativeUnitsMean{iCond}{iArea}] = calc95CI(firingRateNegativeUnits{iCond}{iArea});
      [firingRateNegativeSignificantUnitsCI95{iCond}{iArea}, ~, ~, firingRateNegativeSignificantUnitsMean{iCond}{iArea}] = calc95CI(firingRateNegativeSignificantUnits{iCond}{iArea});
      [firingRateNegativeNonsignificantUnitsCI95{iCond}{iArea}, ~, ~, firingRateNegativeNonsignificantUnitsMean{iCond}{iArea}] = calc95CI(firingRateNegativeNonsignificantUnits{iCond}{iArea});
      [firingRatePositiveMUAsCI95{iCond}{iArea}, ~, ~, firingRatePositiveMUAsMean{iCond}{iArea}] = calc95CI(firingRatePositiveMUAs{iCond}{iArea});
      [firingRatePositiveSignificantMUAsCI95{iCond}{iArea}, ~, ~, firingRatePositiveSignificantMUAsMean{iCond}{iArea}] = calc95CI(firingRatePositiveSignificantMUAs{iCond}{iArea});
      [firingRatePositiveNonsignificantMUAsCI95{iCond}{iArea}, ~, ~, firingRatePositiveNonsignificantMUAsMean{iCond}{iArea}] = calc95CI(firingRatePositiveNonsignificantMUAs{iCond}{iArea});
      [firingRateNegativeMUAsCI95{iCond}{iArea}, ~, ~, firingRateNegativeMUAsMean{iCond}{iArea}] = calc95CI(firingRateNegativeMUAs{iCond}{iArea});
      [firingRateNegativeSignificantMUAsCI95{iCond}{iArea}, ~, ~, firingRateNegativeSignificantMUAsMean{iCond}{iArea}] = calc95CI(firingRateNegativeSignificantMUAs{iCond}{iArea});
      [firingRateNegativeNonsignificantMUAsCI95{iCond}{iArea}, ~, ~, firingRateNegativeNonsignificantMUAsMean{iCond}{iArea}] = calc95CI(firingRateNegativeNonsignificantMUAs{iCond}{iArea});
    end
  end
  
  % Save processed data
  data.firingRatePositiveUnits = firingRatePositiveUnits;
  data.firingRatePositiveSignificantUnits = firingRatePositiveSignificantUnits;
  data.firingRatePositiveNonsignificantUnits = firingRatePositiveNonsignificantUnits;
  data.firingRatePositiveMUAs = firingRatePositiveMUAs;
  data.firingRatePositiveSignificantMUAs = firingRatePositiveSignificantMUAs;
  data.firingRatePositiveNonsignificantMUAs = firingRatePositiveNonsignificantMUAs;
  data.firingRateNegativeUnits = firingRateNegativeUnits;
  data.firingRateNegativeSignificantUnits = firingRateNegativeSignificantUnits;
  data.firingRateNegativeNonsignificantUnits = firingRateNegativeNonsignificantUnits;
  data.firingRateNegativeMUAs = firingRateNegativeMUAs;
  data.firingRateNegativeSignificantMUAs = firingRateNegativeSignificantMUAs;
  data.firingRateNegativeNonsignificantMUAs = firingRateNegativeNonsignificantMUAs;
  data.firingRatePositiveUnitsCI95 = firingRatePositiveUnitsCI95;
  data.firingRatePositiveUnitsMean = firingRatePositiveUnitsMean;
  data.firingRatePositiveSignificantUnitsCI95 = firingRatePositiveSignificantUnitsCI95;
  data.firingRatePositiveSignificantUnitsMean = firingRatePositiveSignificantUnitsMean;
  data.firingRatePositiveNonsignificantUnitsCI95 = firingRatePositiveNonsignificantUnitsCI95;
  data.firingRatePositiveNonsignificantUnitsMean = firingRatePositiveNonsignificantUnitsMean;
  data.firingRateNegativeUnitsCI95 = firingRateNegativeUnitsCI95;
  data.firingRateNegativeUnitsMean = firingRateNegativeUnitsMean;
  data.firingRateNegativeSignificantUnitsCI95 = firingRateNegativeSignificantUnitsCI95;
  data.firingRateNegativeSignificantUnitsMean = firingRateNegativeSignificantUnitsMean;
  data.firingRateNegativeNonsignificantUnitsCI95 = firingRateNegativeNonsignificantUnitsCI95;
  data.firingRateNegativeNonsignificantUnitsMean = firingRateNegativeNonsignificantUnitsMean;
  data.firingRatePositiveMUAsCI95 = firingRatePositiveMUAsCI95;
  data.firingRatePositiveMUAsMean = firingRatePositiveMUAsMean;
  data.firingRatePositiveSignificantMUAsCI95 = firingRatePositiveSignificantMUAsCI95;
  data.firingRatePositiveSignificantMUAsMean = firingRatePositiveSignificantMUAsMean;
  data.firingRatePositiveNonsignificantMUAsCI95 = firingRatePositiveNonsignificantMUAsCI95;
  data.firingRatePositiveNonsignificantMUAsMean = firingRatePositiveNonsignificantMUAsMean;
  data.firingRateNegativeMUAsCI95 = firingRateNegativeMUAsCI95;
  data.firingRateNegativeMUAsMean = firingRateNegativeMUAsMean;
  data.firingRateNegativeSignificantMUAsCI95 = firingRateNegativeSignificantMUAsCI95;
  data.firingRateNegativeSignificantMUAsMean = firingRateNegativeSignificantMUAsMean;
  data.firingRateNegativeNonsignificantMUAsCI95 = firingRateNegativeNonsignificantMUAsCI95;
  data.firingRateNegativeNonsignificantMUAsMean = firingRateNegativeNonsignificantMUAsMean;
  save([mainFolder filesep 'pupilCorrFiringRates.mat'], 'repository','alpha','animals','conditions','areas','data', '-v7.3');
else
  if strcmp(repository, 'allensdk')
    load([mainFolder filesep 'pupilCorrFiringRates.mat']); %#ok<*UNRCH>
  elseif strcmp(repository, 'uol')
    load([mainFolder filesep 'all' filesep 'pupilCorrFiringRates.mat']); %#ok<*UNRCH>
  end
end


if strcmp(repository, 'uol')
  %% Independent-samples t-tests for units (uol):
  yLim = [];
  %yLim = [0 10];
  fH = barPlotUnits(data, 'VB', 'lS1', 'lRSC', 'CA', yLim, 'log');
  
elseif strcmp(repository, 'allensdk')
  %% Independent-samples t-tests for units (allensdk):
  %yLim = [];
  yLim = [0 14];
  fH = barPlotUnits(data, 'VB', 'LGN', 'V1', 'CA', yLim, 'log');
  
end

set(fH(1), 'Name','Firing rate of units positively and negatively correlated to pupil size in different brain areas');
filename = [mainFolder filesep 'pupilCorrFiringRatesUnits_logy'];
savefig(fH(1), filename, 'compact');
print(fH(1), [filename '.png'],'-dpng','-r300');
%close all

figure(fH(2));
yLim = ylim;
if strcmp(repository, 'allensdk')
  ylim([yLim(1) 13]);
else
  ylim([yLim(1) 10]);
end
%xLabel = get(gca, 'xlabel');
yLabel = get(gca, 'ylabel');
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 40, 4/3, 2, [0.005 0], 'out',...
  'on', 'k', {}, xlim, xticks,...
  'on', 'k', {yLabel.String{1}}, ylim, yticks); %#ok<CCAT1>
set(fH(2), 'Name','Firing rate of units positively and negatively correlated to pupil size in different brain areas');

figFileName = [mainFolder filesep 'pupilCorrFiringRatesUnitsRedux_logy'];
figSize = 15;
label = [4.2 1.8];
margin = [0 0.8];
width = ((100-23)/(147-119))*figSize-label(1)-margin(1);
height = figSize-label(2)-margin(2);
paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
savefig(gcf, [figFileName '.fig'], 'compact');
exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);
exportFig(gcf, [figFileName '.eps'],'-depsc','-r1200', paperSize);
%close all


if strcmp(repository, 'uol')
  %% Independent-samples t-tests for units (uol):
  yLim = [];
  %yLim = [0 10];
  [fH, scatterGroups, areaGroups, colourGroups, colorCode, corrGroups] = barPlotUnits(data, 'VB', 'lS1', 'lRSC', 'CA', yLim);
  
elseif strcmp(repository, 'allensdk')
  %% Independent-samples t-tests for units (allensdk):
  %yLim = [];
  yLim = [0 14];
  [fH, scatterGroups, areaGroups, colourGroups, colorCode, corrGroups] = barPlotUnits(data, 'VB', 'LGN', 'V1', 'CA', yLim);
  
end

set(fH(1), 'Name','Firing rate of units positively and negatively correlated to pupil size in different brain areas');
filename = [mainFolder filesep 'pupilCorrFiringRatesUnits'];
savefig(fH(1), filename, 'compact');
print(fH(1), [filename '.png'],'-dpng','-r300');
%close all

set(fH(2), 'Name','Firing rate of units positively and negatively correlated to pupil size in different brain areas');
filename = [mainFolder filesep 'pupilCorrFiringRatesUnitsRedux'];
savefig(fH(2), filename, 'compact');
print(fH(2), [filename '.png'],'-dpng','-r300');
%close all


%% ANOVA for units:
runAnova(scatterGroups, corrGroups, areaGroups, filename)


if strcmp(repository, 'uol')
  %% Independent-samples t-tests for MUAs (uol):
  fH = barPlotMUAs(data, 'VB', 'lS1', 'lRSC', 'CA', yLim, 'log');
  
elseif strcmp(repository, 'allensdk')
  %% Independent-samples t-tests for MUAs (allensdk):
  fH = barPlotMUAs(data, 'VB', 'LGN', 'V1', 'CA', yLim, 'log');
  
end

set(fH(1), 'Name','Firing rate of units and MUAs positively and negatively correlated to pupil size in different brain areas');
filename = [mainFolder filesep 'pupilCorrFiringRatesMUAs_logy'];
savefig(fH(1), filename, 'compact');
print(fH(1), [filename '.png'],'-dpng','-r300');
%close all

figure(fH(2));
yLim = ylim;
if strcmp(repository, 'allensdk')
  ylim([yLim(1) 10]);
else
  ylim([yLim(1) 8]);
end
xLabel = get(gca, 'xlabel');
yLabel = get(gca, 'ylabel');
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 40, 4/3, 2, [0.005 0], 'out',...
  'on', 'k', {}, xlim, xticks,...
  'on', 'k', {yLabel.String{1}}, ylim, yticks); %#ok<CCAT1>
set(fH(2), 'Name','Firing rate of units and MUAs positively and negatively correlated to pupil size in different brain areas');

figFileName = [mainFolder filesep 'pupilCorrFiringRatesMUAsRedux_logy'];
figSize = 15;
label = [4.2 1.8];
margin = [0 0.8];
width = ((100-23)/(147-119))*figSize-label(1)-margin(1);
height = figSize-label(2)-margin(2);
paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
savefig(gcf, [figFileName '.fig'], 'compact');
exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);
exportFig(gcf, [figFileName '.eps'],'-depsc','-r1200', paperSize);
%close all


if strcmp(repository, 'uol')
  %% Independent-samples t-tests for MUAs (uol):
  [fH, scatterGroups, areaGroups, colourGroups, colorCode, corrGroups] = barPlotMUAs(data, 'VB', 'lS1', 'lRSC', 'CA', yLim);
  
elseif strcmp(repository, 'allensdk')
  %% Independent-samples t-tests for MUAs (allensdk):
  [fH, scatterGroups, areaGroups, colourGroups, colorCode, corrGroups] = barPlotMUAs(data, 'VB', 'LGN', 'V1', 'CA', yLim);
  
end

set(fH(1), 'Name','Firing rate of units and MUAs positively and negatively correlated to pupil size in different brain areas');
filename = [mainFolder filesep 'pupilCorrFiringRatesMUAs'];
savefig(fH(1), filename, 'compact');
print(fH(1), [filename '.png'],'-dpng','-r300');
%close all

set(fH(2), 'Name','Firing rate of units and MUAs positively and negatively correlated to pupil size in different brain areas');
filename = [mainFolder filesep 'pupilCorrFiringRatesMUAsRedux'];
savefig(fH(2), filename, 'compact');
print(fH(2), [filename '.png'],'-dpng','-r300');
%close all


%% ANOVA for units:
runAnova(scatterGroups, corrGroups, areaGroups, filename)



%% Local functions
function [dataMean, dataScatter, areaGroup, colourGroup, corrGroup] = barGroup(areaName, colourNumber, barPosition, data)

minDataPoints = 10;
areaCode = determineArea(areaName);
colour = matlabColours(colourNumber);
dataMean = [];
if numel(full(data.firingRatePositiveUnits{1}{areaCode})) >= minDataPoints
  bar(1+barPosition, data.firingRatePositiveUnitsMean{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour, 'FaceAlpha',0.2);
  dataMean = [dataMean full(data.firingRatePositiveUnitsMean{1}{areaCode})];
else
  bar(1+barPosition, 0, 'FaceColor',colour, 'EdgeColor',colour, 'FaceAlpha',0.2);
  dataMean = [dataMean 0];
end
if numel(full(data.firingRatePositiveSignificantUnits{1}{areaCode})) >= minDataPoints
  bar(2+barPosition, data.firingRatePositiveSignificantUnitsMean{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
  dataMean = [dataMean full(data.firingRatePositiveSignificantUnitsMean{1}{areaCode})];
else
  bar(2+barPosition, 0, 'FaceColor',colour, 'EdgeColor',colour);
  dataMean = [dataMean 0];
end
if numel(full(data.firingRateNegativeUnits{1}{areaCode})) >= minDataPoints
  bar(3+barPosition, data.firingRateNegativeUnitsMean{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3), 'FaceAlpha',0.2);
  dataMean = [dataMean full(data.firingRateNegativeUnitsMean{1}{areaCode})];
else
  bar(3+barPosition, 0, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3), 'FaceAlpha',0.2);
  dataMean = [dataMean 0];
end
if numel(full(data.firingRateNegativeSignificantUnits{1}{areaCode})) >= minDataPoints
  bar(4+barPosition, data.firingRateNegativeSignificantUnitsMean{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [dataMean full(data.firingRateNegativeSignificantUnitsMean{1}{areaCode})];
else
  bar(4+barPosition, 0, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [dataMean 0];
end
dataScatter = {full(data.firingRatePositiveUnits{1}{areaCode}),...
  full(data.firingRatePositiveSignificantUnits{1}{areaCode}),...
  full(data.firingRateNegativeUnits{1}{areaCode}),...
  full(data.firingRateNegativeSignificantUnits{1}{areaCode})};
areaGroup = {areaName, areaName, areaName, areaName};
colourGroup = [colour; colour; colour; colour];
corrGroup = {'positive','positive','negative','negative'};
end


function [dataMean, dataScatter, areaGroup, colourGroup, corrGroup] = barGroupRedux(areaName, colourNumber, barPosition, data)

minDataPoints = 10;
areaCode = determineArea(areaName);
colour = matlabColours(colourNumber);
dataMean = [];
if numel(full(data.firingRatePositiveUnits{1}{areaCode})) >= minDataPoints
  bar(1+barPosition, data.firingRatePositiveUnitsMean{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
  dataMean = [dataMean full(data.firingRatePositiveUnitsMean{1}{areaCode})];
else
  bar(1+barPosition, 0, 'FaceColor',colour, 'EdgeColor',colour);
  dataMean = [dataMean 0];
end
if numel(full(data.firingRateNegativeUnits{1}{areaCode})) >= minDataPoints
  bar(2+barPosition, data.firingRateNegativeUnitsMean{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [dataMean full(data.firingRateNegativeUnitsMean{1}{areaCode})];
else
  bar(2+barPosition, 0, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [dataMean 0];
end
dataScatter = {full(data.firingRatePositiveUnits{1}{areaCode}),...
  full(data.firingRateNegativeUnits{1}{areaCode})};
areaGroup = {areaName, areaName};
colourGroup = [colour; colour];
corrGroup = {'positive','negative'};
end


function [dataMean, dataScatter, areaGroup, colourGroup, corrGroup] = barGroupMUAs(areaName, colourNumber, barPosition, data)

minDataPoints = 10;
areaCode = determineArea(areaName);
colour = matlabColours(colourNumber);
dataMean = [];
if numel(full(data.firingRatePositiveMUAs{1}{areaCode})) >= minDataPoints
  bar(1+barPosition, data.firingRatePositiveMUAsMean{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour, 'FaceAlpha',0.2);
  dataMean = [dataMean full(data.firingRatePositiveMUAsMean{1}{areaCode})];
else
  bar(1+barPosition, 0, 'FaceColor',colour, 'EdgeColor',colour, 'FaceAlpha',0.2);
  dataMean = [dataMean 0];
end
if numel(full(data.firingRatePositiveSignificantMUAs{1}{areaCode})) >= minDataPoints
  bar(2+barPosition, data.firingRatePositiveSignificantMUAsMean{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
  dataMean = [dataMean full(data.firingRatePositiveSignificantMUAsMean{1}{areaCode})];
else
  bar(2+barPosition, 0, 'FaceColor',colour, 'EdgeColor',colour);
  dataMean = [dataMean 0];
end
if numel(full(data.firingRateNegativeMUAs{1}{areaCode})) >= minDataPoints
  bar(3+barPosition, data.firingRateNegativeMUAsMean{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3), 'FaceAlpha',0.2);
  dataMean = [dataMean full(data.firingRateNegativeMUAsMean{1}{areaCode})];
else
  bar(3+barPosition, 0, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3), 'FaceAlpha',0.2);
  dataMean = [dataMean 0];
end
if numel(full(data.firingRateNegativeSignificantMUAs{1}{areaCode})) >= minDataPoints
  bar(4+barPosition, data.firingRateNegativeSignificantMUAsMean{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [dataMean full(data.firingRateNegativeSignificantMUAsMean{1}{areaCode})];
else
  bar(4+barPosition, 0, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [dataMean 0];
end
dataScatter = {full(data.firingRatePositiveMUAs{1}{areaCode}),...
  full(data.firingRatePositiveSignificantMUAs{1}{areaCode}),...
  full(data.firingRateNegativeMUAs{1}{areaCode}),...
  full(data.firingRateNegativeSignificantMUAs{1}{areaCode})};
areaGroup = {areaName, areaName, areaName, areaName};
colourGroup = [colour; colour; colour; colour];
corrGroup = {'positive','positive','negative','negative'};
end


function [dataMean, dataScatter, areaGroup, colourGroup, corrGroup] = barGroupMUAsRedux(areaName, colourNumber, barPosition, data)

minDataPoints = 10;
areaCode = determineArea(areaName);
colour = matlabColours(colourNumber);
dataMean = [];
if numel(full(data.firingRatePositiveMUAs{1}{areaCode})) >= minDataPoints
  bar(1+barPosition, data.firingRatePositiveMUAsMean{1}{areaCode}, 'FaceColor',colour, 'EdgeColor',colour);
  dataMean = [dataMean full(data.firingRatePositiveMUAsMean{1}{areaCode})];
else
  bar(1+barPosition, 0, 'FaceColor',colour, 'EdgeColor',colour);
  dataMean = [dataMean 0];
end
if numel(full(data.firingRateNegativeMUAs{1}{areaCode})) >= minDataPoints
  bar(2+barPosition, data.firingRateNegativeMUAsMean{1}{areaCode}, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [dataMean full(data.firingRateNegativeMUAsMean{1}{areaCode})];
else
  bar(2+barPosition, 0, 'FaceColor',colour*(2/3), 'EdgeColor',colour*(2/3));
  dataMean = [dataMean 0];
end
dataScatter = {full(data.firingRatePositiveMUAs{1}{areaCode}),...
  full(data.firingRateNegativeMUAs{1}{areaCode})};
areaGroup = {areaName, areaName};
colourGroup = [colour; colour];
corrGroup = {'positive','negative'};
end


function [p, pSig] = ttestGroup(areaName, data)

areaCode = determineArea(areaName);
[~,p] = ttest2(full(data.firingRatePositiveUnits{1}{areaCode}),...
  full(data.firingRateNegativeUnits{1}{areaCode}));
[~,pSig] = ttest2(full(data.firingRatePositiveSignificantUnits{1}{areaCode}),...
  full(data.firingRateNegativeSignificantUnits{1}{areaCode}));
end


function [p, pSig] = ttestGroupMUAs(areaName, data)

areaCode = determineArea(areaName);
[~,p] = ttest2(full(data.firingRatePositiveMUAs{1}{areaCode}),...
  full(data.firingRateNegativeMUAs{1}{areaCode}));
[~,pSig] = ttest2(full(data.firingRatePositiveSignificantMUAs{1}{areaCode}),...
  full(data.firingRateNegativeSignificantMUAs{1}{areaCode}));
end


function [fH, scatterGroups, areaGroups, colourGroups, colorCode, corrGroups] = barPlotUnits(data, area1, area2, area3, area4, yLim, scaleType)

if nargin < 7
  scaleType = 'normal';
end

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
err = [data.firingRatePositiveUnitsCI95{1}{determineArea(area1)} data.firingRatePositiveSignificantUnitsCI95{1}{determineArea(area1)}...
  data.firingRateNegativeUnitsCI95{1}{determineArea(area1)} data.firingRateNegativeSignificantUnitsCI95{1}{determineArea(area1)}...
  data.firingRatePositiveUnitsCI95{1}{determineArea(area2)} data.firingRatePositiveSignificantUnitsCI95{1}{determineArea(area2)}...
  data.firingRateNegativeUnitsCI95{1}{determineArea(area2)} data.firingRateNegativeSignificantUnitsCI95{1}{determineArea(area2)}...
  data.firingRatePositiveUnitsCI95{1}{determineArea(area3)} data.firingRatePositiveSignificantUnitsCI95{1}{determineArea(area3)}...
  data.firingRateNegativeUnitsCI95{1}{determineArea(area3)} data.firingRateNegativeSignificantUnitsCI95{1}{determineArea(area3)}...
  data.firingRatePositiveUnitsCI95{1}{determineArea(area4)} data.firingRatePositiveSignificantUnitsCI95{1}{determineArea(area4)}...
  data.firingRateNegativeUnitsCI95{1}{determineArea(area4)} data.firingRateNegativeSignificantUnitsCI95{1}{determineArea(area4)}];
err(:, ~dataMean) = 0;

er = errorbar(bars,dataMean,err(2,:),err(1,:));
er.Color = [0 0 0];
er.LineStyle = 'none';


% Scatter
scatterGroups = [scatterGroup1 scatterGroup2 scatterGroup3 scatterGroup4];
areaGroups = [areaGroup1 areaGroup2 areaGroup3 areaGroup4];
colourGroups = [colourGroup1; colourGroup2; colourGroup3; colourGroup4;];
colorCode = [30 30 30]./255;
corrGroups = [corrGroup1 corrGroup2 corrGroup3 corrGroup4];
% for iBar = 1:numel(bars)
%   scatter(bars(iBar)*ones(size(scatterGroups{iBar}))', scatterGroups{iBar}',...
%     'MarkerEdgeColor',colorCode, 'jitter','on'); %colourGroups(iBar,:));
% end


% Graph adjustments
xTickPos = [gaps 20]-2.5;
if nargin < 6 || isempty(yLim)
  yLim = ylim;
end
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Arial', 25, 4/3, 2, [0.003 0.003], 'out', 'on', 'k', {}, [],...
  xTickPos, 'on', 'k', {'Firing rate (APs/s)'}, yLim, yticks);
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
  text(iBar-0.075,-0.01*yLim(2), '+', 'FontSize',20)
  text(iBar+2-0.075,-0.01*yLim(2), '-', 'FontSize',20)
end

if strcmp(scaleType, 'log')
  set(gca,'Yscale','log')
end


%% Draw bar graphs without bars for significant units only
fH2 = figProperties('Bar plot for units', 'normalized', [0, .005, .97, .90], 'w', 'on');
hold on


% Data bars
gap = 1;
nBarsPerGroup = 2;
barGroupRedux1 = barGroupRedux(area1, 1, 0*(nBarsPerGroup+gap), data);
if strcmp(area2,'LGN') || strcmp(area2,'lLGN') || strcmp(area2,'rLGN')
  barGroupRedux2 = barGroupRedux(area2, 2, 1*(nBarsPerGroup+gap), data);
elseif strcmp(area2,'S1') || strcmp(area2,'lS1') || strcmp(area2,'rS1')
  barGroupRedux2 = barGroupRedux(area2, 13, 1*(nBarsPerGroup+gap), data);
end
if strcmp(area3,'V1') || strcmp(area3,'lV1') || strcmp(area3,'rV1')
  barGroupRedux3 = barGroupRedux(area3, 11, 2*(nBarsPerGroup+gap), data);
elseif strcmp(area3,'RSC') || strcmp(area3,'lRSC') || strcmp(area3,'rRSC')
  barGroupRedux3 = barGroupRedux(area3, 3, 2*(nBarsPerGroup+gap), data);
end
barGroupRedux4 = barGroupRedux(area4, 4, 3*(nBarsPerGroup+gap), data);


% Error bars
gaps = 3:3:9;
bars = sort([1 1+gaps 2 2+gaps]);
dataMean = [barGroupRedux1 barGroupRedux2 barGroupRedux3 barGroupRedux4];
err = [data.firingRatePositiveUnitsCI95{1}{determineArea(area1)}...
  data.firingRateNegativeUnitsCI95{1}{determineArea(area1)}...
  data.firingRatePositiveUnitsCI95{1}{determineArea(area2)}...
  data.firingRateNegativeUnitsCI95{1}{determineArea(area2)}...
  data.firingRatePositiveUnitsCI95{1}{determineArea(area3)}...
  data.firingRateNegativeUnitsCI95{1}{determineArea(area3)}...
  data.firingRatePositiveUnitsCI95{1}{determineArea(area4)}...
  data.firingRateNegativeUnitsCI95{1}{determineArea(area4)}];
err(:, ~dataMean) = 0;

er = errorbar(bars,dataMean,err(2,:),err(1,:));
er.Color = [0 0 0];
er.LineStyle = 'none';


% Graph adjustments
xTickPos = [gaps 12]-1.5;
if nargin < 6 || isempty(yLim)
  yLim = ylim;
end
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Arial', 25, 4/3, 2, [0.003 0.003], 'out', 'on', 'k', {}, [],...
  xTickPos, 'on', 'k', {'Firing rate (APs/s)'}, yLim, yticks);
ax1.XTickLabel = xTickLabel;

for iBar = bars([1 3 5 7])
  text(iBar-0.075,-0.01*yLim(2), '+', 'FontSize',20)
  text(iBar+1-0.075,-0.01*yLim(2), '-', 'FontSize',20)
end

if strcmp(scaleType, 'log')
  set(gca,'Yscale','log')
end

fH = [fH1 fH2];
end


function [fH, scatterGroups, areaGroups, colourGroups, colorCode, corrGroups] = barPlotMUAs(data, area1, area2, area3, area4, yLim, scaleType)

if nargin < 7
  scaleType = 'normal';
end

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
err = [data.firingRatePositiveMUAsCI95{1}{determineArea(area1)} data.firingRatePositiveSignificantMUAsCI95{1}{determineArea(area1)}...
  data.firingRateNegativeMUAsCI95{1}{determineArea(area1)} data.firingRateNegativeSignificantMUAsCI95{1}{determineArea(area1)}...
  data.firingRatePositiveMUAsCI95{1}{determineArea(area2)} data.firingRatePositiveSignificantMUAsCI95{1}{determineArea(area2)}...
  data.firingRateNegativeMUAsCI95{1}{determineArea(area2)} data.firingRateNegativeSignificantMUAsCI95{1}{determineArea(area2)}...
  data.firingRatePositiveMUAsCI95{1}{determineArea(area3)} data.firingRatePositiveSignificantMUAsCI95{1}{determineArea(area3)}...
  data.firingRateNegativeMUAsCI95{1}{determineArea(area3)} data.firingRateNegativeSignificantMUAsCI95{1}{determineArea(area3)}...
  data.firingRatePositiveMUAsCI95{1}{determineArea(area4)} data.firingRatePositiveSignificantMUAsCI95{1}{determineArea(area4)}...
  data.firingRateNegativeMUAsCI95{1}{determineArea(area4)} data.firingRateNegativeSignificantMUAsCI95{1}{determineArea(area4)}];
err(:, ~dataMean) = 0;

er = errorbar(bars,dataMean,err(2,:),err(1,:));
er.Color = [0 0 0];
er.LineStyle = 'none';


% Scatter
scatterGroups = [scatterGroup1 scatterGroup2 scatterGroup3 scatterGroup4];
areaGroups = [areaGroup1 areaGroup2 areaGroup3 areaGroup4];
colourGroups = [colourGroup1; colourGroup2; colourGroup3; colourGroup4;];
colorCode = [30 30 30]./255;
corrGroups = [corrGroup1 corrGroup2 corrGroup3 corrGroup4];
% for iBar = 1:numel(bars)
%   scatter(bars(iBar)*ones(size(scatterGroups{iBar}))', scatterGroups{iBar}',...
%     'MarkerEdgeColor',colorCode, 'jitter','on'); %colourGroups(iBar,:));
% end


% Graph adjustments
xTickPos = [gaps 20]-2.5;
if nargin < 6 || isempty(yLim)
  yLim = ylim;
end
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Arial', 25, 4/3, 2, [0.003 0.003], 'out', 'on', 'k', {}, [],...
  xTickPos, 'on', 'k', {'Firing rate (APs/s)'}, yLim, yticks);
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
  text(iBar-0.075,-0.01*yLim(2), '+', 'FontSize',20)
  text(iBar+2-0.075,-0.01*yLim(2), '-', 'FontSize',20)
end

if strcmp(scaleType, 'log')
  set(gca,'Yscale','log')
end


%% Draw bar graphs without bars for significant units+MUAs only
fH2 = figProperties('Bar plot for units', 'normalized', [0, .005, .97, .90], 'w', 'on');
hold on


% Data bars
gap = 1;
nBarsPerGroup = 2;

barGroupRedux1 = barGroupMUAsRedux(area1, 1, 0*(nBarsPerGroup+gap), data);
if strcmp(area2,'LGN') || strcmp(area2,'lLGN') || strcmp(area2,'rLGN')
  barGroupRedux2 = barGroupMUAsRedux(area2, 2, 1*(nBarsPerGroup+gap), data);
elseif strcmp(area2,'S1') || strcmp(area2,'lS1') || strcmp(area2,'rS1')
  barGroupRedux2 = barGroupMUAsRedux(area2, 13, 1*(nBarsPerGroup+gap), data);
end
if strcmp(area3,'V1') || strcmp(area3,'lV1') || strcmp(area3,'rV1')
  barGroupRedux3 = barGroupMUAsRedux(area3, 11, 2*(nBarsPerGroup+gap), data);
elseif strcmp(area3,'RSC') || strcmp(area3,'lRSC') || strcmp(area3,'rRSC')
  barGroupRedux3 = barGroupMUAsRedux(area3, 3, 2*(nBarsPerGroup+gap), data);
end
barGroupRedux4 = barGroupMUAsRedux(area4, 4, 3*(nBarsPerGroup+gap), data);


% Error bars
gaps = 3:3:9;
bars = sort([1 1+gaps 2 2+gaps]);
dataMean = [barGroupRedux1 barGroupRedux2 barGroupRedux3 barGroupRedux4];
err = [data.firingRatePositiveUnitsCI95{1}{determineArea(area1)}...
  data.firingRateNegativeUnitsCI95{1}{determineArea(area1)}...
  data.firingRatePositiveUnitsCI95{1}{determineArea(area2)}...
  data.firingRateNegativeUnitsCI95{1}{determineArea(area2)}...
  data.firingRatePositiveUnitsCI95{1}{determineArea(area3)}...
  data.firingRateNegativeUnitsCI95{1}{determineArea(area3)}...
  data.firingRatePositiveUnitsCI95{1}{determineArea(area4)}...
  data.firingRateNegativeUnitsCI95{1}{determineArea(area4)}];
err(:, ~dataMean) = 0;

er = errorbar(bars,dataMean,err(2,:),err(1,:));
er.Color = [0 0 0];
er.LineStyle = 'none';


% Graph adjustments
xTickPos = [gaps 12]-1.5;
if nargin < 6 || isempty(yLim)
  yLim = ylim;
end
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Arial', 25, 4/3, 2, [0.003 0.003], 'out', 'on', 'k', {}, [],...
  xTickPos, 'on', 'k', {'Firing rate (APs/s)'}, yLim, yticks);
ax1.XTickLabel = xTickLabel;

for iBar = bars([1 3 5 7])
  text(iBar-0.075,-0.01*yLim(2), '+', 'FontSize',20)
  text(iBar+1-0.075,-0.01*yLim(2), '-', 'FontSize',20)
end

if strcmp(scaleType, 'log')
  set(gca,'Yscale','log')
end

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