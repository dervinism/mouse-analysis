% Run this script to produce PCA summary figures.

fclose all;
close all
clear
clc



%% INITIALISE PARAMETERS
dataDir = 'R:\CSN\Shared\Dynamics\Data\mouse_analysis_data';
lists



%% DISPLAY PCA SUMMARY FIGURES FOR ALL ANIMALS AND ALL RECORDINGS
% INITIALISE FIGURES
fH_S1 = initFigsSeries(conditions);
fH_VB = initFigsSeries(conditions);
fH_Th = initFigsSeries(conditions);
fH_DG = initFigsSeries(conditions);
fH_CA1 = initFigsSeries(conditions);
fH_Hp = initFigsSeries(conditions);
fH_RSC = initFigsSeries(conditions);
fH_VB2 = initFigsSeries(conditions);
fH_CA3 = initFigsSeries(conditions);
fH_mPFC = initFigsSeries(conditions);
fH_V1 = initFigsSeries(conditions);
fH_Cx = initFigsSeries(conditions);

% INITIALIASE STORAGE VARIABLES
for iCond = 1:numel(conditions)
  LFR_S1{iCond} = []; %#ok<*SAGROW>
end
LFR_VB = LFR_S1;
LFR_Th = LFR_S1;
LFR_DG = LFR_S1;
LFR_CA1 = LFR_S1;
LFR_Hp = LFR_S1;
LFR_RSC = LFR_S1;
LFR_VB2 = LFR_S1;
LFR_CA3 = LFR_S1;
LFR_mPFC = LFR_S1;
LFR_V1 = LFR_S1;
LFR_Cx = LFR_S1;

% LOOP THROUGH ANIMALS
for a = 1:numel(animals)
  animal = animals{a};
  load([dataDir filesep animal filesep animal '.mat'])
  animalColour = animalColours(animal);
  
  % INITIALISE FIGURE LEGENDS
  if a == 1
    S1Legend = initLegendsSeries(conditions);
    VBLegend = initLegendsSeries(conditions);
    ThLegend = initLegendsSeries(conditions);
    DGLegend = initLegendsSeries(conditions);
    CA1Legend = initLegendsSeries(conditions);
    HpLegend = initLegendsSeries(conditions);
    RSCLegend = initLegendsSeries(conditions);
    VB2Legend = initLegendsSeries(conditions);
    CA3Legend = initLegendsSeries(conditions);
    mPFCLegend = initLegendsSeries(conditions);
    V1Legend = initLegendsSeries(conditions);
    CxLegend = initLegendsSeries(conditions);
  else
    S1Legend = updateLegendsSeries(conditions, S1Legend);
    VBLegend = updateLegendsSeries(conditions, VBLegend);
    ThLegend = updateLegendsSeries(conditions, ThLegend);
    DGLegend = updateLegendsSeries(conditions, DGLegend);
    CA1Legend = updateLegendsSeries(conditions, CA1Legend);
    HpLegend = updateLegendsSeries(conditions, HpLegend);
    RSCLegend = updateLegendsSeries(conditions, RSCLegend);
    VB2Legend = updateLegendsSeries(conditions, VB2Legend);
    CA3Legend = updateLegendsSeries(conditions, CA3Legend);
    mPFCLegend = updateLegendsSeries(conditions, mPFCLegend);
    V1Legend = updateLegendsSeries(conditions, V1Legend);
    CxLegend = updateLegendsSeries(conditions, CxLegend);
  end
  
  % LOOP THROUGH DB ENTRIES
  fnsData = fieldnames(dataStruct.seriesData);
  for dbCount = 1:numel(fnsData)
    dbStruct = dataStruct.seriesData.(fnsData{dbCount});
    
    % TEST FOR EXCEPTIONS
    seriesName = seriesFromEntry(dbStruct.db(dbCount).entryName);
    if exclusionTest(seriesName, except)
      continue
    end
    
    % IDENTIFY CONDITION
    condition = determineCondition(seriesName, awake, anaesthesia);
    if ~condition
      continue
    end
    
    % IDENTIFY AREA
    area = determineAreaFromSeries(seriesName);
    if ~area
      continue
    end
    
    % CALCULATE FIRING RATES
    PR = sum(dbStruct.popData.MUAsAll,1);
    if isempty(PR)
      continue
    end
    params = dbStruct.conf.params;
    [mfr, mfr_1sthalf, mfr_2ndhalf, lfr1] = rateCalc(PR, params.srData);
    
    % DRAW AREA FIGURES AND ASSIGN DATA TO STORAGE VARIABLES
    for i = 1:2 % Repeat twice to deal with non-specific condition assignment
      if i == 2
        condition = 3;
      end
      
      if area == 1
        [fH_S1{condition}, LFR_S1{condition}, S1Legend{condition}] = plotAssignFR(fH_S1{condition},...
          lfr1, animalColour, S1Legend{condition}, animal, LFR_S1{condition});
      elseif area == 2
        [fH_VB{condition}, LFR_VB{condition}, VBLegend{condition}] = plotAssignFR(fH_VB{condition},...
          lfr1, animalColour, VBLegend{condition}, animal, LFR_VB{condition});
      elseif area == 3
        [fH_Th{condition}, LFR_Th{condition}, ThLegend{condition}] = plotAssignFR(fH_Th{condition},...
          lfr1, animalColour, ThLegend{condition}, animal, LFR_Th{condition});
      elseif area == 4
        [fH_DG{condition}, LFR_DG{condition}, DGLegend{condition}] = plotAssignFR(fH_DG{condition},...
          lfr1, animalColour, DGLegend{condition}, animal, LFR_DG{condition});
      elseif area == 5
        [fH_CA1{condition}, LFR_CA1{condition}, CA1Legend{condition}] = plotAssignFR(fH_CA1{condition},...
          lfr1, animalColour, CA1Legend{condition}, animal, LFR_CA1{condition});
      elseif area == 6
        [fH_Hp{condition}, LFR_Hp{condition}, HpLegend{condition}] = plotAssignFR(fH_Hp{condition},...
          lfr1, animalColour, HpLegend{condition}, animal, LFR_Hp{condition});
      elseif area == 7
        [fH_RSC{condition}, LFR_RSC{condition}, RSCLegend{condition}] = plotAssignFR(fH_RSC{condition},...
          lfr1, animalColour, RSCLegend{condition}, animal, LFR_RSC{condition});
      elseif area == 8
        [fH_VB2{condition}, LFR_VB2{condition}, VB2Legend{condition}] = plotAssignFR(fH_VB2{condition},...
          lfr1, animalColour, VB2Legend{condition}, animal, LFR_VB2{condition});
      elseif area == 10
        [fH_CA3{condition}, LFR_CA3{condition}, CA3Legend{condition}] = plotAssignFR(fH_CA3{condition},...
          lfr1, animalColour, CA3Legend{condition}, animal, LFR_CA3{condition});
      elseif area == 11
        [fH_mPFC{condition}, LFR_mPFC{condition}, mPFCLegend{condition}] = plotAssignFR(fH_mPFC{condition},...
          lfr1, animalColour, mPFCLegend{condition}, animal, LFR_mPFC{condition});
      elseif area == 12
        [fH_V1{condition}, LFR_V1{condition}, V1Legend{condition}] = plotAssignFR(fH_V1{condition},...
          lfr1, animalColour, V1Legend{condition}, animal, LFR_V1{condition});
      end
      
      if area == 1 || area == 7 || area == 11 || area == 12
        [fH_Cx{condition}, LFR_Cx{condition}, CxLegend{condition}] = plotAssignFR(fH_Cx{condition},...
          lfr1, animalColour, CxLegend{condition}, animal, LFR_Cx{condition});
      end
    end
  end
end

% UPDATE FIGURES WITH MEAN AND 95% CONFIDENCE INTERVAL
plotUpdate(fH_S1, S1Legend, LFR_S1, conditions, areas{1});
plotUpdate(fH_VB, VBLegend, LFR_VB, conditions, areas{2});
plotUpdate(fH_Th, ThLegend, LFR_Th, conditions, areas{3});
plotUpdate(fH_DG, DGLegend, LFR_DG, conditions, areas{4});
plotUpdate(fH_CA1, CA1Legend, LFR_CA1, conditions, areas{5});
plotUpdate(fH_Hp, HpLegend, LFR_Hp, conditions, areas{6});
plotUpdate(fH_RSC, RSCLegend, LFR_RSC, conditions, areas{7});
plotUpdate(fH_VB2, VB2Legend, LFR_VB2, conditions, areas{8});
plotUpdate(fH_CA3, CA1Legend, LFR_CA3, conditions, areas{10});
plotUpdate(fH_mPFC, mPFCLegend, LFR_mPFC, conditions, areas{11});
plotUpdate(fH_V1, V1Legend, LFR_V1, conditions, areas{12});
plotUpdate(fH_Cx, CxLegend, LFR_Cx, conditions, areas{13});

close all



function breakClause = exclusionTest(seriesName, exclusionsSeries)
% Check if a series should be excluded

for iExclude = 1:numel(exclusionsSeries)
  if strcmpi(seriesName, exclusionsSeries{iExclude})
    breakClause = true;
    break
  else
    breakClause = false;
  end
end
end

function condition = determineCondition(seriesName, awakeSeries, anaesthesiaSeries)
% Determine whether the series belongs to awake or anaesthesia conditions

condition = 0;
for entry = 1:numel(awakeSeries)
  if (numel(seriesName) >= 14 && strcmpi(awakeSeries{entry}, seriesName(1:14))) ||...
      (numel(seriesName) < 14 && strcmpi(awakeSeries{entry}, seriesName))
    condition = 1;
  end
end
for entry = 1:numel(anaesthesiaSeries)
  if numel(seriesName) >= 14 && strcmpi(anaesthesiaSeries{entry}, seriesName(1:14)) ||...
      (numel(seriesName) < 14 && strcmpi(anaesthesiaSeries{entry}, seriesName))
    condition = 2;
  end
end
end

function [fH, areaVar, fLegend] = plotAssignFR(fH, dataVec, animalColour, fLegend, animal, areaVar)
% Storage variable assignment and plotting. A helper functon for globalPCA
% and globalPCA2.

figure(fH); hold on
p = plot(dataVec, 'color',animalColour); hold off
if fLegend.update
  fLegend.lines = [fLegend.lines p];
  if isempty(fLegend.text)
    fLegend.text = {animal};
  else
    fLegend.text{numel(fLegend.text)+1} = animal;
  end
  fLegend.update = false;
end
if isempty(areaVar)
  areaVar = dataVec;
else
  if size(areaVar,2) < numel(dataVec)
    areaVar = [areaVar NaN(size(areaVar,1), numel(dataVec)-size(areaVar,2))];
  elseif size(areaVar,2) > numel(dataVec)
    dataVec = [dataVec NaN(1,size(areaVar,2)-numel(dataVec))];
  end
  areaVar = [areaVar; dataVec];
end
end

function plotUpdate(fH, fLegend, data, conditions, area)
% Figure updating.

for iCond = 1:numel(conditions)
  gLegend.lines = [];
  gLegend.text = {};
  figure(fH{iCond}); hold on
  if ~isempty(data{iCond})
    % Mean
    dataMean = mean(data{iCond},1,'omitnan');
    p = plot(dataMean, 'r:o', 'LineWidth',2, 'MarkerSize',2, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    uistack(p,'bottom');
    fLegend{iCond}.lines = [fLegend{iCond}.lines p];
    fLegend{iCond}.text{numel(fLegend{iCond}.text)+1} = 'Mean';
    
    % 95% confidence interval
    sd = std(data{iCond},1,1, 'omitnan');
    SEM = sd ./ sqrt(size(data{iCond},1));
    CI95 = zeros(2,numel(dataMean));
    meanCI95 = zeros(2,numel(dataMean));
    for i = 1:numel(dataMean)
      CI95(:,i) = (tinv([0.025 0.975], size(data{iCond},1)-1))';
      meanCI95(:,i) = bsxfun(@times, SEM(i), CI95(:,i));
    end
    pC1 = plot(dataMean+meanCI95(1,:), 'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    uistack(pC1,'bottom');
    pC2 = plot(dataMean+meanCI95(2,:), 'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    uistack(pC2,'bottom');
    fLegend{iCond}.lines = [fLegend{iCond}.lines pC2];
    fLegend{iCond}.text{numel(fLegend{iCond}.text)+1} = '95% conf';
  end
  hold off
  
  % Tidy and save the figure
  titleStr = ['firing rate in ' area];
  figName = [conditions{iCond} ': ' titleStr];
  figFileName = [conditions{iCond} '_' strrep(titleStr, ' ', '_')];
  figFileName = strrep(figFileName, ':', '_');
  tidySaveFig(fH{iCond}, 'Firing rate (APs/min)', 'Time (min)', fLegend{iCond}, figName, figFileName, 'northeast')
end
end

function tidySaveFig(fH, yLabel, xLabel, fLegend, figName, figFileName, legendPos)
% Figure tidying and saving

% Label the figure
figure(fH)
xlabel('Principal component');
ylabel(yLabel);
legend(fLegend.lines,fLegend.text, 'Interpreter','None', 'Location',legendPos);
title(figName);
set(gcf, 'Name',figName);

% Tidy the figure
ylimits = ylim;
xlimits = xlim;
set(gcf,'color','w');
legend boxoff
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
  'on', 'k', {xLabel}, xlimits, xticks,...
  'on', 'k', {yLabel}, ylimits, yticks);

% Save the figure in a fig file
hgsave(gcf, figFileName);

% Save the figure in a graphical file
label = [1.8 1.8];
margin = [0.3 0.6];
width = 1*15-label(1)-margin(1);
height = 1*15-label(2)-margin(2);
paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);
end