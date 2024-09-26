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
% [fH_varExpCond, fH_varExpArea] = initFigs(conditions, areas); % variance explained vs PCs
% [fH_cumVarExpCond, fH_cumVarExpArea] = initFigs(conditions, areas); % cumulative variance explained vs PCs
% [fH_corrLFPCond, fH_corrLFPArea] = initFigs(conditions, areas); % correlation coefs with LFP vs PCs
% [fH_varExpLFPCond, fH_varExpLFPArea] = initFigs(conditions, areas); % variance explained by LFP vs PCs
% [fH_corrPRCond, fH_corrPRArea] = initFigs(conditions, areas); % correlation coefs with PR vs PCs
% [fH_varExpPRCond, fH_varExpPRArea] = initFigs(conditions, areas); % variance explained by PR vs PCs
% 
% fH_varExpSeries = initFigsSeries(conditions);
% fH_cumVarExpSeries = initFigsSeries(conditions);
% fH_corrPupilSeries = initFigsSeries(conditions);
% fH_varExpPupilSeries = initFigsSeries(conditions);
% fH_corrMotionSeries = initFigsSeries(conditions);
% fH_varExpMotionSeries = initFigsSeries(conditions);

nPCs2plot = 5;
for iCond = 1:numel(conditions)
  for i = 1:nPCs2plot
    fH_VBvCA1_v{i}{iCond} = figure;
    fH_VBvCA1{i}{iCond} = figure;
  end
end

% INITIALIASE STORAGE VARIABLES
% for iCond = 1:numel(conditions)
%   varExpArea_cond = {};
%   for iArea = 1:numel(areas)
%     varExpArea_cond{iArea} = []; %#ok<*SAGROW>
%   end
%   varExpArea{iCond} = varExpArea_cond; % for individual areas
% end
% cumVarExpArea = varExpArea;
% corrLFPArea = varExpArea;
% varExpLFPArea = varExpArea;
% corrPRArea = varExpArea;
% varExpPRArea = varExpArea;
% 
% for iCond = 1:numel(conditions)
%   varExpSeries{iCond} = [];
% end
% cumVarExpSeries = varExpSeries;
% corrPupilSeries = varExpSeries;
% varExpPupilSeries = varExpSeries;
% corrMotionSeries = varExpSeries;
% varExpMotionSeries = varExpSeries;

% LOOP THROUGH ANIMALS
for a = 1:numel(animals)
  animal = animals{a};
  load([dataDir filesep animal filesep animal '.mat'])
  animalColour = animalColours(animal);
  
  % INITIALISE FIGURE LEGENDS
  if a == 1
%     varExpAreaLegend = initLegends(conditions, areas);
%     cumVarExpAreaLegend = initLegends(conditions, areas);
%     corrLFPAreaLegend = initLegends(conditions, areas);
%     varExpLFPAreaLegend = initLegends(conditions, areas);
%     corrPRAreaLegend = initLegends(conditions, areas);
%     varExpPRAreaLegend = initLegends(conditions, areas);
%     
%     varExpSeriesLegend = initLegendsSeries(conditions);
%     cumVarExpSeriesLegend = initLegendsSeries(conditions);
%     corrPupilSeriesLegend = initLegendsSeries(conditions);
%     varExpPupilSeriesLegend = initLegendsSeries(conditions);
%     corrMotionSeriesLegend = initLegendsSeries(conditions);
%     varExpMotionSeriesLegend = initLegendsSeries(conditions);
    
    for i = 1:nPCs2plot
      VBvCA1Legend_v{i} = initLegendsSeries(conditions); %#ok<*SAGROW>
      VBvCA1Legend{i} = initLegendsSeries(conditions);
    end
  else
%     varExpAreaLegend = updateLegends(conditions, areas, varExpAreaLegend);
%     cumVarExpAreaLegend = updateLegends(conditions, areas, cumVarExpAreaLegend);
%     corrLFPAreaLegend = updateLegends(conditions, areas, corrLFPAreaLegend);
%     varExpLFPAreaLegend = updateLegends(conditions, areas, varExpLFPAreaLegend);
%     corrPRAreaLegend = updateLegends(conditions, areas, corrPRAreaLegend);
%     varExpPRAreaLegend = updateLegends(conditions, areas, varExpPRAreaLegend);
%     
%     varExpSeriesLegend = updateLegendsSeries(conditions, varExpSeriesLegend);
%     cumVarExpSeriesLegend = updateLegendsSeries(conditions, cumVarExpSeriesLegend);
%     corrPupilSeriesLegend = updateLegendsSeries(conditions, corrPupilSeriesLegend);
%     varExpPupilSeriesLegend = updateLegendsSeries(conditions, varExpPupilSeriesLegend);
%     corrMotionSeriesLegend = updateLegendsSeries(conditions, corrMotionSeriesLegend);
%     varExpMotionSeriesLegend = updateLegendsSeries(conditions, varExpMotionSeriesLegend);
    
    for i = 1:nPCs2plot
      VBvCA1Legend_v{i} = updateLegendsSeries(conditions, VBvCA1Legend_v{i});
      VBvCA1Legend{i} = updateLegendsSeries(conditions, VBvCA1Legend{i});
    end
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
    if ~area || area == 1 || area == 3 || area == 4 || area == 6 || area == 7 % these areas contain duplicate channels
      continue
    end
    
    % TEST IF PCA DATA EXISTS
    if ~isfield(dbStruct, 'pcaData2_noS1_noRSC') && ~isfield(dbStruct, 'pcaData2_noS1_noRSC_noVB')
      %error(['No PCA2 data exist for series: ' fnsData{dbCount}])
      continue
    end
    
    % DRAW AREA FIGURES AND ASSIGN DATA TO STORAGE VARIABLES
    for i = 1:2 % Repeat twice to deal with non-specific condition assignment
      if i == 2
        condition = 3;
      end
%       [fH_varExpArea{condition}{area}, varExpArea{condition}{area}, varExpAreaLegend{condition}{area}] = plotAssign(...
%         fH_varExpArea{condition}{area}, torow(dbStruct.pcaData2.explained_area), animalColour,...
%         varExpAreaLegend{condition}{area}, animal, varExpArea{condition}{area});
%       [fH_cumVarExpArea{condition}{area}, cumVarExpArea{condition}{area}, cumVarExpAreaLegend{condition}{area}] = plotAssign(...
%         fH_cumVarExpArea{condition}{area}, cumsum(torow(dbStruct.pcaData2.explained_area)), animalColour,...
%         cumVarExpAreaLegend{condition}{area}, animal, cumVarExpArea{condition}{area});
%       [fH_corrLFPArea{condition}{area}, corrLFPArea{condition}{area}, corrLFPAreaLegend{condition}{area}] = plotAssign(...
%         fH_corrLFPArea{condition}{area}, real(torow(dbStruct.pcaData2.rPCs2lfpSpearman_area)), animalColour,...
%         corrLFPAreaLegend{condition}{area}, animal, corrLFPArea{condition}{area});
%       [fH_varExpLFPArea{condition}{area}, varExpLFPArea{condition}{area}, varExpLFPAreaLegend{condition}{area}] = plotAssign(...
%         fH_varExpLFPArea{condition}{area}, (real(torow(dbStruct.pcaData2.rPCs2lfpSpearman_area))).^2, animalColour,...
%         varExpLFPAreaLegend{condition}{area}, animal, varExpLFPArea{condition}{area});
%       [fH_corrPRArea{condition}{area}, corrPRArea{condition}{area}, corrPRAreaLegend{condition}{area}] = plotAssign(...
%         fH_corrPRArea{condition}{area}, real(torow(dbStruct.pcaData2.rPCs2prSpearman_area)), animalColour,...
%         corrPRAreaLegend{condition}{area}, animal, corrPRArea{condition}{area});
%       [fH_varExpPRArea{condition}{area}, varExpPRArea{condition}{area}, varExpPRAreaLegend{condition}{area}] = plotAssign(...
%         fH_varExpPRArea{condition}{area}, (real(torow(dbStruct.pcaData2.rPCs2prSpearman_area))).^2, animalColour,...
%         varExpPRAreaLegend{condition}{area}, animal, varExpPRArea{condition}{area});
%       if area == 1
%         [fH_varExpSeries{condition}, varExpSeries{condition}, varExpSeriesLegend{condition}] = plotAssign(...
%           fH_varExpSeries{condition}, torow(dbStruct.pcaData2.explained), animalColour,...
%           varExpSeriesLegend{condition}, animal, varExpSeries{condition});
%         [fH_cumVarExpSeries{condition}, cumVarExpSeries{condition}, cumVarExpSeriesLegend{condition}] = plotAssign(...
%           fH_cumVarExpSeries{condition}, cumsum(torow(dbStruct.pcaData2.explained)), animalColour,...
%           cumVarExpSeriesLegend{condition}, animal, cumVarExpSeries{condition});
%         if isfield(dbStruct.pcaData2, 'rPCs2pupilAreaSpearman')
%           [fH_corrPupilSeries{condition}, corrPupilSeries{condition}, corrPupilSeriesLegend{condition}] = plotAssign(...
%             fH_corrPupilSeries{condition}, real(torow(dbStruct.pcaData2.rPCs2pupilAreaSpearman)), animalColour,...
%             corrPupilSeriesLegend{condition}, animal, corrPupilSeries{condition});
%           [fH_varExpPupilSeries{condition}, varExpPupilSeries{condition}, varExpPupilSeriesLegend{condition}] = plotAssign(...
%             fH_varExpPupilSeries{condition}, (real(torow(dbStruct.pcaData2.rPCs2pupilAreaSpearman))).^2, animalColour,...
%             varExpPupilSeriesLegend{condition}, animal, varExpPupilSeries{condition});
%         end
%         if isfield(dbStruct.pcaData2, 'rPCs2motionSpearman')
%           [fH_corrMotionSeries{condition}, corrMotionSeries{condition}, corrMotionSeriesLegend{condition}] = plotAssign(...
%             fH_corrMotionSeries{condition}, real(torow(dbStruct.pcaData2.rPCs2motionSpearman)), animalColour,...
%             corrMotionSeriesLegend{condition}, animal, corrMotionSeries{condition});
%           [fH_varExpMotionSeries{condition}, varExpMotionSeries{condition}, varExpMotionSeriesLegend{condition}] = plotAssign(...
%             fH_varExpMotionSeries{condition}, (real(torow(dbStruct.pcaData2.rPCs2motionSpearman))).^2, animalColour,...
%             varExpMotionSeriesLegend{condition}, animal, varExpMotionSeries{condition});
%         end
%       end
      
      % DRAW DIAGONAL FIGURES
      if area == 2
        if isfield(dbStruct, 'pcaData2_noS1_noRSC')
          VB_PCs_v = torow(dbStruct.pcaData2_noS1_noRSC.explained_area(1:nPCs2plot));
          VB_PCs = cumsum(torow(dbStruct.pcaData2_noS1_noRSC.explained_area(1:nPCs2plot)));
        elseif isfield(dbStruct, 'pcaData2_noS1_noRSC_noVB')
          VB_PCs_v = torow(dbStruct.pcaData2_noS1_noRSC_noVB.explained_area(1:nPCs2plot));
          VB_PCs = cumsum(torow(dbStruct.pcaData2_noS1_noRSC_noVB.explained_area(1:nPCs2plot)));
        else
          error(['No PCA data for ' seriesName '. Please run PCA script for this series of data.']);
        end
      elseif area == 5
        if isfield(dbStruct, 'pcaData2_noS1_noRSC')
          CA1_PCs_v = torow(dbStruct.pcaData2_noS1_noRSC.explained_area(1:nPCs2plot));
          CA1_PCs = cumsum(torow(dbStruct.pcaData2_noS1_noRSC.explained_area(1:nPCs2plot)));
        elseif isfield(dbStruct, 'pcaData2_noS1_noRSC_noVB')
          CA1_PCs_v = torow(dbStruct.pcaData2_noS1_noRSC_noVB.explained_area(1:nPCs2plot));
          CA1_PCs = cumsum(torow(dbStruct.pcaData2_noS1_noRSC_noVB.explained_area(1:nPCs2plot)));
        else
          error(['No PCA data for ' seriesName '. Please run PCA script for this series of data.']);
        end
        for iPC = 1:nPCs2plot
          figure(fH_VBvCA1_v{iPC}{condition}); hold on
          p = plot(VB_PCs_v(iPC),CA1_PCs_v(iPC), '.', 'MarkerSize',10, 'Color',animalColour); hold off
          if VBvCA1Legend_v{iPC}{condition}.update
            VBvCA1Legend_v{iPC}{condition}.lines = [VBvCA1Legend_v{iPC}{condition}.lines p];
            if isempty(VBvCA1Legend_v{iPC}{condition}.text)
              VBvCA1Legend_v{iPC}{condition}.text = {animal};
            else
              VBvCA1Legend_v{iPC}{condition}.text{numel(VBvCA1Legend_v{iPC}{condition}.text)+1} = animal;
            end
            VBvCA1Legend_v{iPC}{condition}.update = false;
          end
          
          figure(fH_VBvCA1{iPC}{condition}); hold on
          p = plot(VB_PCs(iPC),CA1_PCs(iPC), '.', 'MarkerSize',10, 'Color',animalColour); hold off
          if VBvCA1Legend{iPC}{condition}.update
            VBvCA1Legend{iPC}{condition}.lines = [VBvCA1Legend{iPC}{condition}.lines p];
            if isempty(VBvCA1Legend{iPC}{condition}.text)
              VBvCA1Legend{iPC}{condition}.text = {animal};
            else
              VBvCA1Legend{iPC}{condition}.text{numel(VBvCA1Legend{iPC}{condition}.text)+1} = animal;
            end
            VBvCA1Legend{iPC}{condition}.update = false;
          end
        end
      end
    end
  end
end

% DRAW CONDITION FIGURES AND UPDATE AREA FIGURES WITH MEAN AND 95% CONFIDENCE INTERVAL
% plotUpdate(fH_varExpCond, fH_varExpArea, varExpAreaLegend, varExpArea, conditions, areas,...
%   'Variance explained', 'Variance explained', 'northeast');
% plotUpdate(fH_cumVarExpCond, fH_cumVarExpArea, cumVarExpAreaLegend, cumVarExpArea, conditions, areas,...
%   'Cumulative variance explained', 'Cumulative variance explained', 'southeast');
% plotUpdate(fH_corrLFPCond, fH_corrLFPArea, corrLFPAreaLegend, corrLFPArea, conditions, areas,...
%   'Correlation coef', 'Correlation with LFP', 'northeast');
% plotUpdate(fH_varExpLFPCond, fH_varExpLFPArea, varExpLFPAreaLegend, varExpLFPArea, conditions, areas,...
%   'Variance explained', 'Variance with LFP explained', 'northeast');
% plotUpdate(fH_corrPRCond, fH_corrPRArea, corrPRAreaLegend, corrPRArea, conditions, areas,...
%   'Correlation coef', 'Correlation with PR', 'northeast');
% plotUpdate(fH_varExpPRCond, fH_varExpPRArea, varExpPRAreaLegend, varExpPRArea, conditions, areas,...
%   'Variance explained', 'Variance with PR explained', 'northeast');

% UPDATE SERIES FIGURES WITH MEAN AND 95% CONFIDENCE INTERVAL
% plotUpdateSeries(fH_varExpSeries, varExpSeriesLegend, varExpSeries, conditions,...
%   'Variance explained', 'Variance explained', 'northeast');
% plotUpdateSeries(fH_cumVarExpSeries, cumVarExpSeriesLegend, cumVarExpSeries, conditions,...
%   'Cumulative variance explained', 'Cumulative variance explained', 'southeast');
% plotUpdateSeries(fH_corrPupilSeries, corrPupilSeriesLegend, corrPupilSeries, conditions,...
%   'Correlation coef', 'Correlation with pupil area', 'northeast');
% plotUpdateSeries(fH_varExpPupilSeries, varExpPupilSeriesLegend, varExpPupilSeries, conditions,...
%   'Variance explained', 'Variance with pupil area explained', 'northeast');
% plotUpdateSeries(fH_corrMotionSeries, corrMotionSeriesLegend, corrMotionSeries, conditions,...
%   'Correlation coef', 'Correlation with total movement', 'northeast');
% plotUpdateSeries(fH_varExpMotionSeries, varExpMotionSeriesLegend, varExpMotionSeries, conditions,...
%   'Variance explained', 'Variance with total movement explained', 'northeast');

% UPDATE DIAGONAL FIGURES
for iCond = 1:numel(conditions)
  for i = 1:nPCs2plot
    figure(fH_VBvCA1_v{i}{iCond}); hold on
    plot([0 100],[0 100], 'k--'); hold off
    titleStr = ['VBvCA1 explained variance during ' conditions{iCond} ': PC' num2str(i)];
    title(titleStr);
    
    % Tidy and save the figure
    figName = titleStr;
    figFileName = strrep(titleStr, ' ', '_');
    figFileName = strrep(figFileName, ':', '_');
    tidySaveDiagonal(fH_VBvCA1_v{i}{iCond}, 'Explained variance in CA1',...
      'Explained variance in VB', VBvCA1Legend_v{i}{iCond}, figName, figFileName, 'southeast');
    
    %
    figure(fH_VBvCA1{i}{iCond}); hold on
    plot([0 100],[0 100], 'k--'); hold off
    titleStr = ['VBvCA1 cumulative explained variance during ' conditions{iCond} ': PC' num2str(i)];
    title(titleStr);
    
    % Tidy and save the figure
    figName = titleStr;
    figFileName = strrep(titleStr, ' ', '_');
    figFileName = strrep(figFileName, ':', '_');
    tidySaveDiagonal(fH_VBvCA1{i}{iCond}, 'Cumulative explained variance in CA1',...
      'Cumulative explained variance in VB', VBvCA1Legend{i}{iCond}, figName, figFileName, 'southeast');
  end
end

close all



function [condFigs, areaFigs] = initFigs(conditions, areas) %#ok<*DEFNU>
% Initialise figures

for iCond = 1:numel(conditions)
  condFigs{iCond} = figure; %#ok<*AGROW> % for individual conditions
  areaFigs_cond = {};
  for iArea = 1:numel(areas)
    areaFigs_cond{iArea} = figure;
  end
  areaFigs{iCond} = areaFigs_cond; % for individual areas
end
end

function legendParams = initLegends(conditions, areas)
% Initialise figure legends

for iCond = 1:numel(conditions)
  for iArea = 1:numel(areas)
    legendParams{iCond}{iArea}.lines = [];
    legendParams{iCond}{iArea}.text = {};
    legendParams{iCond}{iArea}.update = true;
  end
end
end

function legendParams = updateLegends(conditions, areas, legendParams)
% Update figure legends

for iCond = 1:numel(conditions)
  for iArea = 1:numel(areas)
    legendParams{iCond}{iArea}.update = true;
  end
end
end

function plotUpdate(fH, fH_Area, fLegend, data, conditions, areas, yLabel, titleStr, legendPos)
% Figure updating

for iCond = 1:numel(conditions)
  gLegend.lines = [];
  gLegend.text = {};
  for iArea = 1:numel(areas)
    figure(fH_Area{iCond}{iArea}); hold on
    if ~isempty(data{iCond}{iArea})
      % Mean
      dataMean = sum(data{iCond}{iArea},1)./size(data{iCond}{iArea},1);
      p = plot(dataMean, 'r:o', 'LineWidth',2, 'MarkerSize',2, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
      uistack(p,'bottom');
      fLegend{iCond}{iArea}.lines = [fLegend{iCond}{iArea}.lines p];
      fLegend{iCond}{iArea}.text{numel(fLegend{iCond}{iArea}.text)+1} = 'Mean';
      
      % 95% confidence interval
      sd = std(data{iCond}{iArea},1,1, 'omitnan');
      SEM = sd ./ sqrt(size(data{iCond}{iArea},1));
      CI95 = zeros(2,numel(dataMean));
      meanCI95 = zeros(2,numel(dataMean));
      for i = 1:numel(dataMean)
        CI95(:,i) = (tinv([0.025 0.975], size(data{iCond}{iArea},1)-1))';
        meanCI95(:,i) = bsxfun(@times, SEM(i), CI95(:,i));
      end
      pC1 = plot(dataMean+meanCI95(1,:), 'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
      uistack(pC1,'bottom');
      pC2 = plot(dataMean+meanCI95(2,:), 'r:o', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
      uistack(pC2,'bottom');
      fLegend{iCond}{iArea}.lines = [fLegend{iCond}{iArea}.lines pC2];
      fLegend{iCond}{iArea}.text{numel(fLegend{iCond}{iArea}.text)+1} = '95% conf';
    else
      dataMean = [];
    end
    hold off
    
    % Tidy and save the figure
    figName = [areas{iArea} ' ' conditions{iCond} ': ' titleStr];
    figFileName = [areas{iArea} '_' conditions{iCond} '_' strrep(titleStr, ' ', '_')];
    figFileName = strrep(figFileName, ':', '_');
    tidySaveFig(fH_Area{iCond}{iArea}, yLabel, fLegend{iCond}{iArea}, figName, figFileName, legendPos)
    
    % Draw the condition figure
    if ~isempty(dataMean)
      figure(fH{iCond}); hold on
      areaColour = areaColours(areas{iArea});
      p = plot(dataMean, 'color',areaColour, 'LineWidth',2);
      gLegend.lines = [gLegend.lines p];
      if isempty(gLegend.text)
        gLegend.text = {areas{iArea}}; %#ok<CCAT1>
      else
        gLegend.text{numel(gLegend.text)+1} = areas{iArea};
      end
      pC1 = plot(dataMean+meanCI95(1,:), '--', 'color',areaColour, 'LineWidth',1);
      uistack(pC1,'bottom');
      pC2 = plot(dataMean+meanCI95(2,:), '--', 'color',areaColour, 'LineWidth',1);
      uistack(pC2,'bottom');
      hold off
    end
  end
  
  figName = ['areas ' conditions{iCond} ': ' titleStr];
  figFileName = ['areas ' conditions{iCond} '_' strrep(titleStr, ' ', '_')];
  figFileName = strrep(figFileName, ':', '_');
  tidySaveFig(fH{iCond}, yLabel, gLegend, figName, figFileName, legendPos)
end
end

function tidySaveDiagonal(fH, yLabel, xLabel, fLegend, figName, figFileName, legendPos)
% Figure tidying and saving. A helper function of plotUpdateSeries.

% Label the figure
figure(fH)
xlabel('Principal component');
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
margin = [0.3 0.3];
width = 1*15-label(1)-margin(1);
height = 1*15-label(2)-margin(2);
paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);
end