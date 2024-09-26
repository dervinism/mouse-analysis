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
[fH_varExp, fH_varExpArea] = initFigs(conditions, areas); % variance explained vs PCs
[fH_cumVarExp, fH_cumVarExpArea] = initFigs(conditions, areas); % cumulative variance explained vs PCs
[fH_corrPupil, fH_corrPupilArea] = initFigs(conditions, areas); % correlation coefs with pupil area vs PCs
[fH_varExpPupil, fH_varExpPupilArea] = initFigs(conditions, areas); % variance explained by pupil area vs PCs
[fH_corrMotion, fH_corrMotionArea] = initFigs(conditions, areas); % correlation coefs with total movement vs PCs
[fH_varExpMotion, fH_varExpMotionArea] = initFigs(conditions, areas); % variance explained by total movement vs PCs
[fH_corrLFP, fH_corrLFPArea] = initFigs(conditions, areas); % correlation coefs with LFP vs PCs
[fH_varExpLFP, fH_varExpLFPArea] = initFigs(conditions, areas); % variance explained by LFP vs PCs
[fH_corrPR, fH_corrPRArea] = initFigs(conditions, areas); % correlation coefs with PR vs PCs
[fH_varExpPR, fH_varExpPRArea] = initFigs(conditions, areas); % variance explained by PR vs PCs

% INITIALIASE STORAGE VARIABLES
for iCond = 1:numel(conditions)
  varExpArea_cond = {};
  for iArea = 1:numel(areas)
    varExpArea_cond{iArea} = []; %#ok<*SAGROW>
  end
  varExpArea{iCond} = varExpArea_cond; % for individual areas
end
cumVarExpArea = varExpArea;
corrPupilArea = varExpArea;
varExpPupilArea = varExpArea;
corrMotionArea = varExpArea;
varExpMotionArea = varExpArea;
corrLFPArea = varExpArea;
varExpLFPArea = varExpArea;
corrPRArea = varExpArea;
varExpPRArea = varExpArea;

% LOOP THROUGH ANIMALS
for a = 1:numel(animals)
  animal = animals{a};
  load([dataDir filesep animal filesep animal '.mat'])
  animalColour = animalColours(animal);
  
  % INITIALISE FIGURE LEGENDS
  if a == 1
    varExpAreaLegend = initLegends(conditions, areas);
    cumVarExpAreaLegend = initLegends(conditions, areas);
    corrPupilAreaLegend = initLegends(conditions, areas);
    varExpPupilAreaLegend = initLegends(conditions, areas);
    corrMotionAreaLegend = initLegends(conditions, areas);
    varExpMotionAreaLegend = initLegends(conditions, areas);
    corrLFPAreaLegend = initLegends(conditions, areas);
    varExpLFPAreaLegend = initLegends(conditions, areas);
    corrPRAreaLegend = initLegends(conditions, areas);
    varExpPRAreaLegend = initLegends(conditions, areas);
  else
    varExpAreaLegend = updateLegends(conditions, areas, varExpAreaLegend);
    cumVarExpAreaLegend = updateLegends(conditions, areas, cumVarExpAreaLegend);
    corrPupilAreaLegend = updateLegends(conditions, areas, corrPupilAreaLegend);
    varExpPupilAreaLegend = updateLegends(conditions, areas, varExpPupilAreaLegend);
    corrMotionAreaLegend = updateLegends(conditions, areas, corrMotionAreaLegend);
    varExpMotionAreaLegend = updateLegends(conditions, areas, varExpMotionAreaLegend);
    corrLFPAreaLegend = updateLegends(conditions, areas, corrLFPAreaLegend);
    varExpLFPAreaLegend = updateLegends(conditions, areas, varExpLFPAreaLegend);
    corrPRAreaLegend = updateLegends(conditions, areas, corrPRAreaLegend);
    varExpPRAreaLegend = updateLegends(conditions, areas, varExpPRAreaLegend);
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
    if area == 3 || area == 6 % these areas contain duplicate channels
      continue
    end
    
    % TEST IF PCA DATA EXISTS
    if ~isfield(dbStruct, 'pcaData')
      error(['No PCA data exist for series: ' fnsData{dbCount}])
    end
    
    % PICK UP LFP CHANNEL
    chOI = dbStruct.lfpPowerData.chOIDB;
    iCh = pickChan(area, animal, chOI);
    if isempty(iCh)
      continue
    end
    
    % DRAW AREA FIGURES AND ASSIGN DATA TO STORAGE VARIABLES
    for i = 1:2 % Repeat twice to deal with non-specific condition assignment
      if i == 2
        condition = 3;
      end
      if area == 5 % When area of interest is CA1 (we are interested in ripple rate)
        [fH_varExpArea{condition}{area}, varExpArea{condition}{area}, varExpAreaLegend{condition}{area}] = plotAssign(...
          fH_varExpArea{condition}{area}, dbStruct.pcaData.explainedHp{iCh}, animalColour,...
          varExpAreaLegend{condition}{area}, animal, varExpArea{condition}{area});
        [fH_cumVarExpArea{condition}{area}, cumVarExpArea{condition}{area}, cumVarExpAreaLegend{condition}{area}] = plotAssign(...
          fH_cumVarExpArea{condition}{area}, cumsum(dbStruct.pcaData.explainedHp{iCh}), animalColour,...
          cumVarExpAreaLegend{condition}{area}, animal, cumVarExpArea{condition}{area});
        if isfield(dbStruct.pcaData, 'rPCsHp2pupilAreaSpearman') && ~isempty(dbStruct.pcaData.rPCsHp2pupilAreaSpearman)
          [fH_corrPupilArea{condition}{area}, corrPupilArea{condition}{area}, corrPupilAreaLegend{condition}{area}] = plotAssign(...
            fH_corrPupilArea{condition}{area}, real(dbStruct.pcaData.rPCsHp2pupilAreaSpearman{iCh}), animalColour,...
            corrPupilAreaLegend{condition}{area}, animal, corrPupilArea{condition}{area});
          [fH_varExpPupilArea{condition}{area}, varExpPupilArea{condition}{area}, varExpPupilAreaLegend{condition}{area}] = plotAssign(...
            fH_varExpPupilArea{condition}{area}, (real(dbStruct.pcaData.rPCsHp2pupilAreaSpearman{iCh})).^2, animalColour,...
            varExpPupilAreaLegend{condition}{area}, animal, varExpPupilArea{condition}{area});
        end
        if isfield(dbStruct.pcaData, 'rPCsHp2motionSpearman') && ~isempty(dbStruct.pcaData.rPCsHp2motionSpearman)
          [fH_corrMotionArea{condition}{area}, corrMotionArea{condition}{area}, corrMotionAreaLegend{condition}{area}] = plotAssign(...
            fH_corrMotionArea{condition}{area}, real(dbStruct.pcaData.rPCsHp2motionSpearman{iCh}), animalColour,...
            corrMotionAreaLegend{condition}{area}, animal, corrMotionArea{condition}{area});
          [fH_varExpMotionArea{condition}{area}, varExpMotionArea{condition}{area}, varExpMotionAreaLegend{condition}{area}] = plotAssign(...
            fH_varExpMotionArea{condition}{area}, (real(dbStruct.pcaData.rPCsHp2motionSpearman{iCh})).^2, animalColour,...
            varExpMotionAreaLegend{condition}{area}, animal, varExpMotionArea{condition}{area});
        end
        [fH_corrLFPArea{condition}{area}, corrLFPArea{condition}{area}, corrLFPAreaLegend{condition}{area}] = plotAssign(...
          fH_corrLFPArea{condition}{area}, real(dbStruct.pcaData.rPCsHp2lfpSpearman{iCh}), animalColour,...
          corrLFPAreaLegend{condition}{area}, animal, corrLFPArea{condition}{area});
        [fH_varExpLFPArea{condition}{area}, varExpLFPArea{condition}{area}, varExpLFPAreaLegend{condition}{area}] = plotAssign(...
          fH_varExpLFPArea{condition}{area}, (real(dbStruct.pcaData.rPCsHp2lfpSpearman{iCh})).^2, animalColour,...
          varExpLFPAreaLegend{condition}{area}, animal, varExpLFPArea{condition}{area});
        [fH_corrPRArea{condition}{area}, corrPRArea{condition}{area}, corrPRAreaLegend{condition}{area}] = plotAssign(...
          fH_corrPRArea{condition}{area}, real(dbStruct.pcaData.rPCsHp2prSpearman{iCh}), animalColour,...
          corrPRAreaLegend{condition}{area}, animal, corrPRArea{condition}{area});
        [fH_varExpPRArea{condition}{area}, varExpPRArea{condition}{area}, varExpPRAreaLegend{condition}{area}] = plotAssign(...
          fH_varExpPRArea{condition}{area}, (real(dbStruct.pcaData.rPCsHp2prSpearman{iCh})).^2, animalColour,...
          varExpPRAreaLegend{condition}{area}, animal, varExpPRArea{condition}{area});
      else
        [fH_varExpArea{condition}{area}, varExpArea{condition}{area}, varExpAreaLegend{condition}{area}] = plotAssign(...
          fH_varExpArea{condition}{area}, dbStruct.pcaData.explained{iCh}, animalColour,...
          varExpAreaLegend{condition}{area}, animal, varExpArea{condition}{area});
        [fH_cumVarExpArea{condition}{area}, cumVarExpArea{condition}{area}, cumVarExpAreaLegend{condition}{area}] = plotAssign(...
          fH_cumVarExpArea{condition}{area}, cumsum(dbStruct.pcaData.explained{iCh}), animalColour,...
          cumVarExpAreaLegend{condition}{area}, animal, cumVarExpArea{condition}{area});
        if isfield(dbStruct.pcaData, 'rPCs2pupilAreaSpearman') && ~isempty(dbStruct.pcaData.rPCs2pupilAreaSpearman)
          [fH_corrPupilArea{condition}{area}, corrPupilArea{condition}{area}, corrPupilAreaLegend{condition}{area}] = plotAssign(...
            fH_corrPupilArea{condition}{area}, real(dbStruct.pcaData.rPCs2pupilAreaSpearman{iCh}), animalColour,...
            corrPupilAreaLegend{condition}{area}, animal, corrPupilArea{condition}{area});
          [fH_varExpPupilArea{condition}{area}, varExpPupilArea{condition}{area}, varExpPupilAreaLegend{condition}{area}] = plotAssign(...
            fH_varExpPupilArea{condition}{area}, (real(dbStruct.pcaData.rPCs2pupilAreaSpearman{iCh})).^2, animalColour,...
            varExpPupilAreaLegend{condition}{area}, animal, varExpPupilArea{condition}{area});
        end
        if isfield(dbStruct.pcaData, 'rPCs2motionSpearman') && ~isempty(dbStruct.pcaData.rPCs2motionSpearman)
          [fH_corrMotionArea{condition}{area}, corrMotionArea{condition}{area}, corrMotionAreaLegend{condition}{area}] = plotAssign(...
            fH_corrMotionArea{condition}{area}, real(dbStruct.pcaData.rPCs2motionSpearman{iCh}), animalColour,...
            corrMotionAreaLegend{condition}{area}, animal, corrMotionArea{condition}{area});
          [fH_varExpMotionArea{condition}{area}, varExpMotionArea{condition}{area}, varExpMotionAreaLegend{condition}{area}] = plotAssign(...
            fH_varExpMotionArea{condition}{area}, (real(dbStruct.pcaData.rPCs2motionSpearman{iCh})).^2, animalColour,...
            varExpMotionAreaLegend{condition}{area}, animal, varExpMotionArea{condition}{area});
        end
        [fH_corrLFPArea{condition}{area}, corrLFPArea{condition}{area}, corrLFPAreaLegend{condition}{area}] = plotAssign(...
          fH_corrLFPArea{condition}{area}, real(dbStruct.pcaData.rPCs2lfpSpearman{iCh}), animalColour,...
          corrLFPAreaLegend{condition}{area}, animal, corrLFPArea{condition}{area});
        [fH_varExpLFPArea{condition}{area}, varExpLFPArea{condition}{area}, varExpLFPAreaLegend{condition}{area}] = plotAssign(...
          fH_varExpLFPArea{condition}{area}, (real(dbStruct.pcaData.rPCs2lfpSpearman{iCh})).^2, animalColour,...
          varExpLFPAreaLegend{condition}{area}, animal, varExpLFPArea{condition}{area});
        [fH_corrPRArea{condition}{area}, corrPRArea{condition}{area}, corrPRAreaLegend{condition}{area}] = plotAssign(...
          fH_corrPRArea{condition}{area}, real(dbStruct.pcaData.rPCs2prSpearman{iCh}), animalColour,...
          corrPRAreaLegend{condition}{area}, animal, corrPRArea{condition}{area});
        [fH_varExpPRArea{condition}{area}, varExpPRArea{condition}{area}, varExpPRAreaLegend{condition}{area}] = plotAssign(...
          fH_varExpPRArea{condition}{area}, (real(dbStruct.pcaData.rPCs2prSpearman{iCh})).^2, animalColour,...
          varExpPRAreaLegend{condition}{area}, animal, varExpPRArea{condition}{area});
      end
    end
  end
end

% DRAW CONDITION FIGURES AND UPDATE AREA FIGURES WITH MEAN AND 95% CONFIDENCE INTERVAL
plotUpdate(fH_varExp, fH_varExpArea, varExpAreaLegend, varExpArea, conditions, areas,...
  'Variance explained', 'Variance explained', 'northeast');
plotUpdate(fH_cumVarExp, fH_cumVarExpArea, cumVarExpAreaLegend, cumVarExpArea, conditions, areas,...
  'Cumulative variance explained', 'Cumulative variance explained', 'southeast');
plotUpdate(fH_corrPupil, fH_corrPupilArea, corrPupilAreaLegend, corrPupilArea, conditions, areas,...
  'Correlation coef', 'Correlation with pupil area', 'northeast');
plotUpdate(fH_varExpPupil, fH_varExpPupilArea, varExpPupilAreaLegend, varExpPupilArea, conditions, areas,...
  'Variance explained', 'Variance with pupil area explained', 'northeast');
plotUpdate(fH_corrMotion, fH_corrMotionArea, corrMotionAreaLegend, corrMotionArea, conditions, areas,...
  'Correlation coef', 'Correlation with total movement', 'northeast');
plotUpdate(fH_varExpMotion, fH_varExpMotionArea, varExpMotionAreaLegend, varExpMotionArea, conditions, areas,...
  'Variance explained', 'Variance with total movement explained', 'northeast');
plotUpdate(fH_corrLFP, fH_corrLFPArea, corrLFPAreaLegend, corrLFPArea, conditions, areas,...
  'Correlation coef', 'Correlation with LFP', 'northeast');
plotUpdate(fH_varExpLFP, fH_varExpLFPArea, varExpLFPAreaLegend, varExpLFPArea, conditions, areas,...
  'Variance explained', 'Variance with LFP explained', 'northeast');
plotUpdate(fH_corrPR, fH_corrPRArea, corrPRAreaLegend, corrPRArea, conditions, areas,...
  'Correlation coef', 'Correlation with PR', 'northeast');
plotUpdate(fH_varExpPR, fH_varExpPRArea, varExpPRAreaLegend, varExpPRArea, conditions, areas,...
  'Variance explained', 'Variance with PR explained', 'northeast');

close all



function [condFigs, areaFigs] = initFigs(conditions, areas)
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
% Initialise figure legends

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
  
  figName = [conditions{iCond} ': ' titleStr];
  figFileName = [conditions{iCond} '_' strrep(titleStr, ' ', '_')];
  figFileName = strrep(figFileName, ':', '_');
  tidySaveFig(fH{iCond}, yLabel, gLegend, figName, figFileName, legendPos)
end
end

function tidySaveFig(fH, yLabel, fLegend, figName, figFileName, legendPos)
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
  'on', 'k', {'Principal component'}, xlimits, xticks,...
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