function cohFreqProfilePlotIndividual(areas, conditions, FOI, areaCohFOIindividual, animalIDs, series, figFileName, options)
% cohFreqProfilePlotIndividual(areas, conditions, FOI, areaCohFOIindividual, animalIDs, series, figFileName, options)
%
% Function produces figures containing coherence frequency profiles of
% individual units or individual population rates. The mean coherence and
% 95% confidence intervals are overlayed on top of individual curves.
% Input: areas - area acronyms or their comparisons.
%        conditions - recording conditions corresponding to states of
%                     vigilance.
%        FOI - frequencies of interest. They could be a single vector or a
%              cell array containing multiple vectors.
%        areaCohFOIindividual - a cell array of coherence frequency profiles.
%        animalIDs - a structure with animal IDs in case the figures need
%                    to be generated on a per animal per area basis.
%        series - a structure with animal IDs in case the figures need to
%                 be generated on a per recording per animal per area basis.
%        figFileName - a short figure name (not including the file path)
%                      containing appropriate formating characters
%                      accomodating area, condition, and, if needed, animal
%                      ID and series ID.
%        options - a structure with the following fields:
%          mainFolder
%          coherenceFrequencyProfilesSubfolder
%          figSize
%          figTitle
%          freqLim
%          cohLim
%          iAreasOI.


%% Parse user input
if nargin < 8
  options.mainFolder = pwd;
  options.cohFrequencyProfilesSubfolder = 'coherenceFrequencyProfiles';
  options.figSize = 15;
  options.figTitle = [];
  options.freqLim = [10e-3 - 0.002   30];
  options.cohLim = [-pi pi];
  options.iAreasOI = 1:numel(areas);
else
  if ~isfield(options,'mainFolder')
    options.mainFolder = pwd;
  end
  if ~isfield(options,'cohFrequencyProfilesSubfolder')
    options.cohFrequencyProfilesSubfolder = 'coherenceFrequencyProfiles';
  end
  if ~isfield(options,'figSize')
    options.figSize = 15;
  end
  if ~isfield(options,'figTitle')
    options.figTitle = [];
  end
  if ~isfield(options,'freqLim')
    options.freqLim = [10e-3 - 0.002   30];
  end
  if ~isfield(options,'cohLim')
    options.cohLim = [-pi pi];
  end
  if ~isfield(options,'iAreasOI')
    options.iAreasOI = 1:numel(areas);
  end
end

if ~exist([options.mainFolder filesep options.cohFrequencyProfilesSubfolder], 'file')
  mkdir([options.mainFolder filesep options.cohFrequencyProfilesSubfolder]);
end


%% Display the figures
if isempty(animalIDs) && isempty(series) % Displaying figures on a per area basis
  for iCond = 1:min([2 numel(conditions)]) % Loop through conditions
    for iArea = 1:numel(options.iAreasOI) % Loop through areas
      if ~isempty(areaCohFOIindividual{iCond}{options.iAreasOI(iArea)})
        inds = 1:size(areaCohFOIindividual{iCond}{options.iAreasOI(iArea)},1);
        coreFunc(iCond, options.iAreasOI(iArea), areas, conditions, FOI, areaCohFOIindividual, [], [], inds, figFileName, options);
      end
    end
  end
elseif ~isempty(animalIDs) && isempty(series) % Displaying figures on a per area per animal basis
  for iCond = 1:min([2 numel(conditions)]) % Loop through conditions
    for iArea = 1:numel(options.iAreasOI) % Loop through areas
      if ~isempty(areaCohFOIindividual{iCond}{options.iAreasOI(iArea)})
        uniqueAnimals = unique(animalIDs{iCond}{options.iAreasOI(iArea)});
        nAnimals = numel(uniqueAnimals);
        for iAnimal = 1:nAnimals
          animalInds = find(ismember(animalIDs{iCond}{options.iAreasOI(iArea)}, uniqueAnimals{iAnimal}));
          cohFreqProfilePlotIndividual(iCond, options.iAreasOI(iArea), areas, conditions, FOI, areaCohFOIindividual,...
            uniqueAnimals{iAnimal}, [], animalInds, figFileName, options);
        end
      end
    end
  end
elseif ~isempty(animalIDs) && ~isempty(series) % Displaying figures on a per area per animal per series basis
  for iCond = 1:min([2 numel(conditions)]) % Loop through conditions
    for iArea = 1:numel(options.iAreasOI) % Loop through areas
      if ~isempty(areaCohFOIindividual{iCond}{options.iAreasOI(iArea)})
        uniqueSeries = unique(series{iCond}{options.iAreasOI(iArea)});
        nSeries = numel(uniqueSeries);
        for iSeries = 1:nSeries
          seriesInds = find(ismember(series{iCond}{options.iAreasOI(iArea)}, uniqueSeries{iSeries}));
          cohFreqProfilePlotIndividual(iCond, options.iAreasOI(iArea), areas, conditions, FOI, areaCohFOIindividual,...
            animalIDs{iCond}{options.iAreasOI(iArea)}{seriesInds(1)}, uniqueSeries{iSeries}, seriesInds, figFileName, options);
        end
      end
    end
  end
end



%% Local functions
function coreFunc(iCond, iArea, areas, conditions, FOI, areaCohFOIindividual, animalID, seriesID, inds, figFileName, options)

% Get coherence values
if iscell(FOI)
  areaFOI = FOI{iCond}{iArea};
else
  areaFOI = FOI;
end
[cohMean, cohCI95] = datamean(areaCohFOIindividual{iCond}{iArea});

if sum(~isnan(cohMean))
  % Plot individual curves
  fH = figure();
  txtCondMean = {};
  plotCondMean = [];
  semilogx([10e-4 10e2],[0 0], 'k:');
  hold on
  for u = 1:numel(inds)
    coh = areaCohFOIindividual{iCond}{iArea}(inds(u),:);
    semilogx(areaFOI, coh);
  end
  
  if sum(~isnan(cohMean) & ~isnan(cohCI95(1,:)) & ~isnan(cohCI95(2,:)))
    % Plot the mean curve
    if iCond == 1
      p1 = semilogx(areaFOI(~isnan(cohMean) & ~isnan(cohCI95(1,:)) & ~isnan(cohCI95(2,:))),...
        cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:)) & ~isnan(cohCI95(2,:))),...
        'g:', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Unit mean
      plotCondMean = [plotCondMean p1]; %#ok<*AGROW>
      txtCondMean{numel(txtCondMean)+1} = 'Mean';
    elseif iCond == 2
      p1 = semilogx(areaFOI(~isnan(cohMean) & ~isnan(cohCI95(1,:)) & ~isnan(cohCI95(2,:))),...
        cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:)) & ~isnan(cohCI95(2,:))),...
        'r:', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Unit mean
      plotCondMean = [plotCondMean p1];
      txtCondMean{numel(txtCondMean)+1} = 'Mean';
      uistack(p1,'bottom');
    end
    
    % Plot 95% confidence intervals
    if iCond == 1
      pC1 = semilogx(areaFOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
        cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))) + cohCI95(1,~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
        'g:', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Upper 95% confidence limitof the mean
      uistack(pC1,'bottom');
      pC2 = semilogx(areaFOI(~isnan(cohMean) & ~isnan(cohCI95(2,:))),...
        cohMean(~isnan(cohMean) & ~isnan(cohCI95(2,:))) + cohCI95(2,~isnan(cohMean) & ~isnan(cohCI95(2,:))),...
        'g:', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Lower 95% confidence limitof the mean
      uistack(pC2,'bottom');
    elseif iCond == 2
      pC1 = semilogx(areaFOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
        cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))) + cohCI95(1,~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
        'r:', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Upper 95% confidence limitof the mean
      uistack(pC1,'bottom');
      pC2 = semilogx(areaFOI(~isnan(cohMean) & ~isnan(cohCI95(2,:))),...
        cohMean(~isnan(cohMean) & ~isnan(cohCI95(2,:))) + cohCI95(2,~isnan(cohMean) & ~isnan(cohCI95(2,:))),...
        'r:', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Lower 95% confidence limitof the mean
      uistack(pC2,'bottom');
    end
    plotCondMean = [plotCondMean pC1];
    txtCondMean{numel(txtCondMean)+1} = '95% CI';
  end
  hold off
  
  % Tidy the figure
  yLim = ylim;
  yLim = [options.cohLim(1) min([options.cohLim(2) yLim(2)])];
  if ~isempty(plotCondMean) && ~isempty(txtCondMean) && numel(plotCondMean) == numel(txtCondMean)
    legend(plotCondMean,txtCondMean, 'Interpreter','None');
    legend boxoff
  end
  if isempty(animalID) && isempty(seriesID)
    if isempty(options.figTitle)
      titleStr = sprintf('Individual coherences: %s %s', areas{iArea}, conditions{iCond});
    else
      titleStr = sprintf(options.figTitle, areas{iArea}, conditions{iCond});
    end
  elseif ~isempty(animalID) && isempty(seriesID)
    if isempty(options.figTitle)
      titleStr = sprintf('Individual coherences: %s %s %s', areas{iArea}, conditions{iCond}, animalID);
    else
      titleStr = sprintf(options.figTitle, areas{iArea}, conditions{iCond}, animalID);
    end
  elseif ~isempty(animalID) && ~isempty(seriesID)
    if isempty(options.figTitle)
      titleStr = sprintf('Individual coherences: %s %s %s %s', areas{iArea}, conditions{iCond}, animalID, seriesID);
    else
      titleStr = sprintf(options.figTitle, areas{iArea}, conditions{iCond}, animalID, seriesID);
    end
  end
  titleStr = strrep(titleStr, 'Vs', 'wrt');
  set(fH,'color','w');
  ax1 = axesProperties({titleStr}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
    'on', 'k', {'Frequency (Hz)'}, options.freqLim, [0.001 0.01 0.1 1 10 30 100 1000],...
    'on', 'k', {'Coherence'}, yLim, 0:0.1:1);
  ax1.XTickLabel = {'0.001','0.01','0.1','1','10','30','100','1000'};
  
  % Save the figure
  if isempty(animalID) && isempty(seriesID)
    figFileNameEnd = sprintf(figFileName, areas{iArea}, conditions{iCond});
  elseif ~isempty(animalID) && isempty(seriesID)
    figFileNameEnd = sprintf(figFileName, areas{iArea}, conditions{iCond}, animalID);
  elseif ~isempty(animalID) && ~isempty(seriesID)
    figFileNameEnd = sprintf(figFileName, areas{iArea}, conditions{iCond}, animalID, seriesID);
  end
  figFileName = [options.mainFolder filesep options.cohFrequencyProfilesSubfolder filesep figFileNameEnd];
  hgsave(fH, figFileName);
  label = [2 1.6];
  margin = [0.3 0.55];
  width = 1*options.figSize-label(1)-margin(1);
  height = 1*options.figSize-label(2)-margin(2);
  paperSize = resizeFig(fH, ax1, width, height, label, margin, 0);
  hgsave(fH, [figFileName '.fig']);
  exportFig(fH, [figFileName '.png'],'-dpng','-r300', paperSize);
  close(fH);
end