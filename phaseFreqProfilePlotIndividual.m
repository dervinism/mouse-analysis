function phaseFreqProfilePlotIndividual(areas, conditions, FOI, areaPhaseFOIindividual, animalIDs, series, figFileName, options)
% phaseFreqProfilePlotIndividual(areas, conditions, FOI, areaPhaseFOIindividual, animalIDs, series, figFileName, options)
%
% Function produces figures containing phase frequency profiles of
% individual units or individual population rates. The mean phase and 95%
% confidence intervals are overlayed on top of individual curves.
% Input: areas - area acronyms or their comparisons.
%        conditions - recording conditions corresponding to states of
%                     vigilance.
%        FOI - frequencies of interest. They could be a single vector or a
%              cell array containing multiple vectors.
%        areaPhaseFOIindividual - a cell array of phase frequency profiles.
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
%          phaseFrequencyProfilesSubfolder
%          figSize
%          figTitle
%          freqLim
%          phaseLim
%          cutoffFreq
%          iAreasOI.


%% Parse user input
if nargin < 8
  options.mainFolder = pwd;
  options.phaseFrequencyProfilesSubfolder = 'phaseFrequencyProfiles';
  options.figSize = 15;
  options.figTitle = [];
  options.freqLim = [10e-3 - 0.002   30];
  options.phaseLim = [-pi pi];
  options.cutoffFreq = 30;
  options.iAreasOI = 1:numel(areas);
else
  if ~isfield(options,'mainFolder')
    options.mainFolder = pwd;
  end
  if ~isfield(options,'phaseFrequencyProfilesSubfolder')
    options.phaseFrequencyProfilesSubfolder = 'phaseFrequencyProfiles';
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
  if ~isfield(options,'phaseLim')
    options.phaseLim = [-pi pi];
  end
  if ~isfield(options,'cutoffFreq')
    options.cutoffFreq = 30;
  end
  if ~isfield(options,'iAreasOI')
    options.iAreasOI = 1:numel(areas);
  end
end

if ~exist([options.mainFolder filesep options.phaseFrequencyProfilesSubfolder], 'file')
  mkdir([options.mainFolder filesep options.phaseFrequencyProfilesSubfolder]);
end


%% Display the figures
if isempty(animalIDs) && isempty(series) % Displaying figures on a per area basis
  for iCond = 1:min([2 numel(conditions)]) % Loop through conditions
    for iArea = 1:numel(options.iAreasOI) % Loop through areas
      if ~isempty(areaPhaseFOIindividual{iCond}{options.iAreasOI(iArea)})
        inds = 1:size(areaPhaseFOIindividual{iCond}{options.iAreasOI(iArea)},1);
        coreFunc(iCond, options.iAreasOI(iArea), areas, conditions, FOI, areaPhaseFOIindividual, [], [], inds, figFileName, options);
      end
    end
  end
elseif ~isempty(animalIDs) && isempty(series) % Displaying figures on a per area per animal basis
  for iCond = 1:min([2 numel(conditions)]) % Loop through conditions
    for iArea = 1:numel(options.iAreasOI) % Loop through areas
      if ~isempty(areaPhaseFOIindividual{iCond}{options.iAreasOI(iArea)})
        uniqueAnimals = unique(animalIDs{iCond}{options.iAreasOI(iArea)});
        nAnimals = numel(uniqueAnimals);
        for iAnimal = 1:nAnimals
          animalInds = find(ismember(animalIDs{iCond}{options.iAreasOI(iArea)}, uniqueAnimals{iAnimal}));
          coreFunc(iCond, options.iAreasOI(iArea), areas, conditions, FOI, areaPhaseFOIindividual,...
            uniqueAnimals{iAnimal}, [], animalInds, figFileName, options);
        end
      end
    end
  end
elseif ~isempty(animalIDs) && ~isempty(series) % Displaying figures on a per area per animal per series basis
  for iCond = 1:min([2 numel(conditions)]) % Loop through conditions
    for iArea = 1:numel(options.iAreasOI) % Loop through areas
      if ~isempty(areaPhaseFOIindividual{iCond}{options.iAreasOI(iArea)})
        uniqueSeries = unique(series{iCond}{options.iAreasOI(iArea)});
        nSeries = numel(uniqueSeries);
        for iSeries = 1:nSeries
          seriesInds = find(ismember(series{iCond}{options.iAreasOI(iArea)}, uniqueSeries{iSeries}));
          coreFunc(iCond, options.iAreasOI(iArea), areas, conditions, FOI, areaPhaseFOIindividual,...
            animalIDs{iCond}{options.iAreasOI(iArea)}{seriesInds(1)}, uniqueSeries{iSeries}, seriesInds, figFileName, options);
        end
      end
    end
  end
end



%% Local functions
function coreFunc(iCond, iArea, areas, conditions, FOI, areaPhaseFOIindividual, animalID, seriesID, inds, figFileName, options)

% Get phase values
if iscell(FOI)
  areaFOI = FOI{iCond}{iArea};
else
  areaFOI = FOI;
end
phase = areaPhaseFOIindividual{iCond}{iArea}(inds,:);
[phaseMean, phaseCI95] = datamean(phase, 'circularNP');
% phaseCI95(phaseCI95 == 0) = NaN;
% phaseMean(isnan(phaseCI95(1,:)) | isnan(phaseCI95(2,:))) = NaN;
if size(phase,1) == 1
  phaseCI95 = repmat([-pi; pi], 1, numel(phaseMean));
end

% Adjust the mean phase values
if ~isempty(phaseMean) && sum(~isnan(phaseMean))
  [phaseMean, phaseMeanROI, phaseCI95ROI] = adjustPhase(phaseMean, areaFOI, phaseCI95, options.cutoffFreq);
else
  return
end

if sum(~isnan(phaseMean))
  % Plot individual curves
  fH = figure();
  txtCondMean = {};
  plotCondMean = [];
  semilogx([10e-4 10e2],[0 0], 'k:');
  hold on
  for u = 1:numel(inds)
    unwrappedPhase = recentrePhase(areaPhaseFOIindividual{iCond}{iArea}(inds(u),:),...
      mean(phaseMeanROI(~isnan(phaseMeanROI))));
    if ~isempty(unwrappedPhase) && sum(~isnan(unwrappedPhase))
      semilogx(areaFOI, adjustPhase(unwrappedPhase, areaFOI, [], options.cutoffFreq), 'Color',[200 200 200]./255);
    end
  end
  
  drawMean = false;
  if sum(~isnan(phaseMean) & ~isnan(phaseCI95(1,:)) & ~isnan(phaseCI95(2,:))) && drawMean
    % Plot the mean curve
    if iCond == 1
      p1 = semilogx(areaFOI(~isnan(phaseMean) & ~isnan(phaseCI95(1,:)) & ~isnan(phaseCI95(2,:))),...
        phaseMean(~isnan(phaseMean)  & ~isnan(phaseCI95(1,:)) & ~isnan(phaseCI95(2,:))),...
        'g:', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Unit mean
      plotCondMean = [plotCondMean p1]; %#ok<*AGROW>
      txtCondMean{numel(txtCondMean)+1} = 'Mean';
    elseif iCond == 2
      p1 = semilogx(areaFOI(~isnan(phaseMean)  & ~isnan(phaseCI95(1,:)) & ~isnan(phaseCI95(2,:))),...
        phaseMean(~isnan(phaseMean)  & ~isnan(phaseCI95(1,:)) & ~isnan(phaseCI95(2,:))),...
        'r:', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Unit mean
      plotCondMean = [plotCondMean p1];
      txtCondMean{numel(txtCondMean)+1} = 'Mean';
      uistack(p1,'bottom');
    end
    
    % Plot 95% confidence intervals
    if iCond == 1
      pC1 = semilogx(areaFOI(~isnan(phaseMean)  & ~isnan(phaseCI95(1,:))),...
        phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95(1,:))) + recentrePhase(bestUnwrap(phaseCI95(1,~isnan(phaseMean) & ~isnan(phaseCI95(1,:)))),0),...
        'g:', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Upper 95% confidence limitof the mean
      uistack(pC1,'bottom');
      pC2 = semilogx(areaFOI(~isnan(phaseMean) & ~isnan(phaseCI95(2,:))),...
        phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95(2,:))) + recentrePhase(bestUnwrap(phaseCI95(2,~isnan(phaseMean) & ~isnan(phaseCI95(2,:)))),0),...
        'g:', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Lower 95% confidence limitof the mean
      uistack(pC2,'bottom');
    elseif iCond == 2
      pC1 = semilogx(areaFOI(~isnan(phaseMean) & ~isnan(phaseCI95(1,:))),...
        phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95(1,:))) + recentrePhase(bestUnwrap(phaseCI95(1,~isnan(phaseMean) & ~isnan(phaseCI95(1,:)))),0),...
        'r:', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Upper 95% confidence limitof the mean
      uistack(pC1,'bottom');
      pC2 = semilogx(areaFOI(~isnan(phaseMean) & ~isnan(phaseCI95(2,:))),...
        phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95(2,:))) + recentrePhase(bestUnwrap(phaseCI95(2,~isnan(phaseMean) & ~isnan(phaseCI95(2,:)))),0),...
        'r:', 'LineWidth',1, 'MarkerSize',1, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Lower 95% confidence limitof the mean
      uistack(pC2,'bottom');
    end
    plotCondMean = [plotCondMean pC1];
    txtCondMean{numel(txtCondMean)+1} = '95% CI';
  end
  hold off
  
  % Tidy the figure
  yLim = options.phaseLim;
  if yLim(1) > min([phaseMeanROI phaseMeanROI-phaseCI95ROI])
    yLim(1) = min([phaseMean phaseMeanROI-phaseCI95ROI]);
  end
  if yLim(2) < max([phaseMeanROI phaseMeanROI+phaseCI95ROI])
    yLim(2) = max([phaseMeanROI phaseMeanROI+phaseCI95ROI]);
  end
  if ~isempty(plotCondMean) && ~isempty(txtCondMean) && numel(plotCondMean) == numel(txtCondMean)
    legend(plotCondMean,txtCondMean, 'Interpreter','None');
    legend boxoff
  end
  if isempty(animalID) && isempty(seriesID)
    if isempty(options.figTitle)
      titleStr = sprintf('Individual phases: %s %s', areas{iArea}, conditions{iCond});
    else
      titleStr = sprintf(options.figTitle, areas{iArea}, conditions{iCond});
    end
  elseif ~isempty(animalID) && isempty(seriesID)
    if isempty(options.figTitle)
      titleStr = sprintf('Individual phases: %s %s %s', areas{iArea}, conditions{iCond}, animalID);
    else
      titleStr = sprintf(options.figTitle, areas{iArea}, conditions{iCond}, animalID);
    end
  elseif ~isempty(animalID) && ~isempty(seriesID)
    if isempty(options.figTitle)
      titleStr = sprintf('Individual phases: %s %s %s %s', areas{iArea}, conditions{iCond}, animalID, seriesID);
    else
      titleStr = sprintf(options.figTitle, areas{iArea}, conditions{iCond}, animalID, seriesID);
    end
  end
  titleStr = strrep(titleStr, 'Vs', 'wrt');
  set(fH,'color','w');
  ax1 = axesProperties({titleStr}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
    'on', 'k', {'Frequency (Hz)'}, options.freqLim, [0.001 0.01 0.1 1 10 30 100 1000],...
    'on', 'k', {'Phase (rad)'}, yLim, [-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi]);
  ax1.XTickLabel = {'0.001','0.01','0.1','1','10','30','100','1000'};
  ax1.YTickLabel = {'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'};
  
  % Save the figure
  if isempty(animalID) && isempty(seriesID)
    figFileNameEnd = sprintf(figFileName, areas{iArea}, conditions{iCond});
  elseif ~isempty(animalID) && isempty(seriesID)
    figFileNameEnd = sprintf(figFileName, areas{iArea}, conditions{iCond}, animalID);
  elseif ~isempty(animalID) && ~isempty(seriesID)
    figFileNameEnd = sprintf(figFileName, areas{iArea}, conditions{iCond}, animalID, seriesID);
  end
  figFileName = [options.mainFolder filesep options.phaseFrequencyProfilesSubfolder filesep figFileNameEnd];
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