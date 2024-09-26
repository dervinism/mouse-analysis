function cohFreqProfilePlotMeans(areas, conditions, FOI, areaCohFOIindividual, areaCohFOIindividualPR, figFileName, options)
% cohFreqProfilePlotMeans(areas, conditions, FOI, areaCohFOIindividual, areaCohFOIindividualPR, figFileName, options)
%
% Function produces figures containing mean coherence frequency profiles of
% individual units or individual population rates or both.
% Input: areas - area acronyms or their comparisons.
%        conditions - recording conditions corresponding to states of
%                     vigilance.
%        FOI - frequencies of interest. They could be a single vector or a
%              cell array containing multiple vectors.
%        areaCohFOIindividual - a cell array of coherence frequency profiles.
%        areaCohFOIindividualPR - a cell array of coherence frequency
%                                 profiles for population rates. Can be
%                                 left empty if not needed. Useful when
%                                 displaying unit means to have PR means
%                                 in comparison.
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
if nargin < 7
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
  if ~isfield(options,'coherenceFrequencyProfilesSubfolder')
    options.coherenceFrequencyProfilesSubfolder = 'coherenceFrequencyProfiles';
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


%% Generate figures
% Initialise figures
for iArea = 1:numel(options.iAreasOI)
  fH{options.iAreasOI(iArea)} = figure();
  semilogx([10e-4 10e2],[0 0], 'k:');
  plotCondMean{options.iAreasOI(iArea)} = []; %#ok<*AGROW>
  txtCondMean{options.iAreasOI(iArea)} = {};
end

for iCond = 1:min([2 numel(conditions)]) % Loop through conditions
  for iArea = 1:numel(options.iAreasOI) % Loop through areas
    if ~isempty(areaCohFOIindividual{iCond}{options.iAreasOI(iArea)})
      
      % Get coherence values
      if iscell(FOI)
        areaFOI = FOI{iCond}{options.iAreasOI(iArea)};
      else
        areaFOI = FOI;
      end
      [cohMean, cohCI95] = datamean(areaCohFOIindividual{iCond}{options.iAreasOI(iArea)});
      if ~isempty(areaCohFOIindividualPR)
        [cohMean_PR, cohCI95_PR] = datamean(areaCohFOIindividualPR{iCond}{options.iAreasOI(iArea)});
      end
      
      % Plot mean curves
      figure(fH{options.iAreasOI(iArea)}); hold on
      if iCond == 1 && sum(~isnan(cohMean) & ~isnan(cohCI95(1,:)))
        %           p1 = semilogx(FOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
        %             cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
        %             'g:', 'LineWidth',2, 'MarkerSize',5,...
        %             'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Unit mean (parametric)
        p1 = semilogx(areaFOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
          cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
          'g:', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Unit mean
        txtCondMean{options.iAreasOI(iArea)}{numel(txtCondMean{options.iAreasOI(iArea)})+1} = 'Wakefulness mean';
        if ~isempty(areaCohFOIindividualPR) && sum(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:)))
          p2 = semilogx(areaFOI(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:))),...
            cohMean_PR(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:))),...
            'g:', 'LineWidth',1); % Population mean
          txtCondMean{options.iAreasOI(iArea)}{numel(txtCondMean{options.iAreasOI(iArea)})+1} = 'Wakefulness PR mean';
        end
      elseif iCond == 2 && sum(~isnan(cohMean) & ~isnan(cohCI95(1,:)))
        %           p1 = semilogx(areaFOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
        %             cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
        %             'r:', 'LineWidth',2, 'MarkerSize',5,...
        %             'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Unit mean (parametric)
        p1 = semilogx(areaFOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
          cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
          'r:', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Unit mean
        txtCondMean{options.iAreasOI(iArea)}{numel(txtCondMean{options.iAreasOI(iArea)})+1} = 'Anaesthesia mean';
        uistack(p1,'bottom');
        if ~isempty(areaCohFOIindividualPR) && sum(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:)))
          p2 = semilogx(areaFOI(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:))),...
            cohMean_PR(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:))),...
            'r:', 'LineWidth',1); % Population mean
          txtCondMean{options.iAreasOI(iArea)}{numel(txtCondMean{options.iAreasOI(iArea)})+1} = 'Anaesthesia PR mean';
        end
      end
      if ~isempty(areaCohFOIindividualPR) && sum(~isnan(cohMean) & ~isnan(cohCI95(1,:))) &&...
          sum(~isnan(cohMean_PR) & ~isnan(cohCI95_PR(1,:)))
        plotCondMean{options.iAreasOI(iArea)} = [plotCondMean{options.iAreasOI(iArea)} p1 p2];
        uistack(p2,'bottom');
      elseif sum(~isnan(cohMean) & ~isnan(cohCI95(1,:)))
        plotCondMean{options.iAreasOI(iArea)} = [plotCondMean{options.iAreasOI(iArea)} p1];
      end
      
      % Plot 95% confidence intervals
      if iCond == 1 && sum(~isnan(cohMean) & ~isnan(cohCI95(1,:)))
        %           pC = ciplot(cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))) + cohCI95(1,~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
        %             cohMean(~isnan(cohMean) & ~isnan(cohCI95(2,:))) + cohCI95(2,~isnan(cohMean) & ~isnan(cohCI95(2,:))),...
        %             areaFOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))), 'g', 0.1); % Parametric
        pC = ciplot(cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))) + cohCI95(1,~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
          cohMean(~isnan(cohMean) & ~isnan(cohCI95(2,:))) + cohCI95(2,~isnan(cohMean) & ~isnan(cohCI95(2,:))),...
          areaFOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))), 'g', 0.2); % Non-parametric
        uistack(pC,'bottom');
      elseif iCond == 2 && sum(~isnan(cohMean) & ~isnan(cohCI95(1,:)))
        %           pC = ciplot(cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))) + cohCI95(1,~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
        %             cohMean(~isnan(cohMean) & ~isnan(cohCI95(2,:))) + cohCI95(2,~isnan(cohMean) & ~isnan(cohCI95(2,:))),...
        %             areaFOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))), 'r', 0.05); % Parametric
        pC = ciplot(cohMean(~isnan(cohMean) & ~isnan(cohCI95(1,:))) + cohCI95(1,~isnan(cohMean) & ~isnan(cohCI95(1,:))),...
          cohMean(~isnan(cohMean) & ~isnan(cohCI95(2,:))) + cohCI95(2,~isnan(cohMean) & ~isnan(cohCI95(2,:))),...
          areaFOI(~isnan(cohMean) & ~isnan(cohCI95(1,:))), 'r', 0.1); % Non-parametric
        uistack(pC,'bottom');
      end
      hold off
    end
    
    % Tidy and save the figure
    if (numel(conditions) == 1 && ~isempty(areaCohFOIindividual{1}{options.iAreasOI(iArea)})) ||...
        (iCond == 2 && (~isempty(areaCohFOIindividual{1}{options.iAreasOI(iArea)}) || ~isempty(areaCohFOIindividual{2}{options.iAreasOI(iArea)})))
      yLim = ylim;
      yLim = [options.cohLim(1) min([options.cohLim(2) yLim(2)])];
      if ~isempty(plotCondMean{options.iAreasOI(iArea)}) && ~isempty(txtCondMean{options.iAreasOI(iArea)}) && numel(plotCondMean{options.iAreasOI(iArea)}) == numel(txtCondMean{options.iAreasOI(iArea)})
        legend(plotCondMean{options.iAreasOI(iArea)},txtCondMean{options.iAreasOI(iArea)}, 'Interpreter','None');
        legend boxoff
      end
      if isempty(options.figTitle)
        titleStr = sprintf('Mean coherence comparisons: %s', areas{options.iAreasOI(iArea)});
      else
        titleStr = sprintf(options.figTitle, areas{options.iAreasOI(iArea)});
      end
      titleStr = strrep(titleStr, 'Vs', 'wrt');
      set(fH{options.iAreasOI(iArea)},'color','w');
      ax1 = axesProperties({titleStr}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
        'on', 'k', {'Frequency (Hz)'}, options.freqLim, [0.001 0.01 0.1 1 10 30 100 1000],...
        'on', 'k', {'Coherence'}, yLim, 0:0.1:1);
      ax1.XTickLabel = {'0.001','0.01','0.1','1','10','30','100','1000'};
      
      figFileNameEnd = sprintf(figFileName, areas{options.iAreasOI(iArea)});
      figFileNameFull = [options.mainFolder filesep options.cohFrequencyProfilesSubfolder filesep figFileNameEnd];
      label = [2 1.6];
      margin = [0.3 0.55];
      width = 1*options.figSize-label(1)-margin(1);
      height = 1*options.figSize-label(2)-margin(2);
      paperSize = resizeFig(fH{options.iAreasOI(iArea)}, ax1, width, height, label, margin, 0);
      hgsave(fH{options.iAreasOI(iArea)}, [figFileNameFull '.fig']);
      exportFig(fH{options.iAreasOI(iArea)}, [figFileNameFull '.png'],'-dpng','-r300', paperSize);
    end
  end
end
close all