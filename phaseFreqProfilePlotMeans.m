function phaseFreqProfilePlotMeans(areas, conditions, FOI, areaPhaseFOIindividual, areaPhaseFOIindividualPR, figFileName, options)
% phaseFreqProfilePlotMeans(areas, conditions, FOI, areaPhaseFOIindividual, areaPhaseFOIindividualPR, figFileName, options)
%
% Function produces figures containing mean phase frequency profiles of
% individual units or individual population rates or both.
% Input: areas - area acronyms or their comparisons.
%        conditions - recording conditions corresponding to states of
%                     vigilance.
%        FOI - frequencies of interest. They could be a single vector or a
%              cell array containing multiple vectors.
%        areaPhaseFOIindividual - a cell array of phase frequency profiles.
%        areaPhaseFOIindividualPR - a cell array of phase frequency
%                                   profiles for population rates. Can be
%                                   left empty if not needed. Useful when
%                                   displaying unit means to have PR means
%                                   in comparison.
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
if nargin < 7
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
    if ~isempty(areaPhaseFOIindividual{iCond}{options.iAreasOI(iArea)})
      
      % Get phase values
      if iscell(FOI)
        areaFOI = FOI{iCond}{options.iAreasOI(iArea)};
      else
        areaFOI = FOI;
      end
      [phaseMean, phaseCI95] = datamean(areaPhaseFOIindividual{iCond}{options.iAreasOI(iArea)}, 'circularNP');
%       phaseCI95(phaseCI95 == 0) = NaN;
%       phaseMean(isnan(phaseCI95(1,:)) | isnan(phaseCI95(2,:))) = NaN;
      if size(areaPhaseFOIindividual{iCond}{options.iAreasOI(iArea)},1) == 1
        phaseCI95 = repmat([-pi; pi], 1, numel(phaseMean));
      end
      if ~isempty(areaPhaseFOIindividualPR)
        [phaseMean_PR, phaseCI95_PR] = datamean(areaPhaseFOIindividualPR{iCond}{options.iAreasOI(iArea)}, 'circularNP');
%         phaseCI95_PR(phaseCI95_PR == 0) = NaN;
%         phaseMean_PR(isnan(phaseCI95_PR(1,:)) | isnan(phaseCI95_PR(2,:))) = NaN;
      end
      if ~isempty(areaPhaseFOIindividualPR) && size(areaPhaseFOIindividualPR{iCond}{options.iAreasOI(iArea)},1) == 1
        phaseCI95_PR = repmat([-pi; pi], 1, numel(phaseMean_PR));
      end
      
      % Adjust phase values
      if ~isempty(phaseMean) && sum(~isnan(phaseMean))
        [phaseMean, phaseMeanROI, phaseCI95ROI] = adjustPhase(phaseMean, areaFOI, phaseCI95, options.cutoffFreq);
      else
        continue
      end
      if ~isempty(areaPhaseFOIindividualPR) && ~isempty(phaseMean_PR) && sum(~isnan(phaseMean_PR))
        [phaseMean_PR, phaseMean_PRROI, phaseCI95_PRROI] = adjustPhase(phaseMean_PR, areaFOI, phaseCI95_PR, options.freqLim(2));
      else
        phaseMean_PR = []; phaseMean_PRROI = []; phaseCI95_PRROI = [];
      end
      
      % Plot mean curves
      figure(fH{options.iAreasOI(iArea)}); hold on
      if iCond == 1
        p1 = semilogx(areaFOI(~isnan(phaseMean) & ~isnan(phaseCI95(1,:))),...
          phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95(1,:))),...
          'g:', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','g', 'MarkerFaceColor','g'); % Unit mean
        txtCondMean{options.iAreasOI(iArea)}{numel(txtCondMean{options.iAreasOI(iArea)})+1} = 'Wakefulness mean';
        if ~isempty(areaPhaseFOIindividualPR)
          p2 = semilogx(areaFOI(~isnan(phaseMean_PR) & ~isnan(phaseCI95_PR(1,:))),...
            phaseMean_PR(~isnan(phaseMean_PR) & ~isnan(phaseCI95_PR(1,:))),...
            'g:', 'LineWidth',1); % Population mean
          txtCondMean{options.iAreasOI(iArea)}{numel(txtCondMean{options.iAreasOI(iArea)})+1} = 'Wakefulness PR mean';
        end
      elseif iCond == 2
        p1 = semilogx(areaFOI(~isnan(phaseMean) & ~isnan(phaseCI95(1,:))),...
          phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95(1,:))),...
          'r:', 'LineWidth',2, 'MarkerSize',5, 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); % Unit mean
        txtCondMean{options.iAreasOI(iArea)}{numel(txtCondMean{options.iAreasOI(iArea)})+1} = 'Anaesthesia mean';
        uistack(p1,'bottom');
        if ~isempty(areaPhaseFOIindividualPR)
          p2 = semilogx(areaFOI(~isnan(phaseMean_PR) & ~isnan(phaseCI95_PR(1,:))),...
            phaseMean_PR(~isnan(phaseMean_PR) & ~isnan(phaseCI95_PR(1,:))),...
            'r:', 'LineWidth',1); % Population mean
          txtCondMean{options.iAreasOI(iArea)}{numel(txtCondMean{options.iAreasOI(iArea)})+1} = 'Anaesthesia PR mean';
        end
      end
      if ~isempty(areaPhaseFOIindividualPR)
        plotCondMean{options.iAreasOI(iArea)} = [plotCondMean{options.iAreasOI(iArea)} p1 p2];
        uistack(p2,'bottom');
      else
        plotCondMean{options.iAreasOI(iArea)} = [plotCondMean{options.iAreasOI(iArea)} p1];
      end
      
      % Plot 95% confidence intervals
      if iCond == 1
        pC = ciplot(phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95(1,:))) + phaseCI95(1,~isnan(phaseMean) & ~isnan(phaseCI95(1,:))),...
          phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95(2,:))) + phaseCI95(2,~isnan(phaseMean) & ~isnan(phaseCI95(2,:))),...
          areaFOI(~isnan(phaseMean) & ~isnan(phaseCI95(1,:))), 'g', 0.2); % Non-parametric
      elseif iCond == 2
        pC = ciplot(phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95(1,:))) + phaseCI95(1,~isnan(phaseMean) & ~isnan(phaseCI95(1,:))),...
          phaseMean(~isnan(phaseMean) & ~isnan(phaseCI95(2,:))) + phaseCI95(2,~isnan(phaseMean) & ~isnan(phaseCI95(2,:))),...
          areaFOI(~isnan(phaseMean) & ~isnan(phaseCI95(1,:))), 'r', 0.1); % Non-parametric
      end
      uistack(pC,'bottom');
      hold off
    end
    
    % Tidy and save the figure
    if (numel(conditions) == 1 && ~isempty(areaPhaseFOIindividual{1}{options.iAreasOI(iArea)})) ||...
        (iCond == 2 && (~isempty(areaPhaseFOIindividual{1}{options.iAreasOI(iArea)}) || ~isempty(areaPhaseFOIindividual{2}{options.iAreasOI(iArea)})))
      yLim = options.phaseLim;
      if yLim(1) > min([phaseMeanROI phaseMeanROI-phaseCI95ROI phaseMean_PRROI phaseMean_PRROI-phaseCI95_PRROI])
        yLim(1) = min([phaseMean phaseMeanROI-phaseCI95ROI phaseMean_PRROI phaseMean_PRROI-phaseCI95_PRROI]);
      end
      if yLim(2) < max([phaseMeanROI phaseMeanROI-phaseCI95ROI phaseMean_PRROI phaseMean_PRROI-phaseCI95_PRROI])
        yLim(2) = max([phaseMeanROI phaseMeanROI-phaseCI95ROI phaseMean_PRROI phaseMean_PRROI-phaseCI95_PRROI]);
      end
      if ~isempty(plotCondMean{options.iAreasOI(iArea)}) && ~isempty(txtCondMean{options.iAreasOI(iArea)}) && numel(plotCondMean{options.iAreasOI(iArea)}) == numel(txtCondMean{options.iAreasOI(iArea)})
        legend(plotCondMean{options.iAreasOI(iArea)},txtCondMean{options.iAreasOI(iArea)}, 'Interpreter','None');
        legend boxoff
      end
      if isempty(options.figTitle)
        titleStr = sprintf('Mean phase comparisons: %s', areas{options.iAreasOI(iArea)});
      else
        titleStr = sprintf(options.figTitle, areas{options.iAreasOI(iArea)});
      end
      set(fH{options.iAreasOI(iArea)},'color','w');
      if yLim <= pi/2
        yTicks = -pi/4:pi/8:pi/4;
      elseif yLim <= pi
        yTicks = -pi/2:pi/4:pi/2;
      else
        yTicks = -2*pi:pi/2:2*pi;
      end
      ax1 = axesProperties({titleStr}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
        'on', 'k', {'Frequency (Hz)'}, options.freqLim, [0.001 0.01 0.1 1 10 30 100 1000],...
        'on', 'k', {'Phase (rad)'}, yLim, yTicks);
      ax1.XTickLabel = {'0.001','0.01','0.1','1','10','30','100','1000'};
      if yLim <= pi/2
        ax1.YTickLabel = {'-\pi/4','-\pi/8','0','\pi/8','\pi/4'};
      elseif yLim <= pi
        ax1.YTickLabel = {'-\pi/2','-\pi/4','0','\pi/4','\pi/2'};
      else
        ax1.YTickLabel = {'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'};
      end
      
      figFileNameEnd = sprintf(figFileName, areas{options.iAreasOI(iArea)});
      figFileNameFull = [options.mainFolder filesep options.phaseFrequencyProfilesSubfolder filesep figFileNameEnd];
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