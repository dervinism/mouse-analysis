function [rFOI, pvalFOI, nFOI] = halfPhaseCorr(areas, conditions, FOI, areaPhase1FOIindividual, areaPhase2FOIindividual, options)
% [rFOI, pvalFOI, nFOI] = halfPhaseCorr(areas, conditions, FOI, areaPhase1FOIindividual, areaPhase2FOIindividual, options)
%
% Function runs half-recording phase correlation analyses and displays
% correlation figures.
% Input: areas - area acronyms or their comparisons.
%        conditions - recording conditions corresponding to states of
%                     vigilance.
%        FOI - frequencies of interest. They could be a single vector or a
%              cell array containing multiple vectors.
%        areaPhase1FOIindividual - a cell array of phase frequency profiles
%                                 fot the 1st half of the recording.
%        areaPhase2FOIindividual - a cell array of phase frequency profiles
%                                 fot the 2nd half of the recording.
%        options - a structure with the following fields:
%          figFolder
%          figSize
%          figTitle
%          xLabel
%          yLabel
%          phaseLim
%          individualGraphs
%          summaryGraphs
%          iAreasOI.
% Output: rFOI - circular correlation coefficient for every FOI.
%         pvalFOI - corresponding p-values.
%         nFOI - significant counts structure variable with fields:
%           bothHalvesSignificant
%           oneHalfSignificant
%           total


%% Parse user input
if nargin < 6
  options.figFolder = pwd;
  options.figSize = 15;
  options.figTitle = 'SPIKING';
  options.xLabel = '1st half';
  options.yLabel = '2nd half';
  options.phaseLim = [-pi pi];
  options.individualGraphs = true;
  options.summaryGraphs = true;
  options.iAreasOI = iAreasOI;
else
  if ~isfield(options,'figFolder')
    options.mainFolder = pwd;
  end
  if ~isfield(options,'figSize')
    options.figSize = 15;
  end
  if ~isfield(options,'figTitle')
    options.figTitle = 'SPIKING';
  end
  if ~isfield(options,'xLabel')
    options.xLabel = '1st half';
  end
  if ~isfield(options,'yLabel')
    options.yLabel = '2nd half';
  end
  if ~isfield(options,'phaseLim')
    options.phaseLim = [-pi pi];
  end
  if ~isfield(options,'individualGraphs')
    options.individualGraphs = true;
  end
  if ~isfield(options,'summaryGraphs')
    options.summaryGraphs = true;
  end
  if ~isfield(options,'iAreasOI')
    options.iAreasOI = iAreasOI;
  end
end

if ~exist(options.figFolder, 'file')
  mkdir(options.figFolder);
end


%% Display the figures
for iCond = 1:min([2 numel(conditions)])
  for iArea = 1:numel(options.iAreasOI)
    disp(['Processing half-recording data for ' conditions{iCond} ' ' areas{options.iAreasOI(iArea)}...
      ' (comparison # ' num2str((iCond-1)*numel(areas) + options.iAreasOI(iArea)) '/' num2str(numel(conditions)*numel(areas)) ')']);
    
    if iscell(FOI)
      areaFOI = FOI{iCond}{options.iAreasOI(iArea)};
    else
      areaFOI = FOI;
    end
    phase1 = recentrePhase(areaPhase1FOIindividual{iCond}{options.iAreasOI(iArea)},mean(options.phaseLim));
    phase2 = recentrePhase(areaPhase2FOIindividual{iCond}{options.iAreasOI(iArea)},mean(options.phaseLim));
    if ~isempty(phase1) && ~isempty(phase2)
      [areaFOI, sortInds] = sort(areaFOI, 'descend'); % Sort in descending order as the function was originally coded to work with descending frequencies
      reverseSortInds = 1:numel(areaFOI);
      reverseSortInds = reverseSortInds(sortInds);
      phase1 = phase1(:,sortInds);
      phase2 = phase2(:,sortInds);
      if ~isempty(phase1) && ~isempty(phase2)
        if strcmp(options.figTitle, 'SPIKING')
          suffix = '';
        elseif strcmp(options.figTitle, 'PUPIL')
          suffix = 'VsPupil';
        elseif strcmp(options.figTitle, 'MOTION')
          suffix = 'VsMotion';
        end
        
        % Correlations and individual graphs
        [~, rFOI{iCond}{options.iAreasOI(iArea)}, pvalFOI{iCond}{options.iAreasOI(iArea)}, nFOI{iCond}{options.iAreasOI(iArea)}] = halfPhaseCorrPlot_ca(phase1, phase2,...
          areaFOI, conditions{iCond}, areas{options.iAreasOI(iArea)}, options.figFolder,...
          'Half-recording phase correlations for ', options.individualGraphs, options.figTitle, suffix, options); %#ok<*AGROW>
        close all
        
        % Summary subplots
        if options.summaryGraphs
          sbPhase = halfPhaseCorrSummary(phase1, phase2, areaFOI, rFOI{iCond}{options.iAreasOI(iArea)}, pvalFOI{iCond}{options.iAreasOI(iArea)},...
            nFOI{iCond}{options.iAreasOI(iArea)}, conditions{iCond}, areas{options.iAreasOI(iArea)}, options.figFolder,...
            '_half_recording_phase_correlations_at_p01_p03_p05_p1_p3_p5_1_4_10_Hz', options.figTitle, suffix, options.phaseLim, options);
          close(sbPhase)
        end
        
        % Reverse indexing
        rFOI{iCond}{options.iAreasOI(iArea)} = rFOI{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
        pvalFOI{iCond}{options.iAreasOI(iArea)} = pvalFOI{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
        nFOI{iCond}{options.iAreasOI(iArea)} = nFOI{iCond}{options.iAreasOI(iArea)}(reverseSortInds);
      end
    end
  end
end