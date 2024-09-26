function [figPhase, rFOI, pvalFOI, nFOI] = halfPhaseCorrPlot_ca(phase1, phase2, FOI, condition, area, figPath, nameStr, drawFigures, prefix, suffix, options)
% [figPhase, rFOI, pvalFOI, nFOI] = halfPhaseCorrPlot_ca(phase1, phase2, FOI, condition, area, figPath, nameStr, drawFigures, prefix, suffix, options)
%
% Function displays correlations between recording half phases.
% Input: phase1 - phase sample of the first half of recordings.
%        phase2 - phase sample of the second half of recordings.
%        FOI - a vector of frequencies of interest.
%        condition - a condition string.
%        area - an area string.
%        figPath - a folder for saving figures. If left empty, will save in
%                  the present working directory.
%        nameStr - part of the figure name string that would uniquely
%                  identify figures.
%        drawFigures - if true, will display and save correlation figures.
%                      Otherwise would only carry out the correlation
%                      analyses. Default is true.
%        prefix and suffix of the figure name.
%        options.
% Output: figPhase - figure handles.
%         rFOI - correlation coefficients.
%         pvalFOI - corresponding p-values;
%         nFOI - recording counts.

if nargin < 11
  options.phaseLim = [-pi pi];
  options.xLabel = '1st half';
  options.yLabel = '2nd half';
  options.figSize = 15;
else
  if ~isfield(options, 'phaseLim')
    options.phaseLim = [-pi pi];
  end
  if ~isfield(options, 'xLabel')
    options.xLabel = '1st half';
  end
  if ~isfield(options, 'yLabel')
    options.yLabel = '2nd half';
  end
  if ~isfield(options, 'figSize')
    options.figSize = 15;
  end
end
if nargin < 10
  suffix = '';
end
if nargin < 9
  prefix = 'SPIKING';
end
if nargin < 8 || isempty(drawFigures)
  drawFigures = true;
end
if isempty(figPath)
  figPath = pwd;
end

% Correlation analyses
[rFOI, pvalFOI] = corrMulti(phase1', phase2', 'circular');

% Produce figures
if drawFigures
  figPhase = phaseVphasePlot(phase1, phase2, FOI, options.xLabel, options.yLabel, 'on');
else
  figPhase = [];
end

for f = 1:numel(FOI) % Loop over frequencies of interest
  if FOI(f) == 0.01 || FOI(f) == 0.03 || FOI(f) == 0.05 || FOI(f) == 0.1 || FOI(f) == 0.3 ||...
      FOI(f) == 0.5 || FOI(f) == 1 || FOI(f) == 4 || FOI(f) == 10 || f == numel(FOI)
    
    % Count units
    nFOI(f).bothHalvesSignificant = sum(~isnan(phase1(:,f)) & ~isnan(phase2(:,f))); %#ok<*AGROW>
    nFOI(f).oneHalfSignificant = sum(~isnan(phase1(:,f)) | ~isnan(phase2(:,f)));
    nFOI(f).total = numel(phase1(:,f));
    
    % Tidy the figures
    if drawFigures
      figure(figPhase(f));
      figTitle = [nameStr area ' ' condition ' and ' num2str(FOI(f)) ' Hz: r=' num2str(rFOI(f))...
        ' p=' num2str(pvalFOI(f)) ' n=' num2str(nFOI(f).bothHalvesSignificant) '/' num2str(nFOI(f).oneHalfSignificant) '/'...
        num2str(nFOI(f).total)];
      ax1 = axesProperties({figTitle}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
        'on', 'k', {get(get(gca,'xlabel'),'string')}, options.phaseLim, [-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi],...
        'on', 'k', {get(get(gca,'ylabel'),'string')}, options.phaseLim, [-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi]);
      ax1.XTickLabel = {'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'};
      ax1.YTickLabel = {'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'};
      set(figPhase(f), 'Name',figTitle);
      figFileName = [prefix '_' lower(nameStr(1)) nameStr(2:end) area suffix '_' condition '_' num2str(FOI(f)) '_Hz'];
      figFileName = strrep(figFileName, ' ', '_');
      label = [2 1.6];
      margin = [0.3 0.55];
      width = 1*options.figSize-label(1)-margin(1);
      height = 1*options.figSize-label(2)-margin(2);
      paperSize = resizeFig(figPhase(f), ax1, width, height, label, margin, 0);
      exportFig(figPhase(f), [figPath filesep figFileName '.png'],'-dpng','-r300', paperSize);
      hgsave(figPhase(f), [figPath filesep figFileName '.fig']);
    end
  end
end