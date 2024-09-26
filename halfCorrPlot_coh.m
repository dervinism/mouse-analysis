function [figCoh, rFOI, pvalFOI, nFOI] = halfCorrPlot_coh(coh1, coh2, FOI, condition, area, comp, figPath, nameStr, drawFigures, prefix, suffix, opt)
% [figCoh, rFOI, pvalFOI, nFOI] = halfCorrPlot_ca(coh1, coh2, FOI, condition, area, comp, figPath, nameStr, drawFigures, prefix, suffix, opt)
%
% Function displays correlations between coherences.
% Input: coh1 - the first instance of coherence.
%        coh2 - the second instance of coherence.
%        FOI - a vector of frequencies of interest.
%        condition - a condition string.
%        area - an area string.
%        comp - an area comparison.
%        figPath - a folder for saving figures. If left empty, will save in
%                  the present working directory.
%        nameStr - part of the figure name string that would uniquely
%                  identify figures.
%        drawFigures - if true, will display and save correlation figures.
%                      Otherwise would only carry out the correlation
%                      analyses. Default is true.
%        prefix and suffix of the figure name.
%        opt - structure variable with fields xLabel and yLabel.
% Output: figCoh - figure handles.
%         rFOI - correlation coefficients.
%         pvalFOI - corresponding p-values;
%         nFOI - recording counts.

if nargin < 12
  opt.xLabel = 'Coherence 1st half';
  opt.yLabel = 'Coherence 2nd half';
else
  if ~isfield(opt, 'xLabel')
    opt.xLabel = 'Coherence 1st half';
  end
  if ~isfield(opt, 'yLabel')
    opt.yLabel = 'Coherence 2nd half';
  end
end
if nargin < 11
  suffix = '';
end
if nargin < 10
  prefix = 'SPIKING';
end
if nargin < 9 || isempty(drawFigures)
  drawFigures = true;
end
if isempty(figPath)
  figPath = pwd;
end

% Correlation analyses
[rFOI, pvalFOI] = corrMulti(coh1', coh2', 'Spearman');

% Produce figures
if drawFigures
  figCoh = cohVcohPlot(coh1, coh2, FOI, opt.xLabel, opt.yLabel, 'on');
else
  figCoh = [];
end

for f = 1:numel(FOI) % Loop over frequencies of interest
  
  % Count units
  nFOI(f).bothHalvesSignificant = sum(~isnan(coh1(:,f)) & ~isnan(coh2(:,f))); %#ok<*AGROW>
  nFOI(f).oneHalfSignificant = sum(~isnan(coh1(:,f)) | ~isnan(coh2(:,f)));
  nFOI(f).total = numel(coh1(:,f));
  
  % Tidy the figures
  if drawFigures
    figure(figCoh(f));
    
    xLim = xlim;
    xLim = [xLim(1)-0.001 max([xLim(2)+0.001 0.005])];
    xAxisLength = xLim(2)-xLim(1);
    yLim = ylim;
    yLim = [yLim(1)-0.001 max([yLim(2)+0.001 0.005])];
    yAxisLength = yLim(2)-yLim(1);
    text(xLim(2)-xAxisLength*0.25, yLim(1)+yAxisLength*0.2, ['r=' num2str(rFOI(f))], 'FontSize',16);
    text(xLim(2)-xAxisLength*0.25, yLim(1)+yAxisLength*0.1, ['p=' num2str(pvalFOI(f))], 'FontSize',16);
    
    xTicks = xticks;
    xTicks = xTicks(1:2:end);
    yTicks = yticks;
    yTicks = yTicks(1:2:end);

    figTitle = [nameStr area 'VsLocalPR v ' comp ' ' condition ' and ' num2str(FOI(f)) ' Hz: r=' num2str(rFOI(f))...
      ' p=' num2str(pvalFOI(f)) ' n=' num2str(nFOI(f).bothHalvesSignificant) '/' num2str(nFOI(f).oneHalfSignificant) '/'...
      num2str(nFOI(f).total)];
    axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
      'on', 'k', {get(get(gca,'xlabel'),'string')}, xLim, xTicks,...
      'on', 'k', {get(get(gca,'ylabel'),'string')}, yLim, yTicks);
    title(figTitle);
    set(figCoh(f), 'Name',figTitle);
    figName = [prefix '_' lower(nameStr(1)) nameStr(2:end) area 'VsLocalPR v ' comp...
      '_' suffix '_' condition '_' num2str(FOI(f)) '_Hz'];
    figName = strrep(figName, ' ', '_');
    hgsave(figCoh(f), [figPath filesep figName '.fig']);
  end
end