function [figPhase, rFOI, pvalFOI, nFOI] = halfCorrPlot_ca(phase1, phase2, phase, FOI, condition, area, figPath, nameStr, drawFigures, prefix, suffix)
% [figPhase, rFOI, pvalFOI, nFOI] = halfCorrPlot_ca(phase1, phase2, phase, FOI, condition, area, figPath, nameStr, drawFigures)
%
% Function displays correlations between recording half phases.
% Input: phase1 - phase sample of the first half of recordings.
%        phase2 - phase sample of the second half of recordings.
%        phase - phase sample of full length recordings.
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
% Output: figPhase - figure handles.
%         rFOI - correlation coefficients.
%         pvalFOI - corresponding p-values;
%         nFOI - recording counts.

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
[rFOI, pvalFOI] = corrMulti(phase1', phase2', 'circular');

% Produce figures
if drawFigures
  figPhase = phaseVphasePlot(phase1, phase2, FOI, '1st half', '2nd half', 'on');
else
  figPhase = [];
end

for f = 1:numel(FOI) % Loop over frequencies of interest
  
  % Count units
  nFOI(f).bothHalvesSignificant = sum(~isnan(phase1(:,f)) & ~isnan(phase2(:,f))); %#ok<*AGROW>
  nFOI(f).oneHalfSignificant = sum(~isnan(phase1(:,f)) | ~isnan(phase2(:,f)));
  nFOI(f).fullSignificant = sum(~isnan(phase(:,f)));
  nFOI(f).total = numel(phase(:,f));
  
  % Tidy the figures
  if drawFigures
    figure(figPhase(f));
    figTitle = [nameStr area ' ' condition ' and ' num2str(FOI(f)) ' Hz: r=' num2str(rFOI(f))...
      ' p=' num2str(pvalFOI(f)) ' n=' num2str(nFOI(f).bothHalvesSignificant) '/' num2str(nFOI(f).oneHalfSignificant) '/'...
      num2str(nFOI(f).fullSignificant) '/' num2str(nFOI(f).total)];
    title(figTitle);
    set(gcf, 'Name',figTitle);
    figName = [prefix '_' lower(nameStr(1)) nameStr(2:end) area suffix '_' condition '_' num2str(FOI(f)) '_Hz'];
    figName = strrep(figName, ' ', '_');
    hgsave(gcf, [figPath filesep figName '.fig']);
  end
end