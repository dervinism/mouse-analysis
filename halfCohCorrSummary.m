function sbCoh = halfCohCorrSummary(coh1, coh2, FOI, rFOI, pvalFOI, nFOI, condition, area, figPath, nameStr, prefix, suffix, cohLim, options)
% sbCoh = halfCohCorrSummary(coh1, coh2, FOI, rFOI, pvalFOI, nFOI, condition, area, figPath, nameStr, prefix, suffix, cohLim, options)
%
% Function produces summary subplots of correlations between recording half
% coherences.
% Input: coh1 - coherence sample of the first half of recordings.
%        coh2 - coherence sample of the second half of recordings.
%        FOI - a vector of frequencies of interest.
%        rFOI - correlation coefficients.
%        pvalFOI - corresponding p-values.
%        nFOI - recording counts.
%        condition - a cell array of condition strings.
%        area - a cell array of area strings.
%        figPath - a folder for saving figures. If left empty, will save in
%                  the present working directory.
%        nameStr - part of the figure name string that would uniquely
%                  identify figures.
%        prefix and suffix of the figure name.
%        cohLim.
%        options.
% Output: sbCoh - a figure handle.

if nargin < 13
  cohLim = [0 1];
end
if nargin < 12
  suffix = '';
end
if nargin < 11
  prefix = 'SPIKING';
end
if isempty(figPath)
  figPath = pwd;
end

plotCount = 0;
sbCoh = figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible', 'on');
for f = 1:numel(FOI)
  if plotCount == 9
    break
  end
  if FOI(f) == 0.01 || FOI(f) == 0.03 || FOI(f) == 0.05 || FOI(f) == 0.1 || FOI(f) == 0.3 ||...
      FOI(f) == 0.5 || FOI(f) == 1 || FOI(f) == 4 || FOI(f) == 10
    figure(sbCoh);
    plotCount = plotCount + 1;
    subplot(3, 3, 10 - plotCount);
    if plotCount <= 2
      cohVcohSubplot(coh1(:,f), coh2(:,f), options.xLabel, '', rFOI(f), pvalFOI(f), nFOI(f), cohLim);
      dispFOI(FOI(f));
    elseif plotCount == 3
      cohVcohSubplot(coh1(:,f), coh2(:,f), options.xLabel, '', rFOI(f), pvalFOI(f), nFOI(f), cohLim);
      dispFOI(FOI(f));
    elseif plotCount <= 5
      cohVcohSubplot(coh1(:,f), coh2(:,f), '', '', rFOI(f), pvalFOI(f), nFOI(f), cohLim);
      dispFOI(FOI(f));
    elseif plotCount == 6
      cohVcohSubplot(coh1(:,f), coh2(:,f), '', options.yLabel, rFOI(f), pvalFOI(f), nFOI(f), cohLim);
      dispFOI(FOI(f));
    elseif strcmp(prefix, 'SPIKING') && plotCount <= 8
      cohVcohSubplot(coh1(:,f), coh2(:,f), '', '', rFOI(f), pvalFOI(f), nFOI(f), cohLim)
      dispFOI(FOI(f));
    elseif plotCount == 9 || ((strcmp(prefix, 'PUPIL') || strcmp(prefix, 'MOTION') || strcmp(prefix, 'MIXED')) && plotCount == 7)
      cohVcohSubplot(coh1(:,f), coh2(:,f), '', '', rFOI(f), pvalFOI(f), nFOI(f), cohLim)
      dispFOI(FOI(f));
      set(gcf, 'Color','w');
      print(sbCoh, [figPath filesep prefix '_' area suffix '_' condition nameStr '.png'], '-dpng', '-r300');
      hgsave([figPath filesep prefix '_' area suffix '_' condition nameStr '.fig'])
    end
  end
end