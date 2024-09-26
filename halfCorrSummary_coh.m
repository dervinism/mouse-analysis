function sbCoh = halfCorrSummary_coh(coh1, coh2, FOI, rFOI, pvalFOI, nFOI, condition, area, comparison, figPath, nameStr, prefix, suffix)
% sbCoh = halfCorrSummary_coh(coh1, coh2, FOI, rFOI, pvalFOI, nFOI, condition, area, figPath, nameStr, prefix, suffix, adjustAxes)
%
% Function produces summary subplots of correlations between coherences.
% Input: coh1 - the first instance of coherence.
%        coh2 - the second instance of coherence.
%        FOI - a vector of frequencies of interest.
%        rFOI - correlation coefficients.
%        pvalFOI - corresponding p-values.
%        nFOI - recording counts.
%        condition - a condition string.
%        area - an area string.
%        figPath - a folder for saving figures. If left empty, will save in
%                  the present working directory.
%        nameStr - part of the figure name string that would uniquely
%                  identify figures.
%        prefix and suffix of the figure name.
% Output: sbCoh - a figure handle.

if nargin < 12
  suffix = '';
end
if nargin < 11
  prefix = 'SPIKING';
end
if isempty(figPath)
  figPath = pwd;
end

if ~exist(figPath, 'file')
  mkdir(figPath);
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
      if strcmp(prefix, 'MIXED')
        cohVcohSubplot(coh1(:,f), coh2(:,f), 'Spiking coherence', '', rFOI(f), pvalFOI(f), nFOI(f));
      else
        cohVcohSubplot(coh1(:,f), coh2(:,f), ['Unit coherence with ' area ' PR'], '', rFOI(f), pvalFOI(f), nFOI(f));
      end
      dispFOI(FOI(f));
    elseif plotCount == 3
      if strcmp(prefix, 'MIXED')
        cohVcohSubplot(coh1(:,f), coh2(:,f), 'Spiking coherence', '', rFOI(f), pvalFOI(f), nFOI(f));
      else
        cohVcohSubplot(coh1(:,f), coh2(:,f), ['Unit coherence with ' area ' PR'], '', rFOI(f), pvalFOI(f), nFOI(f));
      end
      dispFOI(FOI(f));
    elseif plotCount <= 5
      cohVcohSubplot(coh1(:,f), coh2(:,f), '', '', rFOI(f), pvalFOI(f), nFOI(f));
      dispFOI(FOI(f));
    elseif plotCount == 6
      if strcmp(prefix, 'MIXED')
        cohVcohSubplot(coh1(:,f), coh2(:,f), '', 'Pupil coherence', rFOI(f), pvalFOI(f), nFOI(f));
      else
        cohVcohSubplot(coh1(:,f), coh2(:,f), '', ['Unit coherence: ' comparison], rFOI(f), pvalFOI(f), nFOI(f));
      end
      dispFOI(FOI(f));
    elseif strcmp(prefix, 'SPIKING') && plotCount <= 8
      cohVcohSubplot(coh1(:,f), coh2(:,f), '', '', rFOI(f), pvalFOI(f), nFOI(f))
      dispFOI(FOI(f));
    elseif plotCount == 9 || ((strcmp(prefix, 'PUPIL') || strcmp(prefix, 'MOTION') || strcmp(prefix, 'MIXED')) && plotCount == 7)
      cohVcohSubplot(coh1(:,f), coh2(:,f), '', '', rFOI(f), pvalFOI(f), nFOI(f))
      dispFOI(FOI(f));
      set(sbCoh, 'Color','w');
      print(sbCoh, [figPath filesep prefix '_' area suffix '_v_' comparison '_' condition nameStr '.png'], '-dpng', '-r300');
      hgsave(sbCoh, [figPath filesep prefix '_' area suffix '_v_' comparison '_' condition nameStr '.fig'])
    end
  end
end