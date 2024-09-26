function sbPhase = halfPhaseCorrSummary(phase1, phase2, FOI, rFOI, pvalFOI, nFOI, condition, area, figPath, nameStr, prefix, suffix, phaseLim, options)
% sbPhase = halfPhaseCorrSummary(phase1, phase2, FOI, rFOI, pvalFOI, nFOI, condition, area, figPath, nameStr, prefix, suffix, phaseLim, options)
%
% Function produces summary subplots of correlations between recording half
% phases.
% Input: phase1 - phase sample of the first half of recordings.
%        phase2 - phase sample of the second half of recordings.
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
%        phaseLim.
%        options.
% Output: sbPhase - a figure handle.

if nargin < 14
  if strcmp(prefix, 'MIXED')
    options.xLabel = 'Spiking phase (rad)';
    options.yLabel = 'Pupil phase (rad)';
  else
    options.xLabel = '1st half phase (rad)';
    options.yLabel = '2nd half phase (rad)';
  end
end
if nargin < 13
  phaseLim = [-pi pi];
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
sbPhase = figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible', 'on');
for f = 1:numel(FOI)
  if plotCount == 9
    break
  end
  if FOI(f) == 0.01 || FOI(f) == 0.03 || FOI(f) == 0.05 || FOI(f) == 0.1 || FOI(f) == 0.3 ||...
      FOI(f) == 0.5 || FOI(f) == 1 || FOI(f) == 4 || FOI(f) == 10
    figure(sbPhase);
    plotCount = plotCount + 1;
    subplot(3, 3, 10 - plotCount);
    if plotCount <= 2
      phaseVphaseSubplot(phase1(:,f), phase2(:,f), options.xLabel, '', rFOI(f), pvalFOI(f), nFOI(f), phaseLim);
      dispFOI(FOI(f));
    elseif plotCount == 3
      phaseVphaseSubplot(phase1(:,f), phase2(:,f), options.xLabel, '', rFOI(f), pvalFOI(f), nFOI(f), phaseLim);
      dispFOI(FOI(f));
    elseif plotCount <= 5
      phaseVphaseSubplot(phase1(:,f), phase2(:,f), '', '', rFOI(f), pvalFOI(f), nFOI(f), phaseLim);
      dispFOI(FOI(f));
    elseif plotCount == 6
      phaseVphaseSubplot(phase1(:,f), phase2(:,f), '', options.yLabel, rFOI(f), pvalFOI(f), nFOI(f), phaseLim);
      dispFOI(FOI(f));
    elseif strcmp(prefix, 'SPIKING') && plotCount <= 8
      phaseVphaseSubplot(phase1(:,f), phase2(:,f), '', '', rFOI(f), pvalFOI(f), nFOI(f), phaseLim)
      dispFOI(FOI(f));
    elseif plotCount == 9 || ((strcmp(prefix, 'PUPIL') || strcmp(prefix, 'MOTION') || strcmp(prefix, 'MIXED')) && plotCount == 7)
      phaseVphaseSubplot(phase1(:,f), phase2(:,f), '', '', rFOI(f), pvalFOI(f), nFOI(f), phaseLim)
      dispFOI(FOI(f));
      set(sbPhase, 'Color','w');
      print(sbPhase, [figPath filesep prefix '_' area suffix '_' condition nameStr '.png'], '-dpng', '-r300');
      hgsave(sbPhase, [figPath filesep prefix '_' area suffix '_' condition nameStr '.fig'])
    end
  end
end