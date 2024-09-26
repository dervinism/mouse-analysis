function sbPhase = phaseLocationProfileSummary(phase, locations, FOI, r, pval, n, area, condition, figPath, prefix, suffix)
% sbPhase = phaseLocationProfileSummary(phase, locations, FOI, area, condition, figPath)
%
% Function produces a subplot summary figure with unit phase location
% profiles for 9 frequencies: 0.01, 0.02, 0.03, 0.1, 0.2, 1, 2, 4, and 10
% Hz.
% Input: phase - a matrix with unit (rows) phase values for every frequency
%                of interest (columns).
%        locations - unit probe location vector.
%        FOI - frequency of interest vector.
%        r - correlation coefficients for every frequency of interest.
%        pval - corresponding p-values.
%        n - count of units with significant phase values count/ total
%            number of units for every frequency of interest.
%        area - an area comparison string.
%        condition - a condition string.
%        figPath - a path to the folder for saving figures; default is
%                  present working directory.
%        prefix and suffix of figure file names.
% Output: sbPhase - the figure handle.

if nargin < 11
  suffix = '';
end

if nargin < 10
  prefix = 'SPIKING';
end

plotCount = 0;
sbPhase = figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible', 'on');
minLoc = min(locations);
maxLoc = max(locations);
yTickLabel = {'-\pi','-\pi/2','0','\pi/2','\pi'};
for f = 1:numel(FOI)
  if plotCount == 9
    break
  end
  [locations, sortOrder] = sort(locations);
  phaseFOI = phase(sortOrder,f);
  if FOI(f) == 0.01 || FOI(f) == 0.02 || FOI(f) == 0.03 || FOI(f) == 0.1 || FOI(f) == 0.2 ||...
      FOI(f) == 1 || FOI(f) == 2 || FOI(f) == 4 || FOI(f) == 10
    figure(sbPhase);
    plotCount = plotCount + 1;
    subplot(3, 3, 10 - plotCount);
    if plotCount <= 2
      xySubplot(locations, phaseFOI, 'Location (\mum)', [], [], [minLoc maxLoc], '', [], yTickLabel, [],...
        r(f), pval(f), n{f}, false, true, 'linear-circular');
      dispFOI(FOI(f));
    elseif plotCount == 3
      xySubplot(locations, phaseFOI, 'Location (\mum)', [], [], [minLoc maxLoc], '', [], yTickLabel, [],...
        r(f), pval(f), n{f}, false, true, 'linear-circular');
      dispFOI(FOI(f));
    elseif plotCount <= 5
      xySubplot(locations, phaseFOI, '', [], [], [minLoc maxLoc], '', [], yTickLabel, [],...
        r(f), pval(f), n{f}, false, true, 'linear-circular');
      dispFOI(FOI(f));
    elseif plotCount == 6
      xySubplot(locations, phaseFOI, '', [], [], [minLoc maxLoc], 'Unit phase (rad)', [], yTickLabel, [],...
        r(f), pval(f), n{f}, false, true, 'linear-circular');
      dispFOI(FOI(f));
    elseif strcmp(prefix, 'SPIKING') && plotCount <= 8
      xySubplot(locations, phaseFOI, '', [], [], [minLoc maxLoc], '', [], yTickLabel, [],...
        r(f), pval(f), n{f}, false, true, 'linear-circular');
      dispFOI(FOI(f));
    elseif plotCount == 9 || ((strcmp(prefix, 'PUPIL') || strcmp(prefix, 'MOTION')) && plotCount == 7)
      xySubplot(locations, phaseFOI, '', [], [], [minLoc maxLoc], '', [], yTickLabel, [],...
        r(f), pval(f), n{f}, false, true, 'linear-circular');
      dispFOI(FOI(f));
      set(gcf, 'Color','w');
      print(sbPhase, [figPath filesep [prefix '_' area suffix '_' condition...
        '_unit_location_and_phase_correlations_at_p01_p02_p03_p1_p3_1_2_4_10_Hz'] '.png'], '-dpng', '-r300');
      hgsave([figPath filesep [prefix '_' area suffix '_' condition...
        '_unit_location_and_phase_correlations_at_p01_p02_p03_p1_p3_1_2_4_10_Hz'] '.fig'])
    end
  end
end