function figH = plotFiringRate(filename, series, period, ylimits)
% figH = figH = plotFiringRate(filename, series, period, ylimits)
% This function plots the pupil area signal for a given time interval
% Input: filename
%        series
%        period
%        ylimits
% Output: figH - figure handle.

if nargin < 4
  ylimits = [];
end
if nargin < 3
  period = [];
end
load(filename); %#ok<*LOAD>

if iscell(series)
  for iSeries = 1:numel(series)
    dbStruct = dataStruct.seriesData.(series{iSeries});
    if iSeries == 1
      PR = sum(dbStruct.popData.MUAsAll,1);
    else
      if numel(PR) > size(dbStruct.popData.MUAsAll,2)
        PR = PR + [sum(dbStruct.popData.MUAsAll,1) zeros(1,numel(PR)-size(dbStruct.popData.MUAsAll,2))];
      elseif numel(PR) < size(dbStruct.popData.MUAsAll,2)
        PR = [PR zeros(1,size(dbStruct.popData.MUAsAll,2)-numel(PR))] + sum(dbStruct.popData.MUAsAll,1);
      else
         PR = PR + sum(dbStruct.popData.MUAsAll,1);
      end
    end
  end
else
  dbStruct = dataStruct.seriesData.(series);
  PR = sum(dbStruct.popData.MUAsAll,1);
end
times = (1:numel(PR)).*(1/dbStruct.conf.samplingParams.srData);

figH = figure;
plot(times,PR, 'k', 'LineWidth',0.5);
if ~isempty(period)
  xlim(period);
else
  period = xlim;
end
if ~isempty(ylimits)
  ylim(ylimits);
else
  ylimits = ylim;
end

ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
  'on', 'k', {'Time (s)'}, period, [],...
  'on', 'k', {'Firing rate (spikes)'}, ylimits, []);
  
label = [0 0];
margin = [0 0];
width = 2*15-label(1)-margin(1);
height = 1*15-label(2)-margin(2);
paperSize = resizeFig(figH, ax1, width, height, label, margin, 0);
if iscell(series)
  [seriesName, animal] = seriesFromEntry(series);
  figFileName = [animal '_s' seriesName(1:14) '_firingRate'];
else
  figFileName = [series '_firingRate'];
end
hgsave(figH, [figFileName '.fig']);
exportFig(figH, [figFileName '.png'],'-dpng','-r300', paperSize);
close(figH);