function h = phaseHistPlotUnadjusted(edges, histPhaseRelation, yLabel, figFileName, visibility, adjustXaxis, figSize)
% h = phaseHistPlotUnadjusted(edges, histPhaseRelation, yLabel, figFileName, visibility, adjustXaxis, figSize)
%
% A helper function to display phase frequency histograms.

if nargin < 7
  figSize = 15;
end
if nargin < 6
  adjustXaxis = false;
end

h = figProperties(figFileName, 'normalized', [0, .005, .97, .90], 'w', visibility);
hold on

if numel(histPhaseRelation) == numel(edges)
  yTick = ceil(max(histPhaseRelation(2:end))/2);
  if histPhaseRelation(1) > max(histPhaseRelation(2:end))
    text(1, 0.2*yTick, num2str(histPhaseRelation(1)), 'FontSize',25)
  else
    bar(2, histPhaseRelation(1), 'k');
  end
  for i = 2:numel(histPhaseRelation)
  	bar(7+i-1, histPhaseRelation(i), 'k');
  end
  ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
  	'on', 'k', {'Phase (rad)'}, [0 24], [2 7.5 11.5 15.5 19.5 23.5],...
    'on', 'k', {yLabel}, [0 yTick*2+1], 0:yTick:yTick*2);
  if adjustXaxis == 0
    ax1.XTickLabel = {'non-signif.','-\pi','-\pi/2','0','\pi/2','\pi'};
  elseif adjustXaxis == pi/2
    ax1.XTickLabel = {'non-signif.','-\pi/2','0','\pi/2','\pi','3\pi/2'};
  elseif adjustXaxis == 3*pi/4
    ax1.XTickLabel = {'non-signif.','-\pi/4','\pi/4','3\pi/4','5\pi/4','7\pi/4'};
  elseif adjustXaxis == pi
    ax1.XTickLabel = {'non-signif.','0','\pi/2','\pi','3\pi/2','2\pi'};
  elseif isnumeric(adjustXaxis)
    error(['The figure is not defined for an x-axis centred on ' num2str(adjustXaxis) ' radians']);
  else
    error('The figure is not defined for an x-axis centred on a non numeric value');
  end
  
elseif numel(histPhaseRelation)+1 == numel(edges)
  yTick = ceil(max(histPhaseRelation)/2);
  for i = 1:numel(histPhaseRelation)
    bar(i, histPhaseRelation(i), 'k');
  end
  ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
  	'on', 'k', {'Phase (rad)'}, [0 17], [7.5 11.5 15.5 19.5 23.5]-7,...
    'on', 'k', {yLabel}, [0 yTick*2+1], 0:yTick:yTick*2);
  if adjustXaxis == 0
    ax1.XTickLabel = {'-\pi','-\pi/2','0','\pi/2','\pi'};
  elseif adjustXaxis == pi/2
    ax1.XTickLabel = {'-\pi/2','0','\pi/2','\pi','3\pi/2'};
  elseif adjustXaxis == 3*pi/4
    ax1.XTickLabel = {'-\pi/4','\pi/4','3\pi/4','5\pi/4','7\pi/4'};
  elseif adjustXaxis == pi
    ax1.XTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
  elseif isnumeric(adjustXaxis)
    error(['The figure is not defined for an x-axis centred on ' num2str(adjustXaxis) ' radians']);
  else
    error('The figure is not defined for an x-axis centred on a non numeric value');
  end
  
else
  error('Supplied histogram edges do not correspond to the data')
end
hold off
    
label = [4 2.65];
margin = [0.7 1];
width = 1.3*figSize-label(1)-margin(1);
height = (1.3*figSize)-label(2)-margin(2);
paperSize = resizeFig(h, ax1, width, height, label, margin, 0);
hgsave(h, [figFileName '.fig']);
exportFig(h, [figFileName '.png'],'-dpng','-r300', paperSize);
close(h);
end