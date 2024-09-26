function ax1 = phaseHistSubplotUnadjusted(edges, phaseHist, yLabel, xLabel, colour, adjustXaxis)
% [ax1, position] = phaseHistSubplotUnadjusted(edges, phaseHist, yLabel, xLabel, colour)
%
% Function produces phase histogram subplots.
% Input: edges - a vector specifying the edges of histogram bins (used in
%                histcounts).
%        phaseHist - a phase histogram vector. The first entry could be the
%                    count of non-significant phase values.
%        yLabel - the y axis label.
%        xLabel - the x axis label.
%        colour - colour code for drawing histogram bars.
%        adjustXaxis.
% Output: ax1 - axes handle.
%         position - position of axes.

if nargin < 6
  adjustXaxis = 0;
end
if nargin < 5
  colour = 'k';
end

hold on
if numel(phaseHist) == numel(edges)
  yTick = ceil(max(phaseHist(2:end))/2);
  if phaseHist(1) > max(phaseHist(2:end))
    text(1, max([0.3 0.2*yTick]), num2str(phaseHist(1)), 'FontSize',25)
  else
    bar(2, phaseHist(1), 'FaceColor',colour, 'EdgeColor',colour);
  end
  for i = 2:numel(phaseHist)
  	bar(7+i-1, phaseHist(i), 'FaceColor',colour, 'EdgeColor',colour);
  end
  ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
  	'on', 'k', {xLabel}, [0 24], [2 7.5 11.5 15.5 19.5 23.5],...
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
  
elseif numel(phaseHist)+1 == numel(edges)
  yTick = ceil(max(phaseHist)/2);
  for i = 1:numel(phaseHist)
    bar(i, phaseHist(i), 'FaceColor',colour, 'EdgeColor',colour);
  end
  ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
  	'on', 'k', {xLabel}, [0 17], [7.5 11.5 15.5 19.5 23.5]-7,...
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

% Position the subplot
figPosition = get(gcf,'OuterPosition');
axesPosition = get(ax1,'Position');
if figPosition(4) <= 0.75
  set(ax1,'Position', [axesPosition(1)-0.05 axesPosition(2)+0.04 axesPosition(3)+0.03 axesPosition(4)+0.03]);
elseif figPosition(4) <= 0.9
  set(ax1,'Position', [axesPosition(1)-0.05 axesPosition(2) axesPosition(3)+0.03 axesPosition(4)+0.03]);
elseif figPosition(4) > 0.9
  if axesPosition(2) >= 2/3
    set(ax1,'Position', [axesPosition(1)-0.05 axesPosition(2)+0.03 axesPosition(3)+0.03 axesPosition(4)+0.03]);
  else
    set(ax1,'Position', [axesPosition(1)-0.05 axesPosition(2) axesPosition(3)+0.03 axesPosition(4)+0.03]);
  end
end
set(ax1, 'Units','centimeters');
axesPosition = get(ax1,'Position');
set(ax1,'Position', [axesPosition(1) axesPosition(2) axesPosition(3) 6.1116]);
set(ax1, 'Units','normalized');