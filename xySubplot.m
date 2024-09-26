function xySubplot(x, y, xStr, xTicks, xTickLabel, xLim, yStr, yTicks, yTickLabel, yLim, r, pval, nFOI, xZero, yZero, bestFitLine)
% xySubplot(x, y, xStr, xTicks, xTickLabel, yStr, yTicks, yTickLabel, r, pval, nFOI)
%
% Function produces 3x3 subplot figures.
% Input: x, y - data vectors to be plotted on x and y axes, respectively.
%        xStr - x-axis label.
%        xTicks - x-axis tick positions. If left empty, then original
%                 values will be used.
%        xTickLabel - x-axis tick labels. If left empty, then original
%                 labels will be used.
%        xLim - x-axis limits.
%        yStr - y-axis label.
%        yTicks - y-axis tick positions. If left empty, then original
%                 values will be used.
%        yTickLabel - y-axis tick labels. If left empty, then original
%                 labels will be used.
%        yLim - y-axis limits.
%        r - data correlation coefficient. If left empty, then it won't
%               appear in the figure.
%        pval - data correlation p-value. If left empty, then it won't
%               appear in the figure.
%        nFOI - a vector with significant recording numbers. Can have
%               several values. If left empty, then these numbers won't
%               appear in the figure.
%        xZero - true if you want to indicate the origin of the x-axis.
%                Default is false.
%        yZero - true if you want to indicate the origin of the y-axis.
%                Default is false.
%        bestFitLine - 'linear-linear', 'linear-circular' or 'none'.
%                      Default is none.

% Remove NaNs
inds = ~isnan(x) & ~isnan(y);
x = x(inds);
y = y(inds);

if nargin < 14
  bestFitLine = 'none';
end
if nargin < 13
  yZero = false;
end
if nargin < 12
  xZero = false;
end

[yFit, slope] = fitLine(x, y, bestFitLine);

% Plot the data
hold on
if sum(sum(~isnan(x))) && sum(sum(~isnan(y)))
  plot(x, y, '.', 'MarkerSize',10);
end

% Adjust and tidy the figure
if isempty(xLim)
  xLim = xlim;
end
xAxisLength = xLim(2)-xLim(1);
if isempty(xTicks)
xTicks = xticks;
xTicks = xTicks(1:2:end);
end
if isempty(yLim)
  yLim = [-pi pi];
end
yAxisLength = yLim(2)-yLim(1);
if isempty(yTicks)
  yTicks = [yLim(1) yLim(1)+0.25*yAxisLength yLim(1)+0.5*yAxisLength yLim(1)+0.75*yAxisLength yLim(2)];
end

if ~isempty(r)
  if slope < 0 && r > 0
    r = -r;
  end
  text(xLim(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.2, ['r=' num2str(r)], 'FontSize',16);
end
if ~isempty(pval)
  text(xLim(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.1, ['p=' num2str(pval)], 'FontSize',16);
end
if ~isempty(nFOI)
  textStr = ['n=' num2str(nFOI(1))];
  if numel(nFOI) > 1
    for iN = 2:numel(nFOI)
      textStr = [textStr '/' num2str(nFOI(iN))]; %#ok<*AGROW>
    end
  end
  text(xLim(1)+xAxisLength*0.025, yLim(2)-yAxisLength*0.05, textStr, 'FontSize',16);
end

ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
  'on', 'k', {xStr}, xLim, xTicks,...
  'on', 'k', {yStr}, yLim, yTicks);
if ~isempty(xTickLabel)
  ax1.XTickLabel = xTickLabel;
end
if ~isempty(yTickLabel)
  ax1.YTickLabel = yTickLabel;
end

% Mark the origins of axes
if xZero
  plot([0 0], ylim, 'k:');
end
if yZero
  plot(xlim, [0 0], 'k:');
end

% Draw the best fit line
if ~strcmpi(bestFitLine, 'none')
  plot(x, yFit, 'r-');
  if strcmpi(bestFitLine, 'linear-circular')
    plot(x, yFit+2*pi, 'r-');
    plot(x, yFit-2*pi, 'r-');
  end
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