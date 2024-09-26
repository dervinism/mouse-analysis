function dispFOI(f, n)
% dispFOI(f, n)
%
% Function displays count and frequency of interest in the upper right-hand
% corner of a figure.

xLim = xlim;
xAxisLength = xLim(2)-xLim(1);
yLim = ylim;
yAxisLength = yLim(2)-yLim(1);

if nargin < 2
  textStr = [num2str(f) 'Hz'];
  text(xLim(2)-xAxisLength*0.15, yLim(2)-yAxisLength*0.05, textStr, 'FontSize',16);
else
  textStr = ['n=' num2str(n) ' ' num2str(f) 'Hz'];
  text(xLim(2)-xAxisLength*0.3, yLim(2)-yAxisLength*0.05, textStr, 'FontSize',16);
end