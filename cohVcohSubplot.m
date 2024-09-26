function cohVcohSubplot(coh1, coh2, xStr, yStr, r, pval, nFOI, cohLim)

if nargin < 8
  cohLim = [0 1];
end

hold on
if sum(sum(~isnan(coh1))) && sum(sum(~isnan(coh2)))
  plot(coh1,coh2, '.', 'MarkerSize',10);
end

xLim = xlim;
xLim = [cohLim(1) min([cohLim(2) xLim(2)])];
xLim = [xLim(1)-0.001 max([xLim(2)+0.001 0.005])];
xlim(xLim);
yLim = ylim;
yLim = [cohLim(1) min([cohLim(2) yLim(2)])];
yLim = [yLim(1)-0.001 max([yLim(2)+0.001 0.005])];
ylim(yLim);
xAxisLength = xLim(2)-xLim(1);
yAxisLength = yLim(2)-yLim(1);
text(xLim(2)-xAxisLength*0.25, yLim(1)+yAxisLength*0.2, ['r=' num2str(r)], 'FontSize',16);
text(xLim(2)-xAxisLength*0.25, yLim(1)+yAxisLength*0.1, ['p=' num2str(pval)], 'FontSize',16);

text(xLim(1)+xAxisLength*0.025, yLim(2)-yAxisLength*0.05, ['n=' num2str(nFOI.bothHalvesSignificant) '/' num2str(nFOI.oneHalfSignificant) '/'...
  num2str(nFOI.total)], 'FontSize',16);

xTicks = xticks;
xTicks = xTicks(1:2:end);
yTicks = yticks;
yTicks = yTicks(1:2:end);

ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 24, 4/3, 2, [0.005 0], 'out',...
  'on', 'k', {xStr}, xLim, xTicks,...
  'on', 'k', {yStr}, yLim, yTicks);
hold off

position = get(ax1,'Position');
set(ax1,'Position', [position(1)-0.05 position(2) position(3)+0.03 position(4)+0.03]);