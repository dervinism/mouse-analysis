function phaseVphaseSubplot(phase1, phase2, xStr, yStr, r, pval, nFOI, phaseLim)

if nargin < 8
  phaseLim = [-pi pi];
end

hold on
if sum(sum(~isnan(phase1))) && sum(sum(~isnan(phase2)))
  plot(phase1,phase2, '.', 'MarkerSize',10);
end

xLim = phaseLim;
yLim = phaseLim;
xAxisLength = xLim(2)-xLim(1);
yAxisLength = yLim(2)-yLim(1);
text(xLim(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.2, ['r=' num2str(r)], 'FontSize',16);
text(xLim(2)-xAxisLength*0.3, yLim(1)+yAxisLength*0.1, ['p=' num2str(pval)], 'FontSize',16);

text(xLim(1)+xAxisLength*0.025, yLim(2)-yAxisLength*0.05, ['n=' num2str(nFOI.bothHalvesSignificant) '/' num2str(nFOI.oneHalfSignificant) '/'...
  num2str(nFOI.total)], 'FontSize',16);

ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 24, 4/3, 2, [0.005 0], 'out',...
  'on', 'k', {xStr}, phaseLim, [-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi],...
  'on', 'k', {yStr}, phaseLim, [-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi]);
ax1.XTickLabel = {'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'};
ax1.YTickLabel = {'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'};
hold off

position = get(ax1,'Position');
set(ax1,'Position', [position(1)-0.05 position(2) position(3)+0.03 position(4)+0.03]);