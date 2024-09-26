function saveFigsStateMFR(initFigName, figMFR)
% A helper function to globalFigs for saving mean firing rate figures.

for i = 1:numel(figMFR)
  h = figure(figMFR(i)); %#ok<*NASGU>
  figFileName = [initFigName '_MFR'];
  figFileName = strrep(figFileName,'.','p');
  set(gcf, 'Name',figFileName);
  hgsave(gcf, figFileName);
  
%   ax1 = findobj(h, 'type','axe');
%   xLabel = get(get(ax1,'xlabel'),'string');
%   xTick = get(ax1,'xtick');
%   yLabel = get(get(ax1,'ylabel'),'string');
%   ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
%   	'on', 'k', {xLabel}, [], xTick,...
%     'on', 'k', {yLabel}, [], xTick);
%   label = [4 2.6];
%   margin = [1 1];
%   width = 1.3*15-label(1)-margin(1);
%   height = (1.3*15)-label(2)-margin(2);
%   paperSize = resizeFig(h, ax1, width, height, label, margin, 0);
%   exportFig(h, [figFileName '.png'],'-dpng','-r300', paperSize);
end