function saveFigsState(initFigName, figPhase, figCoherence, FOI)
% A helper function to compareHalves_figs, compareSeriesMultiple,
% compareSeries_figs, eyeAnalysis_figs, and globalFigs for saving phase and
% coherence figures.

for i = 1:numel(figPhase)
  h = figure(figPhase(i)); %#ok<*NASGU>
  figFileName = [initFigName '_phase_' num2str(round(FOI(i),2)) 'Hz'];
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
  close(h);
  
  h = figure(figCoherence(i));
  figFileName = [initFigName '_coherence_' num2str(round(FOI(i),2)) 'Hz'];
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
  close(h);
end