function h = cohHistPlot(edges,histPhaseRelation,xTickLabel,yLabel,figFileName,visibility)
% A helper function of globalFigs for plotting coherence relation histograms.

h = figProperties(figFileName, 'normalized', [0, .005, .97, .90], 'w', visibility);
hold on

yTick = ceil(max(histPhaseRelation)/2);
if numel(histPhaseRelation) == numel(edges)
  bar(2, histPhaseRelation(1), 'k');
  for i = 2:numel(histPhaseRelation)
  	bar(7+i-1, histPhaseRelation(i), 'k');
  end
  ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
  	'on', 'k', {'Coherence'}, [0 28], [2 7.5 12.5 17.5 22.5 27.5],...
    'on', 'k', {yLabel}, [0 yTick*2+1], 0:yTick:yTick*2);
  ax1.XTickLabel = xTickLabel;
  
elseif numel(histPhaseRelation)+1 == numel(edges)
  for i = 1:numel(histPhaseRelation)
    bar(i, histPhaseRelation(i), 'k');
  end
  ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 25, 4/3, 2, [0.005 0], 'out',...
  	'on', 'k', {'Coherence'}, [0 21], [7.5 12.5 17.5 22.5 27.5]-7,...
    'on', 'k', {yLabel}, [0 yTick*2+1], 0:yTick:yTick*2);
  ax1.XTickLabel = xTickLabel;
  
else
  error('Supplied histogram edges do not correspond to the data')
end
    
label = [4 2.6];
margin = [0 1];
width = 1.3*15-label(1)-margin(1);
height = 1.3*15-label(2)-margin(2);
paperSize = resizeFig(h, ax1, width, height, label, margin, 0);
exportFig(h, [figFileName '.png'],'-dpng','-r300', paperSize);
close(h);
end