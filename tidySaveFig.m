function tidySaveFig(fH, yLabel, fLegend, figName, figFileName, legendPos)
% Figure tidying and saving. A helper function of plotUpdateSeries.

% Label the figure
figure(fH)
xlabel('Principal component');
ylabel(yLabel);
legend(fLegend.lines,fLegend.text, 'Interpreter','None', 'Location',legendPos);
title(figName);
set(gcf, 'Name',figName);

% Tidy the figure
ylimits = ylim;
xlimits = xlim;
set(gcf,'color','w');
legend boxoff
ax1 = axesProperties({}, 1, 'normal', 'off', 'w', 'Calibri', 12, 4/3, 2, [0.005 0], 'out',...
  'on', 'k', {'Principal component'}, xlimits, xticks,...
  'on', 'k', {yLabel}, ylimits, yticks);

% Save the figure in a fig file
hgsave(gcf, figFileName);

% Save the figure in a graphical file
label = [1.8 1.8];
margin = [0.3 0.3];
width = 1*15-label(1)-margin(1);
height = 1*15-label(2)-margin(2);
paperSize = resizeFig(gcf, ax1, width, height, label, margin, 0);
exportFig(gcf, [figFileName '.png'],'-dpng','-r300', paperSize);
end